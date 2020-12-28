import re
import os
import glob
import numpy as np
import pandas as pd
from pandas import HDFStore
from typing import List, Iterable
from LabData.RawData.genome import GenePosition
from LabData.config_global import mb_pipeline_dir
from LabUtils.SeqUtils import translate_codon, codons_inv
from LabData.DataAnalyses.MBSNPs.taxonomy import taxonomy_df
from UseCases.DataAnalyses.mbsnps_annotation import Annotationist
from LabData.DataLoaders.GeneAnnotLoader import GeneAnnotLoader, ANNOTS_03_2020_EXPANDED
from LabData.DataLoaders.MBSNPLoader import get_mbsnp_loader_class, _MISSING_DATA_NOT_IN_GENE


first_columns = ['Contig_with_part', 'MajorAllele', 'MinorAllele', 'MajorCodon', 'MinorCodon', 'MajorAA', 'MinorAA',
                 'GeneID', 'PositionInGene', 'SGB', 'contig', 'start_pos', 'end_pos', 'strand', 'feature', 'gene', 'product']
annotations_file = ANNOTS_03_2020_EXPANDED
annotations = None


def load_annotations():
    global annotations
    if annotations is None:
        annotations = GeneAnnotLoader(annotations_file).get_annot()


def order_columns(snps):
    # order the columns by the predefined first columns
    first_cols = [col for col in first_columns if col in snps.columns]
    other_cols = [col for col in snps.columns if col not in first_columns]

    return snps[first_cols + other_cols]


def find_unique_snps(x_mwas_files_path, y_mwas_files_path, output_path):
    # find unique snps and count comparisons (for p_value calculations)

    x_mwas_files = glob.glob(x_mwas_files_path)
    y_mwas_files = glob.glob(y_mwas_files_path)

    x_species = list(set(['SGB' + file.split('SGB')[-1].split('.')[0] for file in x_mwas_files]))
    y_species = list(set([os.path.basename(file).split('.')[0] for file in y_mwas_files]))

    # pre run
    snps = set()
    pvals_df = pd.DataFrame(index=y_species, columns=x_species)

    # run
    for x_mwas_file in x_mwas_files:
        x_species = 'SGB' + x_mwas_file.split('SGB')[-1].split('.')[0]
        with HDFStore(x_mwas_file, 'r') as x_mwas_df:
            x_mwas_df = x_mwas_df[x_species]

        y_species, y_species_count = np.unique(x_mwas_df.index.get_level_values('Y'), return_counts=True)
        pvals_df.loc[y_species.tolist(), x_species] = y_species_count.tolist()

        snps = snps.union(set(x_mwas_df.index.droplevel('Y').values))

    # post run
    snps = pd.DataFrame(snps)
    snps.columns = ['Species', 'Contig', 'Position']
    snps = snps.set_index(['Species', 'Contig', 'Position'])
    pvals_df.to_csv(os.path.join(output_path, 'pvals_count.csv'))
    snps.to_hdf(os.path.join(output_path, 'snps_unique.h5'), key='snps')

    return order_columns(snps)


def add_alleles(snps, P, output_path):
    # adds major and minor allele to each position

    indices = ['Species', 'Contig', 'Position']
    columns = ['Major', 'Minor']

    alleles = pd.DataFrame(columns=indices + columns).set_index(indices)
    for species in snps.index.get_level_values('Species').unique():
        maf_path = os.path.join(mb_pipeline_dir, 'Analyses', 'MBSNP', P.body_site, 'MAF1', '{}_{}_R{}_S{}.h5'.
                                format('mb_snp_g_maf', species, P.min_reads_per_snp, P.min_subjects_per_snp_cached))
        if os.path.exists(maf_path):
            contigs = snps.loc[species].index.get_level_values('Contig').unique()
            with pd.HDFStore(maf_path, 'r') as hdf:
                for contig in contigs:
                    contig_df = hdf[contig][columns]
                    contig_df['Species'] = species
                    contig_df['Contig'] = contig
                    contig_df = contig_df.reset_index().set_index(indices)
                    alleles = alleles.append(contig_df)

    snps = snps.join(alleles)

    snps = snps.rename(columns={'Major': 'MajorAllele', 'Minor': 'MinorAllele'})

    snps.to_hdf(os.path.join(output_path, 'snps_alleles.h5'), key='snps')

    return order_columns(snps)


def add_contig_without_part(snps, output_path):

    def remove_contig_parts(contigs: Iterable[str]) -> List[str]:
        """Returns the given contig names without the trailing '_P###'."""
        contig_part_re = re.compile('_P\\d+$')
        # Maps old contig name to new contig name, to avoid repetition.
        cache = {}

        def cached_name(contig):
            if contig not in cache:
                cache[contig] = contig_part_re.sub('', contig)
            return cache[contig]

        return [cached_name(x) for x in contigs]

    # Remove the part from the contig key
    idx = snps.index.names
    snps = snps.reset_index()
    snps = snps.rename(columns={'Contig': 'Contig_with_part'})
    snps['Contig'] = remove_contig_parts(snps['Contig_with_part'].values)
    snps = snps.set_index(idx)

    snps.to_hdf(os.path.join(output_path, 'snps_contig.h5'), key='snps')

    return order_columns(snps)


def add_maf_genes_annotations(snps, P, output_path):
    # add annotations from maf_annot file

    def _add_snp_annotations(species, df):
        maf_annot_fname = get_mbsnp_loader_class(P.body_site).\
            get_snp_maf_annot_fname(species, P.min_reads_per_snp, P.min_subjects_per_snp_cached)
        if os.path.exists(maf_annot_fname):
            with HDFStore(maf_annot_fname, 'r') as maf_annot:
                maf_annot = maf_annot[species]
            df = df.join(maf_annot)
        return df

    snps = pd.concat([_add_snp_annotations(species, df) for species, df in snps.groupby('Species')])

    snps.to_hdf(os.path.join(output_path, 'snps_maf_genes_info.h5'), key='snps')

    return order_columns(snps)


def add_surrounding_genes(snps, output_path):
    # add multiple genes

    load_annotations()
    snps = Annotationist(annotations, load_sequences=False).lookup_df(snps)  # current/upstream/downstream, +/- strand
    # TODO: clean up, rename and order annotations columns
    snps.to_hdf(os.path.join(output_path, 'snps_surrounding_genes.h5'), key='snps')

    return order_columns(snps)


def flatten_current_genes(snps, output_path):
    # in case of multiple genes take only the current ones and
    # multiply cases where you have a gene on both the plus and the minus strands
    # TODO: change to stack and include all genes
    snps = snps.drop(['PlusUpstreamGeneID', 'PlusUpstreamGeneDistance',
                      'PlusDownstreamGeneID', 'PlusDownstreamGeneDistance',
                      'MinusUpstreamGeneID', 'MinusUpstreamGeneDistance',
                      'MinusDownstreamGeneID', 'MinusDownstreamGeneDistance'], axis=1)

    # no gene
    no_gene = snps[['PlusCurrentGeneID', 'MinusCurrentGeneID']].isna().sum(axis=1) == 2
    no_gene = snps.loc[no_gene].drop(['MinusCurrentGeneID', 'MinusCurrentGenePos'], axis=1)  # could have been plus
    no_gene = no_gene.rename(columns={'PlusCurrentGeneID': 'GeneID', 'PlusCurrentGenePos': 'PositionInGene'})

    # plus genes
    plus_genes = snps[['PlusCurrentGeneID', 'PlusCurrentGenePos']].dropna().index
    plus_genes = snps.loc[plus_genes].drop(['MinusCurrentGeneID', 'MinusCurrentGenePos'], axis=1)
    plus_genes = plus_genes.rename(columns={'PlusCurrentGeneID': 'GeneID', 'PlusCurrentGenePos': 'PositionInGene'})

    # minus genes
    minus_genes = snps[['MinusCurrentGeneID', 'MinusCurrentGenePos']].dropna().index
    minus_genes = snps.loc[minus_genes].drop(['PlusCurrentGeneID', 'PlusCurrentGenePos'], axis=1)
    minus_genes = minus_genes.rename(columns={'MinusCurrentGeneID': 'GeneID', 'MinusCurrentGenePos': 'PositionInGene'})

    snps = pd.concat([plus_genes, minus_genes, no_gene]).sort_index()

    snps.to_hdf(os.path.join(output_path, 'snps_current_genes.h5'), key='snps')

    return order_columns(snps)


def add_annotations(snps, output_path, on='GeneID', rsuffix=''):
    # add gene annotations

    load_annotations()
    snps = snps.join(annotations, on=on, rsuffix=rsuffix)

    snps.to_hdf(os.path.join(output_path, 'snps_annotations.h5'), key='snps')

    return order_columns(snps)


def add_codons(snps, P, output_path):
    # adds codons

    # match column names to what is used in the MBSNPLoader
    snps = snps.rename(columns={'MajorAllele': 'Major', 'MinorAllele': 'Minor', 'PositionInGene': 'GeneDistance'})

    all_species = []
    for species_name, species_df in snps.groupby('Species'):

        load_annotations()
        species_genes = annotations[annotations['SGB'] == int(species_name.replace('SGB_', ''))]
        gene_id_2_idx = {gene_id: i for i, gene_id in enumerate(species_genes.index.values)}

        gp = [GenePosition(*(list(['SGB_{}'.format(x[0]), 'C_{}'.format(x[1]), x[2] - 1]) + list(x[3:])))
              for x in species_genes[['SGB', 'contig', 'start_pos', 'end_pos', 'strand']].values]
        # in gp, start_pos is already corrected for zero based

        get_mbsnp_loader_class(P.body_site)._add_codon_annot(species_df, gene_id_2_idx, gp)  # updates species_df
        all_species.append(species_df)

    snps = pd.concat(all_species).fillna(_MISSING_DATA_NOT_IN_GENE)

    # revert back column names
    snps = snps.rename(columns={'Major': 'MajorAllele', 'Minor': 'MinorAllele',
                                'MajorCodonMAF': 'MajorCodon', 'MinorCodonMAF': 'MinorCodon',
                                'GeneDistance': 'PositionInGene'})

    is_coding = snps['MajorCodon'] >= 0

    snps.loc[is_coding, 'MajorCodon'] = snps.loc[is_coding, 'MajorCodon'].map(codons_inv)
    snps.loc[is_coding, 'MinorCodon'] = snps.loc[is_coding, 'MinorCodon'].map(codons_inv)

    snps.to_hdf(os.path.join(output_path, 'snps_codons.h5'), key='snps')

    return order_columns(snps)


def add_amino_acids(snps, output_path):
    # adds amino acids

    is_coding = [type(value) == str for value in snps['MajorCodon']]

    snps.loc[is_coding, 'MajorAA'] = snps.loc[is_coding, 'MajorCodon'].apply(translate_codon)
    snps.loc[is_coding, 'MinorAA'] = snps.loc[is_coding, 'MinorCodon'].apply(translate_codon)

    snps.to_hdf(os.path.join(output_path, 'snps_amino_acids.h5'), key='snps')

    return order_columns(snps)


def add_taxonomy(snps, output_path):

    tax_df = taxonomy_df(level_as_numbers=False).set_index('SGB')[
        ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]
    snps = snps.join(tax_df, on='Species')

    snps.to_hdf(os.path.join(output_path, 'snps_taxonomy.h5'), key='snps')

    return order_columns(snps)
