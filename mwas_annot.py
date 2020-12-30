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

indices = ['Species', 'Contig', 'Position']

first_columns = ['SGB', 'contig', 'Contig_with_part',
                 'MajorAllele', 'MinorAllele', 'MajorCodon', 'MinorCodon', 'MajorAA', 'MinorAA',
                 'GeneID', 'GeneRelation', 'GeneDistance', 'strand', 'start_pos', 'end_pos',
                 'feature', 'gene', 'product', 'eggNOG annot']

annotations_file = ANNOTS_03_2020_EXPANDED
annotations = None


def load_annotations():
    global annotations
    if annotations is None:
        annotations = GeneAnnotLoader(annotations_file).get_annot()
        # snps is 0 based while the annotations are 1 based
        annotations['start_pos'] = annotations['start_pos'] - 1
        annotations['end_pos'] = annotations['end_pos'] - 1


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
    y_species = list(set(['SGB' + file.split('SGB')[-1].split('.')[0] for file in y_mwas_files]))

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
    snps.columns = indices
    snps = snps.set_index(indices)

    pvals_df.to_csv(os.path.join(output_path, 'pvals_count.csv'))
    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_path, 'snps_unique.h5'), key='snps')

    return snps


def add_alleles(snps, P, output_path):
    # adds major and minor allele to each position

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
    snps[['MajorAllele', 'MinorAllele']] = snps[['MajorAllele', 'MinorAllele']].astype(int)

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_path, 'snps_alleles.h5'), key='snps')

    return snps


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
    snps = snps.reset_index()
    snps = snps.rename(columns={'Contig': 'Contig_with_part'})
    snps['Contig'] = remove_contig_parts(snps['Contig_with_part'].values)
    snps = snps.set_index(indices)

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_path, 'snps_contig.h5'), key='snps')

    return snps


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

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_path, 'snps_maf_genes.h5'), key='snps')

    return snps


def add_surrounding_genes(snps, output_path):
    # add multiple genes

    load_annotations()
    snps = Annotationist(annotations, load_sequences=False).lookup_df(snps)  # current/upstream/downstream, +/- strand

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_path, 'snps_surrounding_genes.h5'), key='snps')

    return snps


def flatten_surrounding_genes(snps, output_path):
    # multiply cases where you have more than one gene

    snps = snps.rename(columns={'PlusCurrentGenePos': 'PlusCurrentGeneDistance',
                                'MinusCurrentGenePos': 'MinusCurrentGeneDistance'})

    surrounding_columns = ['PlusCurrentGeneID', 'PlusCurrentGeneDistance',
                           'MinusCurrentGeneID', 'MinusCurrentGeneDistance',
                           'PlusUpstreamGeneID', 'PlusUpstreamGeneDistance',
                           'MinusUpstreamGeneID', 'MinusUpstreamGeneDistance',
                           'PlusDownstreamGeneID', 'PlusDownstreamGeneDistance',
                           'MinusDownstreamGeneID', 'MinusDownstreamGeneDistance']

    dfs = []
    for relation in ['Current', 'Upstream', 'Downstream']:

        df = snps.copy()
        df = df.drop([col for col in surrounding_columns if relation not in col], axis=1)

        # no gene
        no_gene = df[[f'Plus{relation}GeneID', f'Minus{relation}GeneDistance']].isna().sum(axis=1) == 2
        no_gene = df.loc[no_gene].drop([f'Minus{relation}GeneID', f'Minus{relation}GeneDistance'], axis=1)  # could have been plus
        no_gene = no_gene.rename(columns={f'Plus{relation}GeneID': 'GeneID', f'Plus{relation}GeneDistance': 'GeneDistance'})

        # plus genes
        plus_genes = df[[f'Plus{relation}GeneID', f'Plus{relation}GeneDistance']].dropna().index
        plus_genes = df.loc[plus_genes].drop([f'Minus{relation}GeneID', f'Minus{relation}GeneDistance'], axis=1)
        plus_genes = plus_genes.rename(columns={f'Plus{relation}GeneID': 'GeneID', f'Plus{relation}GeneDistance': 'GeneDistance'})

        # minus genes
        minus_genes = df[[f'Minus{relation}GeneID', f'Minus{relation}GeneDistance']].dropna().index
        minus_genes = df.loc[minus_genes].drop([f'Plus{relation}GeneID', f'Plus{relation}GeneDistance'], axis=1)
        minus_genes = minus_genes.rename(columns={f'Minus{relation}GeneID': 'GeneID', f'Minus{relation}GeneDistance': 'GeneDistance'})

        df = pd.concat([plus_genes, minus_genes, no_gene])
        df['GeneRelation'] = relation
        dfs.append(df)

    snps = pd.concat(dfs).reset_index()
    snps = snps.drop_duplicates(snps.columns[:-1]).set_index(indices).sort_index()
    snps.loc[snps['GeneID'].isna(), 'GeneRelation'] = ''
    snps['GeneID'] = snps['GeneID'].replace(np.nan, None).astype(int)

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_path, 'snps_surrounding_genes_flat.h5'), key='snps')

    return snps


def add_annotations(snps, output_path, on='GeneID', rsuffix=''):
    # add gene annotations

    original_columns = snps.columns

    load_annotations()
    snps = snps.join(annotations, on=on, rsuffix=rsuffix)  # TODO: take only useful columns, rename unclear ones
    snps = snps.drop(['SGB', 'contig'], axis=1)
    new_columns = list(set(snps.columns) - set(original_columns) - set(['start_pos', 'end_pos']))
    snps[new_columns] = snps[new_columns].astype(str)

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_path, 'snps_annotations.h5'), key='snps')

    return snps


def add_codons(snps, P, output_path):
    # adds codons

    # match column names to what is used in the MBSNPLoader
    snps = snps.rename(columns={'MajorAllele': 'Major', 'MinorAllele': 'Minor'})

    snps['MajorCodonMAF'] = None
    snps['MinorCodonMAF'] = None

    load_annotations()
    for species in snps.index.get_level_values('Species').unique():

        species_genes = annotations[annotations['SGB'] == int(species.replace('SGB_', ''))]
        gene_id_2_idx = {gene_id: i for i, gene_id in enumerate(species_genes.index.values)}

        # in case annotation file is 0 based x[2] -1, x[3] or if is 1 based x[2], x[3] + 1
        gp = [GenePosition(*(list(['SGB_{}'.format(x[0]), 'C_{}'.format(x[1]), x[2], x[3]+1, x[4]])))
              for x in species_genes[['SGB', 'contig', 'start_pos', 'end_pos', 'strand']].values]
        # in gp, start_pos is already corrected for zero based

        condition_species_current = (snps.index.get_level_values('Species') == species) & (snps['GeneRelation'] == 'Current')
        species_current_df = snps[condition_species_current].copy()
        # updates species_current_df
        get_mbsnp_loader_class(P.body_site)._add_codon_annot(species_current_df, gene_id_2_idx, gp, None)
        snps[condition_species_current] = species_current_df

        snps[['MajorCodonMAF', 'MinorCodonMAF']] = snps[['MajorCodonMAF', 'MinorCodonMAF']].fillna(_MISSING_DATA_NOT_IN_GENE)

    # revert back column names
    snps = snps.rename(columns={'Major': 'MajorAllele', 'Minor': 'MinorAllele',
                                'MajorCodonMAF': 'MajorCodon', 'MinorCodonMAF': 'MinorCodon'})

    is_coding = snps['MajorCodon'] >= 0

    snps.loc[is_coding, 'MajorCodon'] = snps.loc[is_coding, 'MajorCodon'].map(codons_inv)
    snps.loc[is_coding, 'MinorCodon'] = snps.loc[is_coding, 'MinorCodon'].map(codons_inv)

    snps[['MajorCodon', 'MinorCodon']] = snps[['MajorCodon', 'MinorCodon']].astype(str)

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_path, 'snps_codons.h5'), key='snps')

    return snps


def add_amino_acids(snps, output_path):
    # adds amino acids

    is_coding = [len(value) == 3 for value in snps['MajorCodon']]

    snps.loc[is_coding, 'MajorAA'] = snps.loc[is_coding, 'MajorCodon'].apply(translate_codon)
    snps.loc[is_coding, 'MinorAA'] = snps.loc[is_coding, 'MinorCodon'].apply(translate_codon)

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_path, 'snps_amino_acids.h5'), key='snps')

    return snps


def add_taxonomy(snps, output_path):
    # adds taxonomy

    tax_df = taxonomy_df(level_as_numbers=False).set_index('SGB')[
        ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]
    snps = snps.join(tax_df, on='Species')

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_path, 'snps_taxonomy.h5'), key='snps')

    return snps


def run(P, output_path, x_mwas_files_path, y_mwas_files_path):

    snps = find_unique_snps(x_mwas_files_path, y_mwas_files_path, output_path)  # TODO: handel the regular mwas case
    snps = add_alleles(snps, P, output_path)  # early because it uses contig with part
    snps = add_contig_without_part(snps, output_path)
    snps = add_surrounding_genes(snps, output_path)
    snps = flatten_surrounding_genes(snps, output_path)
    snps = add_annotations(snps, output_path)
    snps = add_codons(snps, P, output_path)  # late because it depends on column 'feature' from annotations
    snps = add_amino_acids(snps, output_path)

    return snps
