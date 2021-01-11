import re
import os
import glob
import numpy as np
import pandas as pd
from pandas import HDFStore
from typing import List, Iterable
from LabUtils.SeqUtils import translate_codon
from LabData.DataAnalyses.MBSNPs.taxonomy import taxonomy_df
from LabData.DataLoaders.MBSNPLoader import get_mbsnp_loader_class
from LabData.DataAnalyses.MBSNPs.annotationist import Annotationist
from LabData.DataLoaders.GeneAnnotLoader import GeneAnnotLoader, ANNOTS_03_2020_EXPANDED


indices = ['Species', 'Contig', 'Position']

first_columns = ['SGB', 'contig', 'Contig_with_part',
                 'MajorAllele', 'MinorAllele', 'MajorCodon', 'MinorCodon', 'MajorAA', 'MinorAA',
                 'GeneID', 'GeneRelation', 'GeneDistance', 'strand', 'start_pos', 'end_pos',
                 'feature', 'gene', 'product', 'eggNOG annot']

annotations_file = ANNOTS_03_2020_EXPANDED
gene_annotations = None
annotations_list = None


def load_gene_annotations():
    global gene_annotations
    if gene_annotations is None:
        gene_annotations = GeneAnnotLoader(annotations_file).get_annot()
        # snps is 0 based while the annotations are 1 based
        # gene_annotations['start_pos'] = gene_annotations['start_pos'] - 1
        # gene_annotations['end_pos'] = gene_annotations['end_pos'] - 1


def load_annotations_list():
    load_gene_annotations()
    global annotations_list
    if annotations_list is None:
        load_gene_annotations()
        annotations_list = Annotationist(gene_annotations, load_sequences=False)


def order_columns(snps):
    # order the columns by the predefined first columns

    first_cols = [col for col in first_columns if col in snps.columns]
    other_cols = [col for col in snps.columns if col not in first_columns]

    return snps[first_cols + other_cols]


def find_unique_snps(x_mwas_files_path, y_mwas_files_path, output_path, pval_col='Global_Bonferroni', alpha=0.05):
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

        snps = snps.union(set(x_mwas_df[x_mwas_df[pval_col] < alpha].index.droplevel('Y').values))

    # post run
    snps = pd.DataFrame(snps)
    snps.columns = indices
    snps = snps.set_index(indices)

    pvals_df.to_csv(os.path.join(output_path, 'pvals_count.csv'))
    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_path, 'snps_unique.h5'), key='snps')

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

    load_annotations_list()
    snps = annotations_list.lookup_df(snps)  # current/upstream/downstream, +/- strand

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

    load_gene_annotations()
    snps = snps.join(gene_annotations, on=on, rsuffix=rsuffix)  # TODO: take only useful columns, rename unclear ones
    snps = snps.drop(['SGB', 'contig'], axis=1)
    new_columns = list(set(snps.columns) - set(original_columns) - set(['start_pos', 'end_pos']))
    snps[new_columns] = snps[new_columns].astype(str)

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_path, 'snps_annotations.h5'), key='snps')

    return snps


def add_codons(snps, P, output_path):
    # adds codons

    # match column names to what is used in the MBSNPLoader
    # snps = snps.rename(columns={'MajorAllele': 'Major', 'MinorAllele': 'Minor'})

    load_annotations_list()
    snps = annotations_list.assign_codons(snps)

    codon_cols = ['PlusCurrentMajorCodon', 'PlusCurrentMinorCodon', 'MinusCurrentMajorCodon', 'MinusCurrentMinorCodon']
    snps[codon_cols] = snps[codon_cols].astype(str)

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


def run(P, output_path, x_mwas_files_path, y_mwas_files_path=None):

    if y_mwas_files_path:
        snps = find_unique_snps(x_mwas_files_path, y_mwas_files_path, output_path, pval_col='Pval', alpha=0.05)
    else:
        snps = pd.read_hdf(x_mwas_files_path)[[]]
    # snps = add_alleles(snps, P, output_path)  # early because it uses contig with part
    snps = add_contig_without_part(snps, output_path)
    snps = add_surrounding_genes(snps, output_path)
    snps = snps.reset_index()
    snps = snps.rename(columns={'Contig': 'Contig_without_part'})
    snps['Contig'] = snps['Contig_with_part'].values
    snps = snps.set_index(indices)
    snps = add_codons(snps, P, output_path)  # late because it depends on column 'feature' from annotations
    # snps = flatten_surrounding_genes(snps, output_path)
    # snps = add_annotations(snps, output_path)
    # snps = add_amino_acids(snps, output_path)

    return snps
