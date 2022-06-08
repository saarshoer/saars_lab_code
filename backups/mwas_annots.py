import re
import os
import glob
import numpy as np
import pandas as pd
from pandas import HDFStore
from typing import List, Iterable
from LabUtils.SeqUtils import translate_codon
from LabData.DataAnalyses.MBSNPs.taxonomy import taxonomy_df
from LabData.RawData.genome import read_sequences, GenePosition
from LabData.DataAnalyses.MBSNPs.annotationist import Annotationist
from LabData.DataLoaders.GeneAnnotLoader import GeneAnnotLoader
from Bio.Seq import Seq
from LabData import config_global as config

indices = ['Species', 'Contig', 'Position']

first_columns = ['SGB', 'contig', 'Contig_with_part', 'Contig_without_part',
                 'MajorAllele', 'MinorAllele', 'MajorCodon', 'MinorCodon', 'MajorAA', 'MinorAA',
                 'GeneID', 'GeneRelation', 'GeneDistance', 'strand', 'start_pos', 'end_pos',
                 'feature']  # , 'gene', 'product', 'eggNOG annot']

# annotations_file = '/net/mraid08/export/genie/LabData/Data/Annotations/Segata_annots/Segata_annots_2021_02_04_prokka_eggnog.csv'
annotations_file = '/net/mraid08/export/genie/LabData/Data/Annotations/Segal_annots/Segal_annots_2021_07_31_prokka_eggnog.csv'
gene_annotations = None
annotations_list = None
global_body_site = 'Gut'


def load_gene_annotations():
    global gene_annotations
    if gene_annotations is None:
        gene_annotations = GeneAnnotLoader(annotations_file, body_site=global_body_site).get_annot()


def load_annotations_list():
    load_gene_annotations()

    global annotations_list
    if annotations_list is None:
        annotations_list = Annotationist(gene_annotations, body_site=global_body_site, load_sequences=False)


def order_columns(snps):
    # order the columns by the predefined first columns

    first_cols = [col for col in first_columns if col in snps.columns]
    other_cols = [col for col in snps.columns if col not in first_columns]

    return snps[first_cols + other_cols]


def find_unique_snps(x_mwas_files_path, y_mwas_files_path, output_dir, pval_col='Global_Bonferroni', alpha=0.05):
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

    pvals_df.to_csv(os.path.join(output_dir, 'pvals_count.csv'))
    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_dir, 'snps_unique.h5'), key='snps')

    return snps


def choose_contig_type(snps, contig_type):
    def remove_contig_parts(contigs: Iterable[str]) -> List[str]:
        """Returns the given contig names without the trailing '_P###'."""

        # contig_part_re = re.compile('_P\\d+$')
        contig_part_re = re.compile('^C_\\d+')
        # Maps old contig name to new contig name, to avoid repetition.
        cache = {}

        def cached_name(contig):
            if contig not in cache:
                # cache[contig] = contig_part_re.sub('', contig)
                cache[contig] = contig_part_re.findall(contig)[0]
            return cache[contig]

        return [cached_name(x) for x in contigs]

    # first run no matter what
    if 'Contig_with_part' not in snps.columns:
        snps['Contig_with_part'] = snps.index.get_level_values('Contig')
        snps['Contig_with_part'] = snps['Contig_with_part'].astype(str)

    # first time you ask for contig without part
    if 'Contig_without_part' not in snps.columns and contig_type == 'Contig_without_part':
        snps['Contig_without_part'] = remove_contig_parts(snps['Contig_with_part'].values)
        snps['Contig_without_part'] = snps['Contig_without_part'].astype(str)

    # the requested assignment
    snps = snps.reset_index().drop(['Contig'], axis=1)
    snps['Contig'] = snps[contig_type]
    snps = snps.set_index(indices)

    snps = order_columns(snps)

    return snps


def add_surrounding_genes(snps, annotations=None):
    # add multiple genes
    if annotations is None:
        annotations = annotations_list

    snps = annotations.lookup_df(snps)  # current/upstream/downstream, +/- strand

    snps = order_columns(snps)

    return snps


def flatten_surrounding_genes(snps):
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
        no_gene = df[[f'Plus{relation}GeneID', f'Minus{relation}GeneID']].isna().sum(axis=1) == 2
        no_gene = df.loc[no_gene].drop([f'Minus{relation}GeneID', f'Minus{relation}GeneDistance',
                                        'PlusCurrentMajorCodon', 'PlusCurrentMinorCodon',
                                        'MinusCurrentMajorCodon', 'MinusCurrentMinorCodon'], axis=1)
        no_gene = no_gene.rename(columns={f'Plus{relation}GeneID': 'GeneID',
                                          f'Plus{relation}GeneDistance': 'GeneDistance'})
        no_gene['MajorCodon'] = None
        no_gene['MinorCodon'] = None

        # plus genes
        plus_genes = df[[f'Plus{relation}GeneID', f'Plus{relation}GeneDistance']].dropna().index
        plus_genes = df.loc[plus_genes].drop([f'Minus{relation}GeneID', f'Minus{relation}GeneDistance',
                                              'MinusCurrentMajorCodon', 'MinusCurrentMinorCodon'], axis=1)
        plus_genes = plus_genes.rename(columns={f'Plus{relation}GeneID': 'GeneID',
                                                f'Plus{relation}GeneDistance': 'GeneDistance',
                                                'PlusCurrentMajorCodon': 'MajorCodon',
                                                'PlusCurrentMinorCodon': 'MinorCodon'})

        # minus genes
        minus_genes = df[[f'Minus{relation}GeneID', f'Minus{relation}GeneDistance']].dropna().index
        minus_genes = df.loc[minus_genes].drop([f'Plus{relation}GeneID', f'Plus{relation}GeneDistance',
                                                'PlusCurrentMajorCodon', 'PlusCurrentMinorCodon'], axis=1)
        minus_genes = minus_genes.rename(columns={f'Minus{relation}GeneID': 'GeneID',
                                                  f'Minus{relation}GeneDistance': 'GeneDistance',
                                                  'MinusCurrentMajorCodon': 'MajorCodon',
                                                  'MinusCurrentMinorCodon': 'MinorCodon'})

        df = pd.concat([plus_genes, minus_genes, no_gene])
        if relation != 'Current':
            no_gene['MajorCodon'] = None
            no_gene['MinorCodon'] = None
        df['GeneRelation'] = relation
        dfs.append(df)

    snps = pd.concat(dfs)
    snps.loc[snps['GeneID'].isna(), 'GeneRelation'] = None
    snps = snps.drop_duplicates().sort_index()

    # for snps that have both GeneID na and GeneID not na, delete their GeneID na version
    duplicated_indices = snps['GeneID'].isna().index.intersection((~snps['GeneID'].isna()).index)
    snps = snps[~(snps['GeneID'].isna() & snps.index.isin(duplicated_indices))]

    snps = order_columns(snps)

    return snps


def add_gene_annotations(snps, on='GeneID', rsuffix='', annotations=None):
    # add gene annotations
    if annotations is None:
        annotations = gene_annotations

    snps = snps.join(annotations, on=on, rsuffix=rsuffix)
    snps = snps.drop(['SGB', 'contig'], axis=1)

    snps = order_columns(snps)

    return snps


def add_codons(snps, annotations=None):
    # adds codons
    if annotations is None:
        annotations = annotations_list

    snps = annotations.assign_codons(snps)

    snps = snps.rename(columns={'Major': 'MajorAllele', 'Minor': 'MinorAllele'})

    snps = order_columns(snps)

    return snps


def add_amino_acids(snps):
    # adds amino acids

    is_codon = [len(codon) == 3 and codon != 'nan' for codon in snps['MajorCodon'].astype(str)]

    snps.loc[is_codon, 'MajorAA'] = snps.loc[is_codon, 'MajorCodon'].apply(translate_codon)
    snps.loc[is_codon, 'MinorAA'] = snps.loc[is_codon, 'MinorCodon'].apply(translate_codon)

    snps = order_columns(snps)

    return snps

def add_synonymous(snps):
    snps = snps.set_index(['Y', 'GeneRelation', 'strand'], append=True).\
        reorder_levels(['Y', 'Species', 'Contig', 'Position', 'GeneRelation', 'strand'])

    if 'MajorAA' in snps.columns: # so we don't translate twice
        aas = snps[['MajorAA', 'MinorAA']].dropna()
        return snps.join((aas['MajorAA'] == aas['MinorAA']).rename('IsSynonymous'))

    codons = snps[['MajorCodon', 'MinorCodon']].dropna()
    return snps.join((codons['MajorCodon'].apply(lambda x: Seq(x).translate()) == \
                      codons['MinorCodon'].apply(lambda x: Seq(x).translate())).rename('IsSynonymous'))


def add_stop(snps):
    aas = snps[['MajorAA', 'MinorAA', 'IsSynonymous']].dropna()
    return snps.join(( (aas['IsSynonymous'] == False) & ((aas['MajorAA'] == '*') | (aas['MinorAA'] == '*')) ).rename('StopMutation'))


def _get_row_desc(x):
    if config.mb_species_db_basename == 'segata':
        taxa = x.taxa.split('_')
        taxa1 = taxa[2][0] if len(taxa) >= 3 else ''
        taxa2 = taxa[3] if len(taxa) >= 4 else ''
        short_taxa = f'{taxa1}. {taxa2}'
    else:
        taxa = x.taxa.split(' ')
        short_taxa = f'{taxa[0][0]}. {taxa[1]}' if len(taxa) > 1 else taxa[0]

    in_gene = x.name[4] == 'Current'
    gene = ', int' if not in_gene else (', unk' if str(x['gene']) == 'nan' else f', {x["gene"]}')
    syn = '' if not in_gene or np.isnan(x['IsSynonymous']) else (', syn' if x['IsSynonymous'] else ', ns')

    return f'{short_taxa}{gene}{syn}'


def _combine_row_desc(x):
    ret = ''
    if ('Current', '+') in x.index and str(x.loc[('Current', '+')]) != 'nan':
        ret += x.loc[('Current', '+')]
    if ('Current', '-') in x.index and str(x.loc[('Current', '-')]) != 'nan':
        if len(ret) > 0:
            ret += ', '.join(x.loc[('Current', '-')].split(', ')[1:])
        else:
            ret += x.loc[('Current', '-')]
    if len(ret) > 0:
        return ret
    x = x.dropna()
    return x.iloc[0] if len(x) > 0 else ''


def add_snp_descriptions(snps):
    long_desc = snps.apply(_get_row_desc, axis='columns').unstack(level=['GeneRelation', 'strand'])
    long_desc = long_desc.apply(_combine_row_desc, axis='columns').rename('snp_desc')
    return snps.join(long_desc)


def add_taxonomy(snps, mb_species_db_basename='segata', body_site='Gut'):
    # adds taxonomy

    tax_df = taxonomy_df(level_as_numbers=False).set_index('SGB')[['Species']].rename(columns={'Species': 'taxa'})

    # ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    snps = snps.join(tax_df, on='Species')

    snps = order_columns(snps)

    return snps


def add_sequences(snps):
    # add genomic and amino acid sequences to genes

    unique_genes = snps.dropna(subset=['GeneID']).drop_duplicates(['GeneID']).reset_index().set_index('GeneID')

    contig_col = 'Contig_without_part' if 'Contig_without_part' in unique_genes.columns else 'Contig'
    seqs = read_sequences([GenePosition(*[x[0], x[1], int(x[2]), int(x[3] + 1), x[4]]) for x in
                           unique_genes[['Species', contig_col, 'start_pos', 'end_pos', 'strand']].values])

    snps[['nuc_seq', 'aa_seq']] = None, None
    is_protein = unique_genes['feature'] == 'CDS'
    for gene_id, nuc_seq, translate in zip(unique_genes.index, seqs, is_protein):
        aa_seq = nuc_seq.translate() if translate else ''
        snps.loc[snps['GeneID'] == gene_id, ['nuc_seq', 'aa_seq']] = str(nuc_seq), str(aa_seq)

    snps = order_columns(snps)

    return snps


def run(mwas_file_path, output_dir, body_site='Gut'):
    global global_body_site
    global_body_site = body_site

    snps = pd.read_hdf(mwas_file_path)
    load_annotations_list()
    snps = add_surrounding_genes(choose_contig_type(snps, 'Contig_without_part'))
    snps = add_codons(choose_contig_type(snps, 'Contig_with_part'))
    snps = flatten_surrounding_genes(snps)
    snps = add_amino_acids(snps)
    snps = add_taxonomy(snps)
    snps = add_gene_annotations(snps)
    snps.to_hdf(os.path.join(output_dir, 'snps_gene_annotations.h5'), key='snps')
    snps.reset_index().to_excel(os.path.join(output_dir, 'snps_gene_annotations.xlsx'))
    snps = add_sequences(snps)
    snps.to_hdf(os.path.join(output_dir, 'snps_sequences.h5'), key='snps')

    return snps

# if __name__ == '__main__':
#     f = '/Users/eran/Dropbox (Weizmann Institute)/Genie/LabData/Analyses/eran/20210419_175302_mwas_10k_bm/annotations/snps_gene_annotations.h5'
#     df = pd.read_hdf(f)
#     print(add_snp_descriptions(df))