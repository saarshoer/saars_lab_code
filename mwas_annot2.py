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

first_columns = ['SGB', 'contig', 'Contig_with_part', 'Contig_without_part',
                 'MajorAllele', 'MinorAllele', 'MajorCodon', 'MinorCodon', 'MajorAA', 'MinorAA',
                 'GeneID', 'GeneRelation', 'GeneDistance', 'strand', 'start_pos', 'end_pos',
                 'feature', 'gene', 'product', 'eggNOG annot']

annotations_file = ANNOTS_03_2020_EXPANDED
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
        contig_part_re = re.compile('_P\\d+$')
        # Maps old contig name to new contig name, to avoid repetition.
        cache = {}

        def cached_name(contig):
            if contig not in cache:
                cache[contig] = contig_part_re.sub('', contig)
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


def add_surrounding_genes(snps, output_dir):
    # add multiple genes

    load_annotations_list()
    snps = annotations_list.lookup_df(snps)  # current/upstream/downstream, +/- strand

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_dir, 'snps_surrounding_genes.h5'), key='snps')

    return snps


def flatten_surrounding_genes(snps, output_dir):
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
        if relation is not 'Current':
            no_gene['MajorCodon'] = None
            no_gene['MinorCodon'] = None
        df['GeneRelation'] = relation
        dfs.append(df)

    snps = pd.concat(dfs).reset_index()
    snps = snps.drop_duplicates(snps.columns[:-1]).set_index(indices).sort_index()
    snps.loc[snps['GeneID'].isna(), 'GeneRelation'] = None
    snps['GeneID'] = snps['GeneID'].replace(np.nan, None)#.astype(int)
    snps['MajorCodon'] = snps['MajorCodon'].astype(str)
    snps['MinorCodon'] = snps['MinorCodon'].astype(str)

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_dir, 'snps_surrounding_genes_flat.h5'), key='snps')

    return snps


def add_gene_annotations(snps, output_dir, on='GeneID', rsuffix=''):
    # add gene annotations

    original_columns = snps.columns

    load_gene_annotations()
    snps = snps.join(gene_annotations, on=on, rsuffix=rsuffix)  # TODO: take only useful columns, rename unclear ones
    snps = snps.drop(['SGB', 'contig'], axis=1)
    new_columns = list(set(snps.columns) - set(original_columns) - set(['start_pos', 'end_pos']))
    snps[new_columns] = snps[new_columns].astype(str)

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_dir, 'snps_gene_annotations.h5'), key='snps')

    return snps


def add_codons(snps, output_dir):
    # adds codons

    load_annotations_list()
    snps = annotations_list.assign_codons(snps)

    snps = snps.rename(columns={'Major': 'MajorAllele', 'Minor': 'MinorAllele'})
    codon_cols = ['PlusCurrentMajorCodon', 'PlusCurrentMinorCodon', 'MinusCurrentMajorCodon', 'MinusCurrentMinorCodon']
    snps[codon_cols] = snps[codon_cols].astype(str)

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_dir, 'snps_codons.h5'), key='snps')

    return snps


def add_amino_acids(snps, output_dir):
    # adds amino acids

    is_codon = [len(codon) == 3 for codon in snps['MajorCodon'].astype(str)]

    snps.loc[is_codon, 'MajorAA'] = snps.loc[is_codon, 'MajorCodon'].apply(translate_codon)
    snps.loc[is_codon, 'MinorAA'] = snps.loc[is_codon, 'MinorCodon'].apply(translate_codon)

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_dir, 'snps_amino_acids.h5'), key='snps')

    return snps


def add_taxonomy(snps, output_dir):
    # adds taxonomy

    tax_df = taxonomy_df(level_as_numbers=False).set_index('SGB')[
        ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]
    snps = snps.join(tax_df, on='Species')

    snps = order_columns(snps)
    snps.to_hdf(os.path.join(output_dir, 'snps_taxonomy.h5'), key='snps')

    return snps


def run(mwas_file_path, output_dir, body_site='Gut'):

    global global_body_site
    global_body_site = body_site

    snps = pd.read_hdf(mwas_file_path)
    # snps = add_surrounding_genes(choose_contig_type(snps, 'Contig_without_part'), output_dir)
    # snps = add_codons(choose_contig_type(snps, 'Contig_with_part'), output_dir)
    snps = flatten_surrounding_genes(snps, output_dir)
    snps = add_amino_acids(snps, output_dir)
    snps = add_gene_annotations(snps, output_dir)

    return snps

