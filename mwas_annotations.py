import re
import os
import glob
import numpy as np
import pandas as pd
from pandas import HDFStore
from typing import List, Iterable
from LabQueue.qp import qp, fakeqp
from LabUtils.addloglevels import sethandlers
from LabUtils.SeqUtils import translate_codon, codons_inv
from LabData.DataAnalyses.MBSNPs.taxonomy import taxonomy_df
from UseCases.DataAnalyses.mbsnps_annotation import Annotationist
from LabData.DataLoaders.MBSNPLoader import get_mbsnp_loader_class
from LabData.DataLoaders.GeneAnnotLoader import GeneAnnotLoader, ANNOTS_03_2020_EXPANDED


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
    snps.to_hdf(os.path.join(output_path, 'snps.h5'), key='snps')

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
    idx = snps.index.names
    snps = snps.reset_index()
    snps = snps.rename(columns={'Contig': 'Contig_with_part'})
    snps['Contig'] = remove_contig_parts(snps['Contig_with_part'].values)
    snps = snps.set_index(idx)

    snps.to_hdf(os.path.join(output_path, 'snps_contig.h5'), key='snps')

    return snps


def add_single_gene_maf_annotations(snps, P, output_path):
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

    snps.to_hdf(os.path.join(output_path, 'snps_maf_annotations.h5'), key='snps')

    return snps


def map_codon_information(snps, output_path):
    # add codon information

    # remove the codons extracted from the reference genomes, rather than the MAFs
    snps = snps.drop(columns=['MajorCodon', 'MinorCodon'])

    is_coding = snps['MajorCodonMAF'] >= 0

    # codon
    snps.loc[is_coding, 'MajorCodonMAF'] = snps.loc[is_coding, 'MajorCodonMAF'].map(codons_inv)
    snps.loc[is_coding, 'MinorCodonMAF'] = snps.loc[is_coding, 'MinorCodonMAF'].map(codons_inv)

    # amino acid
    snps.loc[is_coding, 'MajorAA'] = snps.loc[is_coding, 'MajorCodonMAF'].apply(translate_codon)
    snps.loc[is_coding, 'MinorAA'] = snps.loc[is_coding, 'MinorCodonMAF'].apply(translate_codon)

    snps.to_hdf(os.path.join(output_path, 'snps_codon.h5'), key='snps')

    return snps


def add_gene_info(snps, output_path, on='GeneID', rsuffix=''):
    # add gene annotations

    genes = GeneAnnotLoader(ANNOTS_03_2020_EXPANDED).get_annot()
    snps = snps.join(genes.drop(columns=['strand', 'feature']), on=on, rsuffix=rsuffix)

    snps.to_csv(os.path.join(output_path, 'snps_gene_annotations.csv'))

    return snps


def add_surrounding_genes(snps, output_path):
    # add multiple genes

    genes = GeneAnnotLoader(ANNOTS_03_2020_EXPANDED).get_annot()
    snps = Annotationist(genes, load_sequences=False).lookup_df(snps)  # upstream/downstream, +/- strand genes

    snps.to_hdf(os.path.join(output_path, 'snps_multiple_genes.h5'), key='snps')

    return snps


def add_taxonomy(snps, output_path):

    tax_df = taxonomy_df(level_as_numbers=False).set_index('SGB')[
        ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]
    snps = snps.join(tax_df, on='Species')

    snps.to_hdf(os.path.join(output_path, 'snps_taxonomy.h5'), key='snps')

    return snps


if __name__ == '__main__':

    def antibiotics():

        # parameters
        x_mwas_files_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_raw/mb_gwas_SGB_*.h5'
        y_mwas_files_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed/*/SGB_*.h5'
        output_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed'
        from antibiotics_mwas import P

        # run
        snps = find_unique_snps(x_mwas_files_path, y_mwas_files_path, output_path)
        snps = add_contig_without_part(snps, output_path)

        add_surrounding_genes(snps, output_path)

        snps = add_single_gene_maf_annotations(snps, P, output_path)
        snps = map_codon_information(snps, output_path)
        snps = add_gene_info(snps, output_path)

        return snps


    jobs_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/jobs'

    os.chdir(jobs_path)
    sethandlers()

    with qp(jobname='annot', _delete_csh_withnoerr=True, q=['himem7.q'], max_r=1, _mem_def='50G') as q:
        q.startpermanentrun()
        q.waitforresult(q.method(antibiotics))
