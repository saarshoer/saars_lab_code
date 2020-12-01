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

# parameters
x_mwas_files_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_raw/mb_gwas_SGB_*.h5'
y_mwas_files_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed/*/SGB_*.h5'
annotations_files_path = '/net/mraid08/export/genie/LabData/Analyses/saarsh/anti_mwas_processed'
jobs_path = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/jobs'
from antibiotics_mwas import P


def do():

    # phase 1
    x_mwas_files = glob.glob(x_mwas_files_path)
    y_mwas_files = glob.glob(y_mwas_files_path)

    x_species = list(set(['SGB' + file.split('SGB')[-1].split('.')[0] for file in x_mwas_files]))
    y_species = list(set([os.path.basename(file).split('.')[0] for file in y_mwas_files]))

    pvals_df = pd.DataFrame(index=y_species, columns=x_species)

    snps = set()

    for x_mwas_file in x_mwas_files:
        species = 'SGB' + x_mwas_file.split('SGB')[-1].split('.')[0]
        with HDFStore(x_mwas_file, 'r') as x_mwas_df:
            x_mwas_df = x_mwas_df[species]

        values, counts = np.unique(x_mwas_df.index.get_level_values('Y'), return_counts=True)
        pvals_df.loc[values.tolist(), species] = counts.tolist()

        snps = snps.union(set(x_mwas_df.index.droplevel('Y').values))

    snps = pd.DataFrame(snps)

    pvals_df.to_csv(os.path.join(annotations_files_path, 'pvals_count.csv'))

    snps.to_hdf(os.path.join(annotations_files_path, 'snps.h5'), key='snps')

    # phase 2
    def _add_snp_annotations(species, df):
        maf_annot_fname = get_mbsnp_loader_class(P.body_site).\
            get_snp_maf_annot_fname(species, P.min_reads_per_snp, P.min_subjects_per_snp_cached)
        if os.path.exists(maf_annot_fname):
            with HDFStore(maf_annot_fname, 'r') as maf_annot:
                maf_annot = maf_annot[species]
            df = df.join(maf_annot)
        return df

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
    snps.columns = ['Species', 'Contig', 'Position']
    snps['Contig'] = remove_contig_parts(snps['Contig'].values)
    snps = snps.set_index(snps.columns.tolist())

    snps = pd.concat([_add_snp_annotations(species, df) for species, df in snps.groupby('Species')])

    snps.to_hdf(os.path.join(annotations_files_path, 'snps_maf_annotations.h5'), key='snps')

    # phase 3
    # remove the codons extracted from the reference genomes, rather than the MAFs
    snps = snps.drop(columns=['MajorCodon', 'MinorCodon'])

    is_coding = snps['MajorCodonMAF'] >= 0

    # codon
    snps.loc[is_coding, 'MajorCodonMAF'] = snps.loc[is_coding, 'MajorCodonMAF'].map(codons_inv)
    snps.loc[is_coding, 'MinorCodonMAF'] = snps.loc[is_coding, 'MinorCodonMAF'].map(codons_inv)

    # amino acid
    snps.loc[is_coding, 'MajorAA'] = snps.loc[is_coding, 'MajorCodonMAF'].apply(translate_codon)
    snps.loc[is_coding, 'MinorAA'] = snps.loc[is_coding, 'MinorCodonMAF'].apply(translate_codon)

    snps.to_hdf(os.path.join(annotations_files_path, 'snps_codon_annotations.h5'), key='snps')

    # phase 4
    genes = GeneAnnotLoader(ANNOTS_03_2020_EXPANDED).get_annot()
    snps = snps.join(genes.drop(columns=['strand', 'feature']), on='GeneID')
    # snps = Annotationist(genes, load_sequences=False).lookup_df(snps)  # upstream/downstream, +/- strand genes

    snps.to_hdf(os.path.join(annotations_files_path, 'snps_gene_annotations.h5'), key='snps')

    # # phase 5
    # tax_df = taxonomy_df(level_as_numbers=False).set_index('SGB')[
    #     ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]
    # snps = snps.join(tax_df, on='Species')
    #
    # snps.to_hdf(os.path.join(annotations_files_path, 'snps_taxonomy_annotations.h5'), key='snps')


os.chdir(jobs_path)
sethandlers()

with qp(jobname='annot', _delete_csh_withnoerr=True, q=['himem7.q'], max_r=1, _mem_def='50G') as q:
    q.startpermanentrun()
    q.waitforresult(q.method(do))
