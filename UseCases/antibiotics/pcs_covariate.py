import os
import glob
import pickle
import pandas as pd
from LabQueue.qp import qp, fakeqp
from sklearn.impute import KNNImputer
from sklearn.decomposition import PCA
from LabUtils.addloglevels import sethandlers
from LabData.DataLoaders.MBSNPLoader import MAF_MISSING_VALUE

study = '10K'
run_type = 'within'  # always

run_dir = f'/net/mraid08/export/jafar/Microbiome/Analyses/saar/antibiotics/{study}/{run_type}'
mafs_dir = f'/net/mraid08/export/mb/MBPipeline/Analyses/MBSNP/Gut/MAFS/{study}'


def do(f):
    species = f'Rep_{f.split("_")[-1].split(".")[0]}'
    samples = pd.read_pickle(os.path.join(os.path.dirname(run_dir), 'data_frames', 'snps.df')).index.tolist()
    positions = pd.read_hdf(f).reset_index().set_index('Contig')['Position']
    contig_parts = positions.index.unique().tolist()

    maf_concat = None

    with pd.HDFStore(os.path.join(mafs_dir, f'mb_snp_maf_{species}_R1_S500.h5'), 'r') as maf:

        for contig_part in contig_parts:
            contig_positions = positions.loc[contig_part].tolist()
            if type(contig_positions) is int:
                contig_positions = [contig_positions]

            contig_maf = maf[f'/{contig_part}']
            if 'Position' in contig_maf.index.names:
                contig_maf = contig_maf.unstack('Position')
            contig_maf = contig_maf.loc[contig_maf.index.intersection(samples), contig_positions]
            contig_maf.columns = contig_maf.columns.astype(str) + f'_{contig_part}'

            maf_concat = maf_concat.join(contig_maf, how='outer') if maf_concat is not None else contig_maf

    maf_concat = maf_concat.fillna(MAF_MISSING_VALUE)
    maf_concat.to_hdf(os.path.join(os.path.dirname(run_dir), 'pcs_covariate', f'mb_snp_maf_{species}.h5'), key='snps', complevel=9)

    maf_concat = pd.read_hdf(os.path.join(os.path.dirname(run_dir), 'pcs_covariate', f'mb_snp_maf_{species}.h5'))

    if maf_concat.shape[1] > 20:
        samples = maf_concat.index.tolist()
        imputer = KNNImputer(missing_values=MAF_MISSING_VALUE, n_neighbors=20, weights='distance')  # n_neighbors with actual value (non-missing)
        maf_concat = imputer.fit_transform(maf_concat)
        pca = PCA(n_components=20, random_state=42)
        maf_concat = pca.fit_transform(maf_concat)
        maf_concat = pd.DataFrame(maf_concat, index=samples)

        maf_concat.to_hdf(os.path.join(os.path.dirname(run_dir), 'pcs_covariate', f'mb_snp_pcs_{species}.h5'), key='snps', complevel=9)

        return pca.explained_variance_ratio_


if __name__ == '__main__':

    jobs_dir = os.path.join(run_dir, 'jobs')
    os.chdir(jobs_dir)
    sethandlers(file_dir=jobs_dir)

    os.makedirs(os.path.join(os.path.dirname(run_dir), 'pcs_covariate'))

    with qp(jobname='pcscov', _tryrerun=True, _delete_csh_withnoerr=True, _mem_def='30G') as q:
        q.startpermanentrun()
        tkttores = {}

        print('start sending jobs')
        tested = glob.glob(os.path.join(run_dir, 'naive_hdfs', 'mb_gwas_Rep_*_Rep_*.h5'))
        for file in tested:
            s = f'Rep_{file.split("_")[-1].split(".")[0]}'
            tkttores[s] = q.method(do, [file], _job_name=f'p{file.split("_")[-1].split(".")[0]}')
        print('finished sending jobs')

        ev = dict()
        with open(os.path.join(os.path.dirname(run_dir), 'pcs_covariate', 'ev'), 'wb') as f:
            print('start waiting for jobs')
            for k, v in tkttores.items():
                ev[k] = q.waitforresult(v)
                if (len(ev) % 10 == 0) | (len(ev) == len(tkttores)):
                    pickle.dump(ev, f)
            print('finished waiting for jobs')

    print('done')
