import os
from LabMBPipeline.mmmbp import get_args, main

input_dir = r'/net/mraid08/export/genie/Data/ExternalSamples/IBD/metagenomics'
output_dir = r'/net/mraid08/export/jafar/Microbiome/Analyses/saar/IBD/mmmbp'

module_seq = 'RAW,QCS,HGF,UNT,URA'  # maybe should have had MSE

df_path = os.path.join(output_dir, 'DFOut')
tmp_dir = os.path.join(output_dir, 'tmp')
email = 'saar.shoer@weizmann.ac.il'
max_r = '550'

args, sm = get_args([

                     df_path,
                     tmp_dir,
                     email,

                     '--module_seq', module_seq,
                     '--use_general_python',
                     '--max_r', max_r,

                     '--raw_input_folders', input_dir

                     ])

main(args, sm)
