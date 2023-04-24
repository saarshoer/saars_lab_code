"""Runs the lab's URA pipeline."""

import os

from LabMBPipeline.mmmbp import get_args, main

# change here for your dir
out_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/saar/PNP3/simulation'
DFOut_param = os.path.join(out_dir, 'DFOut') + ' '
tmp_param = os.path.join(out_dir, 'tmp2') + ' '

Email = ' saar.shoer@weizmann.ac.il '
General_params = ' --max_r 550 --use_general_python '
Modules = ' --module_seq "RAW,UNT,URS" '
RAW_params = f' --raw_read_length 100 --raw_input_folders  {out_dir}/input '
UNT_params = ' --unt_append_to_existing '

# change this parameters
URS_params = (
    ' --urs_num_mapped_to_subsample 5000000 '
    '--urs_min_mapped_to_retain 1000000 --urs_run_type LargeOrNewGenusSGBs ')

args, sm = get_args((DFOut_param + tmp_param + Email + General_params +
                     Modules + RAW_params + UNT_params + URS_params).split())

main(args, sm)
