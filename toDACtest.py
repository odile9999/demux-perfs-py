
import os
import re
import argparse

import xifufwk.process.energyResolutionEstimation as ERes

''' 
For DAC test - march 2020
Usage:
python --src path/to/session

Compute energy resolution estimation in batch mode

'''

#PATTERN YYYYMMDD_HHMMSS
DATE_SELECTER = "^\d{8}\_\d{6}"
OUTPUTDIR = 'OUTPUTS'
XIFUSIM = 'XIFU_INST_XIFUSIM_20191205-164022_PP_PARAM_CONFXIFUSIM.fits'

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--src", 
                    help="Absolute path to session try", 
                    type=str)
    
    return parser.parse_args()

def extract_date_session(src):
    #Extract test date from session folder
    (src_head, src_tail) = os.path.split(src)
    date_test = re.findall(DATE_SELECTER, src_tail)

    if not date_test or date_test is None:
        raise Exception(f'Invalid directory {date_test} - expected format: YYYYMMDD_HHMMSS__freetext')
    else:
        date_test = date_test[0]
        date_test = date_test.replace("_", "-")

    return date_test

def check_output_dir(src):
    #check if output directory exist
    path_ = os.path.join(src, OUTPUTDIR)
    if not os.path.isdir(path_):
        os.makedirs(path_)


if __name__ == "__main__":
    args = parser()
    #date_session = extract_date_session(args.src)
    #check_output_dir(args.src)
    
    #x_path = os.path.join(args.src, 'SCRIPTS_PYTHON', 'demux-perfs-py', XIFUSIM)
    src=os.path. abspath('../..') 
    x_path = os.path.join(src, 'SCRIPTS_PYTHON', 'demux-perfs-py', XIFUSIM)
    print("> src : ", src)
    print("> x_path :", x_path)
    
    #ERes.compute_energy(path_session=args.src, x_path=x_path)
    ERes.compute_energy(path_session=src, x_path=x_path)
    