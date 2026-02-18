#!/usr/bin/env python3
'''
python script that runs the asym_stats.py script in batches
'''

from _utils.logger import logger
log = logger()

def main(args):
    import os
    import numpy as np
    import pandas as pd
    from tqdm import tqdm
    import time
    
    # nroi: HCP = 376, 500sym = 334, aparc = 84, economo = 102, sjh = 1027
    nroi = 376
      
    if args.force:
      fs = ' -f'
    else: fs = ''
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(name = 'pheno',partition = 'icelake',timeout = 15, parallel = 8)
    
    subjs = np.loadtxt(args.subjs,dtype = 'U')
    n_completed = 0
    n_missing = 0
    n_submitted = 0

    for subj in tqdm(subjs, desc = 'Checking progress for subjects'):
      in_fname = args._in.replace('%sub',subj)
      
      # check progress
      if not os.path.isfile(in_fname):
        log.log(f'{subj}: no connectome found')
        n_missing += 1
        continue
      
      skip = True
      if args.force: skip = False
      
      if skip:
        try:
          tmp = np.loadtxt(f'{args.out}/global/{subj}.txt')
          if tmp.size != 17: skip = False
          tmp = np.loadtxt(f'{args.out}/global_asym/{subj}.txt')
          if tmp.size != 17: skip = False
          tmp = pd.read_csv(f'{args.out}/local/{subj}.txt')
          if tmp.shape[0] != nroi or tmp.shape[1] != 7: skip = False
          tmp = pd.read_csv(f'{args.out}/local_asym/{subj}.txt')
          if tmp.shape[0] != nroi/2 or tmp.shape[1] != 21: skip = False
        except: skip = False
      
      if skip: 
        n_completed += 1
      else:
        indir = args._in.replace('%sub', subj)
        n_submitted += 1
        submitter.add(f'python pheno.py {subj} -i {indir} -o {args.out} {fs}')
    
    submitter.submit()
    
    
if __name__ == '__main__':
    # input argument processing
    from _utils.slurm import slurm_parser
    parser = slurm_parser(description='This programme processes the connectome '+
                               ' for one single individual for imaging derived phenotypes')
    parser.add_argument('-i','--in',dest = '_in', help =
        'Target file to screen',
        default = '/home/yh464/rds/rds-rb643-ukbiobank2/Data_Imaging/'+
        'UKB%sub/func/fMRI/parcellations/HCP.fsaverage.aparc_seq/Connectivity_sc2345.txt')
    parser.add_argument('-s', '--subjs', dest = 'subjs', help = 'list of subjs',
        default = '../params/subjlist_rsfmri_hcp.txt')
    parser.add_argument('-o','--out',dest = 'out', help = 'Output directory',
        default = '../pheno/ukb/')
    parser.add_argument('-f','--force', dest = 'force', help = 'Force output',
        default = False,const = True, action = 'store_const')
    args = parser.parse_args()
    import os
    for arg in ['_in','out','subjs']:
        setattr(args, arg, os.path.realpath(getattr(args, arg)))
        
    from _utils import cmdhistory, logger
    logger.splash(args)
    cmdhistory.log()
    try: main(args)
    except: cmdhistory.errlog()