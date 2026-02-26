#!/usr/bin/env python3
'''
python script that runs the asym_stats.py script in batches
'''

from _utils.logger import logger
log = logger()

def qc(input_args):
  subj, nroi, args = input_args
  import pandas as pd
  import os
  in_fname = args._in.replace('%sub',subj)
  if not os.path.isfile(in_fname): return False, False

  skip = True
  try:
    tmp = pd.read_table(f'{args.out}/global_graph/{subj}.txt')
    if tmp.shape[0] != 17: skip = False
    tmp = pd.read_table(f'{args.out}/global_asym/{subj}.txt')
    if tmp.shape[0] != 17: skip = False
    tmp = pd.read_table(f'{args.out}/local/{subj}.txt')
    if tmp.shape[0] != nroi or tmp.shape[1] != 7: skip = False
    tmp = pd.read_table(f'{args.out}/local_asym/{subj}.txt')
    if tmp.shape[0] != nroi/2 or tmp.shape[1] != 21: skip = False
  except: skip = False
  return True, skip

def main(args):
    import os
    import numpy as np
    from tqdm import tqdm
    from multiprocessing import Pool
    
    # nroi: HCP = 376, 500sym = 334, aparc = 84, economo = 102, sjh = 1027
    if 'HCP' in args._in: nroi = 376
    elif os.path.basename(os.path.dirname(args._in)) == 'aparc_seq': nroi = 84
    elif 'economo' in args._in: nroi = 102
    else: raise ValueError('Unknown parcellation')
      
    if args.force:
      fs = ' -f'
    else: fs = ''
    
    # array submitter
    from _utils.slurm import array_submitter
    submitter = array_submitter(name = 'pheno',partition = 'icelake',timeout = 15, parallel = 8)
  
    subjs = np.loadtxt(args.subjs,dtype = 'U')
    pool = Pool(64)

    # check progress for all subjects
    if not args.force:
      progress = list(tqdm(pool.imap(qc, [(subj, nroi, args) for subj in subjs], chunksize = 64), total = len(subjs), desc = 'Checking progress for subjects'))
      completed = [x[1] for x in progress]; missing = [x[0] for x in progress]; to_submit = [x[0] and not x[1] for x in progress]
      n_completed = sum(completed); n_missing = sum(missing); n_to_submit = sum(to_submit)
      log.log(f'{n_completed} subjects completed, {n_missing} subjects missing, {n_to_submit} subjects with incomplete output')
      for subj in subjs[to_submit]: submitter.add(f'python pheno.py {subj} -i {args._in.replace("%sub", subj)} -o {args.out} {fs}')
    else:
      n_missing = 0
      for subj in subjs:
        in_fname = args._in.replace('%sub',subj)
        if not os.path.isfile(in_fname):
          log.log(f'{subj}: no connectome found')
          n_missing += 1
          continue
        submitter.add(f'python pheno.py {subj} -i {in_fname} -o {args.out} {fs}')
      log.log(f'{n_missing} subjects missing, {len(subjs)-n_missing} subjects to submit')
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