from argparse import ArgumentParser
from collections import namedtuple
import os
import random

import setup_robust_wgcna as setup_wgcna

CHARS='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ012345678'
RWGCNA_Args = namedtuple('RWGCNA_Args', ('data', 'out_dir', 'power', 'cut_height', 'min_module_size', 'module_corr', 'loo', 'nboot', 'transposed', 'type', 'fmt_dir'))

QSUB_BASH_WRAPPER="#!/bin/bash\n#$ -cwd -l h_rt={}:00:00,h_data={}G,highp {}\n. /u/local/Modules/default/init/modules.sh\nmodule load R/3.4.0\nmodule load python/2.7.3\nmodule load java/1.8.0_77\n\n{}"

def make_qsub(cmd, rt=12, mem=8, cores=1):
  core_str = '-pe shared {}'.format(cores) if cores > 1 else ''
  return QSUB_BASH_WRAPPER.format(rt, int(float(mem)/cores), core_str, cmd)


def get_args():
  parser = ArgumentParser('run wgcna pipeline')
  parser.add_argument('expression', help='The expression data. Rows are genes, columns are samples.')
  parser.add_argument('outdir', help='The output directory')
  parser.add_argument('--qsub', action='store_true', help='qsub rather than print commands')
  parser.add_argument('--jpfx', help='The job prefix', default=None)
  parser.add_argument('--fmt_dir', help='The directory holding the format files', default='.')
  parser.add_argument('--nboot', help='Number of bootstraps', default=100, type=int)
  parser.add_argument('--memory', help='Memory to use (in G)', default='8', type=str)
  args = parser.parse_args()
  return args

def main():
  args = get_args()
  outdir = os.path.dirname(os.path.abspath(args.outdir)) + '/' + args.outdir
  job_no = 1
  job_pfx = args.jpfx or 'WP{}'.format(random.sample(CHARS, 5))
  commands = []
  # zeroth thing: make the script dir
  os.system('mkdir -p ' + outdir)
  scriptdir = '{}/scripts/'.format(outdir)
  os.system('mkdir -p ' + scriptdir)
  # first estimate the power
  powfile = '{}/{}'.format(outdir, 'power.param')
  powcmd = 'ALLOW_WGCNA_THREADS=4 Rscript {}/estimate_power.R -e {} -o {}'.format(args.fmt_dir, args.expression, outdir)
  scriptfile = '{}/{}_{}.bash'.format(scriptdir, job_pfx, job_no)
  with open(scriptfile, 'w') as out:
    out.write(make_qsub(powcmd, cores=4, mem=args.memory))
  commands.append('qsub -N {}_{} {}'.format(job_pfx, job_no, scriptfile))
  powjob = '{}_{}'.format(job_pfx, job_no)
  job_no += 1
  wgcna_args = RWGCNA_Args(args.expression, outdir, powfile, 0.975, 50, 0.85, False, args.nboot, True, 'signed', args.fmt_dir)
  tom_bash, merge_bash, nsamples = setup_wgcna.make_commands(wgcna_args)
  commands.append('qsub -v INTEGRATED_PIPELINE_DIR={} -v ALLOW_WGCNA_THREADS=5 -t 1:{} -N {}_{} -hold_jid {} {}'.format(args.fmt_dir, nsamples, job_pfx, job_no, powjob, tom_bash))
  tomjob = '{}_{}'.format(job_pfx, job_no)
  job_no += 1
  mergejob = '{}_{}'.format(job_pfx, job_no)
  commands.append('qsub -v INTEGRATED_PIPELINE_DIR={} -v ALLOW_WGCNA_THREADS=5 -N {} -hold_jid {} {}'.format(args.fmt_dir, mergejob, tomjob, merge_bash))

  for cmd in commands:
    print(cmd)
    if args.qsub:
      os.system(cmd)
  

if __name__ == '__main__':
  main()

