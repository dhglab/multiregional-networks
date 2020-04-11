from argparse import ArgumentParser
import os


def get_args():
  parser = ArgumentParser('setup_robust_wgcna')
  parser.add_argument('data', help='the data frame matrix; rows are samples, columns are genes/biomarkers')
  parser.add_argument('out_dir', help='the output directory')
  parser.add_argument('power', help='the WGCNA power', type=int)
  parser.add_argument('cut_height', help='the cut height', type=float)
  parser.add_argument('min_module_size', help='minimum module size', type=int)
  parser.add_argument('module_corr', help='module correlation (for merging)', type=float)
  parser.add_argument('--loo', action='store_true', help='Use leave-one-out as opposed to bootstraps')
  parser.add_argument('--nboot', type=int, default=100, help='How many bootstraps to use')
  parser.add_argument('--transposed', action='store_true', help='The input data is transposed (so transpose it first)')
  parser.add_argument('--type', help='The adjacency type to use', default='signed')
  parser.add_argument('--fmt_dir', help='The location of the directory holding the format files', default='.')
  parser.add_argument('--jsuffix', help='Suffix to add to job names', default='def')

  return parser.parse_args()

def make_commands(args):
  out_dir = os.path.abspath(args.out_dir)  # lazy

  if not os.path.exists(out_dir):
    os.system('mkdir -p {}'.format(out_dir))

  if args.transposed:
    print 'transposing the data'
    data = list()
    nrow = 0
    for line in open(args.data):
      fields = line.strip().split()
      data.append(fields)
      nrow += 1
    ncol = len(data[0]) # note: header is 1 fewer column than body
    outfile = out_dir + '/' + os.path.basename(args.data.strip('txt')) + 'transposed.txt'
    with open(outfile, 'w') as out:
      out.write('\t' + '\t'.join([x[0] for x in data][1:]) + '\n')  # write genes
      for c in xrange(ncol):
        fields = [data[0][c]] # `c`th sample id
        fields += [data[r][c+1] for r in xrange(1, nrow)]
        out.write('\t'.join(fields) + '\n')
    data = os.path.abspath(outfile)
  else:
    data = os.path.abspath(args.data)

  print 'counting samples'
  ns = -1  # header
  for line in open(data):
    ns += 1
  
  ns = ns if args.loo else args.nboot
  out_tom_list = '{}/robust.out.toms.list.txt'.format(out_dir)
  out_tom_dir = '{}/out_toms/'.format(out_dir)
  out_tom_file = '{}/robust.tom'.format(out_tom_dir)
  out_work_dir = '{}/working_dir'.format(out_dir)
  if not os.path.exists(out_dir):
    os.mkdir(out_dir)
  if not os.path.exists(out_tom_dir):
    os.mkdir(out_tom_dir)
  if not os.path.exists(out_work_dir):
    os.mkdir(out_work_dir)
  out_rob_tom = '{}/robust.TOM.txt'.format(out_dir)
  out_tom_bash = '{}/run_LOO_TOM.bash'.format(out_dir) if args.loo else '{}/run_BOOT_TOM.bash'.format(out_dir)
  out_merge_bash = '{}/run_ROBUST_CONSENSUS_WGCNA.bash'.format(out_dir)
  
  tom_fmt = ''.join(open(args.fmt_dir + '/run_LOO_TOM.fmt')) if args.loo else ''.join(open(args.fmt_dir + '/run_BOOT_TOM.fmt'))
  tom_bash = tom_fmt.format(out_work_dir, data, args.power, args.type, out_tom_file)
  with open(out_tom_bash, 'w') as out:
    out.write(tom_bash)

  with open(out_tom_list, 'w') as out:
    for x in range(ns):
      out.write('{}.{}\n'.format(out_tom_file, 1 + x))

  merge_fmt = ''.join(open(args.fmt_dir + '/run_ROBUST_CONSENSUS_WGCNA.fmt'))
  outdir = os.path.abspath(args.out_dir)
  with open(out_merge_bash, 'w') as out:
    out.write(merge_fmt.format(out_work_dir, data, out_tom_list, args.min_module_size, args.cut_height, args.module_corr, out_dir, out_tom_list, out_tom_list))

  return out_tom_bash, out_merge_bash, ns

if __name__ == '__main__':
  args = get_args()
  out_tom_bash, out_merge_bash, ns = make_commands(args)
  print '#run::'
  pipedir = os.getcwd() if args.fmt_dir == '.' else args.fmt_dir
  print 'export INTEGRATED_PIPELINE_DIR={}'.format(pipedir)
  print 'export ALLOW_WGCNA_THREADS=5'
  print '#then::'
  print 'qsub -N mk.rob.{} -t 1:{} {}'.format(args.jsuffix, ns, out_tom_bash)
  print '#then::'
  print 'qsub -hold_jid mk.rob.{} {}'.format(args.jsuffix, out_merge_bash)
