"""\
Script for calculating phi, given gene-core distances.

The gene-core distances are pre-calculated by the R scripts; which in turn are called by 
`data_preprocessing.py` which is also responsible for defining the core gene test sets

"""
import pandas as pd
import numpy as np
from argparse import ArgumentParser

WGCNA_COLORS = ['turquoise', 'blue', 'brown', 'yellow', 'green', 'red',
                'black', 'pink', 'magenta', 'purple', 'greenyellow', 'tan',
                'salmon', 'cyan', 'midnightblue', 'lightcyan', 'grey60',
                'lightgreen', 'lightyellow', 'royalblue', 'darkred', 'darkgreen',
                'darkturquoise', 'darkgrey', 'orange', 'darkorange', 'white',
                'skyblue', 'saddlebrown', 'steelblue', 'paleturquoise', 'violet',
                'darkolivegreen', 'darkmagenta', 'sienna3', 'yellowgreen',
                'skyblue3', 'plum1', 'orangered4', 'mediumpurple3',
                'lightsteelblue1', 'lightcyan1', 'ivory', 'floralwhite',
                'darkorange2', 'brown4', 'bisque4', 'darkslateblue', 'plum2',
                'thistle2', 'thistle1', 'salmon4', 'palevioletred3',
                'navajowhite2', 'maroon', 'lightpink4', 'lavenderblush3',
                'honeydew1', 'darkseagreen4', 'coral1']


def get_args():
    parser = ArgumentParser('calculate_phi.py')
    parser.add_argument('distance_file', help='File listing gene distances, by core set')
    parser.add_argument('rdnv_gene_file', help='File listing RDNV gene sets')
    parser.add_argument('output', help='The output file to which to write')
    parser.add_argument('--core_exclusions', help='optional core exlusion file for phi calculation, '
                        'to exclude genes that were used in the calculation of phi', default=None)

    return parser.parse_args()


def calc_phi(core_genes, symbol_index, distance_vector, method='D1'):
    # normalize the distance vector
    d = (distance_vector - distance_vector.min())/(distance_vector.max()-distance_vector.min())
    # compute the cutoff
    if method is 'D1':
        cutoff = 0.1
    else:
        cutoff = np.percentile(d, 10)  # Q1
    is_Q1 = 1 * (d < cutoff)
    is_core = np.array([1 if x in core_genes else 0 for x in symbol_index])
    return np.dot(is_Q1, is_core), np.sum(is_core), np.sum(is_Q1), np.dot(is_Q1, is_core)/float(np.sum(is_core))


def is_module(str_):
    return any([('M%d' % x) in str_ for x in range(1, 10) ]) or any([x in str_ for x in WGCNA_COLORS])


def ext_module(test_name):
    entr = test_name.split('.')
    test_name = []
    for j in range(len(entr)):
        test_name.append(entr[j])
        if entr[j][0] == 'M' or entr[j] in WGCNA_COLORS:
            break
    return '.'.join(test_name).replace('kME.', '')


def main(args):
    distance_df = pd.read_csv(args.distance_file, sep='\t')
    gene_sets = dict()
    for entry in open(args.rdnv_gene_file):
        key, lst = entry.strip().split('\t')
        gene_sets[key] = lst.split(',')
    known_cores = dict()

    core_excludes = dict()
    if args.core_exclusions:
        for entry in open(args.core_exclusions):
            key, lst = entry.strip().split('\t')
            core_excludes[key] = lst.split(',')

    entries = {'distance': [], 'core': [], 'method': [], 'phi': [], 'overlap': [], 'n_core': [], 'n_D': []}
    for j in range(2, distance_df.shape[1]):
        name = distance_df.columns[j]
        if ':' in name:  # a core set
            csname = name.split('.')[0]
            excludes = core_excludes.get(csname, [])
        else:
            excludes = []
        for core_name, core_set in gene_sets.items():
            core_set = [x for x in core_set if x not in excludes]
            if len(core_set) < 15:
                # too small to really draw inferences from
                continue
            for method in ('D1', 'Q1'):
                ov, ic, id_, phi = calc_phi(core_set, distance_df['symbol'], distance_df.iloc[:,j], method)
                entries['distance'].append(name)
                entries['core'].append(core_name)
                entries['method'].append(method)
                entries['phi'].append(phi)
                entries['overlap'].append(ov)
                entries['n_core'].append(ic)
                entries['n_D'].append(id_)
    results = pd.DataFrame(entries)
    if any([is_module(c) for c in results['distance']]):
        print('pairwise...')
        entries = {'distance': [], 'core': [], 'method': [], 'phi': [], 'overlap': [], 'n_core': [], 'n_D': []}
        for core_name, core_set in gene_sets.items():
            excludes = []
            df_sub = results[results['method'] == 'Q1']
            df_sub = df_sub[df_sub['core'] == core_name]
            df_sub = df_sub.sort_values('phi', ascending=False)
            distance_cols = []
            next_best_idx = 0
            for j in range(5):
                if next_best_idx >= df_sub.shape[0]:
                    break
                # special case co-expression
                test_name = df_sub['distance'].iloc[next_best_idx]
                if is_module(test_name):
                    test_name = ext_module(test_name)
                else:
                    test_name = test_name.split('.')[0].split(':')[0]
                while not is_module(test_name) or any([test_name in prev for prev in distance_cols]):
                    print((j, test_name, is_module(test_name), ext_module(test_name)))
                    next_best_idx += 1
                    if next_best_idx >= df_sub.shape[0]:
                      break
                    test_name = df_sub['distance'].iloc[next_best_idx]
                    if is_module(test_name):
                        test_name = ext_module(test_name)
                    else:
                        test_name = test_name.split('.')[0].split(':')[0]
                distance_cols.append(df_sub['distance'].iloc[next_best_idx])
                if args.core_exclusions:
                    if ':' in distance_cols[-1]:
                        csname = distance_cols[-1].split(':')[0]
                        excludes += core_excludes.get(csname, [])
                core_set = [x for x in core_set if x not in excludes]
                distance_mat = distance_df.loc[:, np.array(distance_cols)]
                distance = np.min(distance_mat, axis=1)
                ov, ic, id_, phi = calc_phi(core_set, distance_df['symbol'], distance, 'Q1')
                if j == 0:  # already have this entry (best alone)
                    continue
                entries['distance'].append('+'.join(distance_cols))
                entries['core'].append(core_name)
                entries['method'].append('Q1')
                entries['phi'].append(phi)
                entries['overlap'].append(ov)
                entries['n_core'].append(ic)
                entries['n_D'].append(id_)
        results_cum = pd.DataFrame(entries)
        results = pd.concat([results, results_cum])

    results = results.sort_values('phi', ascending=False)

    results.to_csv(args.output, sep='\t')


if __name__ == '__main__':
    main(get_args())
