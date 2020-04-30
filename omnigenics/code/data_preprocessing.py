"""\


Networks defined in `network_defs` will be converted to distance files


"""
import os
import pandas as pd

output_dir = 'processed_data'

CORE_GENE_N = [25, 35, 50, 75, 100]
CORE_GENE_FORDISTANCE_N = [10, 20]

R_SCRIPT='Rscript'


syndromic_genes = {
	'ASD': ['FMRP', 'ANK2', 'SYNGAP1', 'CHD8', 'SHANK2', 'SHANK3', 'SCN2A'],
	'SCZ': ['ATP1A3', 'FXYD', 'NRXN1', 'RB1CC1']  # 10.1038/s41380-018-0103-8 + 10.1038/ng.3725
}


# just by happenstance, all of these use HGNC symbols
raw_data_files = {
	'iHart': {
		'file': 'omnigenics/data2/re.TADA.ASC_SSC_iHART.ms1_SandersSmallDel_finalFDR_values.withSim.txt',
		'key': 'BF.total',
		'gene_key': 'gene.id',
		'sep': '\t',
		'log': True
	},
	'extTADA_ASD': {
		'file': 'omnigenics/data2/Stahl.extTADA.ASD.csv',
		'sep': ',',
		'key': 'BF',
		'gene_key': "'geneName'",
		'log': True
	},
	'extTADA_SCZ': {
		'file': 'omnigenics/data2/Stahl.extTADA.SCZ.csv',
		'sep': ',',
		'key': 'BF',
		'gene_key': "'Gene'",
		'log': True
	},
	'extTADA_EPI': {
		'file': 'omnigenics/data2/Stahl.extTADA.EPI.csv',
		'sep': ',',
		'key': 'BF',
		'gene_key': 'geneName',
		'log': True
	},
	'extTADA_DD': {
		'file': 'omnigenics/data2/Stahl.extTADA.DD.csv',
		'sep': ',',
		'key': 'BF',
		'gene_key': 'geneName',
		'log': True
	},
	'extTADA_ID': {
		'file': 'omnigenics/data2/Stahl.extTADA.ID.csv',
		'sep': ',',
		'key': 'BF',
		'gene_key': 'geneName',
		'log': True
	},
	'DuWu_ASD': {
		'file': 'omnigenics/data2/DuWu2019.s41436-019-0610-2.qscores.csv',
		'sep': ',',
		'key': 'Qscore',
		'gene_key': 'Gene_symbol',
		'log': False
	}
}

# because core gene files use HGNC, so do the filter files
filter_files = {
	'protein_coding_noTF_noDNABP_noRNABP': 'omnigenics/data2/human_proteins.noTF.noDNA_BP.noRNA_BP.hgnc.txt'
}

network_defs = {
	'WHOLE_BRAIN': {
		'mod_file': 'omnigenics/data2/WHOLE_BRAIN.modules.mapped.txt',
		'expression_file': 'omnigenics/data2/GTEx.WholeBrain.centered.expr.txt',
		'type': 'coexpression'
		
	},

	'RC_Neuron': {
		'binding_graph': 'omnigenics/data2/neurons.txt.gz',
		'type': 'graph_edges',
		'convert_to_ensembl': True,
	},

	'InWeb': {
		'matrix': 'omnigenics/data2/inweb_brain.ALL.txt',
		'type': 'graph_matrix'
	}
}


def extract_distance(key, info_dict, c_str, Rscript=R_SCRIPT):
	out_file = 'omnigenics/outputs/%s.distances.txt' % key
	if info_dict['type'] == 'coexpression':
		cmd = Rscript + ' omnigenics/code/extract_coexpression_distance.R %s %s %s %s' % (info_dict['mod_file'], info_dict['expression_file'], out_file, c_str)
	elif info_dict['type'] == 'graph_edges':
		cmd = Rscript + ' omnigenics/code/extract_TF_network_distance.R %s %s %s' % (info_dict['binding_graph'], out_file, c_str)
	else:
		cmd = Rscript + ' omnigenics/code/extract_PPI_network_distance.R %s %s %s' % (info_dict['matrix'], out_file, c_str)
	return cmd


def rmquote(s):
	return s.replace("'", "").replace('"','')


def extract_core_genes(file_key, file_info):
	data = pd.read_csv(file_info['file'], sep=file_info['sep'])
	data = data.sort_values(file_info['key'], ascending=False)
	filters = {fkey: {rmquote(x.strip()) for x in open(ffile)} for fkey, ffile in filter_files.items()}
	# first extract the sets without filters
	core_genes = dict()
	for N in CORE_GENE_N:
		core_key = '%s:phi:core_%d' % (file_key, N)
		core_list_raw = [rmquote(x) for x in data.loc[:N, file_info['gene_key']]]
		core_genes[core_key] = core_list_raw
		for fkey, whitelist in filters.items():
			filter_key = '%s:%s' % (core_key, fkey)
			filter_genes = [x for x in core_list_raw if x in whitelist]
			core_genes[filter_key] = filter_genes
	for N in CORE_GENE_FORDISTANCE_N:
		core_key = '%s:known:core_%d' % (file_key, N)
		core_list_raw = [rmquote(x) for x in data.loc[:N, file_info['gene_key']]]
		core_genes[core_key] = core_list_raw
		for fkey, whitelist in filters.items():
			filter_key = '%s:%s' % (core_key, fkey)
			filter_genes = [x for x in core_list_raw if x in whitelist]
			core_genes[filter_key] = filter_genes
	return core_genes

def dictmerge(dict_a, dict_b):
	return {k: v for dct in [dict_a, dict_b] for k, v in dct.items()}


ALL_CORE_GENES = dict()
for datfile in ['iHart', 'extTADA_ASD', 'extTADA_SCZ', 'DuWu_ASD']:
	ALL_CORE_GENES = dictmerge(ALL_CORE_GENES, extract_core_genes(datfile, raw_data_files[datfile]))

core_gene_str = [x for x in ALL_CORE_GENES.keys() if 'known' in x]
core_gene_str = [x + '=' + ','.join(ALL_CORE_GENES[x]) for x in core_gene_str]
core_gene_str = '"' + ';'.join(core_gene_str) + '"'

out1 = open('omnigenics/outputs/core_genes.known.for_distance.txt', 'w')
out2 = open('omnigenics/outputs/core_genes.hypothetical.for_phi.txt', 'w')

for corekey, corelist in ALL_CORE_GENES.items():
	ostr = corekey + '\t' + ','.join(corelist)
	if 'known' in corekey:
		out1.write(ostr + '\n')
	else:
		out2.write(ostr + '\n')
out1.close()
out2.close()

for netkey, netinfo in network_defs.items():
	command = extract_distance(netkey, netinfo, core_gene_str)
	print(command)
	os.system(command)


