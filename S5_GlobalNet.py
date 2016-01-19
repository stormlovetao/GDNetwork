import os
import sys
import shutil
import string
import argparse
import numpy as np
import networkx as nx
import cPickle as pickle


from S1_ExpBinning import ReadGeneExpDic


def ReadDModulesExps(D_modules_dic, D_modules_exp_dic, DModulesExps_path):
	fp = open(DModulesExps_path, 'r')
	tmp_dic = {}
	while True:
		line = fp.readline()
		if not line:
			break
		linesplit = line.strip().split('\t')
		if len(linesplit) == 1 and len(tmp_dic) == 0:
			organelle_name = linesplit[0]
		elif len(linesplit) == 1 and len(tmp_dic) > 0:
			D_modules_dic[organelle_name] = tmp_dic
			organelle_name = linesplit[0]
			tmp_dic = {}
			
		else:
			module_name = linesplit[0]
			tmp_dic[module_name] = linesplit[1:]
			nextline = fp.readline()
			nextline = nextline.strip()
			explist = [string.atof(x) for x in nextline.split('\t')]
			D_modules_exp_dic[module_name] = explist



if __name__ == '__main__':
	cmd_parser = argparse.ArgumentParser()
	cmd_parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	cmd_parser.add_argument('-f', action = 'store', dest = "filename",
							help = "Essential,The name of source gene expression file")
	cmd_parser.add_argument('-g', action = 'store', dest = "global_weight",type = float,
							help = "Essential, Float, global weight")
	result = cmd_parser.parse_args()

	if result.filename == None:
		print "No gene expression source file exit, use argument '-h' for help!"
		sys.exit(1)
	else:
		filename = result.filename
	
	if result.global_weight == None:
		print "-g global_weight, use argument '-h' for help"
		sys.exit()
	else:
		global_weight = string.atof(result.global_weight)


	GlobalNet_dir = './' + filename + '/GlobalNet/'
	if not os.path.exists(GlobalNet_dir):
		os.mkdir(GlobalNet_dir)
	# else:
	# 	shutil.rmtree(GlobalNet_dir)
	# 	os.mkdir(GlobalNet_dir)


	
	D_modules_dic = {}
	D_modules_exp_dic = {}

	#DModulesExps_path = './' + filename + '/DnodeSim/WGCNADModulesExps.txt'
	#ReadDModulesExps(D_modules_dic, D_modules_exp_dic, DModulesExps_path)
	DnodeSim_dir = './' + filename + '/DnodeSim/'
	D_modules_dic_pickle_path = DnodeSim_dir + "D_modules_dic.pickle"
	D_modules_dic = pickle.load(open(D_modules_dic_pickle_path, 'rb'))
	D_modules_exp_dic_pickle_path = DnodeSim_dir + "D_modules_exp_dic.pickle"
	D_modules_exp_dic = pickle.load(open(D_modules_exp_dic_pickle_path, 'rb'))
	

	gene_exp_dic = {}
	gene_exp_path = './SourceData/' + filename
	ReadGeneExpDic(gene_exp_path, gene_exp_dic)

	all_genes = gene_exp_dic.keys()
	genes_in_DModules = []
	for organelle_name in D_modules_dic:
		for module_name in D_modules_dic[organelle_name]:
			genes_in_DModules.extend(D_modules_dic[organelle_name][module_name])
	genes_in_DModules = list(set(genes_in_DModules))

	genes_notin_DModules = []
	for gene in all_genes:
		if gene not in genes_in_DModules:
			genes_notin_DModules.append(gene)

	Dnodes = D_modules_exp_dic.keys()

	all_nodes = genes_notin_DModules + Dnodes 
	
	GD_nodes_exp_array = []
	for gene in genes_notin_DModules:
		GD_nodes_exp_array.append(gene_exp_dic[gene])
	for Dnode in Dnodes:
		GD_nodes_exp_array.append(D_modules_exp_dic[Dnode])

	Cor_Matrix = np.corrcoef(GD_nodes_exp_array)

	GlobalNetwork = nx.Graph()

	for i in range(len(all_nodes)):
		for j in range(i+1, len(all_nodes)):
			node1 = all_nodes[i]
			node2 = all_nodes[j]
			corrcoef = Cor_Matrix[i][j]
			if corrcoef >= global_weight:
				GlobalNetwork.add_weighted_edges_from([(node1, node2, corrcoef)])

	GlobalNetwork.edgesnum = GlobalNetwork.number_of_edges()
	GlobalNetwork.nodesnum = GlobalNetwork.number_of_nodes()
	print "GlobalNetwork:"
	print "GlobalThr:\t", global_weight
	print "Edges#:\t" , GlobalNetwork.edgesnum
	print "Nodes#:\t" , GlobalNetwork.nodesnum
	
	GlobalNetwork.outpath = GlobalNet_dir + "GlobalNetwork_" + str(global_weight) + '.gml'
	nx.write_gml(GlobalNetwork, GlobalNetwork.outpath)
	print "Save GlobalNetwork into:\t" , GlobalNetwork.outpath

	




