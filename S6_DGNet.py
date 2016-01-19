import os
import sys
import argparse
import shutil
import string
import numpy as np
import networkx as nx

def averageList(lst):
	count = len(lst)
	if count == 0:
		return 0
	sum_list = 0.0
	for item in lst:
		sum_list += item
	return float(sum_list)/count

def MaxDepth(core_weight):
	depth = 1
	while True:
		if InfluenceCoefficient(depth, core_weight) <= 0.2:
			break
		else:
			depth += 1
	return depth
def InfluenceCoefficient(d, core_weight):
	return 1/(1+ np.exp(d - 5*core_weight))
	
def Judge(potential_edges, G, influence, inner_nodes, saved_edges, removed_edges):	
	outer_nodes = []
	for edge in potential_edges:
		edge_weight = edge[2]
		if abs(edge_weight) < influence:
			removed_edges[(edge[0], edge[1])] = edge_weight
			G.remove_edge(edge[0], edge[1])
		else:
			saved_edges[(edge[0], edge[1])] = edge_weight
			if edge[0] not in inner_nodes:
				outer_nodes.append(edge[0])
			elif edge[1] not in inner_nodes:
				outer_nodes.append(edge[1])
	return outer_nodes

def CoreInfluence(G, D, core_weight, basic_weight, max_depth):
	distance = 0
	LND = [D]
	saved_edges = {}
	removed_edges = {}
	while G.number_of_edges() != 0:
		potential_edges = []
		distance += 1
		if distance > max_depth:
			break
		influence_coefficient = InfluenceCoefficient(distance, core_weight)
		influence = influence_coefficient * core_weight + (1 - influence_coefficient)*basic_weight
		#print 'LND = ', LND
		for node in LND:
			node_neighbors_dic = G[node]
			for neighbor in node_neighbors_dic:
				weight = G[node][neighbor]['weight']
				if node < neighbor:
					edge = (node, neighbor, weight)
					if edge not in potential_edges:
						potential_edges.append(edge)
				else:
					edge = (neighbor, node, weight)
					if edge not in potential_edges:
						potential_edges.append(edge)
		next_layer_nodes = Judge(potential_edges, G, influence, LND, saved_edges, removed_edges)
		if len(next_layer_nodes) == 0:
			break
		else:
			G.remove_nodes_from(LND)
			LND = next_layer_nodes
	return saved_edges, removed_edges

if __name__ == '__main__':
	
	cmd_parser = argparse.ArgumentParser()
	cmd_parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	cmd_parser.add_argument('-f', action = 'store', dest = "filename",
							help = "Essential,The name of source gene expression file")
	cmd_parser.add_argument('-g', action = 'store', dest = "global_weight",type = float,
							help = "Essential, Float, global weight")
	cmd_parser.add_argument('-u', action = 'store', dest = "user_given_weight",type = float,
							help = "Essential, Float, user_given_weight")

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
	if result.user_given_weight == None:
		print "-u user_given_weight, use argument '-h' for help"
		sys.exit()
	else:
		user_given_weight = string.atof(result.user_given_weight)
	
	DGNet_dir = './' + filename + '/DGNet/'
	if not os.path.exists(DGNet_dir):
		os.mkdir(DGNet_dir)
	# else:
	# 	shutil.rmtree(DGNet_dir)
	# 	os.mkdir(DGNet_dir)

	GlobalNet_dir = './' + filename + '/GlobalNet/'
	DnodeSim_dir = './' + filename + '/DnodeSim/'
	DModulesTerms_path = DnodeSim_dir + 'WGCNADModulesTerms.txt'
	Dnodes = []
	fp = open(DModulesTerms_path, 'r')
	while True:
		line = fp.readline()
		if not line:
			break
		linesplit = line.split('\t')
		Dnodes.append(linesplit[0])
	fp.close()

	
	GlobalNetwork_path = GlobalNet_dir + 'GlobalNetwork_' + str(global_weight) + '.gml'
	ini_network = nx.read_gml(GlobalNetwork_path, relabel = True)
	

	# regard the average weight of edges around D node as D's outer_Core_weight
	outer_D_core_weights = {}
	
	for D_name in Dnodes:
		if ini_network.has_node(D_name):
			weights_around_D = []
			for neighbor_node in ini_network[D_name]:
				weights_around_D.append(ini_network[D_name][neighbor_node]['weight'])
			average_weight = averageList(weights_around_D)
			outer_D_core_weights[D_name] =  average_weight
	print outer_D_core_weights
	sys.exit()
	# For each D node, sparse the ini_network
	
	saved_edges_list = []
	removed_edges_list = []
	for D_name in outer_D_core_weights:
		D_core_network = ini_network.copy()
		for D_name_2 in outer_D_core_weights:
			if (D_name_2 != D_name) and (D_core_network.has_node(D_name_2)):
				D_core_network.remove_node(D_name_2)
		max_depth = MaxDepth(outer_D_core_weights[D_name])
		print "In function CoreInfluence..."
		saved_edges, removed_edges = CoreInfluence( D_core_network, D_name, outer_D_core_weights[D_name], user_given_weight, max_depth)
		saved_edges_list.append(saved_edges)
		removed_edges_list.append(removed_edges)

	# D_index = 0
	# for D_name in outer_D_core_weights:
	# 	saved_edges = saved_edges_list[D_index]
	# 	D_index += 1
	# 	this_D_network = nx.Graph()
	# 	for edge_pair in saved_edges:
	# 		node1 = edge_pair[0]
	# 		node2 = edge_pair[1]
	# 		weight = saved_edges[edge_pair]
	# 		this_D_network.add_weighted_edges_from([(node1, node2, weight)])
	# 	this_D_network_path = organelle_folderpath + '/' + results_label + D_name 
	# 	nx.write_gml(this_D_network, this_D_network_path)



	print "edges_label_dic..."
	edges_label_dic = {}
	for edge in ini_network.edges_iter():
		if edge[0] < edge[1]:
			tmp_edge = edge
		else:
			tmp_edge = (edge[1], edge[0])
		this_edge_label = []
		for i in range(len(outer_D_core_weights)):
			if saved_edges_list[i].has_key(tmp_edge):
				this_edge_label.append(0)
			elif removed_edges_list[i].has_key(tmp_edge):
				this_edge_label.append(1)
			else:
				this_edge_label.append(2)
		edges_label_dic[tmp_edge] = this_edge_label

	# print "write edges_label_dic..."
	# tmp_file_path = organelle_folderpath + '/' + results_label + "ini_networkEdges_saveornot.label"
	# tmp_file = open(tmp_file_path, 'w')
	# for key in edges_label_dic:
	# 	outline = `key` + '\t' + `edges_label_dic[key]` + '\n'
	# 	tmp_file.write(outline)
	# tmp_file.close()

	final_network = ini_network.copy()
	#tmp_file_path = organelle_folderpath + '/' + results_label + "Edges_not_been_judged_" +".txt"
	#tmp_file = open(tmp_file_path, 'w')
	edge_not_labeled = 0
	for edge in edges_label_dic:
		edge_label = edges_label_dic[edge]
		aye_vote = edge_label.count(0)
		against_vote = edge_label.count(1)
		abstention_vote = edge_label.count(2)
		if aye_vote < against_vote:
			final_network.remove_edge(edge[0], edge[1])
		elif aye_vote == 0 and against_vote == 0:
			#print "this edge has never been judged: " , `edge`
			final_network.remove_edge(edge[0], edge[1])
			outline = `edge` + '\n'
			#tmp_file.write(outline)
			edge_not_labeled += 1
		#elif aye_vote == against_vote:
			#final_network.remove_edge(edge[0], edge[1])

	print "Edges not labeled: " + `edge_not_labeled`
	#tmp_file.close()

	# write D-D edges into final_network
	#tmp_file_path = root_folderpath + '/' + "!lookup"
	#tmp_file = open(tmp_file_path, 'w')
	modules_dic = {}
	modules_exp_dic = {}
	from S2_WGCNA import ReadWGCNARawModules
	wgcna_modules_dir = './' + filename + '/WGCNA/RawModules/'
	ReadWGCNARawModules(wgcna_modules_dir, modules_dic, modules_exp_dic)
	
	Dnames = outer_D_core_weights.keys()
	for i in range(len(Dnames)):
		for j in range(i+1, len(Dnames)):
			D1 = Dnames[i]
			D2 = Dnames[j]
			D1_core_weight = outer_D_core_weights[D1]
			D2_core_weight = outer_D_core_weights[D2]
			if D1 < D2:
				D_pair = (D1, D2)
			else:
				D_pair = (D2, D1)
			
			DD_weight = np.corrcoef(modules_exp_dic[D1], modules_exp_dic[D2])[0][1]
			if DD_weight <= 0:
				continue
			D1_influence_coefficient = InfluenceCoefficient(1, D1_core_weight)
			D1_influence = D1_influence_coefficient * D1_core_weight + (1 - D1_influence_coefficient)*user_given_weight
			D2_influence_coefficient = InfluenceCoefficient(1, D2_core_weight)
			D2_influence = D2_influence_coefficient * D2_core_weight + (1 - D2_influence_coefficient)*user_given_weight
			if DD_weight >= D1_influence and DD_weight >= D2_influence:
				final_network.add_weighted_edges_from([(D1, D2, DD_weight)])
			#outline = `D_pair` + '\t' + `DD_weight` + '\t' + `D1_influence` + '\t' + `D2_influence` + '\n'
			#tmp_file.write(outline)
	#tmp_file.close()


	final_network_path = DGNet_dir + "DG_network_" + str(global_weight)+'_'+ str(user_given_weight)  + ".gml"
	nx.write_gml(final_network, final_network_path)
