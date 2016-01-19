# x-axis: maxflow, minimum edge cut, shortest path length between D nodes
# y-axis: GO similarity
import networkx as nx

import string
import shutil
import matplotlib.pyplot as plt
import os
import sys
import argparse

def CalMaxFlows(DG_network, Dnodes, maxFlows):
	"""
	If two nodes are connected in networkx graph G, This function returns maxflow value, shortest path length, minimum edges cut 
	between this two nodes.
	"""
	for i in range(len(Dnodes)):
		for j in range(i+1, len(Dnodes)):

			if nx.has_path(DG_network, Dnodes[i], Dnodes[j]):
				maxflow = nx.maximum_flow_value(DG_network, Dnodes[i], Dnodes[j], capacity = 'weight')
				shortest_path_length = nx.shortest_path_length(DG_network, source = Dnodes[i], target = Dnodes[j])
				min_edges_cut = len(nx.minimum_edge_cut(DG_network, Dnodes[i], Dnodes[j]))
			else:
				continue
			if Dnodes[i] < Dnodes[j]:
				a_path = (Dnodes[i], Dnodes[j])
			else:
				a_path = (Dnodes[j], Dnodes[i])

			if not maxFlows.has_key(a_path):
				maxFlows[a_path] = (maxflow, shortest_path_length, min_edges_cut)
			else:
				print "impossibly!", a_path
				sys.exit(1)

def WriteMaxFlowsSim(maxFlows, sim_dic, out_file_path):
	outfile = open(out_file_path, 'w')
	out1 = []
	out2 = []
	out3 = []
	for DDpair in maxFlows:
		D1 = DDpair[0]
		D2 = DDpair[1]
		
		result = maxFlows[DDpair]
		maxflow = result[0]
		spl = result[1]
		mec = result[2]
		if sim_dic.has_key(DDpair):
			similarity = sim_dic[DDpair]
		else:
			similarity = -1

		if D1.startswith('Nucleus') and D2.startswith('Nucleus'):
			out1.append('\t'.join([D1 , D2, str(maxflow), str(spl), str(mec), str(similarity)]))
		elif D1.startswith('Nucleus') and not D2.startswith('Nucleus'):
			out2.append('\t'.join([D1 , D2, str(maxflow), str(spl), str(mec), str(similarity)]))
		elif not D1.startswith('Nucleus') and D2.startswith('Nucleus'):
			out2.append('\t'.join([D2 , D1, str(maxflow), str(spl), str(mec), str(similarity)]))
		else:
			out3.append('\t'.join([D1 , D2, str(maxflow), str(spl), str(mec), str(similarity)]))
	outfile.write('D\tD\tmax flow\tshortest path length\tminimum edge cut\t similarity\n')
	outfile.write('\n'.join(out1))
	outfile.write('\n\n')
	outfile.write('\n'.join(out2))
	outfile.write('\n\n')
	outfile.write('\n'.join(out3))

def ReadPathsLen(DG_network, Dnodes, paths):
	
	for i in range(len(Dnodes)):
		for j in range(i+1, len(Dnodes)):
			
			if nx.has_path(DG_network, Dnodes[i], Dnodes[j]):
				shortest_path_length = nx.shortest_path_length(DG_network, source = Dnodes[i], target = Dnodes[j])
			else:
				continue

			if Dnodes[i] < Dnodes[j]:
				a_path = (Dnodes[i], Dnodes[j])
			else:
				a_path = (Dnodes[j], Dnodes[i])

			if not paths.has_key(a_path):
				paths[a_path] = shortest_path_length
			else:
				print "impossibly!", a_path
				sys.exit(1)

def SimpleReadSim(sim_file_path, sim_dic):
	sim_file = open(sim_file_path, 'r')
	organelle_pair = ()
	while True:
		line = sim_file.readline()
		if not line:
			break
		if line.startswith('#'):
			continue
		if not line.strip():
			continue
		linesplit = line.strip().split()
		if len(linesplit) == 2:
			D1 = linesplit[0].strip()
			D2 = linesplit[1].strip()
			if D1 < D2 :
				organelle_pair = (D1, D2)
			else:
				organelle_pair = (D2, D1)
		if len(linesplit) == 1:
			similarity = string.atof(linesplit[0])
			sim_dic[organelle_pair] = similarity

def averageList(lst):
	count = len(lst)
	if count == 0:
		return 0
	sum_list = 0.0
	for item in lst:
		sum_list += item
	return float(sum_list)/count

def Median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2

    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def PlotMaxflowSim(maxFlows, sim_dic, title, savepath, whichkind):
	similarity_list = []
	maxflow_list = []
	#shortest_path_length_list = []
	for pair in sim_dic:
		similarity = sim_dic[pair]
		if not maxFlows.has_key(pair):
			continue
		else:
			result = maxFlows[pair]
			if whichkind == 'maxflow':
				maxflow = result[0]
				xlabel = 'maxflow value'
			elif whichkind == 'shortestpath': 
				maxflow = result[1]  # maxflow shortestpathLen  minEdgesCut
				xlabel = 'shortestpathLen'
			elif whichkind == 'minedgecut':
				maxflow = result[2]
				xlabel = 'minEdgesCut'

			similarity_list.append(similarity)
			maxflow_list.append(maxflow)
			#shortest_path_length_list.append(shortest_path_length)
	
	if whichkind == 'maxflow':
		plt.figure()
		plt.plot(maxflow_list, similarity_list, marker = '*', linestyle = 'None')
		plt.xlabel(xlabel)
		plt.ylabel("similarity")
		plt.title(title)
		plt.xlim(xmin = -1)
		plt.ylim(ymax = 1.1)
		plt.savefig(savepath)

	if whichkind == 'shortestpath' or whichkind == 'minedgecut':
		plt.figure()
		plt.plot(maxflow_list, similarity_list, marker = '*', linestyle = 'None')
		
		dic = {}
		for i in range(len(maxflow_list)):
			min_edge_cut = maxflow_list[i]
			similarity = similarity_list[i]
			if not dic.has_key(min_edge_cut):
				dic[min_edge_cut] = [similarity]
			else:
				dic[min_edge_cut].append(similarity)

		dic_key_list = []
		dic_val_ave_list = []
		

		for key in dic:
			ave = averageList(dic[key])
			dic_key_list.append(key)
			dic_val_ave_list.append(ave)

		dic_key_list,dic_val_ave_list = (list(x) for x in zip(*sorted(zip(dic_key_list,dic_val_ave_list))))
		plt.plot(dic_key_list, dic_val_ave_list, marker = "+", label = 'average')

		dic_key_list_1 = []
		dic_val_median_list = []

		for key in dic:
			if len(dic[key]) >= 5:
				median = Median(dic[key])
				dic_key_list_1.append(key)
				dic_val_median_list.append(median)

		dic_key_list_1,dic_val_median_list = (list(x) for x in zip(*sorted(zip(dic_key_list_1,dic_val_median_list))))
		plt.plot(dic_key_list_1, dic_val_median_list, marker = "o", label = 'median')


		plt.xlabel(xlabel)
		plt.ylabel("similarity")
		plt.title(title)
		plt.legend()
		plt.xlim(xmin = -1)
		plt.ylim(ymax = 1.1)
		plt.savefig(savepath)




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
	
	gml_file_path = './' + filename + '/DGNet/DG_network_'+ str(global_weight) + '_' + str(user_given_weight) + '.gml'
	GO_sim_file_path = './' + filename + '/DnodeSim/WGCNADModulesTermsSim.txt'

	DG_network = nx.read_gml(gml_file_path, relabel = True)


	sita1 = 0.001
	sita2 = 0.001
	edges_weights = []
	for edge in DG_network.edges_iter(data = True):
		edge_weight = edge[2]['weight']
		edges_weights.append(edge_weight)
	edges_weights = sorted(edges_weights)


	minimum_edge_weight = edges_weights[0]
	maximum_edge_weight = edges_weights[len(edges_weights)-1]
	if minimum_edge_weight + sita1 >= maximum_edge_weight:
		print 'sita1 too large, exit'
		sys.exit()
	if maximum_edge_weight - sita2 <= minimum_edge_weight:
		print 'sita2 too large, exit'
		sys.exit()

	index_1 = 0
	index_2 = len(edges_weights) - 1
	for i in range(len(edges_weights)):
		if edges_weights[i] > minimum_edge_weight + sita1:
			index_1 = i
			break
		
	for i in range(len(edges_weights)):
		if edges_weights[i] > maximum_edge_weight - sita2:
			index_2 = i - 1
			break

	normalized_num = index_2 - index_1 + 1
	offset = float(1/float(normalized_num + 1))

	normalized_edges_weights_dic = {}
	tmp = 0
	for i in range(len(edges_weights)):
		if i < index_1:
			normalized_edges_weights_dic[edges_weights[i]] = 0
		elif i >= index_1 and i <= index_2:
			tmp += offset
			normalized_edges_weights_dic[edges_weights[i]] = tmp
		elif i > index_2:
			normalized_edges_weights_dic[edges_weights[i]] = 1

	#print sorted(normalized_edges_weights_dic.values())
	
	for edge in DG_network.edges_iter(data = True):
		edge_weight = edge[2]['weight']
		new_edge_weight = normalized_edges_weights_dic[edge_weight]
		edge[2]['weight'] = new_edge_weight
 
	

	# tmp_file = open('./test.txt', 'w')
	# for key in normalized_edges_weights_dic:
	# 	outline = str(key) + '\t' + str(normalized_edges_weights_dic[key]) + '\n'
	# 	tmp_file.write(outline)
	# tmp_file.close()
	# new_gml_file_path = './' + filename + '/DGNet/New_DG_network_'+ str(global_weight) + '_' + str(user_given_weight) + '.gml'
	# nx.write_gml(DG_network, new_gml_file_path)
	# sys.exit()


	Gnodes = []
	Dnodes = []
	for node in DG_network.nodes():
		if node.startswith('AT'):
			Gnodes.append(node)
		else:
			Dnodes.append(node)

	Nucleus_Dnodes = []
	nonNucleus_Dnodes = []
	for Dnode in Dnodes:
		if Dnode.startswith('Nucleus'):
			Nucleus_Dnodes.append(Dnode)
		else:
			nonNucleus_Dnodes.append(Dnode)
		
	maxFlows = {}
	CalMaxFlows(DG_network, Dnodes, maxFlows)

	sim_dic = {}
	SimpleReadSim(GO_sim_file_path, sim_dic)



	maxflow_dir = './' + filename + '/MaxFlow/'
	if not os.path.exists(maxflow_dir):
		os.mkdir(maxflow_dir)
	# else:
	# 	shutil.rmtree(maxflow_dir)
	# 	os.mkdir(maxflow_dir)

	out_file_path = maxflow_dir + filename + '_%s_%s_maxflow_shortestpath_minimumcut.txt' % (str(global_weight), str(user_given_weight))
	WriteMaxFlowsSim(maxFlows, sim_dic, out_file_path)

	plot_save_path = maxflow_dir + filename + '_' +str(global_weight) + '_' + str(user_given_weight) + '.png' 
	PlotMaxflowSim(maxFlows, sim_dic, filename, plot_save_path, 'maxflow' )

	minEdgeCut_dir = './' + filename + '/minEdgeCut/'
	if not os.path.exists(minEdgeCut_dir):
		os.mkdir(minEdgeCut_dir)
	plot_save_path = minEdgeCut_dir + filename + '_' +str(global_weight) + '_' + str(user_given_weight) + '.png' 
	PlotMaxflowSim(maxFlows, sim_dic, filename, plot_save_path, 'minedgecut' )

	shortestpath_dir = './' + filename + '/shortestPathLen/'
	if not os.path.exists(shortestpath_dir):
		os.mkdir(shortestpath_dir)
	plot_save_path = shortestpath_dir + filename + '_' +str(global_weight) + '_' + str(user_given_weight) + '.png' 
	PlotMaxflowSim(maxFlows, sim_dic, filename, plot_save_path, 'shortestpath' )