# -*- coding: utf-8 -*-
"""
Created on 2015-1-8
Binning the expression genes into several Organelle groups.


@author: WangTao
@e-mail: wangtao_bic@hit.edu.cn

"""

#################### PACKAGE NEEDED #######################
import os
import string
import sys
import re
import shutil
import numpy as np

#################### Global  PARAMETERS  ###########################
# The lowest threshold of each value in gene expression profile.
GENE_EXP_VALUE = 1
ORGAN_NUM = 0


organelle_id2go_dic = {}
organelle_go2name_dic = {'GO:0009507':'Chloroplast',
                         'GO:0005783':'Endoplasmic_Reticulum',
                         'GO:0005794':'Golgi_Apparatus',
                         'GO:0005739':'Mitochondria',
                         'GO:0005773':'Vacuole',
                         'GO:0005634':'Nucleus',
                         'GO:0009514':'Glyoxysome',
                         'GO:0005764':'Lysosome',
                         'GO:0005777':'Peroxisome',
                         'GO:0005840':'Ribosome',
                         'GO:0031982':'Vesicle',
                         'GO:0000502':'Proteasome',
                         }
# organelle_id2go_dic = {0:'GO:0009507',
#                        1:'GO:0005783',
#                        2:'GO:0005794',
#                        3:'GO:0005739',
#                        4:'GO:0005773',
#                        5:'GO:0005634',
#                        6:'GO:0009514',
#                        7:'GO:0005764',
#                        8:'GO:0005777',
#                        9:'GO:0005840',
#                        10:'GO:0031982',
#                        11:'GO:0000502',
#                        }

###########################################

def CheckExp(denList):
	flag = False
	for a in denList:
		if (a>=GENE_EXP_VALUE) | (a <= -GENE_EXP_VALUE) :
			flag = True
			break
	return flag

#****************************************************************************************#
#Function Name: ReadGeneExpDic
#Function description: Read the gene expression file into memory.
#Parameters description:
#	--gene_exp_path: string type,the path where gene expression file locates in.
#	--gene_exp_dic: Dictionary type, which contains every gene's expression profile.
#	  e.g. gene_exp_dic = { ('gene1', [0.1, 0.3...]), ('gene2', [0.2, -0.8,...] )... }   
#****************************************************************************************#

def ReadGeneExpDic(gene_exp_path, gene_exp_dic):
	print "Reading the Gene Expression File into memory..."
	try:
		gene_exp = open(gene_exp_path, 'r')
	except:
		print '\nCan not open file :' + gene_exp_path + '\n'
		sys.exit()
	try:
		while True:
			line = gene_exp.readline()
			if not line:
				break				
			linesplit = line.split()
			gene_key = linesplit[0].strip()
			tmpList = [string.atof(exp) for exp in linesplit[1:]]
			if CheckExp(tmpList) == False:
				continue
			if gene_key == 'no_match':
				continue
			gene_exp_dic[gene_key] = tmpList
	except:
		print "File Reading Failure in Function ReadGeneExpDic: "
		print gene_exp_path + " can not be recognized!"
		sys.exit()
	print "Done!"

#****************************************************************************************#
#Function Name: ReadGeneLocation
#Function description: Read the gene location file into memory. Gene location file describes
#					   which organnelle a gene locates.
#Parameters description:
#	--gene_location_file_path: string type,the path where gene location file locates in.
#	--gene_cluster: list type, elements are list type too. Each element means one organelle
#					and contains all the genes locating in that organelle. Gene location.
#	  e.g. gene_cluster = [['genei1','genei2',...], ['genej1','genej2',...]... ]   
#****************************************************************************************#
def ReadGeneLocation(gene_location_file_path, gene_cluster):
	print "Reading the gene location file (organelle annotation file)..."
	try:
		gene_mark = open(gene_location_file_path, 'r')
	except:
		print '\nCan not open file :' + gene_location_file_path + '\n'
		sys.exit()

	organelle_info = gene_mark.readline()
	infosplit = organelle_info.strip().split('\t')
	for item in infosplit:
		if item.startswith('!'):
			continue
		itemsplit = item.split(': ')
		organ_id = string.atoi(itemsplit[0]) - 1
		organ_go = itemsplit[1]
		organelle_id2go_dic[organ_id] = organ_go

	organ_num = len(organelle_id2go_dic)
	for i in range(organ_num):
	    gene_cluster.append([])

	while True:
	    line = gene_mark.readline()
	    if len(line) == 0:
	        break
	    linesplit = line.strip().split('\t')
	    cluster_id = string.atoi(linesplit[2])
	    if cluster_id != 0:
	        if linesplit[0] not in gene_cluster[cluster_id - 1]:
	            gene_cluster[cluster_id - 1].append(linesplit[0])
	print "Done!"

#****************************************************************************************#
#Function Name: ClassifyOrganelleGenes
#Function description: split all the genes into organelle bins
#Parameters description:
#	--gene_exp_dic: Dictionary type, which contains every gene's expression profile.
#	  e.g. gene_exp_dic = { ('gene1', [0.1, 0.3...]), ('gene2', [0.2, -0.8,...] )... }   
#	--gene_cluster: list type, elements are list type too. Each element means one organelle
#					and contains all the genes locating in that organelle. Gene location.
#	  e.g. gene_cluster = [['genei1','genei2',...], ['genej1','genej2',...]... ] 
#	--organ_genes: list type, its elements have the same type with gene_cluster
#****************************************************************************************#
def ClassifyOrganelleGenes(gene_exp_dic, gene_cluster, organ_genes):
	for i in range(len(gene_cluster)):
		organ_genes.append([])
	for gene in gene_exp_dic:
		for i in range(0, len(gene_cluster)):
			if gene in gene_cluster[i]:
				organ_genes[i].append(gene)

if __name__ == '__main__':

	import argparse
	cmd_parser = argparse.ArgumentParser()
	cmd_parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	cmd_parser.add_argument('-f', action = 'store', dest = "filename",
							help = "Essential,The name of source gene expression file")

	result = cmd_parser.parse_args()

	if result.filename == None:
		print "No gene expression source file exit, use argument '-h' for help!"
		sys.exit(1)
	else:
		filename = result.filename
	
	
	
	root_folderpath = './'+ filename + '/'
	if not os.path.exists(root_folderpath):
		os.mkdir(root_folderpath)
	organelle_expression_folderpath = root_folderpath + 'ExpBinning/'
	if not os.path.exists(organelle_expression_folderpath):
		os.mkdir(organelle_expression_folderpath)
	else:
		shutil.rmtree(organelle_expression_folderpath)
		os.mkdir(organelle_expression_folderpath)
	

	# Read the gene expression source file into memory. 
	gene_exp_dic = {}
	gene_exp_path = './SourceData/' + filename
	ReadGeneExpDic(gene_exp_path, gene_exp_dic)

	# Read the gene annotation file into memory, and the file tells which organelle 
	# each gene belongs.
	gene_cluster = []
	gene_location_file_path = './SourceData/gene_organelle_mark.txt'
	ReadGeneLocation(gene_location_file_path, gene_cluster)

	ORGAN_NUM = len(organelle_id2go_dic)
	if ORGAN_NUM == 0:
		print "No organelle is imported! Exit!"
		sys.exit()
	print "Number of Organelles = ", ORGAN_NUM

	# Classify all the genes gene_exp_dic contains to different organelles.
	organ_genes = []
	ClassifyOrganelleGenes(gene_exp_dic, gene_cluster, organ_genes)

	for i in range(len(organ_genes)):
		if len(organ_genes[i]) == 0:
			continue
		organelle_go = organelle_id2go_dic[i]
		organelle_name = organelle_go2name_dic[organelle_go]
		organelle_expression_file_path = organelle_expression_folderpath + organelle_name + "_Expr.txt"
		organelle_expr_file = open(organelle_expression_file_path, 'w')
		for gene in organ_genes[i]:
			outline = gene + '\t' + '\t'.join([str(num) for num in gene_exp_dic[gene]]) + '\n'
			organelle_expr_file.write(outline)
		organelle_expr_file.close()