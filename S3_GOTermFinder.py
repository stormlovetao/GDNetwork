# -*- coding: utf-8 -*-
"""
Created on 2015-1-10
Enrichment analysis using GO-TermFinder


@author: WangTao
@e-mail: wangtao_bic@hit.edu.cn

"""

import os
import argparse
import sys
import shutil
import string



obo_file_path = './SourceData/gene_ontology.1_2.obo'
tair_file_path = './SourceData/gene_association.tair'


def EnrichWGCNAModules(WGCNAModules_dir, pcutoff ):
	
	enrichment_wgcna_dir = enrichment_dir + 'WGCNAModules/'
	if not os.path.exists(enrichment_wgcna_dir):
		os.mkdir(enrichment_wgcna_dir)
	else:
		shutil.rmtree(enrichment_wgcna_dir)
		os.mkdir(enrichment_wgcna_dir)

	wgcna_mfilenames = os.listdir(WGCNAModules_dir)
	wgcna_modulefiles = []
	for fname in wgcna_mfilenames:
		if fname.endswith('MGenes'):
			wgcna_modulefiles.append(fname)
	
	study_file_paths_list = []
	
	for fname in wgcna_modulefiles:
		modules_genes_dic = {}
		module_genes_file_path = WGCNAModules_dir + fname
		fpm = open(module_genes_file_path, 'r')
		while  True:
			line = fpm.readline()
			if not line:
				break
			linesplit = line.strip().split('\t')
			module_name = linesplit[0]
			modules_genes_dic[module_name] = linesplit[1:]
		fpm.close()

		for module_name in modules_genes_dic:
			study_file_path = enrichment_wgcna_dir + module_name
			study_file_paths_list.append(study_file_path)
			study_file_p = open(study_file_path, 'w')
			genes_list = modules_genes_dic[module_name]
			outline = '\n'.join(genes_list)
			study_file_p.write(outline)
			study_file_p.close()
	command = ('perl ./GO-TermFinder-0.86/examples/analyze.pl  %f  ' + tair_file_path + ' ' + backgroundgene_file_path + ' ' +  obo_file_path + ' ' + ' '.join(study_file_paths_list))  % (pcutoff)

 	os.system(command)




if __name__ == '__main__':
	
	cmd_parser = argparse.ArgumentParser()
	cmd_parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	cmd_parser.add_argument('-f', action = 'store', dest = "filename",
							help = "Essential,The name of source gene expression file")
	cmd_parser.add_argument('-p', action = 'store', dest = "pcutoff",
							help = "Optional,default = 0.05, P value cut off")
	result = cmd_parser.parse_args()

	if result.filename == None:
		print "No gene expression source file exit, use argument '-h' for help!"
		sys.exit(1)
	else:
		filename = result.filename

	if result.pcutoff == None:
		pcutoff = 0.05
	else:
		pcutoff = string.atof(result.pcutoff)

	enrichment_dir = './' + filename + '/Enrichment/'
	
	if not os.path.exists(enrichment_dir):
		os.mkdir(enrichment_dir)
	else:
		shutil.rmtree(enrichment_dir)
		os.mkdir(enrichment_dir)

	backgroundgene_file_path = './SourceData/' + filename

	WGCNAModules_dir = './' + filename + '/WGCNA/Modules/'
	EnrichWGCNAModules(WGCNAModules_dir, pcutoff)

	
