import os
import sys
import shutil
import string
import argparse
import cPickle as pickle
from S2_WGCNA import ReadWGCNARawModules

def Modules2DModules(modules_dic, modules_exp_dic, D_modules_dic, D_modules_exp_dic, enrichment_dir):
	for organelle_name in modules_dic:
		tmp_dic = {}
		if len(modules_dic[organelle_name]) == 0:
			continue
		for module_name in modules_dic[organelle_name]:
			if(CheckEnrichment(module_name, enrichment_dir)):
				tmp_dic[module_name] = modules_dic[organelle_name][module_name]
				D_modules_exp_dic[module_name] = modules_exp_dic[module_name]

		D_modules_dic[organelle_name] = tmp_dic

def CheckEnrichment(module_name, enrichment_dir):
	module_enrichment_filepath = enrichment_dir + module_name + '.terms'
	terms_list = GOTermFinderEnrichmentTermsReader(module_enrichment_filepath)
	if len(terms_list) > 0 :
		return True
	else:
		return False

def GOTermFinderEnrichmentTermsReader(enrichment_file_path):
	enrichment_file = open(enrichment_file_path, 'r')
	terms_list = []
	while True:
		line = enrichment_file.readline()
		if not line or line.startswith('Finding terms for C'):
			break
		if line.startswith('GOID'):
			linesplit = line.strip().split()
			term = linesplit[1]
			if term != 'GO:XXXXXXX':
				terms_list.append(term)
	
	return terms_list

def WriteDModulesTerms(Dmodule_names, enrichment_dir, DModulesTerms_path):
	fp = open(DModulesTerms_path, 'w')
	for module_name in Dmodule_names:
		module_enrichment_filepath = enrichment_dir + module_name + '.terms'
		terms_list = GOTermFinderEnrichmentTermsReader(module_enrichment_filepath)
		if len(terms_list) > 0 and  'GO:XXXXXXX' not in terms_list:
			outline = module_name + '\t' + '\t'.join(terms_list) + '\n' 
			fp.write(outline)
	fp.close()

def WriteDModulesExps(D_modules_dic, D_modules_exp_dic, DModulesExps_path):
	fp = open(DModulesExps_path, 'w')
	for organelle_name in D_modules_dic:
		outline = organelle_name + '\n'
		fp.write(outline)

		for module_name in D_modules_dic[organelle_name]:
			outline = module_name + '\t' + '\t'.join(D_modules_dic[organelle_name][module_name]) + '\n'
			fp.write(outline)
			exp = D_modules_exp_dic[module_name]
			outline = '\t'.join([str(x) for x in exp]) + '\n'
			fp.write(outline)

if __name__ == '__main__':
	

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

	DnodeSim_dir = './' + filename + '/DnodeSim/'
	if not os.path.exists(DnodeSim_dir):
		os.mkdir(DnodeSim_dir)
	else:
		shutil.rmtree(DnodeSim_dir)
		os.mkdir(DnodeSim_dir)

	modules_dic = {}
	modules_exp_dic = {}
	D_modules_dic = {}
	D_modules_exp_dic = {}
	
	wgcna_modules_dir = './' + filename + '/WGCNA/RawModules/'
	wgcna_mdoules_enrichment_dir = './' + filename + '/Enrichment/WGCNAModules/'
	ReadWGCNARawModules(wgcna_modules_dir, modules_dic, modules_exp_dic)
	modules_enrichment_dir =  wgcna_mdoules_enrichment_dir

	Modules2DModules(modules_dic, modules_exp_dic, D_modules_dic, D_modules_exp_dic, modules_enrichment_dir)

	print "Organelles modules number and enriched modules number"
	for organelle_name in D_modules_dic:
		print organelle_name 
		print len(modules_dic[organelle_name]) , len(D_modules_dic[organelle_name])

	Dnodes = D_modules_exp_dic.keys()
	
	DModulesTerms_path = DnodeSim_dir + "WGCNADModulesTerms.txt"
	WriteDModulesTerms(Dnodes, modules_enrichment_dir, DModulesTerms_path)

	D_modules_dic_pickle_path = DnodeSim_dir + "D_modules_dic.pickle"
	pickle.dump(D_modules_dic, open(D_modules_dic_pickle_path, 'wb'), True)

	D_modules_exp_dic_pickle_path = DnodeSim_dir + "D_modules_exp_dic.pickle"
	pickle.dump(D_modules_exp_dic, open(D_modules_exp_dic_pickle_path, 'wb'), True)

	#DModulesExps_path =  DnodeSim_dir + "WGCNADModulesExps.txt"
	#WriteDModulesExps(D_modules_dic, D_modules_exp_dic, DModulesExps_path)

	command = 'python DnodeSim/DnodeSim.py -f %s' % filename
	os.system(command)
