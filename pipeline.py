import sys
import os
import argparse
import numpy as np

minModuleSize=5
softPower=6
pval_cutoff=0.05

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

	global_weight_list = np.arange(0.96, 0.961, 0.01)
	
	for gw in global_weight_list:
		for uw in np.arange(gw+0.01, 1.001, 0.01):
			print gw, uw
			command = 'python S1_ExpBinning.py   -f %s' % filename
			#os.system(command)
			command = 'python S2_WGCNA.py        -f %s -s %d -sp %d ' % (filename, minModuleSize, softPower)
			#os.system(command)
			command = 'python S3_GOTermFinder.py -f %s -p %f ' % (filename, pval_cutoff)
			#os.system(command)
			command = 'python S4_DnodeSim.py     -f %s' % filename
			#os.system(command)
			command = 'python S5_GlobalNet.py    -f %s -g %f' % (filename, gw)
			#os.system(command)
			command = 'python S6_DGNet.py        -f %s -g %f -u %f' % (filename, gw, uw)
			os.system(command)
			command = 'python E1_MaxFlow.py      -f %s -g %f -u %f' % (filename, gw, uw)
			os.system(command)

	