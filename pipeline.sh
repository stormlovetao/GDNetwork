#! /bin/bash
filename="x_final.txt"
echo $filename
global_weight=0.96
user_weight=0.98
minModuleSize=5
softPower=1
pval_cutoff=0.05
python S1_ExpBinning.py   -f $filename
python S2_WGCNA.py        -f $filename -s $minModuleSize -sp $softPower               # -s minModuleSize   -p softPower
python S3_GOTermFinder.py -f $filename -p $pval_cutoff                               # -p pval_cutoff
python S4_DnodeSim.py     -f $filename
python S5_GlobalNet.py    -f $filename -g $global_weight                             # -g global_weight
python S6_DGNet.py        -f $filename -g $global_weight -u $user_weight

python E1_MaxFlow.py      -f $filename -g $global_weight -u $user_weight
