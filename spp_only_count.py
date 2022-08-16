import sys, os
from ete3 import Tree

orthofinder_dir = sys.argv[1]
if orthofinder_dir[-1] != "/" :
	orthofinder_dir = orthofinder_dir + "/"

tree_dir = orthofinder_dir+"Resolved_Gene_Trees/"
spp = sys.argv[2]
A = 0
for i in os.listdir(tree_dir):
	B = True
	tree = Tree(tree_dir+i, format=1)
	for n in Tree.get_leaf_names(tree):
		if n.split("_")[0] != spp and n.find(spp) == -1: B = False
	if B == True: 
		A = A + 1


print(spp+" only:\t", A)
