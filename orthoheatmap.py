
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from ete3 import Tree
import scipy
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', required=True, help="Fasta file containing the assembly to analyze")
	parser.add_argument('-r', '--relation', choices=["1to1", "1toM", "MtoM"], default="1to1")
	
args = parser.parse_args()


orthofinder_dir = args.input
if orthofinder_dir[-1] != "/" :
	orthofinder_dir = orthofinder_dir + "/"
orthodir = orthofinder_dir	

gene_numbers_list = []

data1 = open(orthodir+"Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv")
if args.relation == "1to1":
	data = open(orthodir+"Comparative_Genomics_Statistics/OrthologuesStats_one-to-one.tsv")
	N = 10
if args.relation == "1toM":
	data = open(orthodir+"Comparative_Genomics_Statistics/OrthologuesStats_one-to-many.tsv")
	N = 0
if args.relation == "MtoM":
	data = open(orthodir+"Comparative_Genomics_Statistics/OrthologuesStats_many-to-many.tsv")
	N = 0

t = Tree(orthodir+"Species_Tree/SpeciesTree_rooted_node_labels.txt", format=1)
t.convert_to_ultrametric()


switch = False
for line in data1:
	chunk = line[:-1].split("\t")
	if switch == False:
		if chunk[0] == "Number of genes":
			gene_numbers_list.append(chunk[1:])
			switch = True
		else:
			gene_numbers_list.append(chunk[1:])
	else:
		break

	
gene_numbers = {}
for i in range(len(gene_numbers_list[0])):
	gene_numbers[gene_numbers_list[0][i]] = int(gene_numbers_list[1][i])

taxalist = []
for node in t.get_leaves():
	taxalist.append(node.name)

taxalist.sort()

distmatrix = []
for x in taxalist:
	distrow = []
	for y in taxalist:
		if x == y: distrow.append(0.0)
		else:
			distrow.append(t.get_distance(x,y))
	distmatrix.append(distrow)

datamatrix = []
for line in data:
	datarow = []
	chunk = line.split("\t")
	datarow = []
	for i in chunk[1:]:
		if i.find("\n") > -1: i = i[:i.find("\n")]
		if i in taxalist:
			datarow.append(i)
		else:
			datarow.append(float(i))
	datamatrix.append(datarow)

newmatrix = []

for i in range(1,len(datamatrix)):
	newrow = []
	for n in datamatrix[i]:
		newrow.append( n * 100 / gene_numbers[datamatrix[0][i-1]])
	newmatrix.append(newrow)

matrixframe = pd.DataFrame(newmatrix, columns=datamatrix[0], index=datamatrix[0])
print (matrixframe)
mask = np.zeros_like(matrixframe, dtype=np.bool)
mask[np.triu_indices_from(mask)] = True

d = scipy.cluster.hierarchy.linkage(distmatrix)

sns.set(font_scale=0.5)
sns.clustermap(matrixframe, linewidths=0, cmap="RdBu_r", vmin=N, annot=True, annot_kws={"size": 6}, row_linkage=d, col_linkage=d)
plt.savefig(args.relation+".svg")
plt.show()
