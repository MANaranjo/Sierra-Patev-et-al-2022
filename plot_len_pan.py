import sys, os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import scipy
from Bio import SeqIO
import sys
		
orthofinder_dir = sys.argv[1]
if orthofinder_dir[-1] == "/":
	orthofinder_dir = orthofinder_dir[:-1]
sppfile, [] = sys.argv[2], []
accfile, accdict, coredict, acclist = sys.argv[3], {}, {}, []

for sp in open(sppfile):
	accdict[sp[:-1]], coredict[sp[:-1]] = [], []

for og in open(accfile):
	acclist.append(og[:-1])

for fasta in os.listdir(orthofinder_dir+'/Orthogroup_Sequences/'):
	ogfasta = SeqIO.index(orthofinder_dir+'/Orthogroup_Sequences/'+fasta, "fasta")
	acc = False
	if fasta.split(".")[0] in acclist:
		acc = True
	for seq in ogfasta:
		for sp in accdict:
			if seq.split('|')[1].find(sp) > -1:
				if acc == True: accdict[sp].append(len(ogfasta[seq].seq))
				else: coredict[sp].append(len(ogfasta[seq].seq))

alldataset = {"length":[], "spp":[], "core/acc":[]}

for i in coredict:
	for e in coredict[i]:
		alldataset["length"].append(int(e))
		alldataset["spp"].append(i)
		alldataset["core/acc"].append("Core")
for i in accdict:
	for e in accdict[i]:
		alldataset["length"].append(int(e))
		alldataset["spp"].append(i)
		alldataset["core/acc"].append("Acc")

alldataset = pd.DataFrame(alldataset)

fig, ax = plt.subplots(figsize=(36, 10))

corearray, accarray, namearray = [], [], []
for i in coredict:
	corearray.append(coredict[i])
	namearray.append(i)
corearray = np.asarray(corearray)

for i in accdict:
	accarray.append(accdict[i])
accarray = np.asarray(accarray)
ax = sns.violinplot(x="spp", y="length", hue="core/acc", data=alldataset, palette="Set2", split=False, scale="count", inner="quartile")

ax.yaxis.grid(False)
ax.set_xticks([y for y in range(len(corearray))])

plt.setp(ax, xticks=[y for y in range(len(corearray)+1)],
         xticklabels=namearray)
plt.xticks(rotation=45)
plt.ylim(0,1200)

plt.savefig("lenplot.png")
