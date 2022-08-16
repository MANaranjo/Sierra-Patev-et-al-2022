from ete3 import Tree
import sys, os
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats
from ete3 import Tree, TreeStyle, TextFace

orthofinder_dir = sys.argv[1]
if orthofinder_dir[-1] != "/" :
	orthofinder_dir = orthofinder_dir + "/"
maintree = Tree(orthofinder_dir+"Species_Tree/SpeciesTree_rooted_node_labels.txt", format=1)
dupfile = open(orthofinder_dir+"Gene_Duplication_Events/Duplications.tsv")
orthocount = open(orthofinder_dir+"Orthogroups/Orthogroups.GeneCount.tsv")

def get_duplicates(dupfile):
	orthodict = {}
	for line in dupfile:
		chunk = line.split()
		if chunk[0] == "Orthogroup":
			continue
		else:
			if chunk[0] not in orthodict:
				orthodict[chunk[0]] = {}
			if chunk[1] not in orthodict[chunk[0]]: orthodict[chunk[0]][chunk[1]] = 1
			else: orthodict[chunk[0]][chunk[1]] = orthodict[chunk[0]][chunk[1]] + 1
	return (orthodict)				
		

def get_presabs(orthocount):
	spps = {}
	spplist = []
	for line in orthocount:
		chunk = line[:-1].split("\t")
		if chunk[0] == "Orthogroup":
			for i in chunk[1:]:
				spplist.append(i)
		else:
			spps[chunk[0]] = {}
			for i in range(len(chunk[1:])):
				spps[chunk[0]][spplist[i]] = chunk[i+1]
	return (spps)

def get_zeroes(sppsentry, maintree):
	presdict = {}
	for leaf in maintree.get_leaves():
		if int(sppsentry[leaf.name]) > 0:
			presdict[leaf.name] = True
		else:
			presdict[leaf.name] = False
	for node in maintree.traverse():
		plist = []
		childs = node.get_children()
		for i in childs:
			pa = "False"
			for l in i.get_leaves():
				if int(sppsentry[l.name]) > 0:
					pa = "True"
			plist.append(pa)
		if len(plist) == 0: continue
		elif "False" not in plist:
			presdict[node.name] = True
		else: 
			presdict[node.name] = False
	return (presdict)
		

def get_duplrate(dupfile, orthocount, maintree):
	nodedict = {}
	for node in maintree.traverse():
		nodedict[node.name] = [[],[]]
	orthodict = get_duplicates(dupfile)
	for i in orthodict:
		for e in orthodict[i]:
			n = int(orthodict[i][e])
			nodedict[e][0].append(n)
			if e != "N0":
				nodedict[e][1].append(int(n)/max(maintree.search_nodes(name=e)[0].dist, 0.001))
			else:
				nodedict[e][1].append(int(n) * 0)
	spps = get_presabs(orthocount)
	for i in spps:
		presdict = get_zeroes(spps[i], maintree)
		for e in presdict:
			if presdict[e] == True:
				if i not in orthodict:
					orthodict[i] = {}
				if e not in orthodict[i]:
					nodedict[e][0].append(0)
					nodedict[e][1].append(0)

	statsdict = {}
	for i in nodedict:
		statsdict[i] = [np.mean(nodedict[i][0]), np.mean(nodedict[i][1]), scipy.stats.variation(nodedict[i][0]), scipy.stats.skew(nodedict[i][0])]
	return (statsdict, nodedict)

statsdict, nodedict = get_duplrate(dupfile, orthocount, maintree)

#print(nodedict)

for t in maintree.traverse():
	meantext = TextFace(str(round(statsdict[t.name][0], 3))+" / "+str(round(statsdict[t.name][1],3)))
	meantext.ftype = "Verdana"
	meantext.fsize = 6
	meantext.fgcolor = "DarkRed"
	numtext = TextFace(str(len(nodedict[t.name][1])-nodedict[t.name][0].count(0))+" / "+str(len(nodedict[t.name][1])))
	numtext.ftype = "Verdana"
	numtext.fsize = 6
	numtext.fgcolor = "Purple"
	devtext = TextFace(str(round(statsdict[t.name][2],3))+" / "+str(round(statsdict[t.name][3],3)))
	devtext.ftype = "Verdana"
	devtext.fsize = 6
	devtext.fgcolor = "DarkBlue"
	branchtext = TextFace(str(t.dist)+" / "+str(t.support))
	branchtext.ftype = "Verdana"
	branchtext.fsize = 6
	branchtext.fgcolor = "DarkGreen"
	t.add_face(meantext, column=0, position = "branch-top")
	t.add_face(numtext, column=0, position = "branch-top")
	t.add_face(devtext, column=1, position = "branch-bottom")
	t.add_face(branchtext, column=1, position = "branch-bottom")

maintree.show()
maintree.render("LentinulaAgaricales.pdf")
