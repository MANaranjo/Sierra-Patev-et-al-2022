import sys, os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

orthofinder_dir = sys.argv[1]
if orthofinder_dir[-1] != "/" :
	orthofinder_dir = orthofinder_dir + "/"
import seaborn as sns


ogfile = open(orthofinder_dir+"/Orthogroups/Orthogroups.GeneCount.tsv")
ortho = open(orthofinder_dir+"/Orthogroups/Orthogroups.tsv")
sppfile = sys.argv[2]
num = 3
if len(sys.argv) == 4:
	num = int(sys.argv[3])

spplist, spplist2 = [], []
for line in open(sppfile): spplist.append(line[:-1])

acclist, univlist, commonlist, Bdict= [], [], [], {}

def Byoungnam_criteria(c2, spplist2):
	A = str(spplist2)
	nw = A.count("_1'")+A.count("_2'")+A.count("_3'")+A.count("_4'")+A.count("_5'")
	ow = A.count("_6'")+A.count("_8'")+A.count("_9'")+A.count("_10'")+A.count("_11'")
	veredict = 'Other'
	if c2[1:].count("0") == len(c2[1:])-1:
		for e in c2[:-1]:
			if e != "0":
						strain = spplist2[c2.index(e)]
		veredict = strain+"-specific"
	if c2.count("0") == 0:
		veredict = "Core"
	ITSlist = []
	for i in range(1, len(c2)):
		if c2[i-1] != "0":
			ITSlist.append(spplist2[i-1].split("_")[-1])		
	if len(ITSlist[0]) < 3 and len(set(ITSlist)) == 1 and len(ITSlist) == str(spplist2).count("_"+ITSlist[0]+"'"):
		veredict = "ITS"+ITSlist[0]+"-specific"
	if str(ITSlist).count("'1'")+str(ITSlist).count("'2'")+str(ITSlist).count("'3'")+str(ITSlist).count("'4'")+str(ITSlist).count("'5'") == ow:
		veredict = "OldWorld-specific"
	if str(ITSlist).count("'6'")+str(ITSlist).count("'8'")+str(ITSlist).count("'9'")+str(ITSlist).count("'10'")+str(ITSlist).count("'11'") == nw:
		veredict = "NewWorld-specific"
	return veredict


for line in ogfile:
	chunk = line.split("\t")
	if chunk[0] == "Orthogroup":
		for s in chunk[1:]:
			spplist2.append(s)
	else:
		c1, c2 = [], []
		for i in range(len(chunk[1:])):
			c2.append(chunk[i])
			if spplist2[i] in spplist:
				c1.append(chunk[i])
		Bdict[chunk[0]] = Byoungnam_criteria(c2, spplist2)
		if c2.count("0") == 0 and c2.count("1") < len(c2):
			univlist.append(chunk[0])
		elif c1.count("0") == 0:
			commonlist.append(chunk[0])
		elif c1.count("0") >= num and c1.count("0") <= len(c1)-num:
			acclist.append(chunk[0])
		else: continue


data = open(orthofinder_dir+"phylodata.tsv", "w")
dataB = open(orthofinder_dir+"phylodataB.tsv", "w")
header = "gene;species;orthogroup;type\n"
data.write(header)
dataB.write(header)


spplist3 = []
for line in ortho:
	line = line[:-1]
	if line.find("Orthogroup") == 0:
		spplist3 = line.split("\t")
	else:
		chank = line.split("\t")
		for n in range(1, len(chank)):
			for m in chank[n].split(", "):
				if len(m) != 0:
					if chank[0] in acclist:
						data.write(m+";"+spplist3[n]+";"+chank[0]+";acc\n")
					if chank[0] in commonlist:
						data.write(m+";"+spplist3[n]+";"+chunk[0]+";common\n")
					if chank[0] in univlist:
						data.write(m+";"+spplist3[n]+";"+chank[0]+";univ\n")
					else:
						data.write(m+";"+spplist3[n]+";"+chank[0]+";none\n")
					dataB.write(m+";"+spplist3[n]+";"+chank[0]+";"+Bdict[chank[0]]+"\n")
data.close()

df = pd.read_csv(orthofinder_dir+"phylodataB.tsv", sep=";")
print (df)
