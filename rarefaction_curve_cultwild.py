import sys
import matplotlib.pyplot as plt
import random
from statistics import median, mean

orthofinder_dir = sys.argv[1]
if orthofinder_dir[-1] == "/":
	orthofinder_dir = orthofinder_dir[:-1]
orthofile = orthofinder_dir + "/Orthogroups/Orthogroups.GeneCount.tsv"

sppfile1, sppfile2, spplist = sys.argv[2], sys.argv[3], []
figname = sys.argv[-1]
spppos = {}
firstlist, secondlist = [], []

for i in open(sppfile1):
	spplist.append(i[:-1])
	firstlist.append(i[:-1])
for i in open(sppfile2):
	spplist.append(i[:-1])
	secondlist.append(i[:-1])

setdict = {}
for line in open(orthofile):
	chunk = line[:-1].split("\t")
	if chunk[0] == "Orthogroup":
		for i in range(len(chunk)):
			if chunk[i] in spplist:
				spppos[i] = chunk[i] 
				setdict[chunk[i]] = set()
	else:
		for i in range(len(chunk)):
			if i in spppos and int(chunk[i]) > 0:
				setdict[spppos[i]].add(chunk[0])
				
rand_pan, rand_core = [], []
	
for i in spplist:
	rand_pan.append([])
	rand_core.append([])

while len(rand_pan[0]) < 1000:
	names, pan, core = [], [] , []
	Aset, Bset = set(), set()
	random.shuffle(firstlist)
	random.shuffle(secondlist)
	newlist = firstlist+secondlist
	for spp in newlist:
		names.append(spp)
		if len(pan) == 0:
			pan.append(len(setdict[spp]))
			core.append(len(setdict[spp]))
			Aset, Bset = setdict[spp], setdict[spp]
		else:
			Aset, Bset = Aset | setdict[spp], Bset & setdict[spp]
			pan.append(len(Aset))
			core.append(len(Bset))
	for i in range(0, len(pan)):
		rand_pan[i].append(pan[i])
		rand_core[i].append(core[i])

med_pan, med_core = [], []
for i in rand_pan:
	med_pan.append(mean(i))
for i in rand_core:
	med_core.append(mean(i))

fig = plt.figure()
plt.axhline(y=med_pan[-1], color='#CCCCCC', linestyle='-.')
plt.axhline(y=med_core[-1], color='#CCCCCC', linestyle='-.')
plt.axvline(x=len(firstlist), color='#CCCCCC', linestyle='-.')
p1 = plt.violinplot(rand_pan[1:], positions=range(2,len(spplist)+1))
p2 = plt.violinplot(rand_core)
p3 = plt.plot(range(1,len(spplist)+1), med_pan, "b-")
p4 = plt.plot(range(1,len(spplist)+1), med_core, "r-")

plt.text(len(firstlist)-5, med_pan[-1]+20, str("Cultivars"), fontsize=12)
plt.text(len(firstlist)+2, med_pan[-1]+20, str("Wild"), fontsize=12)

plt.xrange=[1,len(spplist)]
plt.text(len(spplist), med_pan[-1]+20, str(int(med_pan[-1])), fontsize=12)
plt.text(len(spplist), med_core[-1]+20, str(int(med_core[-1])), fontsize=12)
plt.text(len(firstlist)-3, med_pan[len(firstlist)-1]+100, str(int(med_pan[len(firstlist)-1])), fontsize=12)
plt.text(len(firstlist), med_core[len(firstlist)-1]+100, str(int(med_core[len(firstlist)-1])), fontsize=12)
plt.title(figname)

plt.ylabel('Cumulative number of orthogroups', fontsize=12)

fig.set_size_inches(18.5, 10.5)
plt.savefig("rarefactioncurve_"+figname+".svg")
plt.savefig("rarefactioncurve_"+figname+".png")
plt.show()

