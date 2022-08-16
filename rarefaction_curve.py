import sys
import matplotlib.pyplot as plt
import random
from statistics import median

orthofinder_dir = sys.argv[1]
if orthofinder_dir[-1] == "/":
	orthofinder_dir = orthofinder_dir[:-1]
orthofile = orthofinder_dir + "/Orthogroups/Orthogroups.GeneCount.tsv"

sppfile, spplist = sys.argv[2], []
figname = sys.argv[-1]
spppos = {}

for i in open(sppfile):
	spplist.append(i[:-1])

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
	random.shuffle(spplist)
	for spp in spplist:
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
	med_pan.append(median(i))
for i in rand_core:
	med_core.append(median(i))

fig = plt.figure()
p1 = plt.violinplot(rand_pan[1:], positions=range(2,len(spplist)+1))
p2 = plt.violinplot(rand_core)
p3 = plt.plot(range(1,len(spplist)+1), med_pan, "b-")
p4 = plt.plot(range(1,len(spplist)+1), med_core, "r-")
plt.axhline(y=med_pan[-1], color='#CCCCCC', linestyle='-.')
plt.axhline(y=med_core[-1], color='#CCCCCC', linestyle='-.')
plt.xrange=[1,len(spplist)]
plt.text(len(spplist), med_pan[-1]+20, str(int(med_pan[-1])), fontsize=12)
plt.text(len(spplist), med_core[-1]+20, str(int(med_core[-1])), fontsize=12)
plt.title(figname)

#p2 = plt.bar(x=range(1,len(names)+1), height=core, color="tab:blue")
#plt.xticks(rotation = 90)
plt.ylabel('Cumulative number of orthogroups', fontsize=12)

fig.set_size_inches(18.5, 10.5)
plt.savefig("rarefactioncurve_"+figname+".svg")
plt.savefig("rarefactioncurve_"+figname+".png")
plt.show()

