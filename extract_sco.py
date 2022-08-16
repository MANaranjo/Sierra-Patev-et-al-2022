import sys, os
from Bio import AlignIO

orthofinder_dir = sys.argv[1]
if orthofinder_dir[-1] == "/":
	orthofinder_dir = orthofinder_dir[:-1]
orthofile = orthofinder_dir + "/Orthogroups/Orthogroups.GeneCount.tsv"

ortholist = []

for line in open(orthofile):
	chunk = line.split()
	if chunk[0] == "Orthogroup":
		continue
	else:
		subchunk = []
		for i in chunk[int(sys.argv[2]):int(sys.argv[3])]:
			subchunk.append(int(i))
		if sum(subchunk) == len(subchunk) and 0 not in subchunk:
			ortholist.append(chunk[0])

for i in ortholist:
	aln = AlignIO.read(open(orthofinder_dir+"/MultipleSequenceAlignments/"+i+".fa"), "fasta")
	aln_out = open(i+".fa", "w")
	for seq in aln:
		if seq.id.find("Ledo") > -1:
			aln_out.write(seq.format("fasta"))
			aln_out.write("\n")


