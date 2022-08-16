
import sys, argparse, random
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle, TextFace

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_dir', default=False, help='OrthoFinder output directory. Unless specified otherwise, it will search the necessary input files within that directory')
	parser.add_argument('-p', '--pfam', required=True, help='Output from a HMMrsearch run with tableoutput. Required.')
	parser.add_argument('-s', '--search', nargs='+', default=False, help='Terms to search in the database. It will search for Pfam terms and Orthogroupsby ID')
	parser.add_argument('-f', '--search_file', default=False, help='You can provide a file for the search. In that case, it will assume that each line of the file contains a single search term')
	parser.add_argument('-e', '--exclude', nargs='+', default=False, help='These terms will be excluded from the search')
	parser.add_argument('-x', '--exclude_file', default=False, help='You can provide a file for excluding it. In that case, it will assume that each line of the file contains a single search term')
	parser.add_argument('-O', '--orthologues',default=False, help='File relating the identificator of orthogroups to the gene sequences included in that orthogroup. By default it will search for them within the input directory as Orthogroups/Orthogroups_total.tsv.')
	parser.add_argument('-D','--duplications',default=False, help='File indicating duplication events. By default it will search for them within the input directory.')
	parser.add_argument('-R','--ref_tree',default=False, help='File indicating the species tree created by OrthoFinder. By default it will search for them within the input directory.')
	parser.add_argument('--svg', default=False,  help='Render the tree into an svg file.')

args = parser.parse_args()

if args.input_dir == False: inputdir = './'
else: inputdir = args.input_dir
if inputdir[-1] == "/": inputdir = inputdir[:-1]

pfamfile = open(args.pfam)

if args.orthologues == False:
	orthofile = open(inputdir+"/Orthogroups/Orthogroups.tsv")
else: orthofile = open(args.orthologues)

if args.duplications == False:
	dupfile = open(inputdir+"/Gene_Duplication_Events/Duplications.tsv")
else: dupfile = open(args.duplications)

if args.ref_tree == False:
	reftree = Tree(inputdir+"/Species_Tree/SpeciesTree_rooted_node_labels.txt", format=1)
else: reftree = Tree(args.ref_tree, format=1)

searchlist, excludelist = [], []

taxaset = set()

if args.exclude != False:
	for i in args.exclude: excludelist.append(i)

if args.exclude_file != False:
	for line in open(args.exclude_file):
		excludelist.append(line[:-1])

if args.search != False:
	for i in args.search: searchlist.append(i)

if args.search_file != False:
	for line in open(args.search_file):
		searchlist.append(line[:-1])

def pfam_dicts(pfamfile):
	pfam_id_dict, pfam_name_dict, pfam_descr_dict, gene_pfam_dict, pfam_gene_dict = {}, {}, {}, {}, {}
	pfam_id_set = set()
	for line in pfamfile:
		if line[0] == "#": continue
		else:
			chunk = line.split()
			descriptor = ''
			for word in chunk[18:]:
				descriptor = descriptor + word + " "
			pfam_id_set.add((chunk[1], chunk[0], descriptor[:-1]))
			if chunk[2] in gene_pfam_dict:
				gene_pfam_dict[chunk[2]].append((chunk[1], chunk[4]))
			else:
				if len(chunk) < 5: continue
				gene_pfam_dict[chunk[2]] = [(chunk[1], chunk[4])]
			if chunk[1]+"/"+chunk[0]+"/"+descriptor in pfam_gene_dict:	
				pfam_gene_dict[chunk[1]+"/"+chunk[0]+"/"+descriptor].append(chunk[2])
			else:	
				pfam_gene_dict[chunk[1]+"/"+chunk[0]+"/"+descriptor] = [chunk[2]]
	pfamfile.seek(0)
	for i in pfam_id_set:
		pfam_id_dict[i[0]] = [i[1], i[2]]
		pfam_name_dict[i[1]] = [i[0], i[2]]
		pfam_descr_dict[i[2]] = [i[0], i[1]]
	return pfam_id_dict, pfam_name_dict, pfam_descr_dict, gene_pfam_dict, pfam_gene_dict

def ortho_dict(orthogroups):
	ortho_dict, gene_dict = {}, {}
	orthogroups.seek(0)
	spp_name_list = []
	gene2spp_dict = {}
	for line in orthogroups:
		orthogroup = line.split("\t")[0]
		chunk = line[:-1].split("\t")
		genelist = []
		spset = set()
		if orthogroup == "Orthogroup":
			for i in range(1,len(chunk)):
				spp_name_list.append(chunk[i])
				taxaset.add(chunk[i])
		else:
			for e in range(1, len(chunk)):
				for i in chunk[e].split(","):
					i = i.replace(" ", "")
					if i == "": continue
					else:
						genelist.append(i)
						sp = spp_name_list[e-1]
						spset.add(sp)
						gene2spp_dict[i] = sp
		ortho_dict[orthogroup] = [genelist, spset, taxaset-spset]
		for gene in genelist:
			gene_dict[gene] = orthogroup
	return ortho_dict, gene_dict, gene2spp_dict

def duplication_dict(dupfile):
	duplidict = {}
	for line in dupfile:
		chunk = line.split()
		if chunk[0] in duplidict:
			duplidict[chunk[0]] = duplidict[chunk[0]] + "\t" + chunk[1]
		else:
			duplidict[chunk[0]] = chunk[1]
	return duplidict

def retrieve_data_ortho(ortholist, pfamfile, orthofile, pfam_id_dict, pfam_name_dict, pfam_descr_dict, gene_pfam_dict, pfam_gene_dict, orthodict, genedict):
	pfam_ortho_dict = {}
	for ortho in ortholist:
		pfam_ortho_dict[ortho] = []
		for i in orthodict[ortho][0]:
			if i in gene_pfam_dict:
				pfam_ortho_dict[ortho].append(gene_pfam_dict[i])
				for pfam in gene_pfam_dict[i]: print("@", i, "\t", pfam[0], "\t", pfam_id_dict[pfam[0]][0], "\t", pfam_id_dict[pfam[0]][1])
			else:
				print ("@", i, "no pfam domains identified")

def get_dict_for_tree(ortholist, orthodict, dictfortree, gene2spp_dict):
	for ortho in ortholist:
		dictfortree[ortho] = {}
		for taxa in taxaset:
			dictfortree[ortho][taxa] = 0
		for i in orthodict[ortho][0]:
			dictfortree[ortho][gene2spp_dict[i]] = dictfortree[ortho][gene2spp_dict[i]] + 1
		
			
def search_pfam(pfamlist, pfamfile, orthofile, dupfile, reftree, excludelist):
	pfam_id_dict, pfam_name_dict, pfam_descr_dict, gene_pfam_dict, pfam_gene_dict = pfam_dicts(pfamfile)
	orthodict, genedict, gene2spp_dict = ortho_dict(orthofile)
	duplidict = duplication_dict(dupfile)
	nodeduplidict, dictfortree, orthoentrydict = {}, {}, {}
	excludeorthos, includeorthos = [], []
	for ex in excludelist:
		if ex in orthodict: excludeorthos.append(ex)
	for inc in pfamlist:
		if inc in orthodict: includeorthos.append(inc)
	for pfam in pfamlist:
		for entry in pfam_gene_dict:
			toprintdict = {}
			if entry.find(pfam) > -1:
				for ex in excludelist:
					if entry.find(ex) > -1: break
				else:
					print("###\t"+entry+"\t###")
					for gene in pfam_gene_dict[entry]:
						if gene not in genedict: continue
						if genedict[gene] not in excludeorthos and gene in genedict:	
							ortho = genedict[gene]
							orthoentrydict[ortho] = entry
							toprintdict[ortho] = (orthodict[ortho][1], orthodict[ortho][2])
						else:
							toprintdict[gene] = "???", "???"
			else:
				if pfam in includeorthos and pfam not in excludeorthos and pfam not in orthoentrydict:
					print("###\t"+pfam+"\t###")
					if pfam in orthodict:
						orthoentrydict[pfam] = " + / + / + "
						toprintdict[pfam] = (orthodict[pfam][1], orthodict[pfam][2])
			for i in toprintdict:
				if i in genedict:
					if genedict[i] in excludeorthos: continue
				absentstring, presentstring = "", ""
				for string in toprintdict[i][1]:
					absentstring = absentstring + " " + string
				for string in toprintdict[i][0]:
					presentstring = presentstring + " " + string
				print ("%", i, "\tpresent in", len(toprintdict[i][0]), ":", presentstring, "\tabsent in", len(toprintdict[i][1]), absentstring)
				retrieve_data_ortho([i], pfamfile, orthofile, pfam_id_dict, pfam_name_dict, pfam_descr_dict, gene_pfam_dict, pfam_gene_dict, orthodict, genedict)
				get_dict_for_tree([i], orthodict, dictfortree, gene2spp_dict)
				if i in duplidict: 
					print ("{", i, entry, ":", duplidict[i])
					for e in duplidict[i].split():
						if e in nodeduplidict:
							nodeduplidict[e] = nodeduplidict[e] + 1
						else:
							nodeduplidict[e] = 1
		print ("\n")
	plot_tree(reftree, nodeduplidict, dictfortree, orthoentrydict, includeorthos)
				

def plot_tree(mytree, nodeduplidict, dictfortree, orthoentrydict, includeorthos):
	maxperspdict = {}
	currentry = ('', '')
	orderlist = []
	for i in dictfortree:
		orderlist.append(orthoentrydict[i]+"@"+i)
	orderlist.sort()
	for i in orderlist:
		for e in dictfortree[i.split("@")[1]]:
			if e in maxperspdict:
				maxperspdict[e] = maxperspdict[e] + dictfortree[i.split("@")[1]][e]
			else:
				maxperspdict[e] = dictfortree[i.split("@")[1]][e]
	nodestyledict = {}
	ts = TreeStyle()
	standardstyle = NodeStyle()
	standardstyle["vt_line_width"] = 1
	standardstyle["hz_line_width"] = 1
	standardstyle["vt_line_color"] = "#aaaaaa"
	standardstyle["hz_line_color"] = "#aaaaaa"
	standardstyle["fgcolor"] = "SteelBlue"
	standardstyle["size"] = 3
	
	fakestyle = NodeStyle()
	fakestyle["vt_line_width"] = 0
	fakestyle["vt_line_color"] = "#ffffff"
	fakestyle["hz_line_color"] = "#ffffff"
	fakestyle["size"] = 0
	for i in nodeduplidict:
		nodestyledict[i] = NodeStyle()
		nodestyledict[i]["vt_line_width"] = 1
		nodestyledict[i]["hz_line_width"] = nodeduplidict[i]
		nodestyledict[i]["vt_line_color"] = "#aaaaaa"
		nodestyledict[i]["hz_line_color"] = "#000000"
		nodestyledict[i]["fgcolor"] = "FireBrick"
		nodestyledict[i]["size"] = 3
		
	for node in mytree.traverse():
		if node.name in nodestyledict:
			node.img_style = nodestyledict[node.name]
		else:
			node.img_style = standardstyle
	for leaf in mytree.iter_leaves():
		count = 0
		for ortho in orderlist:
			e = ortho.split("@")[1]
			currentface = TextFace(str(dictfortree[e][leaf.name]))
			currentface.border.width = 1
			currentface.margin_top = 10
			currentface.margin_right = 10
			currentface.margin_left = 10
			currentface.margin_bottom = 10
			currentface.margin_bottom = 10
			currentface.fsize = 10
			if dictfortree[e][leaf.name] == 0:
				currentface.fgcolor = "Black"
				currentface.background.color = "LightGray"
				if ortho.split("@")[-1] in includeorthos:
					currentface.border.width = 3
			else:
				if ortho.split("@")[-1] in includeorthos:
					r, g, b = int(255-255*(dictfortree[e][leaf.name]/maxperspdict[leaf.name])), int(255-255*(dictfortree[e][leaf.name]/maxperspdict[leaf.name])), 255
					currentcolor = "#"+hex(r).lstrip("0x").rstrip("L")+hex(g).lstrip("0x").rstrip("L")+hex(b).lstrip("0x").rstrip("L")
					if  g < 10 or b < 10:
						currentcolor = "#"+hex(r).lstrip("0x").rstrip("L")+"0"+hex(g).lstrip("0x").rstrip("L")+"0"+hex(b).lstrip("0x").rstrip("L")
					if g == 0 or b == 0:
						currentcolor = "#0000ff"
					currentface.border.width = 3
				else:
					r, g, b = 255, int(255-255*(dictfortree[e][leaf.name]/maxperspdict[leaf.name])), int(255-255*(dictfortree[e][leaf.name]/maxperspdict[leaf.name]))
					currentcolor = "#"+hex(r).lstrip("0x").rstrip("L")+hex(g).lstrip("0x").rstrip("L")+hex(b).lstrip("0x").rstrip("L")
					if  g < 10 or b < 10:
						currentcolor = "#"+hex(r).lstrip("0x").rstrip("L")+"0"+hex(g).lstrip("0x").rstrip("L")+"0"+hex(b).lstrip("0x").rstrip("L")
					if g == 0 or b == 0:
						currentcolor = "#ff0000"
				currentface.background.color = currentcolor
			leaf.add_face(currentface, column=count, position="aligned")
			count = count + 1
	faketree = Tree('(FAKENODE);')
	fakenode = faketree.search_nodes(name="FAKENODE")[0]
	fakeroot = faketree.get_tree_root()
	fakenode.add_sister(mytree)
	fakenode.img_style = fakestyle
	fakeroot.img_style = fakestyle
	fakenode.name = ''
	fakecount = 0
	for i in orderlist:
		ortho = i.split("@")[1]
		if orthoentrydict[ortho] != currentry[0]:
			entry = orthoentrydict[ortho]
			entryface = TextFace(entry.split("/")[1])
			entryface.fsize = 10
			entryface.margin_left = 5
			entryface.margin_right = 5
			entryface.margin_bottom = 10
			r, g, b = random.random(), random.random(), random.random()
			entryface.background.color = "#"+hex(127 + int(127*r)).lstrip("0x").rstrip("L")+hex(127 + int(127*g)).lstrip("0x").rstrip("L")+hex(127 + int(127*b)).lstrip("0x").rstrip("L")
			fakenode.add_face(entryface, column=fakecount, position="aligned")
			currentry = (entry, [r,g,b])
		else:
			emptyface = TextFace("")
			emptyface.background.color = entryface.background.color = "#"+hex(127 + int(127*r)).lstrip("0x").rstrip("L")+hex(127 + int(127*g)).lstrip("0x").rstrip("L")+hex(127 + int(127*b)).lstrip("0x").rstrip("L")
			emptyface.margin_bottom = 10
			fakenode.add_face(emptyface, column=fakecount, position="aligned")
		currentfakeface = TextFace(ortho)
		currentfakeface.border.width = 1
		currentfakeface.margin_top = 10
		currentfakeface.margin_right = 5
		currentfakeface.margin_left = 5
		currentfakeface.margin_bottom = 10
		currentfakeface.background.color = "#"+hex(90+int(127*r)).lstrip("0x").rstrip("L")+hex(90+int(127*g)).lstrip("0x").rstrip("L")+hex(90+int(127*b)).lstrip("0x").rstrip("L")
		currentfakeface.fgcolor = "Black"
		currentfakeface.fsize = 6
		currentfakeface.tight_text = True
		fakenode.add_face(currentfakeface, column=fakecount, position="aligned")

		fakecount = fakecount + 1
	searchstring, excludestring = '', ''
	for i in searchlist:
		searchstring = searchstring + i + " "
	for i in excludelist:
		excludestring = excludestring + i + " "
	searchface = TextFace("Search: " + searchstring)	
	searchface.margin_right = 20
	fakenode.add_face(searchface, column=0, position="branch-top")
	if len(excludelist) > 0:
		excludeface = TextFace("Exclude: " + excludestring)
		excludeface.margin_right = 20
		fakenode.add_face(excludeface, column=0, position="branch-top")
	if args.svg == False:
		faketree.show(tree_style=ts)
	else:
		faketree.render(args.svg+".svg", h=700, tree_style=ts)

search_pfam(searchlist, pfamfile, orthofile, dupfile, reftree, excludelist)
