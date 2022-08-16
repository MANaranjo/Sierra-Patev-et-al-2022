# Patev-et-al-2022-scripts
Repository for: Sean Sierra-Patev, Byoungnam Min, Miguel Naranjo-Ortiz, Brian Looney, Zachary Konkel, Jason C. Slot, Yuichi Sakamoto, Jacob Steenwyk, Antonis Rokas, Juan Carro, Susana Camarero, Patricia Ferreira, Gonzalo Molpeceres, Francisco J. Ruiz-Dueñas, Ana Serrano, Bernard Henrissat, Karen W. Hughes, Juan L. Mata, Noemia Kazue Ishikawa, Ruby Vargas-Isla, Shuji Ushijima, Chris A. Smith, Steven Ahrendt, Willian Andreopoulos, Guifen He, Kurt LaButti, Anna Lipzen, Vivian Ng, Robert Riley, Laura Sandor, Kerrie Barry, Angel Martinez, Yang Xiao, John G. Gibbons, Kazuhisa Terashima, Igor V. Grigoriev, and David Hibbett. A Global Phylogenomic Analysis of the Shiitake Genus Lentinula.

Orthoheatmap.py
This script generates a heatmap representing the percent of genes having a 1 to 1, 1  to many and many to many relationships between pairs of species. It only requires as input an OrthoFinder run directory.


Map_dupl2tree.py
The script maps duplication and some duplication-related stats to the nodes of the phylogenetic tree generated by OrthoFinder. It requires only an OrthoFinder run directory as input. The numbers represent, for each node:
Red: Duplication rate (Duplications mapped to that node / number of orthogroups represented in that node) / Relative duplication rate (Duplication rate / branch length).
Purple: Number of duplications / Number of orthogroups mapped to that node.
Blue: Variation coefficient (standard deviation / mean) / Skewness
Green: Branch length / branch support


Pfam_interrogate.py
This script is basically a browser. As input It requires an OrthoFinder run directory and a HMMerscan tabular output file with all proteomes. The HMMer file might contain more sequences that the ones used for OrthoFinder, and I think it won’t crash if used with less, it’s just that it won’t be able to find anything.

The script requires a list of search terms that can be either orthogroup numbers (for that particular OrthoFinder run), Pfam identifiers or any word. Words will be used to search matches in the description and name of the Pfam families, and it’s case sensitive. It is possible to introduce terms to exclude, as well.

The script returns a representation of the phylogenetic tree generated by OrthoFinder where branch thickness in pixels equals the number of duplications mapped to all the orthogroups found by the search criteria. Next to the tree it produces a table with all orthogroups found where at least one member has a Pfam domain that matches the search (or, if orthogroup numbers are used as query, that orthogroup in particular), and the number of members that orthgroup has per species.
