# Draw BEAST tree

library(phytools)
library(phyloch)
library(strap)
library(coda)
#library(OutbreakTools)

setwd("~/Desktop")
tree.file = "Lentinula.tree"
root.time = 131.277  # From BEAST results TreeHeight median

###########################
#### PLOT ALL BRANCHES ####
###########################

# Draw background tree
t = read.beast(tree.file)
t$root.time = root.time
geoscalePhylo(
  tree=ladderize(t, right=FALSE), units=c("Period", "Epoch"),
  boxes="Epoch", cex.tip=.5, cex.age=.7, cex.ts=0.5, label.offset=1,
  lwd=3, width=2
)

lastPP = get("last_plot.phylo", envir=.PlotPhyloEnv)
names_list = t$tip
nids <- vector()
pos <- 1
len_nl <- length(names_list)
for(n in names_list){
  for(nn in names_list[pos:len_nl]){
    if(n != nn){
      m <- getMRCA(t, c(n, nn))
      if(m %in% nids == FALSE){
        nids <- c(nids, m)
      }
    }
  }
  pos = pos + 1
}

num_taxa <- length(t$tip.label)
for (nv in nids) {
  bar_xx_a <- c(
    lastPP$xx[nv] + (t$height[nv-num_taxa] - t$"height_95%_HPD_MIN"[nv-num_taxa]),
    lastPP$xx[nv] - (t$"height_95%_HPD_MAX"[nv-num_taxa] - t$height[nv-num_taxa])
  )
  lines(
    bar_xx_a, c(lastPP$yy[nv], lastPP$yy[nv]), col=rgb(0, 0, 1, alpha=0.4), lwd=8
  )
}

t$node.label = t$posterior
p = character(length(t$node.label))
p[t$node.label >= 0.95] <- "black"
p[t$node.label < 0.95 & t$node.label >= 0.75] <- "gray"
p[t$node.label < 0.75] <- "white"
nodelabels(pch=21, cex=.6, bg=p)
