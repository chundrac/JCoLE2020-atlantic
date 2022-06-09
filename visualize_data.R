require(ggtree)
require(phytools)
require(phangorn)
require(ggplot2)
require(tikzDevice)

atlantic.data <- read.csv('../data_raw/Atlantic.csv')
rownames(atlantic.data) <- atlantic.data$Glottocode
atlantic.data <- atlantic.data[,c('ART','DEM','ADJ','CARD','PRO','VERB','PREF')]

trees <- read.nexus('trees.nex')

tree <- maxCladeCred(trees)

atlantic.data[atlantic.data==1]='Marker present'
atlantic.data[atlantic.data==0]='Marker absent'

g <- ggtree(tree) + geom_tiplab(size=2.5,align=TRUE)

tikz('data_and_tree')
gheatmap(g,apply(atlantic.data,2,as.factor), offset=1500, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0) + 
  scale_fill_manual(breaks=c("Marker absent","Marker present"), 
                    values=c("#00BFC4","#F8766D"),name="marker")
dev.off()