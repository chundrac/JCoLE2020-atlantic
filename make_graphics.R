library(reshape2)
library(rstan)
library(ggplot2)
library(tikzDevice)

load('all_rates_hierarchical.Rdata')

ordered <- order(colMeans(extract(fit)$pi))
ordered.names <- colnames(atlantic.data)[ordered]

diff.matrix <- matrix(nrow=length(ordered),ncol=length(ordered))
colnames(diff.matrix) <- ordered.names
rownames(diff.matrix) <- ordered.names

for (i in 1:length(ordered)) { 
  for (j in 1:length(ordered)) { 
    if (ordered[i] != ordered[j]) {
      print(c(i,j));
      diff.matrix[i,j]<-length(which(extract(fit.full)$contrasts_pi[,ordered[i],ordered[j]]>0))/length(extract(fit.full)$contrasts_pi[,ordered[i],ordered[j]])
    }
  }
}

diff.melt <- melt(diff.matrix, na.rm = TRUE)

diff.melt$clf.pref <- rep(NA,nrow(diff.melt))

diff.melt$clf.pref[diff.melt$value >= .75] = "greater in $> 75$\\%"
diff.melt$clf.pref[diff.melt$value >= .85] = "greater in $> 85$\\%"
diff.melt$clf.pref[diff.melt$value >= .95] = "greater in $> 95$\\%"

diff.melt <- na.omit(diff.melt)

tikz(file='diff-pref.tex',width=4,height=4)
ggplot(data = diff.melt, aes(Var2, Var1)) + geom_tile(aes(fill = clf.pref)) + 
  geom_text(aes(label = round(value, 3))) + 
  scale_y_discrete(limits=rev(levels(diff.melt$Var1))) + 
  #scale_fill_manual(values = c("#F8766D1A","#F8766D4C","#F8766DE6","#00BFC41A","#00BFC44C","#00BFC4E6"))
  scale_fill_manual(name='Marker preference',values = c("#F8766D1A","#F8766D99","#F8766DE6","#00BFC41A","#00BFC499","#00BFC4E6")) + 
  theme_classic() + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position=c(.82,.9))

dev.off()

###
###

ordered <- order(colMeans(extract(fit)$r))
ordered.names <- colnames(atlantic.data)[ordered]

diff.matrix <- matrix(nrow=length(ordered),ncol=length(ordered))
colnames(diff.matrix) <- ordered.names
rownames(diff.matrix) <- ordered.names

for (i in 1:length(ordered)) { 
  for (j in 1:length(ordered)) { 
    if (ordered[i] != ordered[j]) {
      print(c(i,j));
      diff.matrix[i,j]<-length(which(extract(fit.full)$contrasts_r[,ordered[i],ordered[j]]>0))/length(extract(fit.full)$contrasts_r[,ordered[i],ordered[j]])
    }
  }
}

diff.melt <- melt(diff.matrix, na.rm = TRUE)

diff.melt$speed <- rep(NA,nrow(diff.melt))

diff.melt$speed[diff.melt$value >= .75] = "greater in $> 75$\\%"
diff.melt$speed[diff.melt$value >= .85] = "greater in $> 85$\\%"
diff.melt$speed[diff.melt$value >= .95] = "greater in $> 95$\\%"

diff.melt <- na.omit(diff.melt)

tikz(file='diff-speed.tex',width=4,height=4)
ggplot(data = diff.melt, aes(Var2, Var1)) + geom_tile(aes(fill = speed)) + 
  geom_text(aes(label = round(value, 3))) + 
  scale_y_discrete(limits=rev(levels(diff.melt$Var1))) + 
  #scale_fill_manual(values = c("#F8766D1A","#F8766D4C","#F8766DE6","#00BFC41A","#00BFC44C","#00BFC4E6"))
  scale_fill_manual(name='Speed',values = c("#F8766D1A","#F8766D99","#F8766DE6","#00BFC41A","#00BFC499","#00BFC4E6")) + 
  theme_classic() + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position=c(.8,.9))

dev.off()

load('all_deltas.Rdata')

ordered <- order(colMeans(extract(fit)$b/extract(fit)$a))
ordered.names <- colnames(atlantic.data)[ordered]

diff.matrix <- matrix(nrow=length(ordered),ncol=length(ordered))
colnames(diff.matrix) <- ordered.names
rownames(diff.matrix) <- ordered.names

for (i in 1:length(ordered)) { 
  for (j in 1:length(ordered)) { 
    if (ordered[i] != ordered[j]) {
      print(c(i,j));
      diff.matrix[i,j]<-length(which(extract(fit.full)$delta_diff[,ordered[i],ordered[j]]>0))/length(extract(fit.full)$delta_diff[,ordered[i],ordered[j]])
    }
  }
}

diff.melt <- melt(diff.matrix, na.rm = TRUE)

diff.melt$stability <- rep(NA,nrow(diff.melt))

diff.melt$stability[diff.melt$value >= .75] = "greater in $> 75$\\%"
diff.melt$stability[diff.melt$value >= .85] = "greater in $> 85$\\%"
diff.melt$stability[diff.melt$value >= .95] = "greater in $> 95$\\%"

diff.melt <- na.omit(diff.melt)

tikz(file='diff-stab.tex',width=4,height=4)
ggplot(data = diff.melt, aes(Var2, Var1)) + geom_tile(aes(fill = stability)) + 
  geom_text(aes(label = round(value, 3))) + 
  scale_y_discrete(limits=rev(levels(diff.melt$Var1))) + 
  #scale_fill_manual(values = c("#F8766D1A","#F8766D4C","#F8766DE6","#00BFC41A","#00BFC44C","#00BFC4E6"))
  scale_fill_manual(name='Stability',values = c("#F8766D1A","#F8766D99","#F8766DE6","#00BFC41A","#00BFC499","#00BFC4E6")) + 
  theme_classic() + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position=c(.82,.9))

dev.off()
