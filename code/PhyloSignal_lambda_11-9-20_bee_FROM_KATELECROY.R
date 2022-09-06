library(picante)
library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)
library(epicalc)

zap()

mytree=read.tree(file.choose())
 mytree
 traits.dat.Deltasf
traits.dat <- read.csv(file.choose(), header=TRUE)
traits.dat.Deltas
Deltas
treeData <- data.frame(traits.dat$Deltas)
treeData
row.names(treeData) <- traits.dat[,1] 
treeData2 <-treeData
attach(treeData)
treeData
#testing if attached:
traits.dat$Deltas
adding traits w/ rownames:
names(b1) <- rownames(treeData)
names(b2) <- rownames(treeData)
names(b3) <- rownames(treeData)
names(s1r) <- rownames(treeData)
names(s1g) <- rownames(treeData)
names(s1b) <- rownames(treeData)
names(s1u) <- rownames(treeData)
names(s1y) <- rownames(treeData)
names(s1v) <- rownames(treeData)
names(s2) <- rownames(treeData)
names(s3) <- rownames(treeData)
names(s5a) <- rownames(treeData)
names(s5b) <- rownames(treeData)
names(s5c) <- rownames(treeData)
names(s6) <- rownames(treeData)
names(s7) <- rownames(treeData)
names(s8) <- rownames(treeData)
names(s9) <- rownames(treeData)
names(h1) <- rownames(treeData)
names(h3) <- rownames(treeData)
names(h4a) <- rownames(treeData)
names(h4b) <- rownames(treeData)
names(h4c) <- rownames(treeData)
names(traits.dat.Deltas) <- rownames(treeData)
names(DeltaQ) <- rownames(treeData)
names(PCA1) <- rownames(treeData)
names(PCA2) <- rownames(treeData)
names(PCA3) <- rownames(treeData)
names(IP) <- rownames(treeData)
names(STD) <- rownames(treeData)
names(pca1l) <- rownames(treeData)
names(pca1ll) <- rownames(treeData)
names(traits.data.Deltasf)
treeData
b1label <- character(length(mytree$tip.label))
names(b1label) <- names(b1)
b2label <- character(length(mytree$tip.label))
names(b2label) <- names(b2)
b3label <- character(length(mytree$tip.label))
names(b3label) <- names(b3)
s1rlabel <- character(length(mytree$tip.label))
names(s1rlabel) <- names(s1r)
s1glabel <- character(length(mytree$tip.label))
names(s1glabel) <- names(s1g)
s1blabel <- character(length(mytree$tip.label))
names(s1blabel) <- names(s1b)
s1ulabel <- character(length(mytree$tip.label))
names(s1ulabel) <- names(s1u)
s1ylabel <- character(length(mytree$tip.label))
names(s1ylabel) <- names(s1y)
s1vlabel <- character(length(mytree$tip.label))
names(s1vlabel) <- names(s1v)
s2label <- character(length(mytree$tip.label))
names(s2label) <- names(s2)
s3label <- character(length(mytree$tip.label))
names(s3label) <- names(s3)
s5alabel <- character(length(mytree$tip.label))
names(s5alabel) <- names(s5a)
s5blabel <- character(length(mytree$tip.label))
names(s5blabel) <- names(s5b)
s5clabel <- character(length(mytree$tip.label))
names(s5clabel) <- names(s5c)
s6label <- character(length(mytree$tip.label))
names(s6label) <- names(s6)
s7label <- character(length(mytree$tip.label))
names(s7label) <- names(s7)
s8label <- character(length(mytree$tip.label))
names(s8label) <- names(s8)
s9label <- character(length(mytree$tip.label))
names(s9label) <- names(s9)
h1label <- character(length(mytree$tip.label))
names(h1label) <- names(h1)
h3label <- character(length(mytree$tip.label))
names(h3label) <- names(h3)
h4alabel <- character(length(mytree$tip.label))
names(h4alabel) <- names(h4a)
h4blabel <- character(length(mytree$tip.label))
names(h4blabel) <- names(h4b)
h4clabel <- character(length(mytree$tip.label))
names(h4clabel) <- names(h4c)
pca1label <- character(length(mytree$tip.label))
names(pca1label) <- names(PCA1)
pca2label <- character(length(mytree$tip.label))
names(pca2label) <- names(PCA2)
pca3label <- character(length(mytree$tip.label))
names(pca3label) <- names(PCA3)
IPlabel <- character(length(mytree$tip.label))
names(IPlabel) <- names(IP)
Deltaslabel <- character(length(mytree$tip.label))
names(Deltaslabel) <- names(traits.dat.Deltas)
Deltasflabel
Deltaqlabel <- character(length(mytree$tip.label))
names(Deltaqlabel) <- names(DeltaQ)
STDlabel <- character(length(mytree$tip.label))
names(STDlabel) <- names(STD)
pca1llabel <- character(length(mytree$tip.label))
names(pca1llabel) <- names(pca1l)
pca1lllabel <- character(length(mytree$tip.label))
names(pca1lllabel) <- names(pca1ll)



#

mytree
treeDataf
phylotraits <- phylo4d(mytree, treeData)
phylotraits

moran.test <- abouheif.moran(phylotraits,method="Abouheif")
moran.test
plot(moran.test)
traits.dat.Deltasf

phylosigDS <- phylosig(mytree, traits.data.Deltasf[mytree$tip.label], method="lambda", test=TRUE)
phylosigDS

phylosigDS <- phylosig(mytree, traits.data.Deltasf[mytree$tip.label], method="lambda", test=TRUE)
phylosigDS

phylosigDS <- phylosig(mytree, traits.data.Deltasf[mytree$tip.label], method="K", test=TRUE)
phylosigDS

phylosigDS <- phylosig(mytree, traits.dat.Deltas[mytree$tip.label], method="K", test=TRUE)
phylosigDS



Optional:
plot(mytree)
plot(mytree, cex=0.5, show.node.label=false)
mytree$node.label
