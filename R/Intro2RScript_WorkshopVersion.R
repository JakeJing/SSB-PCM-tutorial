## # Introduction to R in phylogenetic comparative methods
## ## Prepare the environment
install.packages("geiger")
install.packages("devtools")
library(devtools)
library(geiger)

install_github("richfitz/sowsear")
install_github("arborworkflows/aRbor")
## Note if you get an error in this like: "(error 403) Forbidden"  you likely need to update R and reinstall the packages. This is
## an error because the older versions of devtools that are installed on older versions of R have the wrong url (github changed
## the base URL for downloads)

## ## Directory management
getwd()
setwd("~/repos/SSB/Intro/R")

## Read the tree
tree <- read.tree("../data/schenk.tre")
tree
plot(tree)

class(tree)
str(tree)

tree['tip.label']
length(tree['tip.label'])
tree[['tip.label']]
length(tree[['tip.label']])

tree$tip.label

tree[[3]]

tree$edge
tree$edge.length

tree
drop.tip(tree, c("Eliomys_quercinus", "Sciurus_sp."))

dist <- cophenetic(tree)
dist[1:5, 1:5]
## ## Matching a tree with a dataset
pantheria <- read.table("../data/rodents.csv", sep=",")
pantheria <- read.csv("../data/rodents.csv")

head(pantheria)
rownames(pantheria) <- gsub(" ", "_", pantheria$MSW05_Binomial)
rownames(pantheria)

td <- treedata(tree, pantheria)
td$data <- as.data.frame(td$data)

require(aRbor)
td <- make.treedata(tree, pantheria)
str(td)
length(td)
names(td)
td$dat==-999
td$dat[td$dat==-999] <- NA

td <- mutate(td, lnMass=log(X5.1_AdultBodyMass_g), lnBMR=log(X5.2_BasalMetRateMass_g), lifespan=X17.1_MaxLongevity_m, desert=ifelse(X28.1_Precip_Mean_mm <21, 1, 0))
td <- select(td, lnMass, lnBMR, lifespan, desert)

td <- filter(td, !is.na(lnMass), !is.na(lnBMR), !is.na(lifespan), !is.na(desert))

tree <- td$phy
dat <- td$dat

## ## Visualization
plot(tree)
plot(tree, cex=0.2)
plot
?plot
?plot.phylo

## Plotting a discrete character
plot(tree, show.tip.label=FALSE)
## Typing ?par is very helpful
?par
dat$desert
tiplabels(pch=21, bg=dat$desert+1)
nodelabels(cex=0.75, bg="white")
nodelabels(pch=21, bg="purple")

require(phytools)
phenogram(tree, dat$lnMass)
phenogram(tree, setNames(dat$lnMass, tree$tip.label), spread.labels=FALSE, ftype="off")

fancyDat <- as.matrix(dat[,1:3])
rownames(fancyDat) <- tree$tip.label
fancyTree(tree, type="scattergram", X=fancyDat, fsize=0.3)

## ## V. Simple analyses of traits
## ### Discrete trait models
trait <- setNames(dat$desert, tree$tip.label)+1
mER <- fitDiscrete(tree, trait, model="ER")
mARD <- fitDiscrete(tree, trait, model="ARD")
mARD
str(mARD)
aiccs <- c(mER$opt$aicc, mARD$opt$aicc)

## ## Continuous trait models
## Brownian motion
trait <- setNames(dat$lnMass, tree$tip.label)
mBM <- fitContinuous(tree, trait, model="BM")
mBM
## Ornstein Uhlenbeck
mOU <- fitContinuous(tree, trait, model="OU")
mOU

## ## Phylogenetic signal
log(2)/mOU$opt$alpha

attributes(dat)

rownames(dat) == tree$tip.label

treeply(td, drop.tip, "Jaculus_jaculus")





