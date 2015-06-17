## # An introduction to phylogenetic comparative methods in R
## This tutorial is intended to get you familiar with the R environment for conducting
## phylogenetic comparative analyses.

## ## I. Preparing the R environment
##+echo=FALSE
options(max.print=50)

## ### Installing packages
## To install a package directly from CRAN, use:
##+eval=FALSE
install.packages("geiger")

## Load the package using 'library'
library(geiger)

## These days, many packages may not be on CRAN, but instead hosted on github, bitbucket etc.
## To install these packages, it's useful to use the `devtools` package:
##+eval=FALSE
install_github("richfitz/sowsear")
install_github("arborworkflows/aRbor")
##+eval=TRUE
library(sowsear)
library(devtools)
library(aRbor)

## For phylogenetic comparative methods, there are a number of packages available. One place
## to find them all is the 'Phylogenetics Task View' run by Brian O'meara: (http://cran.r-project.org/web/views/Phylogenetics.html)
## All of these packages can be installed at once using the following lines (however, please note that
## this process can take a considerable amount of time)

##+eval=FALSE
install.packages("ctv")
install.views("Phylogenetics")

## ### Directory management
## R looks for files in the 'working directory'. 
getwd()
##+eval=FALSE
dir.create("~/repos/SSBphylogenetics")
dir.create("~/repos/SSBphylogenetics/R")
dir.create("~/repos/SSBphylogenetics/data")
dir.create("~/repos/SSBphylogenetics/output")
##+eval=TRUE
setwd("~/repos/SSBphylogenetics/R/")

## Save this script into the folder './R' and place your data files into the directory './data'
## Now we are ready to read in data and tree files.

## ### II. Reading in a tree & tree data structures
## There are a number of functions to read in phylogenetic tree data into R. We are going
## to use as our example tree phylogeny of Muroid rodents (Schenk, Rowe & Steppan, 2013; Systematic Biology).
## The tree can be downloaded from my github site, or it be accessed directly from treebase (submission 12303).
tree <- read.tree("../data/schenk.tre")
tree
## We can plot the tree:
plot(tree)
plot(tree, type="fan", cex=0.25)
## You may notice that the help file for `plot` is not very helpful for options for phylogenies: 
?plot
## This is because `plot` is a "method" that behaves differently on objects of different classes. 
plot
methods(plot)
class(tree)
?plot.phylo

## How is a tree stored? How can it be manipulated and modified? 
## A tree object in R is a data structure called a "list", and given a species class, called a "phylo" object.
## It will be worth your time getting familiar with [data structures in R](http://adv-r.had.co.nz/Data-structures.html).

## `str` is a useful commands for determining what is in an R object:
str(tree)
## Each element of the list can be accessed by multiple means:
tree['tip.label']
tree[['tip.label']]
tree[[3]]
tree$tip.label
## The structure of the tree is stored in the edge matrix:
tree$edge
## And branch lengths are stored in the list element `edge.length`
tree$edge.length
## Tips can be dropped with the function `drop.tip`:
#tree before:
tree
#tree after dropping two species:
drop.tip(tree, c("Eliomys_quercinus", "Sicista_tianshanica"))
## We can get a distance matrix from the tree as follows:
dist <- cophenetic(tree)
dist[1:5,1:5]

## ### III. Matching a tree with a dataset
pantheria <- read.table("http://www.esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR05_Aug2008.txt",header=TRUE,sep="\t")
rownames(pantheria) <- gsub(" ", "_", pantheria$MSW05_Binomial)
head(pantheria)

td <- treedata(tree, pantheria)
td$data <- as.data.frame(td$data)

## This part is a bit experimental, so forgive me if it causes problems. It utilizes the dplyr package, which is a powerful
## package for manipulating data frames written by Hadley Wickham. The aRbor package simply wraps the functions from dplyr 
## into versions that work on a list with a tree and data frame, so that the tree and data frame always match each other.
td <- make.treedata(tree, pantheria)
td$dat[td$dat==-999] <- NA
td <- mutate(td, lnMass = log(X5.1_AdultBodyMass_g), lnBMR= log(X5.2_BasalMetRateMass_g), lifespan=X17.1_MaxLongevity_m, desert=ifelse(X28.1_Precip_Mean_mm < 21, 1, 0))
td <- select(td, lnMass, lnBMR, lifespan, desert)
td <- filter(td, !is.na(lnMass), !is.na(lnBMR), !is.na(lifespan), !is.na(desert))

tree <- td$phy
dat <- td$dat

## ## IV. Visualization
## Ape plotting:
plot(tree)
plot(tree, type="fan", cex=0.5)
plot(tree, show.tip.label=FALSE)
tiplabels(pch=21, bg=factor(dat$desert))
nodelabels(cex=0.75,bg="white" )

## Plot a traitgram to visualize a continuous trait:
require(phytools)
phenogram(tree, setNames(dat$lnMass, tree$tip.label), spread.labels=FALSE)

## Another useful tool is to plot the contrasts vs. the node height of the tree:
picMass <- pic(dat$lnMass, tree)
plot(tree, cex=0.5)
nodelabels(pch=21, bg=topo.colors(100)[round((picMass-min(picMass))/diff(range(picMass))*100,0)])

## **Challenge:** Make the size of the points at the nodes also vary with the size of the contrast

## Phytools' scattergram can be useful for exploring correlations between traits across the phylogeny:
fancyDat <- as.matrix(dat[,1:3])
rownames(fancyDat) <- tree$tip.label
fancyTree(tree, type="scattergram", X=fancyDat, fsize=0.3)

## ## V. Simple analyses of traits

## ### Discrete trait models
## We generally use "continuous-time Markov models" to model discrete traits on the phylogeny. These models are
## used to model transitions between discrete character states, whether they are the presence/absence of phenotypic trait
## or nucleotide a specific site. We can fit such a model to a trait and a phylogeny using the `fitDiscrete` function in
## in `geiger`. 
tdDiscrete <- filter(td, !is.na(desert))
trait <- setNames(tdDiscrete$dat$desert, tdDiscrete$phy$tip.label)+1
mER <- fitDiscrete(tdDiscrete$phy, trait, model="ER")
mARD <- fitDiscrete(tdDiscrete$phy, trait, model="ARD")
mER
mARD

## ### Continuous trait models
## For continuous traits, we generally use "Gaussian models" that result in multivariate normal distributions.
## The most basic and familiar of these is Brownian motion. We can fit a Brownian motion model using the function 
## `fitContinuous` from the `geiger` package. 
trait <- setNames(dat$lnMass, tree$tip.label)
mBM <- fitContinuous(tree, trait, model="BM")
mBM
## An more general model than the Brownian motion model is the Ornstein-Uhlenbeck model, which has one additional 
## parameter describing the central tendency of traits to evolve toward an intermediate optimum (think of a rubber 
## band pulling back towards an optimum trait value).
mOU <- fitContinuous(tree, trait, model="OU")
mBM
## Another model is one in which the rate of evolution is initially rapid, but slows over time. This is called
## the early-burst model (EB). Let's fit all models simultaneously using first 1) a loop and then 2) the lapply function.
## This time, we'll use the BMR data instead of the mass data:
trait <- setNames(dat$lnBMR, tree$tip.label)
models <- c("BM", "OU", "EB")
mFits <- list()
for(i in 1:length(models)){
  mFits[[i]] <- fitContinuous(tree, trait, model=models[i])
}

#Above can also be accomplished by using lapply:
mFits <- lapply(1:length(models), function(i) fitContinuous(tree, trait, model=models[i]))

## Model selection can be conducted by comparing AIC values:
aiccs <- lapply(mFits, function(x) x$opt$aicc)
names(aiccs) <- models
aiccs

## ### Phylogenetic signal
## Phylogenetic signal is covariation in traits between related species in the phylogeny.
## There are several methods for testing for phylogenetic signal. Remember to be cautious drawing
## strong conclusions about biological processes from a simple summary statistic such as phylogenetic signal.
## However, in general you can think of two sides of the continuum-- Strong phylogenetic signal indicates that 
## closely related species are more similar to each other than other species. Whereas weak phylogenetic signal 
## indicates that the variation among species is essentially random with respect to phylogeny.

## #### Pagel's Lambda test
## Pagel's lambda test compares the fit of a BM model on the original tree to a BM model with a lambda rescaled
## tree. A lambda equal to 0 produces the following tree:
l0tree <- rescale(tree, "lambda", 0)
plot(l0tree)

## A lambda equal to 1 is the original tree. We can also estimate the value of lambda that best explains the data.
trait <- setNames(dat$lnMass, tree$tip.label) #Make a named data vector
m1 <- fitContinuous(l0tree, trait, model = "BM") #Fit the BM model to the lambda=0 tree
m2 <- fitContinuous(tree, trait, model = "lambda") #Fit the BM model while simultaneously estimating lambda
lnlValues <- c(m1$opt$lnL, m2$opt$lnL) #Create a vector of the resulting log-likelihood values 
names(lnlValues) <- c("Lambda fixed at zero", "Lambda estimated") #Give the vector names
lambdaValue <- m2$opt$lambda #Pull out the estimated lambda value
chisqTestStat <- 2 * (m2$opt$lnL - m1$opt$lnL) #Calculate the Chi square test statistic
chisqPVal <- pchisq(chisqTestStat, 1, lower.tail = F) #Calculate the probability of the Chi-square statistic under the null
aiccScores <- c(m1$opt$aicc, m2$opt$aicc) #Pull out the aicc values
names(aiccScores) <- c("Lambda fixed at zero", "Lambda estimated")#Name them
res <- list(lnlValues = lnlValues, chisqTestStat = chisqTestStat, 
            chisqPVal = chisqPVal, aiccScores = aiccScores, lambdaValue = lambdaValue) #Create a list of results
print(res) #Print out the results

## Make this into a function:
physigLambda <- function(tree, trait){
  m1 <- fitContinuous(l0tree, trait, model = "BM")
  m2 <- fitContinuous(tree, trait, model = "lambda")
  lnlValues <- c(m1$opt$lnL, m2$opt$lnL)
  names(lnlValues) <- c("Lambda fixed at zero", "Lambda estimated")
  lambdaValue <- m2$opt$lambda
  chisqTestStat <- 2 * (m2$opt$lnL - m1$opt$lnL)
  chisqPVal <- pchisq(chisqTestStat, 1, lower.tail = F)
  aiccScores <- c(m1$opt$aicc, m2$opt$aicc)
  names(aiccScores) <- c("Lambda fixed at zero", "Lambda estimated")
  res <- list(lnlValues = lnlValues, chisqTestStat = chisqTestStat, 
              chisqPVal = chisqPVal, aiccScores = aiccScores, lambdaValue = lambdaValue)
  return(res)
}

physigLambda(tree, setNames(dat$lnBMR, tree$tip.label))

## **Challenge:** How could you design a version that measures phylogenetic signal for a discrete trait that
## also uses a lambda transformation?

## Note that some have criticized Pagel's Lambda for being hard to interpret, sensitive to the taxa included
## and be absent an evolutionary model. A particularly nice explanation can be found [here](http://www.carlboettiger.info/2013/10/11/is-it-time-to-retire-pagels-lambda.html) 

## ### Phylogenetic half-life
## Phylogenetic half-life is a measure of phylogenetic signal proposed by Hansen et al. (2008) which
## relies on an Ornstein-Uhlenbeck model of evolution. We will talk more about OU models later for studying 
## adaptive evolution. For interpreting in terms of phylogenetic signal, it helps to consider 
## in terms of the covariance between traits expected under a given phylogeny. For Brownian motion
## the expected covariance is simply equal to the branch length of the tree. 
{
vcvBM <- vcvPhylo(tree, anc.nodes=FALSE)
dend <- as.dendrogram(as.hclust(tree))
heatmap(vcvBM, Rowv = dend, Colv = dend)

par(mfrow=c(3,3), mar=c(0,0,0,0))
image(vcvBM, xaxt="n", yaxt="n")
halflifeseq <- seq(50, 1, length.out=7)
for(i in 1:7){
  vcv <- vcvPhylo(rescale(tree, model="OU", log(2)/halflifeseq[i]), anc.nodes=FALSE)
  image(vcv, xaxt="n", yaxt="n")
}
iidData <- matrix(rnorm(200*length(tree$tip.label), mean(dat$lnMass), 2), ncol=length(tree$tip.label))
vcv <- var(iidData)
image(vcv, xaxt="n", yaxt="n")
}

mOU <- fitContinuous(tree, setNames(dat$lnMass, tree$tip.label), model="OU")
str(mOU)
halflife <- log(2)/mOU$opt$alpha
halflife/max(branching.times(tree))

## A (very) rough guide for interpreting: 
## half-life/tree length = 0-0.1 -> Very strong; = 0.1-0.3 -> Strong; 0.3-0.5 -> Moderate; 0.5-1 -> Weak; >1 Probably not distinct from BM. 

## #### PGLS - Determining relationship between traits
## Phylogenetic generalized least squares is a method for assessing whether phenotypic traits are related to one 
## another. This is a general method of the original independent contrasts methods proposed by Felsenstein (1985).
corBM <- corBrownian(1, phy=tree)
corOU <- corMartins(0.001, phy=tree)

pglsBM <- gls(lnBMR~lnMass, correlation=corBM, data=dat)
pglsOU <- gls(lnBMR~lnMass, correlation=corOU, data=dat)

summary(pglsBM)
summary(pglsOU)

## ## Questions & exercises:

## 1. If you have phylogenetic signal in the trait of interest, that means you MUST USE  a phylogenetic
## regression method (e.g. PGLS, PICs) to analyze the trait. Explain (Hint: Look up the assumptions of 
## Ordinary Least Squares regression).

## 2. You conduct model selection using fitContinuous and find that the OU model is a better fit to trait
## evolution over BM. You should prefer an OU correlation structure for any PGLS model you fit. Why or why not?

## 3. Each of the following lines of code produces an error. Fix the error to accomplish the task requested.
## a. Fit a BM model using fitContinuous
ex3a <- readRDS("../data/ex3a.rds") 
phy <- ex3a$tree
trait <- ex3a$dat
try(fitContinuous(phy, trait, model="BM"))
## b. Plot a phenogram of the data
ex3b <- readRDS("../data/ex3b.rds")
phy <- ex3b$tree
trait <- ex3b$dat
try(phenogram(phy, trait[,1], spread.labels=FALSE))

## 4. Below we have a dataset of 23 species of Icterus (New World Orioles) and a posterior distribution of phylogenies from 
## www.birdtree.org based on genetic sequence data (Hackett backbone). Fit a Brownian motion model to each tree for Tarsus length, 
## save all the results as a list, and extract the estimate of the sigma^2 parameter. Plot a histogram of that parameter.
icterus <-  read.table("../data/icterus.dat", header=TRUE, sep="\t")
icterus.trees <- read.nexus("../data/icterus.tre")

icterus.trees

BMfitter <- function(i){
  ## Fill this in with code to A) match tree and data B) select Tarsus as the dataset C) perform a fitContinuous fit
  ## and D) return the fitted object
}

## Use a loop or an lapply call to call the `BMfitter` function for all values of i from 1 to 100

## Use sapply to extract all the estimates of the parameter `sigsq`

## Plot the resulting vector of sigma^2 estimates.

