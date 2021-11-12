#### load Packages ##########
library(Matching);library(rmatio)
library(Zelig);library(Metrics)
library(CBPS);library(pcalg)
library(MatchIt);library(cobalt)
library(fastmatch)
setwd("~/DAVS")
source("./R/DAVS.R")
source("./R/Create.Ck.R")
source("./R/Estimator.R")
source("./R/Graphical.functions.R")
alpha =0.05
expdat <- read.mat("./datasets/barley_1w_05.mat")

dat <- expdat$data
dat<-data.frame(dat)
expdat <-dat
Tr <- expdat$X17
Y <- expdat$X47
expdat<-expdat[,-c(17,47)]
expdat <- cbind(expdat,Tr,Y)
ACE_NULL <- mean(expdat$Y[expdat$Tr == 1])-mean(expdat$Y[expdat$Tr == 0])



backdoor.rest<-est_reg_bin(c(4,8),expdat) 
backdoor.ACE<-backdoor.rest[1]
backdoor.sd<-backdoor.rest[2] 

##Removing serval variables as hidden variables
expdat <- expdat[,-c(8,36,21)]
starttime<-proc.time()
suffStat_bin<-list(dm=expdat,adaptDF=FALSE)
fci.est <- rfci(suffStat_bin, indepTest = binCItest, p = ncol(expdat), alpha=alpha)
fci.est@amat[length(expdat)-1,length(expdat)]<-2
fci.est@amat[length(expdat),length(expdat)-1]<-3
pag <- fci.est@amat
runtime<-proc.time()-starttime
rfci_Runtime<-runtime[1]
# plot(fci.est)

Q <- 2
starttime<-proc.time()
est_bin <- Davs.bin.causaleffect(expdat,Q,pag,alpha=alpha,possdsep="small",models="regression")
runtime<-proc.time()-starttime
DASE_Q_Runtime<-runtime[1]+rfci_Runtime
Bias_dase <- Cal_bias(est_bin$ACE_DCE,backdoor.ACE)
