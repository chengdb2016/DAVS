# DAVS-Q: Towards unique and unbiased causal effect estimation from data with hidden variables
1. Some pakeages are necessary to run the code:

install.packages("Matching");library("rmatio")

install.packages("Zelig");install.packages("Metrics")

install.packages("CBPS");install.packages("pcalg")

install.packages("MatchIt");install.packages("cobalt")

install.packages("fastmatch");install.packages("stdReg")

-----------------------------------------------
2. Changing the path into your local envionment 

setwd("~/DAVS")

source("./R/DAVS.R")

source("./R/Create.Ck.R")

source("./R/Estimator.R")

source("./R/Graphical.functions.R")

expdat <- read.mat("./datasets/barley_1w_05.mat")

