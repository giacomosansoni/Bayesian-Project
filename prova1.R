library(gRbase)
library(GLSME)
library (BCDAG)
library(gss)
library(graph)
library(Rcpp)

setwd("C:/Users/Jack/Documents/Scuola/Universita/Magistrale/V/Bayesian Statistics/Project/R")
source("propose_DAG.R")
source("operation.R")
source("acceptreject_DAG.R")
source("new_bcdag.R")
source("update_DAG.R")

q <- 4
DAG <- matrix(c(0,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0), nrow = q)
DAG
outDL <- rDAGWishart(n = 1, DAG = DAG, a = q, U = diag(1, q))
outDL$D
outDL$L

q <- 8
set.seed(123)
DAG <- rDAG(q = q, w = 0.2)
DAG
L <- matrix(runif(n = q*(q-1), min = 0.1, max = 1), q, q)*DAG; 
diag(L) <- 1
L
D <- diag(1, q)
D
Omega <- L%*%solve(D)%*%t(L)
Omega

#we generate n = 200 multivariate zero-mean Gaussian data
X <- mvtnorm::rmvnorm(n = 200, sigma = solve(Omega))
X

#We run function learn DAG to approximate the posterior over DAGs and DAG-parameters
#(collapse = FALSE) by fixing the number of final MCMC iterations and burn-in period as
#S = 5000, B = 1000, while prior hyperparameters as a = q, U = Iq, w = 0.1. We implement the
#approximate MCMC proposal (Section 4.2) by setting fast = TRUE:
out_mcmc <- learn_DAG(S = 5000, burn = 1000, data = X, a = q, U = diag(1,q),
                      w = 0.1, fast = TRUE, save.memory = FALSE, collapse = FALSE)
out_mcmc$Graphs #null, why?
out_mcmc$L
out_mcmc$D
get_edgeprobs(learnDAG_output = out_mcmc) #diverso
get_MAPdag(learnDAG_output = out_mcmc) #uguale
get_MPMdag(learnDAG_output = out_mcmc) #uguale

out_mcmc <- update_DAG(DAG, X, a=q, U=diag(1,q), w=0.1, fast=TRUE, collapse=FALSE)
