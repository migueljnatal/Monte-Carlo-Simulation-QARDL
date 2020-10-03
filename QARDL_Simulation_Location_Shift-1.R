

setwd("C:/Users/MIGUEL/OneDrive/Mestrado PPGEst/DISSERTAÇÃO")

#install.packages("knitr")
#install.packages("jtools")
#install.packages("conquer")
library(conquer)
library(jtools)
library(knitr)
library(quantreg)  
source('sqr_function.r') 


# Packages for QADL and Table visualization @github

path =  paste("C://Users//MIGUEL///OneDrive//Mestrado PPGEst//DISSERTAÇÃO//")

knitr::opts_knit$set(root.dir = normalizePath(path))
ipak <- function(pkg){
  new_pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new_pkg)) 
    install.packages(new_pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
#Create list of packages, call function and remove list and function
packages <- c("quantmod", 'quantreg', "forecast", "readxl", "tidyverse", "pastecs", "ggpubr", "ggthemes","tinytex")
ipak(packages)
rm(packages, ipak, new_pkg)

# The variables

x = list()
nrep = 5000                                         # number of replications
u = list()
Y = list()
ones = list()
l.x = list()
l.Y = list()
X = list()
X1 = list()
taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
A = array(0, dim = c(length(taus), 4, nrep))           # arrays for storing the simulated coef; A stands for rq, B & C for conquer and D for sqr functions
B = array(0, dim = c(length(taus), 4, nrep))
C = array(0, dim = c(length(taus), 4, nrep))
D = array(0, dim = c(length(taus), 4, nrep))
h_scale = 1
z = list()
iqr_z = list() 
sigma_z = list()
h = list()

# Running the simulation

for(j in 1:nrep){
  x[[j]] = rnorm(n = 200, 0, 1)
  u[[j]] = rnorm(n = 200, 0, 1)
  Y[[j]] = 0 
  
  # The model to be estimated  
  
  for(t in 2:200){
    
    Y[[j]][1] = 0
    
    Y[[j]][t] = 0.5*Y[[j]][t-1] + 0.5*x[[j]][t] + 0.5*x[[j]][t-1] + u[[j]][t]
    
    
  }
  ones[[j]] = rep(1, 200)
  
  l.x[[j]] <- lag(x[[j]], n=1L)
  l.x[[j]][1] = 0                 # inserting a 0 for the 1st observation in the lagged x variable (instead of removing the NA value)
  
  l.Y[[j]] <- lag(Y[[j]], n=1L)
  l.Y[[j]][1] = 0                 # same thing for Y  
  
  X[[j]] = cbind(l.Y[[j]], x[[j]], l.x[[j]]) 
  X1[[j]] = cbind(ones[[j]],l.Y[[j]], x[[j]], l.x[[j]]) 
  
  
  
  for (i in 1:length(taus)){
    qadl_rq <- rq(Y[[j]]~ l.Y[[j]] + x[[j]] + l.x[[j]] , tau=taus[i])  
    #print(summary(qadl_rq))
    A[i,,j] = qadl_rq$coefficients # stardard QR 
    
    
    qadl_conquer <- conquer(X[[j]], Y[[j]], tau=taus[i], kernel = "Gaussian", h = 0)
    B[i,,j] = conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i])$coef                   #1st conquer estimation               # qadl_conquer$coef 
    # Getting the optimal bandwidth
    z[[j]] = Y[[j]] - cbind(rep(0,length(Y[[j]])),as.matrix(X[[j]]))%*%B[i,,j]                                             # A[i,,j]    #*B[i,,j]   #%*% B[, ,j]
    iqr_z[[j]] = quantile(z[[j]], .75) - quantile(z[[j]], .25)
    
    sigma_z[[j]] = min(sd(z[[j]]), ((iqr_z[[j]])/1.34898))
    h[[j]] = (1.06*sigma_z[[j]])/(length(Y[[j]])^(1/5))
    
    qadl_conquer_2 <- conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i], h = h[[j]])   # 2nd conquer estimation (optimal bandwidth)
    C[i,,j] = qadl_conquer_2$coeff
    
    
    #qadl_sqr <- sqr(X1[[j]], Y[[j]], tau = taus[i], h = 'rule-of-thumb', estimate_sd = TRUE, return_Hessian = TRUE) 
    qadl_sqr = sqr(h='rule-of-thumb', h_scale=h_scale, tau=taus[i],x=X1[[j]],y=Y[[j]], initial_value='standardQR', estimate_sd = TRUE, order=1, return_Hessian = TRUE)
    D[i,,j] = qadl_sqr$sqrfit[,1]  
    
  } 
}    


