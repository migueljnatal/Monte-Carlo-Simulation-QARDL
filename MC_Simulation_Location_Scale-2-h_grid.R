

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
nrep = 5000                 # number of replications
u = list()
Y = list()
ones = list()
l.x = list()
l.Y = list()
X = list()
X1 = list()
taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
A = array(0, dim = c(length(taus), 4, nrep))
B = array(0, dim = c(length(taus), 4, nrep))
h_grid = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)
C = array(0, dim = c(length(taus), 4, nrep, length(h_grid)))
D = array(0, dim = c(length(taus), 4, nrep))
h_scale = 1
z = list()
iqr_z = list() 
sigma_z = list()
h = list()

h_grid = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1) 

# Running the simulation

for(j in 1:nrep){
  x[[j]] = rchisq(n = 200, 3)
  u[[j]] = rt(n = 200, 3)
  Y[[j]] = 0 
  
  # The model to be estimated  
  
  for(t in 2:200){
    
    Y[[j]][1] = 0
    
    Y[[j]][t] = 0.5*Y[[j]][t-1] + 0.5*x[[j]][t] + 0.5*x[[j]][t-1] + (0.2*x[[j]][t])*u[[j]][t]
    
    
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
    A[i,,j] = qadl_rq$coefficients # standard QR 
    
    
    qadl_conquer <- conquer(X[[j]], Y[[j]], tau=taus[i], kernel = "Gaussian", h = 0)
    B[i,,j] = conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i])$coef                   # 1st conquer estimation               #qadl_conquer$coef 
    for (k in 1:length(h_grid)){
      
      qadl_conquer_2 <- conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i], h = h_grid[[k]])
      C[i,,j,k] = qadl_conquer_2$coeff
    }
    
    # Getting the optimal bandwidth
    #z[[j]] = Y[[j]] - cbind(rep(0,length(Y[[j]])),as.matrix(X[[j]]))%*%B[i,,j]                                             # A[i,,j]    #*B[i,,j]   #%*% B[, ,j]
    #iqr_z[[j]] = quantile(z[[j]], .75) - quantile(z[[j]], .25)
    
    #sigma_z[[j]] = min(sd(z[[j]]), ((iqr_z[[j]])/1.34898))
    #h[[j]] = (1.06*sigma_z[[j]])/(length(Y[[j]])^(1/5))
    
    #qadl_conquer_2 <- conquer(as.matrix(X[[j]]), Y[[j]], tau=taus[i], h = h[[j]])   # 2nd conquer estimation (optimal bandwidth)
    #C[i,,j] = qadl_conquer_2$coeff
    
    
    #qadl_sqr <- sqr(X1[[j]], Y[[j]], tau = taus[i], h = 'rule-of-thumb', estimate_sd = TRUE, return_Hessian = TRUE) 
    #qadl_sqr = sqr(h='rule-of-thumb', h_scale=h_scale, tau=taus[i],x=X1[[j]],y=Y[[j]], initial_value='standardQR', estimate_sd = TRUE, order=1, return_Hessian = TRUE)
    #D[i,,j] = qadl_sqr$sqrfit[,1]  
    
  }
}   


# Evaluating the estimators

# rq - standard QR


# 1st quantile: 0.1

bias_alfa_0.1_rq = mean(A[1,2,]) - 0.5 
rmse_alfa_0.1_rq = sqrt(mean((A[1,2,]-0.5)^2))

bias_beta1_0.1_rq = mean(A[1,3,]) - 0.5 
rmse_beta1_0.1_rq = sqrt(mean((A[1,3,]-0.5)^2))

bias_beta2_0.1_rq = mean(A[1,4,]) - 0.5 
rmse_beta2_0.1_rq = sqrt(mean((A[1,4,]-0.5)^2))

# 2nd quantile: 0.25

bias_alfa_0.25_rq = mean(A[2,2,]) - 0.5 
rmse_alfa_0.25_rq = sqrt(mean((A[2,2,]-0.5)^2))

bias_beta1_0.25_rq = mean(A[2,3,]) - 0.5 
rmse_beta1_0.25_rq = sqrt(mean((A[2,3,]-0.5)^2))

bias_beta2_0.25_rq = mean(A[2,4,]) - 0.5 
rmse_beta2_0.25_rq = sqrt(mean((A[2,4,]-0.5)^2))

# 3rd quantile: 0.5

bias_alfa_0.5_rq = mean(A[3,2,]) - 0.5 
rmse_alfa_0.5_rq = sqrt(mean((A[3,2,]-0.5)^2))

bias_beta1_0.5_rq = mean(A[3,3,]) - 0.5 
rmse_beta1_0.5_rq = sqrt(mean((A[3,3,]-0.5)^2))

bias_beta2_0.5_rq = mean(A[3,4,]) - 0.5 
rmse_beta2_0.5_rq = sqrt(mean((A[3,4,]-0.5)^2))

# 4th quantile: 0.75

bias_alfa_0.75_rq = mean(A[4,2,]) - 0.5 
rmse_alfa_0.75_rq = sqrt(mean((A[4,2,]-0.5)^2))

bias_beta1_0.75_rq = mean(A[4,3,]) - 0.5 
rmse_beta1_0.75_rq = sqrt(mean((A[4,3,]-0.5)^2))

bias_beta2_0.75_rq = mean(A[4,4,]) - 0.5 
rmse_beta2_0.75_rq = sqrt(mean((A[4,4,]-0.5)^2))

# 5th quantile: 0.9

bias_alfa_0.9_rq = mean(A[5,2,]) - 0.5 
rmse_alfa_0.9_rq = sqrt(mean((A[5,2,]-0.5)^2))

bias_beta1_0.9_rq = mean(A[5,3,]) - 0.5 
rmse_beta1_0.9_rq = sqrt(mean((A[5,3,]-0.5)^2))

bias_beta2_0.9_rq = mean(A[5,4,]) - 0.5 
rmse_beta2_0.9_rq = sqrt(mean((A[5,4,]-0.5)^2))


# conquer - smoothed QR

# alpha

bias_alfa_0.1_conquer = list()
bias_alfa_0.25_conquer = list()
bias_alfa_0.5_conquer = list()
bias_alfa_0.75_conquer = list()
bias_alfa_0.9_conquer = list()

rmse_alfa_0.1_conquer = list()
rmse_alfa_0.25_conquer = list()
rmse_alfa_0.5_conquer = list()
rmse_alfa_0.75_conquer = list()
rmse_alfa_0.9_conquer = list()

# beta1

bias_beta1_0.1_conquer = list()
bias_beta1_0.25_conquer = list()
bias_beta1_0.5_conquer = list()
bias_beta1_0.75_conquer = list()
bias_beta1_0.9_conquer = list()

rmse_beta1_0.1_conquer = list()
rmse_beta1_0.25_conquer = list()
rmse_beta1_0.5_conquer = list()
rmse_beta1_0.75_conquer = list()
rmse_beta1_0.9_conquer = list()

# beta2

bias_beta2_0.1_conquer = list()
bias_beta2_0.25_conquer = list()
bias_beta2_0.5_conquer = list()
bias_beta2_0.75_conquer = list()
bias_beta2_0.9_conquer = list() 

rmse_beta2_0.1_conquer = list()
rmse_beta2_0.25_conquer = list()
rmse_beta2_0.5_conquer = list()
rmse_beta2_0.75_conquer = list()
rmse_beta2_0.9_conquer = list()




for (r in 1:length(h_grid)){
  
  # 1st quantile: 0.1 
  
  bias_alfa_0.1_conquer[r] = mean(C[1,2,,r]) - 0.5 
  rmse_alfa_0.1_conquer[r] = sqrt(mean((C[1,2,,r]-0.5)^2))
  
  bias_beta1_0.1_conquer[r] = mean(C[1,3,,r]) - 0.5 
  rmse_beta1_0.1_conquer[r] = sqrt(mean((C[1,3,,r]-0.5)^2))
  
  bias_beta2_0.1_conquer[r] = mean(C[1,4,,r]) - 0.5 
  rmse_beta2_0.1_conquer[r] = sqrt(mean((C[1,4,,r]-0.5)^2))
  
  # 2nd quantile: 0.25
  
  bias_alfa_0.25_conquer[r] = mean(C[2,2,,r]) - 0.5 
  rmse_alfa_0.25_conquer[r] = sqrt(mean((C[2,2,,]-0.5)^2))
  
  bias_beta1_0.25_conquer[r] = mean(C[2,3,,r]) - 0.5 
  rmse_beta1_0.25_conquer[r] = sqrt(mean((C[2,3,,r]-0.5)^2))
  
  bias_beta2_0.25_conquer[r] = mean(C[2,4,,r]) - 0.5 
  rmse_beta2_0.25_conquer[r] = sqrt(mean((C[2,4,,r]-0.5)^2))
  
  # 3rd quantile: 0.5
  
  bias_alfa_0.5_conquer[r] = mean(C[3,2,,r]) - 0.5 
  rmse_alfa_0.5_conquer[r] = sqrt(mean((C[3,2,,r]-0.5)^2))
  
  bias_beta1_0.5_conquer[r] = mean(C[3,3,,r]) - 0.5 
  rmse_beta1_0.5_conquer[r] = sqrt(mean((C[3,3,,r]-0.5)^2))
  
  bias_beta2_0.5_conquer[r] = mean(C[3,4,,r]) - 0.5 
  rmse_beta2_0.5_conquer[r] = sqrt(mean((C[3,4,,r]-0.5)^2))
  
  # 4th quantile: 0.75
  
  bias_alfa_0.75_conquer[r] = mean(C[4,2,,r]) - 0.5 
  rmse_alfa_0.75_conquer[r] = sqrt(mean((C[4,2,,r]-0.5)^2))
  
  bias_beta1_0.75_conquer[r] = mean(C[4,3,,r]) - 0.5 
  rmse_beta1_0.75_conquer[r] = sqrt(mean((C[4,3,,r]-0.5)^2))
  
  bias_beta2_0.75_conquer[r] = mean(C[4,4,,r]) - 0.5 
  rmse_beta2_0.75_conquer[r] = sqrt(mean((C[4,4,,r]-0.5)^2))
  
  # 5th quantile: 0.9
  
  bias_alfa_0.9_conquer[r] = mean(C[5,2,,r]) - 0.5 
  rmse_alfa_0.9_conquer[r] = sqrt(mean((C[5,2,,r]-0.5)^2))
  
  bias_beta1_0.9_conquer[r] = mean(C[5,3,,r]) - 0.5 
  rmse_beta1_0.9_conquer[r] = sqrt(mean((C[5,3,,r]-0.5)^2))
  
  bias_beta2_0.9_conquer[r] = mean(C[5,4,,r]) - 0.5 
  rmse_beta2_0.9_conquer[r] = sqrt(mean((C[5,4,,r]-0.5)^2))
  
}

# Bias and RMSE results

# alfa = 0.5 

bias_alfa_rq = list(bias_alfa_0.1_rq, bias_alfa_0.25_rq, bias_alfa_0.5_rq, bias_alfa_0.75_rq, bias_alfa_0.9_rq)
bias_alfa_conquer = list(bias_alfa_0.1_conquer, bias_alfa_0.25_conquer, bias_alfa_0.5_conquer, bias_alfa_0.75_conquer, bias_alfa_0.9_conquer)

rmse_alfa_rq = list(rmse_alfa_0.1_rq, rmse_alfa_0.25_rq, rmse_alfa_0.5_rq, rmse_alfa_0.75_rq, rmse_alfa_0.9_rq) 
rmse_alfa_conquer = list(rmse_alfa_0.1_conquer, rmse_alfa_0.25_conquer, rmse_alfa_0.5_conquer, rmse_alfa_0.75_conquer, rmse_alfa_0.9_conquer)

# beta1 = 0.5

bias_beta1_rq = list(bias_beta1_0.1_rq, bias_beta1_0.25_rq, bias_beta1_0.5_rq, bias_beta1_0.75_rq, bias_beta1_0.9_rq)
bias_beta1_conquer = list(bias_beta1_0.1_conquer, bias_beta1_0.25_conquer, bias_beta1_0.5_conquer, bias_beta1_0.75_conquer, bias_beta1_0.9_conquer)

rmse_beta1_rq = list(rmse_beta1_0.1_rq, rmse_beta1_0.25_rq, rmse_beta1_0.5_rq, rmse_beta1_0.75_rq, rmse_beta1_0.9_rq) 
rmse_beta1_conquer = list(rmse_beta1_0.1_conquer, rmse_beta1_0.25_conquer, rmse_beta1_0.5_conquer, rmse_beta1_0.75_conquer, rmse_beta1_0.9_conquer)

 
# beta2 = 0.5

bias_beta2_rq = list(bias_beta2_0.1_rq, bias_beta2_0.25_rq, bias_beta2_0.5_rq, bias_beta2_0.75_rq, bias_beta2_0.9_rq)
bias_beta2_conquer = list(bias_beta2_0.1_conquer, bias_beta2_0.25_conquer, bias_beta2_0.5_conquer, bias_beta2_0.75_conquer, bias_beta2_0.9_conquer)

rmse_beta2_rq = list(rmse_beta2_0.1_rq, rmse_beta2_0.25_rq, rmse_beta2_0.5_rq, rmse_beta2_0.75_rq, rmse_beta2_0.9_rq) 
rmse_beta2_conquer = list(rmse_beta2_0.1_conquer, rmse_beta2_0.25_conquer, rmse_beta2_0.5_conquer, rmse_beta2_0.75_conquer, rmse_beta2_0.9_conquer)

# Saving main results 

save(bias_alfa_rq, bias_alfa_conquer, file = "bias_alfa_lsc2.RData") 
save(bias_beta1_rq, bias_beta1_conquer, file = "bias_beta1_lsc2.RData")
save(bias_beta2_rq, bias_beta2_conquer, file = "bias_beta2_lsc2.RData")

save(rmse_alfa_rq, rmse_alfa_conquer,  file = "rmse_alfa_lsc2.RData")
save(rmse_beta1_rq, rmse_beta1_conquer,  file = "rmse_beta1_lsc2.RData")
save(rmse_beta2_rq, rmse_beta2_conquer, file = "rmse_beta2_lsc2.RData")
