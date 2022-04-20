#Fitting script for sampling frame model with free theta to aggregate
library(tidyr)
library(ggplot2)
library(tictoc)
library(dplyr)

source("models/sampling_frames_model_free_theta.R")
exp1 <- read.csv("data/exp1.csv")
exp2 <- read.csv("data/exp2.csv")
##Define exp to be fit here.
experiment <-  "E2"
if (experiment=="E1"){
  data <-  exp1
  samplesizes <-  c(2,8,20)
} else if(experiment=="E2"){
  data <-  exp2
  samplesizes <-  c(8,20)
}


#Specify hypotheses space and conditions
x_grid = 1:6                # locations of the test items
y_grid <- seq(.1, .9, .2)   # possible values for the probability

# construct the hypothesis space: each row specifies the
# value of the unknown function at each of the 6 test points. note that
# this is constructed on the *raw* (i.e., probability) scale. not
# the transformed (i.e., logit) scale. this means that the likelihood
# functions can use these values transparently, but the prior needs to
# logit-transform them
hypotheses <- tibble::as_tibble(expand.grid(
  test1 = y_grid,
  test2 = y_grid,
  test3 = y_grid,
  test4 = y_grid,
  test5 = y_grid,
  test6 = y_grid
))

# conditions in the experiment
conditions <- tibble::as_tibble(expand.grid(
  n = samplesizes,
  #frame = c("category", "property"),
  stringsAsFactors = FALSE
))


#RMSD function (to be minimised)
minRMSD <- function(parms,data,conditions,hypotheses) {
  #tau,rho,sigma,mu,thetaC,thetaP
  parms <-  parmxform(parms,direction=-1)
  tau <- parms[1]
  rho <- parms[2]
  sigma <- parms[3]
  mu <- parms[4]
  thetaP <- parms[5]
  thetaC <- parms[6]

  modelP <- sampling_frames_model(tau,rho,sigma,mu,thetaP,
                                  x_grid=x_grid,conditions=conditions,hypotheses=hypotheses)
  modelC <- sampling_frames_model(tau,rho,sigma,mu,thetaC,
                                  x_grid=x_grid,conditions=conditions,hypotheses=hypotheses)


  modelAll<- c(modelC$response,modelP$response)
  n <- length(data)
  rmsd <- sqrt(sum((modelAll - data)^2)/n)
  return(rmsd)
}

#logit/logistic function
logit <- function(x,direction=1){
  if (direction==1){
    out <- log(x/(1-x))
  } else if (direction==-1){
    out <- 1/(1+exp(-x))
  }
  return(out)
}

#parm transformation function for easy fitting
parmxform <- function(parms,direction=1){
  if (direction==1){
    parms[1:3] <- log(parms[1:3])
    parms[4:6] <- logit(parms[4:6],direction)
  } else if (direction==-1) {
    parms[1:3] <- exp(parms[1:3])
    parms[4:6] <- logit(parms[4:6],direction)
  }
  return(parms)
}

initp <- c(3,1,.5,.05,.5,.5)
initp <-  parmxform(initp,direction=1)

allfits = c()

#Define grid of starting points
taus <-  c(2,8)
rhos <- c(1,3)
sigmas <- c(.2,.8)
mus <- c(.2,.8)
thetaPs <- c(.2,.8)
thetaCs <- c(.2,.8)
grid <- expand.grid(taus,rhos,sigmas,mus,thetaPs,thetaCs)
nStart <- nrow(grid)

bestfits <- tibble(rmsd=NA , tau=NA, rho=NA, sigma=NA, mu=NA, thetaP=NA, thetaC=NA)


#Prepare aggregate
pdataC <-  data %>%
  dplyr::filter(sampling_frame=='category') %>%
  group_by(sample_size,test_item) %>%
  summarise(response = mean(response)) %>%
  ungroup()
pdataP <-  data %>%
  dplyr::filter(sampling_frame=='property') %>%
  group_by(sample_size,test_item) %>%
  summarise(response = mean(response)) %>%
  ungroup()
pResponse <- c(pdataC$response,pdataP$response)/10

for (rowi in seq(1,nStart)){
  initp <- unlist(grid[rowi,])
  cat('Fitting startpoint ', rowi,'/',nStart,'\n')
  cat('Start parms: ',initp,"\n")
  initpr <-  parmxform(initp,direction=1)

  tic()

  #print(pResponse)
  fits <- optim(initpr,minRMSD,gr = NULL, pResponse, conditions, hypotheses)
  fits$par  <-  parmxform(fits$par,direction=-1)

  toc()
  cat("\t Final parms: ",fits$par,"\n")
  if (bestfits$rmsd > fits$value || is.na(bestfits$rmsd)) {
    cat("\t   Updating previous best fit (",bestfits$rmsd, ") with current (", fits$value,"). \n")
    #Update best fit if rmsd is lower
    bestfits[1] = fits$value
    bestfits[2:7] = as.list(fits$par)
    #Save
    save(bestfits,file=paste("bestfits",experiment,"Agg.RData",sep=""))
  }
}


