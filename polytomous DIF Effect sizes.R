rm(list=ls())   # clear out old junk then read in data

#first load some needed libraries
library("psych")
#library("lessR")
library("mirt")

# difficulty (b) = -d/a 

################################################################################################
#                                              Functions
################################################################################################


make.data <- function(N){
  set.seed(1234)
  a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
  d <- matrix(rnorm(15,0,.7),ncol=1)
  d1 <- d2 <- cbind(d, d-1, d-2)  # b parameters for both groups
  d2[13:15, ] <- d1[13:15, ] + 1  # here is the DIF
  itemtype <- rep('graded', nrow(a))
  dataset1 <- simdata(a, d1, N, itemtype, mu=1.0)
  dataset2 <- simdata(a, d2, N, itemtype)
  dat <- rbind(dataset1, dataset2)
  return(dat)
}

################################################# END FUNCTIONS ########################


N <- 1000
dat <- make.data(N)
group <- c(rep('Ref', N), rep('Foc', N))



#foc.group <- 'Foc'
itemnames <- colnames(dat)
anc.items.names <- itemnames[c(2,8,9,10,12)]
#test.items <- c(1,2:7,11,13:15)
model_anchor <- multipleGroup(dat, model = 1, group = group,
  invariance = c(anc.items.names, 'free_means', 'free_var'))  # sets mean of group 1 to 0 so ref as 1

# get factor means from model
# coef <- coef(model_anchor,simplify = TRUE)
# theta.mean.1 <- coef[[1]][[2]]
# theta.mean.2 <- coef[[2]][[2]]
# theta.mean.dif <- theta.mean.2 - theta.mean.1 
# 
# 
# ### focal thetas these will be input into the function, check their mean
theta.obs <- fscores(model_anchor, full.scores = TRUE)
focal.theta.obs <- theta.obs[1:1000,]
focal.theta.mean.obs <- mean(focal.theta.obs)
# 
# theta.mean.dif
# focal.theta.mean.obs
# #set mean difference
# if(round(focal.theta.mean.obs,1) == 0){
#   use.focal.theta.mean = theta.mean.dif
# }else{
#   use.focal.theta.mean = focal.theta.mean.obs
# }
# use.focal.theta.mean



#### quadrature nodes - just using base R for this ####
theta.normal   <- seq(-4,4,length=30)
theta.den   <- dnorm(theta.normal,mean=focal.theta.mean.obs, sd=1)
theta.density <- theta.den / sum(theta.den)
rm(theta.den)


list.item_ES_foc.obs <- list()
list.item_ES_ref.obs <- list()
list.item_ES_foc.nrm <- list()
list.item_ES_ref.nrm <- list()

###### compute the expected scores (ES)
for(i in 1:ncol(dat)){
  foc.extract<-extract.item(model_anchor,i,group=1)
  ref.extract<-extract.item(model_anchor,i,group=2)
  foc.ES.obs <- expected.item(foc.extract,focal.theta.obs)
  ref.ES.obs <- expected.item(ref.extract,focal.theta.obs)
  foc.ES.nrm <- expected.item(foc.extract,theta.normal)
  ref.ES.nrm <- expected.item(ref.extract,theta.normal)
  
  list.item_ES_foc.obs[[i]] <- foc.ES.obs
  list.item_ES_ref.obs[[i]] <- ref.ES.obs
  list.item_ES_foc.nrm[[i]] <- foc.ES.nrm
  list.item_ES_ref.nrm[[i]] <- ref.ES.nrm
}
rm(foc.extract,ref.extract,foc.ES.obs,ref.ES.obs,foc.ES.nrm,ref.ES.nrm)
# put these in dataframe
df.ref.obs <- do.call("cbind",list.item_ES_ref.obs) 
df.foc.obs <- do.call("cbind",list.item_ES_foc.obs) 
df.ref.nrm <- do.call("cbind",list.item_ES_ref.nrm) 
df.foc.nrm <- do.call("cbind",list.item_ES_foc.nrm) 


head(df.ref.obs)
head(df.foc.obs)
df.foc.nrm
df.ref.nrm

#### means for each item in dataframe
#mean.ES.foc <- colMeans(df.foc.obs)
#mean.ES.ref <- colMeans(df.ref.obs)

### dataframe of observed difference scores
df.dif.obs <- df.ref.obs - df.foc.obs
df.abs.dif.obs <- abs(df.dif.obs)
SIDS <- colMeans(df.dif.obs)
UIDS <- colMeans(df.abs.dif.obs)
SIDS
UIDS

################## MAX D section ################
get.max.D <- function(f.d){
  f.d.abs <- abs(f.d)
  f.max.D.location <- which.max(f.d.abs) #which.max returns location
  f.max.D <- f.d[f.max.D.location]
  f.max.D.theta <-  focal.theta.obs[f.max.D.location] 
  f.df.d.max <- c(f.max.D.theta,f.max.D)
  return(f.df.d.max)
}
list.max.D <- apply(df.dif.obs,2,get.max.D)
list.max.D



############# Cohen D section ###########
#pooled SD
# v.sd.foc <- apply(df.foc.obs,2,sd) # get SD per sample
# v.sd.ref <- apply(df.ref.obs,2,sd)
# 
# f2 <- function(f2.d, f2.n){
#     f2.return <- (f2.d^2)*(f2.n-1)
#     return(f2.return)
# } 
# f.numerator.1 <- sapply(v.sd.foc,f2,N)
# f.numerator.2 <- sapply(v.sd.ref,f2,N)
# f.numerator <- f.numerator.1 + f.numerator.2
# f.denom <- 2*N - 2
# f.pooled.SD <- sqrt(f.numerator/f.denom)
# 
# v.mean.ES.foc.obs <- colMeans(df.foc.obs)
# v.mean.ES.ref.obs <- colMeans(df.ref.obs)
# v.ESSD <- (v.mean.ES.foc.obs - v.mean.ES.ref.obs) / f.pooled.SD  
# rm(f.numerator.1,f.numerator.2,f.numerator,f.denom,f.pooled.SD)

f.cohen.d <- function (vector.1,vector.2){
  f.mean.v1 <- mean(vector.1)
  f.mean.v2 <- mean(vector.2)
  f.sd.v1 <- sd(vector.1)
  f.sd.v2 <- sd(vector.2)
  f.n.1 <- length(vector.1)
  f.n.2 <- length(vector.2)
  f.numerator <- (f.sd.v1^2)*(f.n.1 - 1) + (f.sd.v2^2)*(f.n.2 - 1)
  f.denom <- f.n.1 + f.n.2 - 2
  f.sd.pooled <- f.numerator / f.denom
  f.return <- (f.mean.v2 - f.mean.v1) / f.sd.pooled
  return(f.return)
}


list.item_CohenD.obs <- list()
for(i in 1:ncol(dat)){
  list.item_CohenD.obs[i] <- f.cohen.d(df.foc.obs[,i],df.ref.obs[,i])
}
#ESSD <- do.call("cbind",list.item_CohenD.obs) #dataframe 
ESSD <- unlist(list.item_CohenD.obs) #vector


###### Test level observed  #######
STDS <- sum(SIDS)
UTDS <- sum(UIDS)
test.ES.foc.obs <- rowSums(df.foc.obs)  ### avg across items
test.ES.ref.obs <- rowSums(df.ref.obs)
test.ES.dif.obs <- test.ES.foc.obs - test.ES.ref.obs
test.ES.dif.abs.obs <- abs(test.ES.dif.obs)
UETSDS <- mean(test.ES.dif.abs.obs)

ETSSD <- f.cohen.d(test.ES.ref.obs, test.ES.foc.obs) # cohen's D
ETSSD

test.Dmax <- get.max.D(test.ES.dif.obs)





############################### normal distribution ####################
#### dataframe of normal dist difference scocres
df.dif.nrm <- df.foc.nrm - df.ref.nrm  # dif in ES at each level of theta
df.abs.dif.nrm <- abs(df.dif.nrm)

# weight normal by density
weighted.dif.nrm <- apply(df.dif.nrm,2, function(x) x*theta.density)
weighted.dif.abs.nrm <- apply(df.abs.dif.nrm,2, function(x) x*theta.density)
SIDN <- colSums(weighted.dif.nrm)           #now sum these
UIDN <- colSums(weighted.dif.abs.nrm)
SIDN
UIDN

# test level
test.ES.foc.nrm <- rowSums(df.foc.nrm)  ### sum across items at each theta
test.ES.ref.nrm <- rowSums(df.ref.nrm)  ### sum across items at each theta
test.ES.dif.nrm <- rowSums(df.dif.nrm)  ### dif at each theta

#test.dif.nrm <- rowSums(df.dif.nrm)  # sum across items

Starks.DTF <- sum(test.ES.dif.nrm*theta.density)

#test.dif.abs.nrm <- rowMeans(df.abs.dif.nrm) #avg across items
test.ES.abs.dif.nrm <- rowSums(df.abs.dif.nrm)
UDTFR <- sum(test.ES.abs.dif.nrm*theta.density)
UETSDN <- sum(test.ES.dif.nrm*theta.density)
