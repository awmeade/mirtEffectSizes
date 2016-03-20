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


# 
# ### focal thetas these will be input into the function, check their mean
theta.obs <- fscores(model_anchor, full.scores = TRUE)
focal.theta.obs <- theta.obs[1:1000,]
focal.theta.mean.obs <- mean(focal.theta.obs)
# 



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
ESSD <- unlist(list.item_CohenD.obs) #vector


###### Test level observed  #######
STDS <- sum(SIDS)
UTDS <- sum(UIDS)

##### UETSDS
ETS.foc.obs <- rowSums(df.foc.obs)  ### test ES (avg across items)
ETS.ref.obs <- rowSums(df.ref.obs)  ### test ES (avg across items)
ETS.dif.obs <- ETS.foc.obs - ETS.ref.obs  ### df in test scores
UETSDS <- mean(abs(ETS.dif.obs))    ### no cancel across theta, yes across items
##### 

#cohen D and D max
ETSSD <- f.cohen.d(test.ES.ref.obs, test.ES.foc.obs) # cohen's D
test.Dmax <- get.max.D(test.ES.dif.obs)


############################### normal distribution ####################
#### dataframe of normal dist difference scocres
df.dif.nrm <- df.foc.nrm - df.ref.nrm        # DF in ES at each level of theta
weighted.dif.nrm <- apply(df.dif.nrm,2, function(x) x*theta.density)
SIDN <- colSums(weighted.dif.nrm)           #now sum these

df.abs.dif.nrm <- abs(df.dif.nrm)            # abs(DF in ES at each level of theta) 
weighted.dif.abs.nrm <- apply(df.abs.dif.nrm,2, function(x) x*theta.density)
UIDN <- colSums(weighted.dif.abs.nrm)

############# test level
ETS.abs.dif.nrm <- rowSums(df.abs.dif.nrm) # [abs(DF in ES at each level of theta)] summed across items
UDTFR <- sum(ETS.abs.dif.nrm*theta.density)
UDTFR.b <- sum(UIDN)

ETS.dif.nrm <- rowSums(df.dif.nrm)  ### ETS dif at each theta, summed across items
Starks.DTFR <- sum(ETS.dif.nrm*theta.density)
Starks.DTFR.b <- sum(SIDN)

ETS.foc.nrm <- rowSums(df.foc.nrm)  ### ETS at each theta (df cancels across items)
ETS.ref.nrm <- rowSums(df.ref.nrm)  ### ETS at each theta (df cancels across items)
ETS.dif.nrm <- ETS.foc.nrm - ETS.ref.nrm   ### DF in ETS at each theta
ETS.abs.dif.nrm <- abs(ETS.dif.nrm)
UETSDN <- sum(ETS.abs.dif.nrm*theta.density)
UETSDN
UDTFR.b


plot.df <- data.frame(focal.theta.obs,ETS.foc.obs,ETS.ref.obs)
plot.df <- plot.df[order(focal.theta.obs),]
plot(plot.df$focal.theta.obs,plot.df$ETS.foc.obs, type="o")
points(plot.df$focal.theta.obs,plot.df$ETS.ref.obs,col="red")
title(main="ETS", font.main=4)
