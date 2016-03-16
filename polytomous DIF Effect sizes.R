rm(list=ls())   # clear out old junk then read in data

#first load some needed libraries
library("psych")
library("lessR")
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
  dataset1 <- simdata(a, d1, N, itemtype)
  dataset2 <- simdata(a, d2, N, itemtype)
  dat <- rbind(dataset1, dataset2)
  return(dat)
}

################################################# END FUNCTIONS ########################


N <- 1000
dat <- make.data(N)
group <- c(rep('D1', N), rep('D2', N))




itemnames <- colnames(dat)
anc.items.names <- itemnames[c(2,8,9,10,12)]
#test.items <- c(1,2:7,11,13:15)
model_anchor <- multipleGroup(dat, model = 1, group = group,
  invariance = c(anc.items.names, 'free_means', 'free_var'))
coef(model_anchor,simplify = TRUE)

theta <- fscores(model_anchor, full.scores = TRUE)
focal.theta <- theta[1:N,]

list.item_ES <- list()
for(i in 1:ncol(dat)){
  foc.extract<-extract.item(model_anchor,i,group=1)
  ref.extract<-extract.item(model_anchor,i,group=2)
  
  foc.ES <- expected.item(foc.extract,focal.theta)
  ref.ES <- expected.item(ref.extract,focal.theta)
  
  list.item_ES[[i]] <- data.frame(foc.ES,ref.ES)
}
list.item_ES


lapply(list.item_ES,cor)

