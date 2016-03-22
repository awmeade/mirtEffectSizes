rm(list=ls())   # clear out old junk then read in data

#first load some needed libraries
library("psych")
library("mirt")


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
# ### focal thetas these will be input into the function, check their mean

source("empirical_ES.R")
effect.sizes.out <- empirical_ES(mod = model_anchor)
effect.sizes.out$test.level.results
effect.sizes.out$item.level.results
