#' Empirical effect sizes based on latent trait estimates
#'
#' Description FIXME!
#'
#' @param mod a multipleGroup object which estimated only 2 groups
#' @param focal_items a numeric vector indicating which items to include the tests. The
#'   default uses all of the items. Selecting fewer items will result in tests of
#'   'differential bundle functioning' when \code{DIF = FALSE}
#' @param npts number of points to use in the integration. Default is 61
#' @param theta_lim lower and upper limits of the latent trait (theta) to be evaluated, and is
#'   used in conjunction with \code{npts}
#' @param Theta.focal an optional matrix of Theta values from the focal group to be evaluated. If not supplied
#'   the default values to \code{\link{fscores}} will be used in conjunction with the \code{...}
#'   arguments passed
#' @param DIF logical; return a data.frame of item-level imputation properties? If \code{FALSE},
#'   DBF and DTF statistics will be computed
#' @param ref.group either 1 or 2 to indicate which group is considered the 'reference' group. Default
#'   is 1
#' @param plot logical; plot the effects? FIXME!
#' @param ... additional arguments to be passed to \code{\link{fscores}}
#'
#' @author Adam Meade and Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Meade, A. (2010)...... #FIXME
#'
#' @export empirical_ES
#' @examples
#' \dontrun{
#'
#' #no DIF
#' set.seed(12345)
#' a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' itemtype <- rep('dich', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('Ref', N), rep('Focal', N))
#'
#' mod <- multipleGroup(dat, 1, group = group,
#'    invariance = c(colnames(dat)[1:5], 'free_means', 'free_var'))
#' coef(mod, simplify=TRUE)
#'
#' empirical_ES(mod)
#' empirical_ES(mod, DIF=FALSE)
#'
#' empirical_ES(mod, plot=TRUE)
#' empirical_ES(mod, plot=TRUE, DIF=FALSE)
#'
#' ###---------------------------------------------
#DIF
set.seed(12345)
a1 <- a2 <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d1 <- d2 <- matrix(rnorm(15,0,.7),ncol=1)
a2[10:15,] <- a2[10:15,] + rnorm(6, 0, .3)
d2[10:15,] <- d2[10:15,] + rnorm(6, 0, .3)
itemtype <- rep('dich', nrow(a1))
N <- 1000
dataset1 <- simdata(a1, d1, N, itemtype)
dataset2 <- simdata(a2, d2, N, itemtype, mu = .1, sigma = matrix(1.5))
dat <- rbind(dataset1, dataset2)
group <- c(rep('Ref', N), rep('Focal', N))

mod <- multipleGroup(dat, 1, group = group,
   invariance = c(colnames(dat)[1:5], 'free_means', 'free_var'))
coef(mod, simplify=TRUE)

empirical_ES(mod)
empirical_ES(mod, DIF=FALSE)

empirical_ES(mod, plot=TRUE)
empirical_ES(mod, plot=TRUE, DIF=FALSE)

#' }
empirical_ES <- function(mod, Theta.focal = NULL, focal_items = 1L:extract.mirt(mod, 'nitems'),
                 DIF = TRUE, npts = 61, theta_lim=c(-6,6), ref.group = 1, plot=TRUE){
    stopifnot(extract.mirt(mod, 'nfact') == 1L)
    stopifnot(extract.mirt(mod, 'ngroups') == 2L)
    ref <- extract.group(mod, ref.group)
    focal <- extract.group(mod, ifelse(ref.group == 1, 2, 1))
    focal_select <- extract.mirt(mod, 'group') != extract.mirt(mod, 'groupNames')[ref.group]
    if(is.null(Theta.focal)){
        Theta <- fscores(mod, full.scores = TRUE, full.scores.SE = FALSE)
        Theta.focal <- Theta[focal_select, , drop = FALSE]
    } else Theta.focal <- as.matrix(Theta.focal)
    if(sum(focal_select) != nrow(Theta.focal))  
        stop('Theta elements do not match the number of individuals in the focal group')

    # add the rest here
    
    ############# helper function -  Cohen D ###########
    f.cohen.d <- function (vector.1,vector.2){
      f.mean.v1 <- mean(vector.1)
      f.mean.v2 <- mean(vector.2)
      f.sd.v1 <- sd(vector.1)
      f.sd.v2 <- sd(vector.2)
      f.n.1 <- length(vector.1)
      f.n.2 <- length(vector.2)
      f.numerator <- (f.sd.v1^2)*(f.n.1 - 1) + (f.sd.v2^2)*(f.n.2 - 1)
      f.denom <- f.n.1 + f.n.2 - 2
      f.sd.pooled <- sqrt(f.numerator / f.denom)
      f.return <- (f.mean.v1 - f.mean.v2) / f.sd.pooled
      return(f.return)
    }
    
    ############# helper function -  Max D ###########
    get.max.D <- function(f.d){
      f.d.abs <- abs(f.d)
      f.max.D.location <- which.max(f.d.abs) #which.max returns location
      f.max.D <- f.d[f.max.D.location]
      f.max.D.theta <-  focal.theta.obs[f.max.D.location] 
      f.df.d.max <- c(f.max.D.theta,f.max.D)
      return(f.df.d.max)
    }

  # item level stuff needed for test
    focal.theta.mean.obs <- mean(Theta.focal)
    
    
    #### quadrature nodes - just using base R for this ####
    theta.normal   <- seq(theta_lim[1],theta_lim[2],length=npts)
    theta.den   <- dnorm(theta.normal,mean=focal.theta.mean.obs, sd=1)
    theta.density <- theta.den / sum(theta.den)
    rm(theta.den)
    
    nitems <- length(focal_items)
    list.item_ES_foc.obs <- list()
    list.item_ES_ref.obs <- list()
    list.item_ES_foc.nrm <- list()
    list.item_ES_ref.nrm <- list()
    
    ###### compute the expected scores (ES)
    for(i in 1:nitems){
      foc.extract<-extract.item(model,i,group=2)
      ref.extract<-extract.item(model,i,group=1)
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
    
        
#DTF goes here. output df in list [1]

    if(plot){  # DTF plot
      ## return plot object. mirt imports lattice so those functions can be used if you are comfortable
      plot <- lattice::xyplot()
      return(plot)
    }
    
    
    
        
    if(DIF){
#Dif
      
      #### means for each item in dataframe
      mean.ES.foc <- colMeans(df.foc.obs)
      mean.ES.ref <- colMeans(df.ref.obs)
      
      ### dataframe of observed difference scores
      df.dif.obs <- df.foc.obs - df.ref.obs
      df.abs.dif.obs <- abs(df.dif.obs)
      SIDS <- colMeans(df.dif.obs)
      UIDS <- colMeans(df.abs.dif.obs)
      SIDS
      UIDS
      
      ################## MAX D section ################
      item.max.D <- apply(df.dif.obs,2,get.max.D)
      mat.item.max.d <- t(item.max.D)
      colnames(mat.item.max.d) <- c('theta of max D',"max D")
      mat.item.max.d
      
      
      
      list.item_CohenD.obs <- list()
      for(i in 1:nitems){
        list.item_CohenD.obs[i] <- f.cohen.d(df.foc.obs[,i],df.ref.obs[,i])
      }
      ESSD <- unlist(list.item_CohenD.obs) #vector
      
      
      ############# normal distribution #########
      df.dif.nrm <- df.foc.nrm - df.ref.nrm        # DF in ES at each level of theta
      weighted.dif.nrm <- apply(df.dif.nrm,2, function(x) x*theta.density)
      SIDN <- colSums(weighted.dif.nrm)           #now sum these
      
      df.abs.dif.nrm <- abs(df.dif.nrm)            # abs(DF in ES at each level of theta) 
      weighted.dif.abs.nrm <- apply(df.abs.dif.nrm,2, function(x) x*theta.density)
      UIDN <- colSums(weighted.dif.abs.nrm)
      
      df.item.output <- round(data.frame(SIDS,UIDS,SIDN,UIDN,ESSD,mat.item.max.d,mean.ES.foc,mean.ES.ref),3)
      row.names(df.item.output)<-paste0("item.",1:nrow(df.item.output))
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      if(plot){
        ## return plot object. mirt imports lattice so those functions can be used if you are comfortable
        plot <- lattice::xyplot()
        return(plot)
      }
      
      
    } 
    
    return(ret)
}
