



############# Function will start here. 
f.effect.sizes <- local(function(focal.theta.obs, model){
  
  ############# helper function -  Cohen D section ###########
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
  
  get.max.D <- function(f.d){
    f.d.abs <- abs(f.d)
    f.max.D.location <- which.max(f.d.abs) #which.max returns location
    f.max.D <- f.d[f.max.D.location]
    f.max.D.theta <-  focal.theta.obs[f.max.D.location] 
    f.df.d.max <- c(f.max.D.theta,f.max.D)
    return(f.df.d.max)
  }
  
  
  require("mirt")
  focal.theta.mean.obs <- mean(focal.theta.obs)
  
  
  #### quadrature nodes - just using base R for this ####
  theta.normal   <- seq(-4,4,length=30)
  theta.den   <- dnorm(theta.normal,mean=focal.theta.mean.obs, sd=1)
  theta.density <- theta.den / sum(theta.den)
  rm(theta.den)
  
  nitems <- nrow(coef(model,simplify=TRUE)[[1]][[1]])
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
  
  
  
  ############################## TEST LEVEL ##################################
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
  ETSSD <- f.cohen.d(ETS.foc.obs, ETS.ref.obs) # cohen's D
  test.Dmax <- get.max.D(ETS.dif.obs)
  
  
  
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
  
  
  out.test.stats <- round(c(STDS,UTDS,UETSDS,ETSSD,Starks.DTFR,UDTFR,UETSDN,test.Dmax),3)
  out.test.names <- c("STDS","UTDS","UETSDS","ETSSD","Starks.DTFR","UDTFR","UETSDN","theta.of.max.test.D","Test.Dmax")
  df.test.output <- data.frame(out.test.names,out.test.stats)
  df.test.output
  
  
  
  
  ############# plots #############################
#   list.plots <- list()
   plot.df <- data.frame(focal.theta.obs,ETS.foc.obs,ETS.ref.obs)
   plot.df <- plot.df[order(focal.theta.obs),]

     plot(plot.df$focal.theta.obs,plot.df$ETS.ref.obs, 
          type="o",xlab="Focal Group Theta",
          ylab="Expected Test Score",
          main="ETS using Foc and Ref Parameters")
   points(plot.df$focal.theta.obs,plot.df$ETS.foc.obs,col="red")
   legend("topleft",pch=c(1,1),legend=c("ref","foc"),col=c("black","red"))
#   

   for(i in 1:nitems){
     #  head(df.ref.obs)
     #head(df.foc.obs)
     plot.df <- data.frame(focal.theta.obs,df.foc.obs[,i],df.ref.obs[,i])
     plot.df <- plot.df[order(focal.theta.obs),]
     the.name <- paste0("Item ",i)
     plot(plot.df$focal.theta.obs,plot.df[,3], 
          type="o",xlab="Focal Group Theta",
          ylab="Expected Score",
          main=the.name)
     points(plot.df$focal.theta.obs,plot.df[,2],col="red")
     legend("topleft",pch=c(1,1),legend=c("ref","foc"),col=c("black","red"))
   }
   
  r.list <- list(item.level.results = df.item.output, test.level.results = df.test.output)
  return(r.list)
}, as.environment(2))

