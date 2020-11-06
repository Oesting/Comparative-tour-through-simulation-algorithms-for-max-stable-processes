############################################################################
########## Approach from DEO (2016) paper (incl. Dieker-Mikosch) ###########
############################################################################

## internal functions: do not use any of them directly!
simu_px_brownresnick <- function(no.simu=1, idx,  N, trend, chol.mat) {
  res <- t(chol.mat)%*%matrix(rnorm(N*no.simu), ncol=no.simu)
  if (!is.matrix(trend)) {
    stopifnot(length(unique(idx))==1)
    res <- exp(t(res - trend))
  } else {
    stopifnot(length(idx)==no.simu)
    res <- exp(t(res - trend[,idx]))   
  }
  return(res/res[cbind(1:no.simu,idx)])
}

simu_px_extremalt <- function(no.simu=1, idx, N, dof, mu, chol.mat) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  if (!is.list(chol.mat)) {  
    res <- t(chol.mat)%*%matrix(rnorm(N*no.simu), ncol=no.simu)
    res <- sqrt((dof+1)/rgamma(no.simu, shape=0.5*(dof+1), scale=2))*(t(res) - res[idx,])  
    res <- pmax(matrix(mu, nrow=no.simu, ncol=N, byrow=TRUE)+res,0)^dof
  } else {
    res <- matrix(NA, nrow=no.simu, ncol=N) 
    for (i in 1:no.simu) {
      res[i,] <- as.vector(t(chol.mat[[idx[i]]])%*%matrix(rnorm(N), ncol=1))
      res[i,] <- sqrt((dof+1)/rgamma(1, shape=0.5*(dof+1), scale=2))*(res[i,] - res[i,idx[i]])
      res[i,] <- pmax(mu[[idx[i]]]+res[i,],0)^dof
    }  
  }
  return(res)
}       


###########################################################################
################## Dieker--Mikosch (sum normalized) #######################
###########################################################################

simu_sumnormal <- function(no.simu=1, coord, vario, corr, dof, type,
                           rel.thresh=1, calc.err=TRUE) {
  
  cat("Dieker-Mikosch simulation...\n")
  if (!is.matrix(coord)) {
    coord <- matrix(coord, ncol=1)
  }
  d <- ncol(coord)    
  N <- nrow(coord)
  
  stopifnot(rel.thresh > 0)
  stopifnot(rel.thresh <= 1)
  thresh <- rel.thresh*N
  if (calc.err) {
    max.thresh <- N
  } else {
    max.thresh <- thresh
  }
  
  if (type == "brownresnick") {
    stopifnot(is.function(vario))
    stopifnot(length(vario(coord[1,]))==1)  
    trend <- sapply(1:N, function(i) sapply(1:N, function(j) vario(coord[i,]-coord[j,])))
    cov.mat <- sapply(1:N, function(i) sapply(1:N, function(j) 
                                       vario(coord[i,]) + vario(coord[j,]) - vario(coord[i,]-coord[j,])))
    cov.mat <- cov.mat + 1e-8 #add constant random effect to avoid numerical problems            
    chol.mat <- chol(cov.mat)
  } else if (type == "extremalt") {
    stopifnot(is.function(corr))
    stopifnot(length(corr(coord[1,]))==1)
    stopifnot(abs(corr(rep(0,times=d))-1)<1e-8)
    stopifnot(dof>1e-12)
    diff.vector <- if (d==1) {
                     matrix(as.vector(outer(coord[,1],coord[,1],'-')),ncol=1)  
                   } else if (d==2) {
                     cbind(as.vector(outer(coord[,1],coord[,1],'-')),
                           as.vector(outer(coord[,2],coord[,2],'-')))  
                   }
    cov.mat.tmp <- matrix(apply(diff.vector, 1, function(x) corr(x)), ncol=N)   
    chol.mat <- vector("list", length=N)
    mu <- vector("list", length=N)
    for (k in 1:N) {
      cov.vec  <- apply(coord, 1, function(x) corr(x-coord[k,]))
      cov.mat  <- (cov.mat.tmp - outer(cov.vec, cov.vec, '*'))/(dof+1) + 1e-6
      chol.mat[[k]] <- chol(cov.mat)
      mu[[k]]       <- apply(coord, 1, function(x) corr(coord[k,]-x))   
    }    
  } else {
    stopifnot(FALSE)
  }
  
  ## initializations 
  res      <- matrix(0, nrow=no.simu, ncol=N) ## results
  res.full <- matrix(0, nrow=no.simu, ncol=N)
  poi <- rexp(no.simu)                 ## poisson points
  time.count      <- rep(0,no.simu)    ## number of Gaussian processes
  time.count.full <- rep(0,no.simu)
  
  while (any(1/poi * max.thresh > apply(res.full,1,min)))  {
    ind.full <- which(1/poi * max.thresh > apply(res,1,min))  
    ## indices of processes where potential updates might take place
    time.count.full[ind.full] <- time.count.full[ind.full] + 1      
    N.ind <- length(ind.full)
    cat(N.ind, " ")
    S <- sample(1:N,N.ind,replace=T,prob=rep(1/N,times=N))
    if (type == "brownresnick") {
      W <- simu_px_brownresnick(no.simu=N.ind, idx=S, N=N, trend=trend, 
                                chol.mat=chol.mat)
    } else if (type == "extremalt") {
      W <- simu_px_extremalt(no.simu=N.ind, idx=S, N=N, dof=dof, mu=mu, 
                             chol.mat=chol.mat)
    } else {
      stopifnot(FALSE)
    }
    stopifnot(nrow(W)==N.ind)
    W <- W/rowSums(W)    
    res.full[ind.full,] <- pmax(1/poi[ind.full]*N*W, 
                                res.full[ind.full,]) ## update proposal
    ind <- which(1/poi * thresh > apply(res,1,min))
    if (length(ind) > 0) {
      res[ind,] <- res.full[ind,]
      time.count[ind] <- time.count.full[ind] 
    }
    poi[ind.full] <- poi[ind.full] + rexp(N.ind) 
  }
  cat("\n")
  
  if (calc.err) {
    err <- apply(res.full - res, 1, function(x) any(x > 1e-8))
    return(list(res=res, time.count=time.count, err=err,
                res.full=res.full, time.count.full= time.count.full))  
  } else {
    return(list(res=res, time.count=time.count))
  }
}

###########################################################################
############ Extremal functions approach ##################################
###########################################################################

simu_extrfcts <- function(no.simu=1, coord, vario, corr, dof, type, 
                          rel.subgrid=1, calc.err=TRUE) {
 
  cat("EF simulation...\n")
  stopifnot((no.simu==round(no.simu)) & (no.simu>=1))
  
  if (!is.matrix(coord)) {
    coord <- matrix(coord, ncol=1)
  }
  d <- ncol(coord)    
  N <- nrow(coord)
  
  stopifnot(rel.subgrid > 0)
  stopifnot(rel.subgrid <= 1)
  n.sub <- round(rel.subgrid*N)
  if (calc.err) {
    max.N <- N
  } else {
    max.N <- n.sub
  }
  subindices <- round(seq(1, N, length=n.sub))  
  
  stopifnot(type %in% c("brownresnick","extremalt"))
  if (type == "brownresnick") {
    stopifnot(is.function(vario))
    stopifnot(length(vario(coord[1,]))==1)
    cov.mat <- sapply(1:N, function(i) sapply(1:N, function(j) 
      vario(coord[i,]) + vario(coord[j,]) - vario(coord[i,]-coord[j,])))
    cov.mat <- cov.mat + 1e-8 #add constant random effect to avoid numerical problems            
    chol.mat <- chol(cov.mat)
  } else if (type == "extremalt") {
    stopifnot(is.function(corr))
    stopifnot(length(corr(coord[1,]))==1)
    stopifnot(abs(corr(rep(0,times=d))-1)<1e-8)
    stopifnot(dof>1e-12)
    diff.vector <- if (d==1) {
      matrix(as.vector(outer(coord[,1],coord[,1],'-')),ncol=1)  
    } else if (d==2) {
      cbind(as.vector(outer(coord[,1],coord[,1],'-')),
            as.vector(outer(coord[,2],coord[,2],'-')))  
    }
    cov.mat.tmp <- matrix(apply(diff.vector, 1, function(x) corr(x)), ncol=N)
  } else {
    stopifnot(FALSE)
  }
  
  ## initializations 
  res      <- matrix(0, nrow=no.simu, ncol=N) ## results
  res.full <- matrix(0, nrow=no.simu, ncol=N)
  poi <- rexp(no.simu)                 ## poisson points
  time.count      <- rep(0,no.simu)    ## number of Gaussian processes
  time.count.full <- rep(0,no.simu)
     
  res <- matrix(0, nrow=no.simu, ncol=N)
  counter <- rep(0, times=no.simu)
  extrfcts.counter <- rep(0, times=no.simu)
  order <- c(subindices,(1:N)[-subindices])
  stopifnot(length(unique(order))==N)
  order <- matrix(order, nrow=no.simu, ncol=N, byrow=TRUE)
 
  for (k in 1:max.N) {
     cat(k, " ")
     stopifnot(all(order[,k]==order[1,k]))
     stopifnot(all(!is.na(order[,k])))
     poisson <- rexp(no.simu)
     if (type == "brownresnick") {
       trend <- sapply(1:N, function(j) vario(coord[j,]-coord[order[1,k],]))
     } else if (type == "extremalt") {
       cov.vec  <- apply(coord, 1, function(x) corr(x-coord[order[1,k],]))
       cov.mat  <- (cov.mat.tmp - outer(cov.vec, cov.vec, '*'))/(dof+1) + 1e-6
       chol.mat <- chol(cov.mat)
       mu <- apply(coord, 1, function(x) corr(coord[order[1,k],]-x))    
     } else {
       stopifnot(FALSE)
     }
     bound <- sapply(1:no.simu, function(i) res.full[i,order[i,k]])
     while (any(1/poisson > bound)) {
       ind <- (1/poisson > bound)
       n.ind <- sum(ind)
       idx <- (1:no.simu)[ind]
       time.count.full[ind] <- time.count.full[ind] + 1
       if (type == "brownresnick") {       
         proc <- simu_px_brownresnick(no.simu=n.ind, idx=order[idx,k], N=N, trend=trend, chol.mat=chol.mat)
       } else if (type == "extremalt") {
         proc <- simu_px_extremalt(no.simu=n.ind, idx=order[1,k], N=N, dof=dof, mu=mu, chol.mat=chol.mat)
       } else {
         stopifnot(FALSE)
       }   
       stopifnot(dim(proc)==c(n.ind, N))
       if (k==1) {
         ind.upd <- rep(TRUE, times=n.ind)
       } else {
         ind.upd <- sapply(1:n.ind, function(i) 
                                    all(1/poisson[idx[i]] * proc[i,order[idx[i],1:(k-1)]] <= res.full[idx[i],order[idx[i],1:(k-1)]]))
       }
       if (any(ind.upd)) {
         idx.upd <- idx[ind.upd]
         res.full[idx.upd,] <- pmax(res.full[idx.upd,], 
                                    1/poisson[idx.upd]*proc[ind.upd,])
       }
       poisson[ind] <- poisson[ind] + rexp(n.ind)
       bound <- sapply(1:no.simu, function(i) res.full[i,order[i,k]])
    } 
    if (k==n.sub) { 
       time.count <- time.count.full
       res <- res.full
    }
  }
  cat("\n")
  
  cat("\n")
  
  if (calc.err) {
    err <- apply(res.full - res, 1, function(x) any(x > 1e-8))
    return(list(res=res, time.count=time.count, err=err,
                res.full=res.full, time.count.full= time.count.full))  
  } else {
    return(list(res=res, time.count=time.count))
  }   
}

############################################################################################
################## sup-normalized (via rejection sampling) #################################
############################################################################################

simu_supnormal <- function(no.simu=1, coord, vario, corr, dof, type, 
                           rel.thresh=1, calc.err=TRUE) { 
    
    cat("sup-normalized simulation...\n")
    if (!is.matrix(coord)) {
      coord <- matrix(coord, ncol=1)
    }
    d <- ncol(coord)    
    N <- nrow(coord)
    
    stopifnot(rel.thresh > 0)
    stopifnot(rel.thresh <= 1)
    if (calc.err) {
      max.thresh <- 1
    } else {
      max.thresh <- rel.thresh
    }
    
    stopifnot(type %in% c("brownresnick","extremalt"))
    if (type == "brownresnick") {
      stopifnot(is.function(vario))
      stopifnot(length(vario(coord[1,]))==1)  
      trend <- sapply(1:N, function(i) sapply(1:N, function(j) vario(coord[i,]-coord[j,])))
      cov.mat <- sapply(1:N, function(i) sapply(1:N, function(j) 
        vario(coord[i,]) + vario(coord[j,]) - vario(coord[i,]-coord[j,])))
      cov.mat <- cov.mat + 1e-8 #add constant random effect to avoid numerical problems            
      chol.mat <- chol(cov.mat)
    } else if (type == "extremalt") {
      stopifnot(is.function(corr))
      stopifnot(length(corr(coord[1,]))==1)
      stopifnot(abs(corr(rep(0,times=d))-1)<1e-8)
      stopifnot(dof>1e-12)
      diff.vector <- if (d==1) {
        matrix(as.vector(outer(coord[,1],coord[,1],'-')),ncol=1)  
      } else if (d==2) {
        cbind(as.vector(outer(coord[,1],coord[,1],'-')),
              as.vector(outer(coord[,2],coord[,2],'-')))  
      }
      cov.mat.tmp <- matrix(apply(diff.vector, 1, function(x) corr(x)), ncol=N)   
      chol.mat <- vector("list", length=N)
      mu <- vector("list", length=N)
      for (k in 1:N) {
        cov.vec  <- apply(coord, 1, function(x) corr(x-coord[k,]))
        cov.mat  <- (cov.mat.tmp - outer(cov.vec, cov.vec, '*'))/(dof+1) + 1e-6
        chol.mat[[k]] <- chol(cov.mat)
        mu[[k]]       <- apply(coord, 1, function(x) corr(coord[k,]-x))   
      }    
    } else {
      stopifnot(FALSE)
    }

    ## initializations 
    res <- matrix(0, nrow=no.simu, ncol=N) ## results
    res.full <- matrix(0, nrow=no.simu, ncol=N)
    poi <- rexp(no.simu)                 ## poisson points
    time.count      <- rep(0,no.simu)    ## number of spectral processes
    time.count.full <- rep(0,no.simu)
    Gauss.time.count      <- rep(0,no.simu) ## number of Gaussian processes 
    Gauss.time.count.full <- rep(0,no.simu)
    rej.counter <- 0
    extr.coeff  <- 0
        
    while (any(1/poi * max.thresh > apply(res.full,1,min)))  {
      ind.full <- which(1/poi * max.thresh > apply(res.full,1,min))  ## indices of processes where potential updates might take place according to thresh.slow
      time.count.full[ind.full] <- time.count.full[ind.full] + 1      
      N.ind <- length(ind.full)
      cat(N.ind, "")
      W <- matrix(NA, nrow=N.ind, ncol=N)
      while (any(is.na(W))) {
        prop.ind <- which(is.na(W[,1]))
        S <- sample(1:N, length(prop.ind),replace=T,prob=rep(1/N,times=N))
        if (type == "brownresnick") {
          W.prop <- simu_px_brownresnick(no.simu=length(prop.ind), idx=S, N=N, trend=trend, chol.mat=chol.mat)
        } else if (type == "extremalt") {
          W.prop <- simu_px_extremalt(no.simu=length(prop.ind), idx=S, N=N, dof=dof, mu=mu, chol.mat=chol.mat)
        } else {
          stopifnot(FALSE)
        }
        Gauss.time.count.full[ind.full][prop.ind] <- 
               Gauss.time.count.full[ind.full][prop.ind] + 1
        extr.coeff <- extr.coeff + sum(apply(W.prop, 1, max))
        acc.probab <- apply(W.prop,1,max)/apply(W.prop,1,sum)
        coins <- runif(length(prop.ind))
        W[prop.ind[coins < acc.probab],] <- W.prop[coins < acc.probab,]           
      }
      W <- W/apply(W, 1, max)     
      res.full[ind.full,] <- pmax(1/poi[ind.full]*W, 
                                  res.full[ind.full,])  ## update proposal
      ind <- which(1/poi * rel.thresh > apply(res,1,min))
      if (length(ind) > 0) {
        res[ind,] <- res.full[ind,]
        time.count[ind] <- time.count.full[ind] 
        Gauss.time.count[ind] <- Gauss.time.count.full[ind]
      }
      poi[ind.full] <- poi[ind.full] + rexp(N.ind) 
    }
    cat("\n")
    
    extr.coeff  <- extr.coeff / sum(Gauss.time.count.full)
    res.full <- res.full * extr.coeff
    res      <- res      * extr.coeff
    
    if (calc.err) {
      err <- apply(res.full - res, 1, function(x) any(x > 1e-8))
      return(list(res=res, spec.proc.count=time.count, 
                  Gauss.proc.count=Gauss.time.count, err=err,
                  res.full=res.full, 
                  spec.proc.count.full= time.count.full,
                  Gauss.proc.count.full=Gauss.time.count.full))  
    } else {
      return(list(res=res, spec.proc.count=time.count, 
                  Gauss.proc.count=Gauss.time.count))
    }
}
