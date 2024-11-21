opMWW.FS <- function(alpha, beta, p0, p1, v1, v2, n1, v3, n2, nsim){
  out1=pts <- c()
  for(sim in 1:nsim){
    n.interim <- c(2*n1, 2*n2)
    J <- length(p1)
    ne <- sample(1:J, size=n2, replace=TRUE, prob=p0)
    nestg1 <- ne[sample(1:n2, size=n1, replace=F)]
    
    group11 <- sum(nestg1 == 1)
    group12 <-  sum(ne == 1)
    for(j in 2:J){
      aux1 <- sum(nestg1 == j)
      group11 <- c(group11, aux1)
      aux2 <- sum(ne == j)
      group12 <- c(group12, aux2)
    }
    
    nc <- sample(1:J, size=n2, replace=TRUE, prob=p1)
    ncstg1 <-nc[sample(1:n2, size=n1, replace=F)] 
    
    group21 <- sum(ncstg1 == 1)
    group22 <-  sum(nc == 1)
    for(j in 2:J){
      aux1 <- sum(ncstg1 == j)
      group21 <- c(group21, aux1)
      aux2 <- sum(nc == j)
      group22 <- c(group22, aux2)
    }
    
    num <- sum(group11[-1]*cumsum(group21[-J])) - sum(group21[-1]*cumsum(group11[-J]))
    denom <- sum(group11)*sum(group21)
    wk <- num/denom 
    ndotjk <- group11+group21
    wkdenom <- (sum(ndotjk)+1)*(1-sum(ndotjk^3-ndotjk)/((sum(ndotjk))^3-sum(ndotjk)))/(3*sum(group11)*sum(group21))
    W1 <- wk
    
    num <- sum(group12[-1]*cumsum(group22[-J])) -sum(group22[-1]*cumsum(group12[-J]))
    denom <- sum(group12)*sum(group22)
    wk <- num/denom 
    ndotjk <- group12+group22
    wkdenom <- (sum(ndotjk)+1)*(1-sum(ndotjk^3-ndotjk)/((sum(ndotjk))^3-sum(ndotjk)))/(3*sum(group12)*sum(group22))
    W2 <- wk
    
    if(W1 <= v1){
      out1[sim] <-  "Early Stop"
      pts[sim] <- n.interim[1]
    }else if(W1>v2){
      out1[sim] <-  "Reject all"
      pts[sim] <- n.interim[1]
    }else{
      if(W2>v3){
        out1[sim] <- "Reject all"
        pts[sim] <- n.interim[2]
      }else{
        out1[sim] <- "Fail stage 2"
        pts[sim] <- n.interim[2]
      }
    }
    
  }
  phat <- length(which(out1=="Reject all"))/nsim
  mpts <- mean(pts)
  res <-list(phat, mpts)
  names(res) <- c("phat", "mpts")
  return(res)
  
}
