ruleF <- function(alpha, beta, p0, p1, criteria){
  
  Q0 = QR_fun(p0,p0)
  R0 = QR_fun(p0,p0)
  Q1 = QR_fun(p0,p1)
  R1 = QR_fun(p1,p0)
  D1 = P_plus(p0,p1)-P_minus(p0,p1)
  
  ### optimal design
  alpha1 <- seq(alpha, 0.99, by=0.02)
  beta1 <- seq(beta/2, beta, by=0.01)
  lalpha1 <- length(alpha1)
  lbeta1 <- length(beta1)
  EN0 <- matrix(rep(0, lalpha1*lbeta1), ncol=lbeta1)
  
  for(i in 1:lalpha1){
    for (j in 1:lbeta1){
      za1 <- qnorm(alpha1[i],lower.tail = FALSE)
      zb1 <- qnorm(beta1[j],lower.tail = FALSE)
      za <- qnorm(alpha,lower.tail = FALSE)
      zb2 <- qnorm(beta-beta1[j],lower.tail = FALSE)
      n1 <- (za1*sqrt(Q0+R0)+zb1*sqrt(Q1+R1))^2/D1^2
      t1 <- za1*sqrt((Q0+R0)/n1)
      n2 <- (za*sqrt(Q0+R0)+zb2*sqrt(Q1+R1))^2/D1^2
      t2 <- za*sqrt((Q0+R0)/n2)
      
      if(criteria==1){
        EN0[i,j] <- 2*n1 + 2*(n2-n1)*alpha1[i]
      }else if(criteria==2){
        EN0[i,j] <- 2*n1 + 2*(n2-n1)*(1-beta1[j])
      }else if(criteria==3){
        EN0[i,j] <- 2*n1 + (n2-n1)*(alpha1[i]+1-beta1[j])
      }
      
    }
  }
  
  ind <- which(EN0 == min(EN0, na.rm=T), arr.ind = TRUE)
  
  i=ind[1,1]; j=ind[1,2]
  za1 <- qnorm(alpha1[i],lower.tail = FALSE)
  zb1 <- qnorm(beta1[j],lower.tail = FALSE)
  za <- qnorm(alpha,lower.tail = FALSE)
  zb2 <- qnorm(beta-beta1[j],lower.tail = FALSE)
  n1 <- (za1*sqrt(Q0+R0)+zb1*sqrt(Q1+R1))^2/D1^2
  t1 <- za1*sqrt((Q0+R0)/n1)
  n2 <- (za*sqrt(Q0+R0)+zb2*sqrt(Q1+R1))^2/D1^2
  t2 <- za*sqrt((Q0+R0)/n2)
  
  res <-list(t1, n1, t2, n2, alpha1[i], beta1[j])
  names(res) <- c("t1", "n1", "t2", "n2", "alpha1", "beta1")
  return(res)
  
}
