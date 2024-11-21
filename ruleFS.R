ruleFS <- function(alpha, beta, p0, p1, criteria){
  
  Q0 = QR_fun(p0,p0)
  R0 = QR_fun(p0,p0)
  Q1 = QR_fun(p0,p1)
  R1 = QR_fun(p1,p0)
  D1 = P_plus(p0,p1)-P_minus(p0,p1)
  
  min_result <- Inf
  best_params <- c(0, 0, 0)
  out_n1=out_n2=out_v1=out_v2=out_v3 <- 0
  alpha1values <- seq(alpha, 0.999, by=0.001)
  beta1values <- seq(0.01, beta-0.01, by=0.01)
  
  for (alpha1 in alpha1values){
    for(beta1 in beta1values){
      
      alpha2values <- seq(0.001, alpha-0.001, by=0.001)
      for(alpha2 in alpha2values){
        
        za1.2 <- qnorm(alpha1,lower.tail = FALSE)
        za2.2 <- qnorm(alpha2,lower.tail = FALSE)
        zb1.2 <- qnorm(1-beta1,lower.tail = FALSE)
        
        beta3 <- beta-beta1
        zb2.2 <- qnorm(1-beta3,lower.tail = FALSE)
        za3.2 <- qnorm(alpha-alpha2,lower.tail = FALSE)
        
        n1 <- (za1.2*sqrt(Q0+R0)-zb1.2*sqrt(Q1+R1))^2/D1^2
        v1 <- za1.2*sqrt((Q0+R0)/n1)
        v2 <- za2.2*sqrt((Q0+R0)/n1)
        
        n2 <- (za3.2*sqrt(Q0+R0)-zb2.2*sqrt(Q1+R1))^2/D1^2
        v3 <-  za3.2*sqrt((Q0+R0)/n2)
        
        beta2 <- pnorm(sqrt(n1)*(v2-D1)/sqrt(Q1+R1))
        
        if(criteria==1){
          result <- 2*n1 + 2*(n2-n1)*(alpha1-alpha2)
        }else if(criteria==2){
          result <- 2*n1 + 2*(n2-n1)*(beta2-beta1)
        }else if(criteria==3){
          result <- 2*n1 + (n2-n1)*(alpha1-alpha2+beta2-beta1)
        }
        
        if(n1>n2) result <- Inf
        
        if(result < min_result){
          min_result <- result
          best_params <- c(alpha1, alpha2, beta1)
          out_n1 <- n1
          out_n2 <- n2
          out_v1 <- v1
          out_v2 <- v2
          out_v3 <- v3
          paraout <- c(alpha1, alpha2, beta1, beta2)
        }
      }
    }
  }
  return(list(paraout=paraout, 
              v1=out_v1, v2=out_v2, n1=out_n1,
              v3=out_v3,  n2=out_n2))
  
}
