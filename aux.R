P_plus <- function(p1, p2){
  # This is the function to compute P+ in the Mann-Whitney-Wilcoxon test, 
  # Input:
  # p1: a vector contains the probabilities of an outcome falling in each level
  # p2: a vector contains the probabilities of an outcome falling in each level for the second treatment
  
  # Note: the length of p1 and p2 should be equal.
  
  # Output:
  # The value of P+
  
  k = length(p1)
  
  p1_rev = rev(p1[-1])
  p1_sum = rev(cumsum(p1_rev))
  
  p2_short = p2[-k]
  res = sum(p2_short*p1_sum)
  
  return(res)
}

P_minus <- function(p1, p2){
  # This is the function to compute P- in the Mann-Whitney-Wilcoxon test, 
  # Input:
  # p1: a vector contains the probabilities of an outcome falling in each level
  # p2: a vector contains the probabilities of an outcome falling in each level for the second treatment
  
  # Note: the length of p1 and p2 should be equal.
  
  # Output:
  # The value of P-
  
  k = length(p1)
  
  p2_rev = rev(p2[-1])
  p2_sum = rev(cumsum(p2_rev))
  
  p1_short = p1[-k]
  res = sum(p1_short*p2_sum)
  
  return(res)
}


QR_fun <- function(p1, p2){
  # This is the function to compute Q or R in the Mann-Whitney-Wilcoxon test, 
  # Input:
  # p1: a vector contains the probabilities of an outcome falling in each level
  # p2: a vector contains the probabilities of an outcome falling in each level for the second treatment
  
  # Note: the length of p1 and p2 should be equal. Q = QR_fun(p1,p2), R = QR_fun(p2,p1)
  
  # Output:
  # The value of Q
  
  k = length(p1)
  
  DD2 = (P_plus(p1,p2) - P_minus(p1,p2))^2
  
  # compute the first part
  p2_rev = rev(p2[-1])
  p2_sum = rev(cumsum(p2_rev))
  
  p1_short = p1[-k]
  part1 = sum(p1_short*(p2_sum^2))
  
  # compute the second part
  p1_short2 = p1[-1]
  p2_short = p2[-k]
  
  p2_sum2 = cumsum(p2_short)
  part2 = sum(p1_short2*(p2_sum2^2))
  
  # compute the third part
  p1_short3 =  p1_short[-1]
  p2_short2 =  p2_short[-(k-1)]
  p2_rev2 = rev(p2[-c(1,2)])
  
  p2_sum3 = cumsum(p2_short2)
  p2_sum4 = rev(cumsum(p2_rev2))
  part3 = 2*sum(p1_short3*p2_sum3*p2_sum4)
  
  res = part1+part2-part3-DD2
  
  return(res)
}
