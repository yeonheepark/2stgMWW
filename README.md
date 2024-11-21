# 2stgMWW

Still in development. References will be posted when published.

## Example
alpha <- 0.05

beta <- 0.2

p0 <- c(1/3, 1/3, 1/3)

p1 <- c(1/2, 1/3, 1/6)

ruleF(alpha, beta, p0, p1, criteria=1)

t1 <- 0.060;
n1 <- 32;
t2 <- 0.125;
n2 <- 104;

set.seed(1234)

opMWW.F(alpha, beta, p0, p0, t1, n1, t2, n2, nsim=10000)

opMWW.F(alpha, beta, p0, p1, t1, n1, t2, n2, nsim=10000)
