n <- 25, n <- 100, n <- 1000         # number of observations
p <- 0.7
y <- rbinom(n,1,p)                   # data

N <- 1000                            # number of values of 
theta <- numeric(N)                  # simulated p

p <- function(zeta,y){

      abs( sum(zeta)-sum(y) ) / n

   }

e_0 <- 0                             # tolerance threshold
for (i in 1:N) {

    d <- e_0+1
    while ( d>e_0 ) {

         th <- rbeta(1,1,1)
         zeta <- rbinom(n,1,th)
         d <- p(zeta,y)
    }
    theta[i]<-th
}