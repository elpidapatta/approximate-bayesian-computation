n <- 25 ; n <- 100 ; n <- 1000             # number of observations
y <- rnorm(n,0,1)                          # data

N <- 1000
m <- numeric(N)                            # number of values of
                                           # simulated mean
p1 <- function(zeta,y){

    abs( mean(zeta)-mean(y) )

  }

e_0 <- 0.005 ; e_0 <- 0.001                # tolerance threshold

for (i in 1:N) {

    d <- e_0 + 1
    while ( d>e_0 ) {
         theta <- rnorm(1,0,10)
         zeta <- rnorm(n,theta,1)
         d <- p1(zeta,y)
    }

    m[i] <- theta
}