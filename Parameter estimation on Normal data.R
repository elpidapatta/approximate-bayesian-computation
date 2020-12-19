n <- 100 ; n <- 500                  # number of observations

real_data <- rnorm(n,10,2)           # data

N <- 1000                            # number of values of simulated mean 
results <- matrix(nrow=N,ncol=3)     # and variance
sortedy <- sort(y)                                          

p1 <- function(zeta,y) {   
      zeta <- sort(zeta)
      sqrt( sum( (  (zeta-y)/ sd(zeta) )^2 ) ) / n
  } 

e_0 <- 0.02 ; e_0 <- 0.01            # tolerance threshold

for ( i in 1:N ) {

     d <- e_0+1
     while ( d>e_0 ){
        var_sim <- 1/rgamma(1,shape=3,rate=3)
        mean_sim <- rnorm(1,5,sqrt(var_sim))
        simulated_data <- rnorm(n,mean_sim,sqrt(var_sim))
        d <- p1(zeta,sortedy)
     }
     results[i,1] <- mean_sim
     results[i,2] <- var_sim
     results[i,3] <- d
 
 }
