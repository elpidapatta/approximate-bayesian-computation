n <- 100 ; n <- 1000 
p <- 0.7 ; y <- rbinom(n,1,p)

nsim <- 2000                                 # store the simulated values
results <- matrix(ncol=2,nrow=nsim)          # and the corresponding model

i <- 1 ; e_0 <- 0

while (i <= nsim){

     choose.model <- sample( c('model1','model2'),1 )
     if (choose.model == 'model1'){
             th <- 0.6
     }
     else{

     th <- rbeta(1,1,1)
     }

     zeta <- rbinom(n,1,th)
     d <- p(zeta,y)
     if (d <= e_0){
        results[i,1] <- th
        results[i,2] <- choose.model
        i <- i+1
     }
}

