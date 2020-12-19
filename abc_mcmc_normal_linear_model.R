library('bayesianRegression')
library('mvtnorm')
library('MASS')

sim_data <- simulate_data(n = 50, seed = 1,true_beta = c(1,0,-2,3), true_sigma2 = 3)

n <- dim(sim_data)[1]   
p <- dim(sim_data)[2]
constant <- rep(1, length = n)
X <- cbind(constant, sim_data[,-1])
V <- n*solve(t(X)%*%X)
     
y_mean <- function(b,Y) {
     # calculate the mean value of data
     b <- matrix(t(b), nrow = p, byrow = FALSE)
     return(Y%*%b)
} 

y <- sim_data[,1]                                                   # distance 
suf_stat <- c(t(X)%*%y/n, var(y))
p1 <- function(zeta,y) {    
	suf_stat_sim <- c(t(X)%*%zeta/n, var(zeta))
	return(dist(rbind(suf_stat, suf_stat_sim)))
} 

a <- b <- 0.1                                              # hyperparameters

# set initial values by least squares estimates
fit <- lm(y~. , data = as.data.frame(sim_data))
sf <- summary(fit)
chain <- c(coefficients(fit),(sf$sigma)^2)
 
mcmcIterations <- 10000000
iter <- 0
acceptRate <- 0
proposal <- numeric(p+1)
thin <- 10000	                                       
savedTheta <- matrix(NA, nrow = mcmcIterations/thin, ncol = p+1)

randomWalkSD <- 0.1	# sd of random walk proposal
e_0 <- 0.5		# distance threshold

par(mfrow = c(2, p+1), mar =c(2,2,3,1))
for(i in 2:mcmcIterations){
	
	proposal[1:p] <- rnorm(p, mean = chain[1:p], sd= rep(randomWalkSD,p) )
	proposal[p+1] <- rlnorm(1, meanlog = log(chain[p+1]), sdlog = randomWalkSD)
        zeta <- mvrnorm(n = 1, y_mean(proposal[1:p],X), proposal[p+1]*diag(n), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

	if ( p1(zeta, y) < e_0 ){
		logPriorRatio <- dmvnorm(proposal[1:p], mean = rep(0, p), sigma = V*proposal[p+1], log = TRUE) - 
					dmvnorm(chain[1:p], mean = rep(0, p), sigma = V*chain[p+1], log = TRUE) + 
					(a+1)*(log(chain[p+1]) - log(proposal[p+1])) - b*(1/proposal[p+1] - 1/chain[p+1])	
		logProposalRatio <- log(proposal[p+1]) - log(chain[p+1])
		logAcceptanceRatio <- logPriorRatio + logProposalRatio
		if(log(runif(1)) < logAcceptanceRatio){
			chain <- proposal
			acceptRate <- acceptRate+1
		}
	}
	if((i %% thin) == 0){
		iter <- iter + 1
		savedTheta[iter, ] <- chain
		cat(paste0('[iteration ',i, '    ] acceptance rate of ABC-MCMC sampler: ', 100*acceptRate/i, '%.','\n'))
		if(iter > 1){
			for(j in 1:p){
				plot(savedTheta[1:iter,j],type='l', main = bquote('trace '*beta[.(j-1)]), xlab='', ylab='')
			}
				plot(savedTheta[1:iter,p+1],type='l', main = bquote('trace '*sigma^2), xlab='', ylab='')
			for(j in 1:p){
				hist(savedTheta[1:iter,j], main = bquote('hist '*beta[.(j-1)]), xlab='', ylab='')
			}
				hist(savedTheta[1:iter,p+1], main = bquote('hist '*sigma^2), xlab='', ylab='')
		}
	} 	

}