set.seed(99)
sim_data <- simulate_data(n = 100, seed = 1,true_beta = c(1,0,-2,3), true_sigma2 = 3)

n <- dim(sim_data)[1]  
p <- dim(sim_data)[2]
constant <- rep(1, length = n)
X <- cbind(constant, sim_data[,-1])  
V <- n*solve(t(X)%*%X)

y <- sim_data[,1]
p1 <- function(zeta, y, activeCols) {    
	suf_stat <- c(t(X[,activeCols])%*%y/n, var(y))
	suf_stat_sim <- c(t(X[,activeCols])%*%zeta/n, var(zeta))
	return(dist(rbind(suf_stat, suf_stat_sim)))
} 

a <- b <- 0.1
fit <- lm(y~. , data = as.data.frame(sim_data))
sf <- summary(fit)
chain <- c(coefficients(fit),(sf$sigma)^2)
s <- rep(1, p)	# full model
 
mcmcIterations <- 10000000   
thin <- 10000
iter <- 0
acceptRate <- 0  
acceptRateRJMOVE <- 0
proposal <- NA*numeric(p+1)
	
savedTheta <- matrix(NA, nrow = mcmcIterations/thin, ncol = p+1)
savedModel <- matrix(NA, nrow = mcmcIterations/thin, ncol = p)

randomWalkSD <- 0.3   
epsilon <- 0.75	 
betaProposalSD <- 2	
betaBirth <- function(sd){	
	beta <- rnorm(1, mean=0,sd=sd)
	log_pdf <- dnorm(beta, mean=0, sd=sd, log=TRUE)
	return(c(beta, log_pdf))
}

birthProb <- function(n1){   
            
	if(n1 < 1){stop('n1 cannot be < 0')}
	if(n1 > p){stop(paste0('n1 cannot be > ',p))}
	if(n1 == 1){
		return(1)
	}else{
		if(n1 < p){
			return(0.5)
		}else{
			return(0)
		}
	}
}


for(i in 2:mcmcIterations){

   n1 <- sum(s)
   indexOf1s <- which(s==1)
   indexOf0s <- (1:p)[-indexOf1s]  
   if(runif(1) < birthProb(n1)){             #	BIRTH MOVE

	if(p-n1>1){
		h <- sample(indexOf0s, 1)    
	}else{
		h <- indexOf0s
	}
	indexOf0sProposed <- indexOf0s[-which(indexOf0s==h)]   
	indexOf1sProposed <- sort(c(h, indexOf1s))             
	newBeta <- betaBirth(betaProposalSD)		
	chainProposed <- chain                                 
	chainProposed[h] <- newBeta[1]                         
       	zeta <- mvrnorm(n = 1, y_mean(chainProposed[indexOf1sProposed],X[,indexOf1sProposed]), chain[p+1]*diag(n),tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

	if(p1(zeta, y) < epsilon){

		priorMeanBetaProposed <- rep(0, n1 + 1)            
		priorMeanBeta <- rep(0, n1)                         
		priorVarBetaProposed <- n*solve(t(X[,indexOf1sProposed]) %*% X[,indexOf1sProposed])*chain[p+1]
		priorVarBeta <- n*solve(t(X[,indexOf1s]) %*% X[,indexOf1s])*chain[p+1]
		jacobian <- 0                                                                             
		priorRatio <- dmvnorm(chainProposed[indexOf1sProposed], priorMeanBetaProposed, priorVarBetaProposed, log=T) - 
				dmvnorm(chain[indexOf1s], priorMeanBeta, priorVarBeta, log=T)
		proposalRatio <- log(1 - birthProb(n1+1)) - log(birthProb(n1)) + log(p - n1) - log( n1 + 1 ) - newBeta[2] 	
		log_acceptance_ratio <- priorRatio + proposalRatio + jacobian

		if(log(runif(1)) < log_acceptance_ratio){
			s[h] <- 1
			chain[h] <- newBeta[1]
			acceptRateRJMOVE <- acceptRateRJMOVE + 1
		}
	}

   }else{                  # DEATH MOVE

	      h <- sample(indexOf1s, 1)
	      indexOf1sProposed <- indexOf1s[-which(indexOf1s==h)]
	      deletedBeta <- chain[h]
	      chainProposed <- chain[c(indexOf1sProposed,p+1)]
      	      zeta <- mvrnorm(n = 1, y_mean(chain[indexOf1sProposed],X[,indexOf1sProposed]), chain[p+1]*diag(n),tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

	      if(p1(zeta, y) < epsilon){
		priorMeanBetaProposed <- rep(0, n1 - 1)
		priorMeanBeta <- rep(0, n1)
		priorVarBetaProposed <- n*solve(t(X[,indexOf1sProposed]) %*% X[,indexOf1sProposed])*chain[p+1]
		priorVarBeta <- n*solve(t(X[,indexOf1s]) %*% X[,indexOf1s])*chain[p+1]
		priorRatio <- dmvnorm(chain[indexOf1sProposed], priorMeanBetaProposed, priorVarBetaProposed, log=T) - dmvnorm(chain[indexOf1s], priorMeanBeta, priorVarBeta, log=T)
		proposalRatio <- log(birthProb(n1-1)) - log(1-birthProb(n1)) + log(n1) - log( p + n1 - 1 ) + dnorm(deletedBeta, mean = 0, sd = betaProposalSD, log=TRUE)
			
		log_acceptance_ratio <- priorRatio + proposalRatio + jacobian
		if(log(runif(1)) < log_acceptance_ratio){
			s[h] <- 0
			chain[h] <- NA
			acceptRateRJMOVE <- acceptRateRJMOVE + 1
		}
	     }
   }


   n1 <- sum(s)
   indexOf1s <- which(s==1)
   indexOf0s <- (1:p)[-indexOf1s]

   proposal <- NA*numeric(p+1)
   proposal[indexOf1s] <- rnorm(n1, mean = chain[indexOf1s], sd= rep(randomWalkSD, n1) )
   proposal[p+1] <- rlnorm(1, meanlog = log(chain[p+1]), sdlog = randomWalkSD)
   zeta <- mvrnorm(n = 1, y_mean(proposal[indexOf1s],X[,indexOf1s]), proposal[p+1]*diag(n),tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

	if ( p1(zeta, y) < epsilon ){
		V <- n*solve(t(X[,indexOf1s]) %*% X[,indexOf1s])
		SigmaProp <- V*proposal[p+1]
		Sigma <- V*chain[p+1]
		logPriorRatio <- dmvnorm(proposal[indexOf1s], mean = rep(0, n1), sigma = SigmaProp, log = TRUE) - 
					 dmvnorm(chain[indexOf1s], mean = rep(0, n1), sigma = Sigma, log = TRUE) + 
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
		savedModel[iter, ] <- s
                cat(paste0('[iteration ',i, '    ] acceptance rate of ABC-MCMC sampler: ', 100*acceptRate/i, '% (parameters updating), ', 100*acceptRateRJMOVE/i, '% (model updating).','\n'))
		if(iter > 1){
			barplot(colMeans(savedModel[1:iter, ]), ylab = 'inclusion probability', xlab = 'variable', names.arg = colnames(X))

		}
	} 	

}
	