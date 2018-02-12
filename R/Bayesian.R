#' Auxiliary function for \code{sp_msfa}.
#'
#' Set parameters for posterior sampling and prior hyperparameters for the  Bayesian MSFA model for sparse setting.
#' The notation follows closely the paper by De Vito et al. (2018).
#'
#' @param nrun Number of posterior simulations. Default is 30000.
#' @param burn Burn-in trials. Default is 20000.
#' @param thin Thinning of posterior samples. Default is 1 (no thinning).
#' @param nu Parameter entering the gamma distribution assumed for omega. Default is 3.
#' @param nus Parameter entering the gamma distribution assumed for omegas. Default is 3.
#' @param a1 Shape parameter for the gamma distribution assumed for delta1. Default is 1.1.
#' @param b1 Scale parameter for the gamma distribution assumed for delta1. Default is 1.
#' @param a2 Shape parameter for the gamma distribution assumed for deltal. Default is 1.1.
#' @param b2 Scale parameter for the gamma distribution assumed for deltal. Default is 1.
#' @param a1s Shape parameter for the gamma distribution assumed for delta1_s. Default is 1.1.
#' @param b1s Scale parameter for the gamma distribution assumed for delta1_s. Default is 1.
#' @param a2s Shape parameter for the gamma distribution assumed for deltal_s. Default is 2.1.
#' @param b2s Scale parameter for the gamma distribution assumed for deltal_s. Default is 1.
#' @param apsi Shape parameter for the gamma distribution assumed for 1/psi. Default is 1.
#' @param bpsi Scale parameter for the gamma distribution assumed for 1/psi. Default is 0.3.
#' @export
#' @references De Vito, R., Bellio, R., Trippa, L. and Parmigiani, G. (2018).
#' Bayesian Multi-study Factor Analysis for High-throughput Biological Data.
#' Manuscript in preparation.
sp_msfa_control <- function(nrun = 30000, burn = 20000, thin = 1,
                          nu = 3, nus = 3,
                          a1 = 1.1, b1 = 1,
                          a2 = 1.1, b2 = 1,
                          a1s = 1.1 , b1s = 1,
                          a2s = 2.1, b2s = 1,
                          apsi = 1, bpsi = 0.3)
{
  return(list(nrun = nrun, burn = burn, thin = thin,
              nu = nu, nus = nus, a1 = a1, b1 = b1, a2 = a2, b2 = b2,
              a1s = a1s, b1s = b1s, a2s = a2s, b2s = b2s, apsi = apsi, bpsi = bpsi))
}



#' Posterior sampling for the  Bayesian MSFA model for sparse setting.
#'
#' The function implements the Gibbs sampler described in De Vito et al. (2018). The code
#' is suitable for small to moderate-size data, and therefore can readily reproduce some of the
#' results of the paper, but not those for large data which would require larger computational times.
#' The \code{outputlevel} argument has an important role for practical usage. A value
#' \code{outputlevel = 1} (the default) will save all the MCMC chains, and this would create
#' a rather bulky output. The option \code{outputlevel = 2} will save only the chains for the loading
#' matrix of common factors, whereas  option \code{outputlevel = 3} will not save any chain, reporting
#'in the output also the posterior means of the crossproduct of the loading matrices.
#'
#'
#'
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
#' Standardization is carried out by the function.
#' @param k Number of common factors.
#' @param j_s Number of study-specific factors. A vector of positive integers of length \eqn{S}{S}.
#' @param trace If \code{TRUE} then trace information is being printed every \code{nprint} iterations of the Gibbs sampling.
#' Default is \code{TRUE}.
#' @param nprint Frequency of tracing information. Default is every 1000 iterations.
#' @param outputlevel Detailed level of output data. See Details. Default is 1.
#' @param control A list of hyperparameters for the prior distributions and for controlling the Gibbs sampling.
#' See \code{sp_msfa_control}.
#' @param ... Arguments to be used to form the default \code{control} argument if it is not supplied directly.
#' @return A list  containing the  posterior samples for the model parameters.  If \code{outputlevel = 1}, the components of the list are:
#' \item{\code{Phi}}{Common factor loadings. An array of dimension \eqn{p \times k \times (nrun - burn)/thin}{p x k x (nrun - burn)/thin}.}
#' \item{\code{Lambda}}{Study-specific factor loadings. A list of arrays of dimension \eqn{p \times j_s[s] \times (nrun - code)/thin}{p x j_s[s] x (nrun - burn)/thin}.}
#' \item{\code{psi}}{Study-specific uniquenesses. A list of arrays of dimension  \eqn{p \times 1 \times (nrun - burn)/thin}{p x 1 x (nrun - burn)/thin}.}
#' \item{\code{f_s}}{Study-specific latent factors associated to common factor loadings.
#' A list of arrays of dimension \eqn{nrow(X_s[[s]]) \times k \times (nrun - burn)/thin}{nrow(X_s[[s]]) x k x (nrun - burn)/thin}.}
#' \item{\code{l_s}}{Study-specific latent factors associated to study-specific factor loadings. A list of arrays of dimension
#' \eqn{nrow(X_s[[s]]) \times j_s[s] \times (nrun - burn)/thin}{nrow(X_s[[s]]) x j_s[s] x (nrun - burn)/thin}.}
#' When instead \code{outputlevel > 1}, the arrays are replaced by posterior means. The matrices \code{SigmaPhi} or the list of
#' matrices \code{SigmaLambda}, containing the posterior means of the these quantities, will be returned when \code{outputlevel}
#' is different from 1.
#' @import statmod matlab
#' @importFrom stats rgamma rnorm
#' @export
#' @references De Vito, R., Bellio, R., Trippa, L. and Parmigiani, G. (2018).
#' Bayesian Multi-study Factor Analysis for High-throughput Biological Data.
#' Manuscript in preparation.
sp_msfa <- function(X_s,  k,  j_s, trace = TRUE, nprint = 1000,
                    outputlevel = 1, control = list(...), ...)
{
  #### read data ####
  S <- length(X_s)
  p <- dim(X_s[[1]])[2]                    #variables
  Y_s <- list()                         #centered data
  n_s <- c()                            #sample size

  ################ setting up priors and initialize
  control <- do.call("sp_msfa_control", control)
  nrun <- control$nrun
  thin <- control$thin
  burn <- control$burn
  sp  <- (nrun - burn) / thin
  a_psi_s <- b_psi_s <-  c()    ####gamma hyperparameters for Psi_s
  nu_s <- c()                # gamma hyperpameters for omega_s
  a1_s <- b1_s <- c()      # gamma hyperparameters for delta_1
  a2_s <- b2_s <- c()      #gamma hyperparameters for delta_l^1

  #### priors
  apsi <- control$apsi
  bpsi<- control$bpsi
  nu <- control$nu
  nus <- control$nus
  a1 <- control$a1
  b1 <- control$b1
  a2 <- control$a2
  b2 <-  control$b2
  a1s <- control$a1s
  b1s <- control$b1s
  a2s <- control$a2s
  b2s <-  control$b2s

  ####initial values
  psi_s <- Psi_s <- Lambda_s <- l_s <- list()
  Meta_s <-  Veta_s <- f_s <- MetaP_s <- VetaP_s <- list()

  ####prior setup
  omegajh_s <- delta_s <- tauh_s <- Plam_s <- list()
  #### output

  Lambdaout <- psiout <- l_sout <- f_sout  <- SigmaLambda <-  list()

   for (s in 1:S){
  	n_s[s] <- nrow(X_s[[s]])
  	Y_s[[s]] <- scale(X_s[[s]])
  	a_psi_s[s] <- apsi
  	b_psi_s[s] <- bpsi
  	nu_s[s] <- nus
  	a1_s[s] <- a1s
  	b1_s[s] <- b1s
  	a2_s[s] <- a2s
  	b2_s[s] <- b2s

  	#### initial values ####
  	###for covariance error terms
  	psi_s[[s]] <- rgamma(p, shape = a_psi_s[s], scale = 1 / b_psi_s[s])
  	Psi_s[[s]] <- diag(1 / psi_s[[s]])
  	###for study-specific f.l.
  	Lambda_s[[s]] <- zeros(p, j_s[s])
  	l_s[[s]] <- matrix(rnorm(n_s[s] * j_s[s]), n_s[s], j_s[s])
  	Meta_s[[s]] <- zeros(n_s[s], j_s[s])
  	Veta_s[[s]] <- eye(j_s[s])
  	###for common f.l.
  	f_s[[s]] <- matrix(rnorm(n_s[s]*k), n_s[s], k)
  	MetaP_s[[s]] <- zeros(n_s[s], k)
  	VetaP_s[[s]] <- eye(k)
  	###prior for study-specific f.l
  	omegajh_s[[s]] <- matrix(rgamma(p * j_s[s], shape = nu_s[s] / 2,
  	                                scale = 2 / nu_s[s]), p, j_s[s])
  	delta_s[[s]] <- rgamma(j_s[s], shape=c(a1_s[s], rep(a2_s[s], j_s[s] - 1)),
  	                       scale = c(b1_s[s], rep(b2_s[s], j_s[s] - 1)))
  	tauh_s[[s]] <- cumprod(delta_s[[s]])
  	Plam_s[[s]] <- matvec(omegajh_s[[s]], tauh_s[[s]])
  	#### save lambda and Psi####
    if(outputlevel == 1) {
  	  Lambdaout[[s]] <- zeros(p, j_s[s], sp)
     	psiout[[s]] <- zeros(p, 1, sp)
  	  l_sout[[s]] <- zeros(n_s[s], j_s[s], sp)
  	  f_sout[[s]] <- zeros(n_s[s], k, sp)
    }
  	else {
  	  Lambdaout[[s]] <- zeros(p, j_s[s])
  	  psiout[[s]] <- zeros(p, 1)
  	  l_sout[[s]] <- zeros(n_s[s], j_s[s])
  	  f_sout[[s]] <- zeros(n_s[s], k)
  	  SigmaLambda[[s]] <- zeros(p, p)
  	}
  	if(outputlevel == 3) SigmaPhi <- zeros(p, p)

  }
  #Phi
  Phi <- zeros(p, k)
  omegajh <- matrix(rgamma(p * k, shape = nu / 2, scale = 2 / nu), p, k)
  delta <- rgamma(k, shape = c(a1, rep(a2, k-1)), scale=c(b1, rep(b2, k - 1)))
  tauh <- cumprod(delta)
  Pphi <- matvec(omegajh, tauh)
  if(outputlevel<3) Phiout <- zeros(p, k, sp)
  else Phiout <- zeros(p, k)
  #### start posterior sampling
  for(r in 1:nrun)
  {
    f_s2 <- l_s2 <- list()
    for (s in 1:S){
      ###Step 1: common latent factors
   		LmsgP1 <- vecmat(psi_s[[s]], Phi)
    	VetaP1 <- diag(k) + t(LmsgP1) %*% Phi
    	TP1 <- chol(VetaP1)
    	qrTP1 <- qr(TP1)
    	QP1 <- qr.Q(qrTP1)
    	RP1 <- qr.R(qrTP1)
    	SP1 <- solve(RP1)
    	VetaP11 <- tcrossprod(SP1)
    	MetaP1 <- (Y_s[[s]] %*% LmsgP1 - l_s[[s]] %*% t(Lambda_s[[s]]) %*% LmsgP1) %*% VetaP11
    	xP1 <- matrix(rnorm(n_s[s] * k), nrow = n_s[s], ncol = k)
    	f_s[[s]] <- MetaP1 + xP1 %*% t(SP1)
      f_s2[[s]] <- crossprod(f_s[[s]])

      #Step 2: study-specific latent factors
   		Lmsg1 <- vecmat(psi_s[[s]], Lambda_s[[s]])
    	Veta1 <- diag(j_s[s]) + t(Lmsg1) %*% Lambda_s[[s]]
    	T1 <- chol(Veta1)
    	qrT1 <- qr(T1)
    	Q1 <- qr.Q(qrT1)
    	R1 <- qr.R(qrT1)
    	S1 <- solve(R1)
    	Veta11 <- tcrossprod(S1)
    	Meta1 <- (Y_s[[s]] %*% Lmsg1 -  f_s[[s]] %*% t(Phi) %*% Lmsg1 ) %*% Veta11
    	x1 <- matrix(rnorm(n_s[s] * j_s[s]), nrow = n_s[s], ncol = j_s[s])
    	l_s[[s]] <- Meta1 + x1 %*% t(S1)
      l_s2[[s]] <- crossprod(l_s[[s]])
   }

    ### Step 3: common factor loadings
    for(i in 1:p){
        q1 <- b1list <- list()
     	  for(s in 1:S){
     		  q1[[s]] <- psi_s[[s]][i] * f_s2[[s]]
     		  b1list[[s]] <- psi_s[[s]][i]*(t(f_s[[s]]) %*% Y_s[[s]][, i]) - psi_s[[s]][i] * (t(f_s[[s]]) %*% l_s[[s]] %*% Lambda_s[[s]][i,])
     		}
		   q1_sum <- Reduce('+', q1)
    	 bphi <- Reduce('+', b1list)
    	 Qphi <- diag(Pphi[i,], k) + q1_sum
    	 Lphi <-  t(chol(Qphi))
    	 zphi <- rnorm(k)
    	 vphi <- forwardsolve(Lphi, bphi)
    	 mphi <-  backsolve(t(Lphi), vphi)
    	 yphi <- backsolve(t(Lphi), zphi)
       Phi[i,] <- t(yphi + mphi)
    	}

    ### Step 4: specific factor loadings, with constraints
    for(s in 1:S){
      for(i in 1:p){
       	Qlam <- diag(Plam_s[[s]][i,], j_s[s]) + psi_s[[s]][i] * l_s2[[s]]
       	blam <- psi_s[[s]][i] * (t(l_s[[s]]) %*% Y_s[[s]][, i]) - psi_s[[s]][i] * (t(l_s[[s]]) %*% f_s[[s]] %*% Phi[i,])
       	Llam <- t(chol(Qlam))
       	zlam <- rnorm(j_s[s])
       	vlam <- forwardsolve(Llam, blam)
       	mlam <- backsolve(t(Llam), vlam)
       	ylam <- backsolve(t(Llam), zlam)
       	zlam4 <- zlam
       	Lambda_s[[s]][i,] <- t(ylam + mlam)
       	}
     }

   ### Step 5: omegajh for phi
    for(h in 1:k){
      omegajh[, h] <- rgamma(p, shape = (nu + 1) / 2, rate = (nu + tauh[h] * Phi[,h]^2) / 2)
      }
    mat <- omegajh * Phi^2
    ad <- a1 + 0.5 * p * k
    bd <- b1 + 0.5 * (1 / delta[1]) * sum(tauh * colSums(mat))
    delta[1] <- rgamma(1, shape = ad, scale = 1 / bd)
    tauh <- cumprod(delta)
    if (k>1){
      for(h in 2:k)
        {
         ad <- a2 + 0.5 * p * (k-h+1)
         bd <- b2 + 0.5 * (1 / delta[h]) * sum(tauh[h:k] * colSums(as.matrix(mat[, h:k])))
         delta[h] <- rgamma(1, shape = ad, scale = 1 / bd)
         tauh <- cumprod(delta)
        }
      }
    Pphi <- matvec(omegajh, tauh)
    ###################### omegajh for Lambda_s
    for(s in 1:S){
	    for(h in 1:j_s[s]){
      		omegajh_s[[s]][, h] <- rgamma(p, shape= (nu_s[[s]] + 1) / 2,
      		                              rate = (nu_s[[s]] + tauh_s[[s]][h] * Lambda_s[[s]][,h]^2) / 2)
    		}
    	mat_s <- omegajh_s[[s]] * Lambda_s[[s]]^2
    	ad <- a1_s[s] + 0.5 * p*j_s[s]
    	bd <- b1_s[s] + 0.5 * (1 / delta_s[[s]][1]) * sum(tauh_s[[s]] * colSums(mat_s))
    	delta_s[[s]][1] <- rgamma(1, shape = ad, scale = 1 / bd)
	    tauh_s[[s]] <- cumprod(delta_s[[s]])
    	if (j_s[s]>1){
    		for(h in 2:j_s[s]){
       		ad <- a2_s[s] + 0.5 * p * (j_s[s] - h + 1)
       		bd <- b2_s[s] + 0.5 * (1 / delta_s[[s]][h]) * sum(tauh_s[[s]][h:j_s[s]] * colSums(as.matrix(mat_s[, h:j_s[s]])))
       		delta_s[[s]][h] <- rgamma(1, shape = ad, scale = 1 / bd)
       		tauh_s[[s]] <- cumprod(delta_s[[s]])
      	 }
      }
    	Plam_s[[s]] <- matvec(omegajh_s[[s]], tauh_s[[s]])
    }

    ### Step 6: Psi_s
    for (s in 1:S){
    	Ytil_s <- Y_s[[s]] - (l_s[[s]] %*% t(Lambda_s[[s]]) + f_s[[s]] %*% t(Phi))
    	psi_s[[s]] <- rgamma(p, shape = a_psi_s[s] + 0.5 * n_s[s], rate = b_psi_s[s] + 0.5 * colSums(Ytil_s^2))
    	Psi_s[[s]] <- diag(1 / psi_s[[s]])
    }

     if(r > burn){
     neff <- (r - burn) / thin
     teff <- (nrun - burn) / thin
      if(outputlevel < 3) Phiout[, , neff] <- Phi
     else {
           Phiout <- Phiout + Phi / teff
           SigmaPhi <- SigmaPhi + tcrossprod(Phi) / teff
     }
     if(outputlevel==1)
     {
      for(s in 1:S){
        Lambdaout[[s]][, , neff] <- Lambda_s[[s]]
      	psiout[[s]][, , neff] <- 1 / psi_s[[s]]
      	l_sout[[s]][, , neff] <- l_s[[s]]
      	f_sout[[s]][, , neff] <- f_s[[s]]
       }
     }
     else{
       for(s in 1:S){
         Lambdaout[[s]] <- Lambdaout[[s]] + Lambda_s[[s]]  / teff
         psiout[[s]] <- psiout[[s]] + (1 / psi_s[[s]]) / teff
         l_sout[[s]] <- l_sout[[s]] + l_s[[s]] / teff
         f_sout[[s]] <- f_sout[[s]] + f_s[[s]] / teff
         SigmaLambda[[s]] <- SigmaLambda[[s]] + tcrossprod(Lambda_s[[s]]) / teff
       }
     }
   }
  if (trace & r %% nprint == 0) cat("r=",r,"\n")
  }
  ##### Save and exit
  if(outputlevel < 3) SigmaPhi <- NULL
  if(outputlevel==1) SigmaLambda <- NULL
  out <- list(Phi = Phiout, Lambda = Lambdaout, psi = psiout, l_s = l_sout, f_s = f_sout,
              SigmaPhi = SigmaPhi, SigmaLambda = SigmaLambda)
 return(structure(out,  class="sp_msfa"))
}


#' Normalized eigenvalues of a symmetric matrix.
#'
#' This is a rather simple function, performing the eigenvalue computation for a
#' covariance matrix. Used for the selection of common latent factor dimensions.
#'
#' @param SigPhi Symmetric matrix, of size \eqn{p \times p}{p x p}.
#' @return The normalized eigenvalues of \code{SigPhi}, a vector of length \eqn{p}{p}.
#' @export
sp_eigen <- function(SigPhi){
	#svdPhi <- svd(SigPhi)  #Since SigPhi is squared, eigen or svd are equivalent
	eigenPhi <- eigen(SigPhi)
  d_value <- eigenPhi$values #svdPhi$d
	choiceK <- d_value / sum(d_value)
	return(choiceK)
}




#' OP algorithm for loading matrix estimation
#'
#' Iterative algorithm for estimating a loading matrix from the posterior samples, resolving
#' the orthogonality indeterminacy.
#'
#' @param Phi Posterior samples for a loading matrix, a 3-dimensional array of
#'  dimension \eqn{p \times k \times nchain}{p x k x nchain}.
#' @param tol Tolerance  for the OP algorithm. Default is 0.001.
#' @param itermax Maximum number of iterations. Default is 10.
#' @return A list containing the estimated loading matrix, plus the number of iterations
#' performed and the Frobenius norm of the difference between the last two iterations
#' of the algorithm.
#' @export
sp_OP <- function(Phi, tol = 10^-3, itermax = 10){
	nchain <- dim(Phi)[3]
	p <- dim(Phi)[1]
	k <- dim(Phi)[2]
	S <- array(0, dim = c(p, p, nchain))
	U_phi <- array(0, dim = c(p, p, nchain))
	V_phi <- array(0, dim = c(p, p, nchain))
	D_phi <- array(0, dim = c(p, p, nchain))
	Phi_final <- array(0, dim = c(p, k, nchain))
	Phi_star0 <- Phi[,, nchain]
	iter <- 0
	repeat
	{
	  iter <- iter + 1
	  for (i in 1:nchain){
		S[,,i] <- tcrossprod(Phi[,, i], Phi_star0)
		SVD_phi <- svd(S[,, i])
		U_phi[,, i] <-  SVD_phi$u
		V_phi[,, i] <-  SVD_phi$v
		D_phi[,, i] <- tcrossprod(U_phi[,, i], V_phi[,, i])
		Phi_final[,, i] <-  D_phi[,, i] %*% Phi[,, i]
		}
	  Phi_star1 <- apply(Phi_final, c(1,2), mean)
	  normd <- norm(Phi_star0-Phi_star1, type = "F")
	  Phi_star0 <- Phi_star1
	  if(normd < tol) break
	  if(iter==itermax) break
 	  }
	  return(list(Phi = Phi_star1, iter = iter, normd = normd))
}


#' Auxiliary function for \code{sp_fa}.
#'
#' Set parameters for posterior sampling and prior hyperparameters for the sparse Bayesian infinite factor model. For the latter, the notation
#' follows closely the paper by Bhattacharya and Dunson (2011).
#'
#' @param nrun Number of posterior simulations. Default is 30000.
#' @param burn Burn-in trials. Default is 20000.
#' @param thin Thinning of posterior samples. Default is 1 (no thinning).
#' @param nu Parameter entering the gamma distribution assumed for phi. Default is 3.
#' @param asigma Shape parameter for the gamma distribution assumed for 1/sigma^2. Default is 1.
#' @param bsigma Scale parameter for the gamma distribution assumed for 1/sigma^2. Default is 0.3.
#' @param a1 Shape parameter for the gamma distribution assumed for delta1. Default is 2.1.
#' @param b1 Scale parameter for the gamma distribution assumed for delta1. Default is 1.
#' @param a2 Shape parameter for the gamma distribution assumed for deltal. Default is 2.1.
#' @param b2 Scale parameter for the gamma distribution assumed for deltal. Default is 1.
#' @export
#' @references Bhattacharya, A. and Dunson, D.B. (2011). Sparse Bayesian infinite factor models. Biometrika,
#' 98, p. 291-306.
sp_fa_control <- function(nrun = 30000, burn = 20000, thin = 1,
                          nu = 3,
                          asigma = 1, bsigma = 0.3,
                          a1 = 2.1, b1 = 1, a2 = 2.1, b2 = 1)
{
  return(list(nrun = nrun, burn = burn, thin = thin,
              a1 = a1, a2 = a2, b1 = b1, b2 = b2, nu = nu,
              asigma = asigma, bsigma = bsigma))
}



#' Posterior sampling for the sparse Bayesian infinite factor model.
#'
#' Ported to R from matlab code supporting the paper by Bhattacharya and Dunson (2011).
#' Courtesy of A. Bhattacharya.
#'
#' @param dat Data matrix, of size \eqn{n \times p}{n x p}.
#' @param k Number of latent factor.
#' @param trace If \code{TRUE} then trace information is being printed every \code{nprint} iterations of the Gibbs sampling.
#' Default is \code{TRUE}.
#' @param nprint Frequency of tracing information. Default is every 1000 iterations.
#' @param control A list of hyperparameters for the prior distributions and for controlling the Gibbs sampling.
#' See \code{sp_fa_control}.
#' @param ... Arguments to be used to form the default \code{control} argument if it is not supplied directly.
#' @return  A list, with components \code{Sigma} and \code{Lambda}, containing the posterior samples for the
#' covariance matrix of the model,
#' with dimension  \eqn{p \times p \times (nrun - code)/thin}{p x p x (nrun - burn)/thin}, and
#' of the loading matrix, with dimension  \eqn{p \times k \times (nrun - code)/thin}{p x k x (nrun - burn)/thin},
#' respectively.
#' @import statmod matlab
#' @importFrom stats rgamma rnorm
#' @export
#' @references Bhattacharya, A. and Dunson, D.B. (2011). Sparse Bayesian infinite factor models. Biometrika,
#' 98, p. 291-306.
sp_fa <- function(dat, k, trace = TRUE, nprint = 1000, control = list(...), ...)
{
  ##### read data
  p <- ncol(dat)
  n <- nrow(dat)
  Y <- scale(dat)

  #### setting up priors and initialize
  control <- do.call("sp_fa_control", control)
  nrun <- control$nrun
  thin <- control$thin
  burn <- control$burn
  sp  <- (nrun - burn) / thin

  asigma <- control$asigma
  bsigma <- control$bsigma
  nu <- control$nu
  a1 <- control$a1
  a2 <- control$a2
  b1 <- control$b1
  b2 <- control$b2

  sigma <- rgamma(p, shape = asigma, scale = 1 / bsigma)
  Sigma <- diag(1 / sigma)
  Lambda <- zeros(p, k)
  eta <-  matrix(rnorm(n * k), n, k)
  meta <- zeros(n, k)
  veta <- eye(k)
  phijh <- matrix(rgamma(p * k, shape = nu / 2, scale = 2 / nu), p, k)
  delta <- rgamma(k, shape=c(a1, rep(a2, k-1)), scale=c(b1, rep(b2, k-1)))
  tauh <- cumprod(delta)
  Plam <- matvec(phijh, tauh)
  nofout <- zeros(nrun + 1, 1)
  Omegaout <- zeros(p, p, sp)
  Lambdaout <- zeros(p, k, sp)

  #### start posterior sampling
  for(i in 1:nrun)
  {
    ### Step 3
    Lmsg <- vecmat(sigma, Lambda)
    Veta1 <- diag(k) + t(Lmsg) %*% Lambda
    T <- chol(Veta1)
    qrT <- qr(T)
    Q <- qr.Q(qrT)
    R <- qr.R(qrT)
    S <- solve(R)
    Veta <- tcrossprod(S)
    Meta <- Y %*% Lmsg %*% Veta
    x <- matrix(rnorm(n * k), nrow = n, ncol = k)
    eta <- Meta + x %*% t(S)

    ### Step 1
    eta2 <- crossprod(eta)
    for(j in 1:p)
    {
      Qlam <- diag(Plam[j,]) + sigma[j] * eta2
      blam <- sigma[j] * (t(eta) %*% Y[,j])
      Llam <- t(chol(Qlam))
      zlam <- rnorm(k)
      vlam <- forwardsolve(Llam, blam)
      mlam <- backsolve(t(Llam), vlam)
      ylam <- backsolve(t(Llam), zlam)
      Lambda[j,] <- t(ylam + mlam)
    }

    ### Step 4
    for(h in 1:k)
      phijh[, h] <- rgamma(p, shape= (nu + 1) / 2, rate = (nu + tauh[h] * Lambda[,h]^2) / 2)

    ### Step 5
    mat <- phijh * Lambda^2
    ad <- a1 + 0.5 * p * k
    bd <- b1 + 0.5 * (1 / delta[1])  * sum(tauh * colSums(mat))
    delta[1] <- rgamma(1, shape = ad, scale = 1 / bd)
    tauh <- cumprod(delta)
    for(h in 2:k)
    {
      ad <- a2 + 0.5 * p * (k - h + 1)
      bd <- b2 + 0.5 * (1 / delta[h]) * sum(tauh[h:k] * colSums(as.matrix(mat[, h:k])))
      delta[h] <- rgamma(1, shape = ad, scale = 1 / bd)
      tauh <- cumprod(delta)
    }

    ### Step 2
    Ytil <- Y - eta %*% t(Lambda)
    sigma <- rgamma(p, shape = asigma + 0.5 * n, rate = bsigma + 0.5 * colSums(Ytil^2))
    Sigma <- diag(1 / sigma)
    Plam <- matvec(phijh, tauh)

    ###Check burn in and then conclude
    if( i > burn)
    {
      Omega <- tcrossprod(Lambda) + Sigma
      Omegaout[,, (i - burn) / thin] <- Omega
      Lambdaout[,, (i - burn) / thin] <- Lambda
    }
    if(trace & i%%nprint == 0) print(i)
  }
  out <- list(Sigma = Omegaout, Lambda = Lambdaout)
  return(structure(out,  class="sp_fa"))
}



