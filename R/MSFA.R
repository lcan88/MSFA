#' @keywords internal
tr <- function(A) sum(diag(A))


#' @import statmod
#' @keywords internal
exp_values <- function(Phi, Lambda_s, Psi_s, Psi_s1, cov_s, getdet = FALSE)
{
   k <- dim(Phi)[2]
   I_k <- diag(1, k)
   S <- length(Lambda_s)

   ###defining objects
   j_s <- numeric(S)
   I_j <- list()
   Sig_s <- list()
   ds_s <- list()
   I_tot <- list()
   LambTOT <- list()
   Sig_s1 <- list()
   delta_Lambda <- list()
   delta_Phi <- list()
   Delta_Lambda <- list()
   Delta_Phi <- list()
   Covfcfs <- list()
   Txsfs <- list()
   Txsfcs <- list()
   Tfsfs <- list()
   Tfcsfcs <- list()
   Tfcsfs <- list()

  for (s in 1:S){
   	ds_s[[s]] <- NULL
   	j_s[s] <- c(dim(Lambda_s[[s]])[[2]])
    I_j[[s]] <- diag(1, j_s[s])
    Sig_s[[s]] <- Phi %*% t(Phi) + Lambda_s[[s]] %*% t(Lambda_s[[s]]) + Psi_s[[s]]
    if (getdet)  {ds_s[[s]] <- det(Sig_s[[s]])}
    I_tot[[s]] <- diag(1, k + j_s[s])
    LambTOT[[s]] <- cbind(Phi, Lambda_s[[s]])
    Sig_s1[[s]] <- Psi_s1[[s]] - (statmod::vecmat(diag(Psi_s1[[s]]), LambTOT[[s]]) %*%
                   solve(I_tot[[s]] + (t(LambTOT[[s]]) %*% statmod::vecmat(diag(Psi_s1[[s]]),
                  LambTOT[[s]]))) %*% statmod::matvec(t(LambTOT[[s]]), diag(Psi_s1[[s]])))
    delta_Lambda[[s]] <- t(Lambda_s[[s]]) %*% Sig_s1[[s]]
    delta_Phi[[s]] <- t(Phi) %*% Sig_s1[[s]]
    Delta_Lambda[[s]] <- I_j[[s]] - (t(Lambda_s[[s]]) %*% Sig_s1[[s]] %*% Lambda_s[[s]])
    Delta_Phi[[s]] <- I_k - (t(Phi) %*% Sig_s1[[s]] %*% Phi)
    Covfcfs[[s]] <- -t(Phi) %*% Sig_s1[[s]] %*% Lambda_s[[s]]
    Txsfs[[s]] <- cov_s[[s]] %*% t(delta_Lambda[[s]])
    Txsfcs[[s]] <- cov_s[[s]] %*% t(delta_Phi[[s]])
    Tfsfs[[s]] <- delta_Lambda[[s]] %*% cov_s[[s]] %*% t(delta_Lambda[[s]]) + Delta_Lambda[[s]]
    Tfcsfcs[[s]] <- delta_Phi[[s]] %*% cov_s[[s]] %*% t(delta_Phi[[s]]) + Delta_Phi[[s]]
    Tfcsfs[[s]] <- delta_Phi[[s]] %*% cov_s[[s]] %*% t(delta_Lambda[[s]]) + Covfcfs[[s]]
   }
 return(list(Txsfs = Txsfs, Txsfcs = Txsfcs, Tfsfs = Tfsfs,
             Tfcsfcs =  Tfcsfcs, Tfcsfs = Tfcsfs, ds_s=ds_s,  Sig_s1 = Sig_s1))
}




#' Provides some starting values for the parameters of a MSFA model
#'
#' This is a supporting function for \code{ecm_msfa}. The method employed is documented in the reference.
#'
#' The upper-triangular zero constraint is adopted to achieve identification,
#' as detailed in the reference, though the function can also be run without such constraint.
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
#' No standardization is carried out by the function.
#' @param k Number of common factors.
#' @param j_s Number of study-specific factors. A vector of positive integers of length \eqn{S}{S}.
#' @param block_lower Should the upper-triangular zero constraint be enforced? Default is \code{TRUE}.
#' @param method Which method should be used to find the starting values? The two possibilities are \code{"adhoc"} for
#' the method described in De Vito et al. (2016), and \code{"fa"} for averaging over separate study-specific FA models.
#' Default is \code{"adhoc"}.
#' @param robust If \code{TRUE}, robust covariance matrix is used in place of the sample covariance. Default
#' is \code{FALSE}.
#' @param corr If \code{TRUE}, the analysis will employ the correlation matrix instead of the covariance matrix.
#' @param mcd If \code{TRUE}, the robust estimator used for the covariance is the same proposed in Pison et al. (2003),
#' otherwise the default value of the function \code{CovRob} of the \code{robust} library is employed. Default is
#' \code{FALSE}.
#' @return A list  containing  \code{Phi},\code{Lambda_s} and  \code{Psi_s}, starting values for the model matrices.
#' @import psych
#' @export
#' @references De Vito, R., Bellio, R., Parmigiani, G. and Trippa, L. (2016). Multi-study Factor Analysis,
#' arXiv:1611.06350. \url{https://arxiv.org/abs/1611.06350}
start_msfa <- function(X_s, k, j_s, block_lower = TRUE, method = "adhoc", robust = FALSE, corr = FALSE, mcd = FALSE)
{
  X_used_s <- X_s
  S <- length(X_s)
  if(corr & !robust)
    for(s in 1:S)  X_used_s[[s]] <- scale(X_s[[s]])
  if(robust & corr & method=="adhoc"){
    for(s in 1:S){
      ogg_s <- if(mcd) covRob(X_s[[s]], estim = "mcd", quan = .75, ntrial = 1000) else covRob(X_s[[s]])
    }
  X_used_s[[s]] <- scale(X_s[[s]], center = ogg_s$center, scale = sqrt(diag(ogg_s$cov)))
  }
  p <- dim(X_s[[1]])[2]
  Phi <- matrix(0, nrow=p, ncol=k)
  Lambda_s <- Psi_s <- list()
  if(method=="adhoc"){
     X <- Reduce(rbind, X_used_s)
     X.pcr <- prcomp(X)
     Phi <- matrix(X.pcr$rotation[,1:k], nrow=p, ncol=k, byrow=FALSE)
     #totvar <- sum(X.pcr$sdev[1:k])
     if (block_lower) Phi[upper.tri(Phi)] <- 0
     for(s in 1:S){
     #   Xscaled.s <- statmod::matvec(X_s[[s]], rep(1 / totvar, p))
     iniLS <- array(prcomp(X_used_s[[s]])$rotation, dim=c(p, j_s[s]))
     iniTot <- cbind(Phi, iniLS)
     if (block_lower) iniTot[upper.tri(iniTot)] <- 0
     Lambda_s[[s]] <-  matrix(iniTot[,(k+1):(k+j_s[s])], p , j_s[s])
     Psi_s[[s]] <- diag(fa(X_used_s[[s]], nfactors = k+j_s[s])$uniq)
     }
  }
 #### method "fa" is still temporary, since the loading sign should be post-processed before taking the mean
 if(method=="fa"){
      est <- ecm_fa(X_s, tot_s = k + j_s, robust = robust, mcd = mcd, corr = corr, tol = 10^-5, nIt = 5000)
      Phi <- est$Omega_s[[1]][,1:k] / S
      Lambda_s[[1]] <-  est$Omega_s[[1]][,(k+1):(k+j_s[1])]
      Psi_s[[1]] <- est$Psi_s[[1]]
      for(s in 2:S){
        Phi <- Phi + est$Omega_s[[s]][,1:k] / S * sign(Phi) * sign(est$Omega_s[[s]][,1:k]) ###to avoid sign changes
        Lambda_s[[s]] <-  est$Omega_s[[s]][,(k+1):(k+j_s[s])]
        Psi_s[[s]] <- est$Psi_s[[s]]
       }
  }
  out <- list(Phi=Phi, Lambda_s=Lambda_s, Psi_s=Psi_s)
  return(out)
}



#' @keywords internal
loglik_ecm <- function(Sig_s1,  ds_s, n_s, cov_s)
{
   S <- length(n_s)
   #####log likelihood value for each study
   val_s <- c()
   for(s in 1:S){
	    val_s[s] <- - (n_s[s]/2) * log(ds_s[[s]]) - (n_s[s]/2) * tr(Sig_s1[[s]] %*% cov_s[[s]])
	    }
   #####sum of each study-likelihood
   val_tot <- sum(val_s)
   return(val_tot)
}


#' @keywords internal
#### loglikelihood function re-expressed as a function of the model parameters
#### theta: c(Phi, Lambda_1,..,Lambda_S,Psi_1,..,Psi_S)
loglik_int <- function(theta, n_s, cov_s, k, j_s)
{
  S <- length(n_s)
  p <- ncol(cov_s[[1]])
  nP <- k * p - k * ( k - 1) / 2
  nL <- j_s * (p - k)  - j_s *  (j_s - 1) / 2
  phi_vals <- theta[1:nP]
  out <- 0
  for(s in 1:S){
    nL_s  <- if(s==1) 0 else sum(nL[1:(s-1)])
    ind <-  (nP + nL_s + 1):(nP + nL_s + nL[s])
    omega_vals_s <- c(phi_vals, theta[ind])
    Omega_s <- matrix(0, p, k + j_s[s])
    Omega_s[lower.tri(Omega_s, diag = TRUE)] <- omega_vals_s
    ind_s <- (nP + sum(nL) + p * (s-1) + 1):(nP + sum(nL) + p * s)
    psi_vals_s <- theta[ind_s]
    #Sigma_s <- tcrossprod(Omega_s) + diag(psi_vals_s)
    Psi_s1 <- diag(1/psi_vals_s)
    D1L_s <- statmod::vecmat(1/sqrt(psi_vals_s), Omega_s)
    LDL_s <- crossprod(D1L_s)
    A <- diag(k + j_s[s]) + LDL_s
    A1 <- chol2inv(chol(A))
    D2L_s <- statmod::vecmat(1/psi_vals_s, Omega_s)
    Sig_s1 <- Psi_s1 - D2L_s  %*% A1 %*% t(D2L_s)
    #Sig_s1 <- chol2inv(chol(Sigma_s))
    log_ds_s <-  log(det(A)) + sum(log(psi_vals_s))
    out  <- out - (n_s[s]/2) * log_ds_s - (n_s[s]/2) * tr(Sig_s1 %*% cov_s[[s]])
    }
  return(out)
}


#' Variance matrix of MLE estimates for a MSFA model
#'
#' Computes the inverse observed information for a MSFA model
#'
#'
#' Numerical differentiation is employed to obtain the observed information matrix at a
#' given parameter values, so that when the parameter values equals the MLE the function
#' returns the estimated variance matrix of the fitted model. The method is rather inefficient, and
#' it may lead to long computations, though the function is designed to be called only once after the
#' estimation has been carried out. However, it would be relatively straightforward to employ analytical
#' differentiation at least for the log-likelihood gradient, and this may be implemented in future
#' releases of the code.
#'
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all
#' the studies.
#' @param param A list with the same slots as that returned by \code{ecm_msfa}.
#' @param getgrad Should the function return also the gradient at \code{param}? Default is \code{FALSE}.
#' @return A list with exactly the same structure of the three slots \code{Phi}, \code{Lambda_s} and
#' \code{Psi_s} of \code{param}, but containing the standard errors rather than the point estimates.
#' Furthemore, slots for the hessian matrix and the gradient at \code{param} are included, the latter
#' is not NULL when \code{getgrad}  is \code{TRUE}.
#' @export
#' @import statmod
#' @importFrom pracma grad
#' @importFrom pracma hessian
vcov_msfa <- function(X_s, param, getgrad = TRUE)
{
  p <- ncol(X_s[[1]])
  phi_vals <- as.vector(param$Phi[lower.tri(param$Phi, diag = TRUE)])
  S <- length(X_s)
  k <- ncol(param$Phi)
  lambda_vals <- psi_vals <- j_s <- c()
  for(s in 1:S){
    Lambda_s <- param$Lambda_s[[s]][-(1:k),]
    lambda_vals <- c(lambda_vals, as.vector(Lambda_s[lower.tri(Lambda_s, diag = TRUE)]))
    psi_vals <- c(psi_vals, diag(param$Psi_s[[s]]))
    j_s[[s]] <- ncol(param$Lambda_s[[s]])
  }
  theta <- c(phi_vals, lambda_vals, psi_vals)
  gout <- NULL
  if(getgrad) gout <- pracma::grad(loglik_int, x = theta, n_s = param$n_s, cov_s = param$cov_s, k = k, j_s = j_s)
  hout <- pracma::hessian(loglik_int, x = theta, n_s = param$n_s, cov_s = param$cov_s, k = k, j_s = j_s)
  seout <- sqrt(diag(chol2inv(chol(-hout))))
  ###### re-arrange the ses like in param
  nP <- k * p - k * ( k - 1) / 2
  nL <- j_s * (p - k)  - j_s *  (j_s - 1) / 2
  phi_vals <- seout[1:nP]
  Phi <- matrix(0, nrow=p, ncol=k)
  Phi[lower.tri(Phi, diag = TRUE)] <- phi_vals
  Lambda_s <- param$Lambda_s
  Psi_s <- param$Psi_s
  for(s in 1:S){
    nL_s  <- if(s==1) 0 else sum(nL[1:(s-1)])
    ind <-  (nP + nL_s + 1):(nP + nL_s + nL[s])
    lambda_vals_s <-  seout[ind]
    Lambda_s[[s]][Lambda_s[[s]]!=0] <- lambda_vals_s
    ind_s <- (nP + sum(nL) + p * (s-1) + 1):(nP + sum(nL) + p * s)
    psi_vals_s <- seout[ind_s]
    Psi_s[[s]] <-  diag(psi_vals_s)
  }
  out <- list(Phi = Phi, Lambda_s = Lambda_s, Psi_s = Psi_s,  grad = gout, hessian = hout)
  return(out)
}




#' Estimates the parameters of a MSFA model
#'
#' Maximum likelihood estimation of the MSFA model parameters via the ECM
#' algorithm.
#'
#' The upper-triangular zero constraint is adopted to achieve identification,
#' as detailed in the reference, though the function can also be run without such constraint.
#' No checking is done on the starting value for the various model matrices,
#' since a suitable value for them  is produced by the function \code{start_msfa}.
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
#' @param start A list containing the slots \code{Phi}, \code{Lambda_s} and \code{Psi_s}, containing the starting
#' values for the matrix  \code{Phi} of common factor loadings, of size \eqn{P \times K}{P x K}, for
#' the matrices \code{Lambda_s} of study-specific factor loadings, a list of size \eqn{S}{S}  where each element
#' contains a matrix with \eqn{P \times J_s}{P x J_s}, and finally for the study-specific matrices of uniquenesses,
#' a list of size \eqn{S}{S}, where each element contains a diagonal matrix with \eqn{P \times P}{P x P}.
#' Note that a suitable list of this kind is produced by \code{start_msfa}.
#' @param nIt Maximum number of iterations for the ECM algorithm. Default is 20000.
#' @param tol Tolerance for declaring convergence of the ECM algorithm. Default is 10^-7.
#' @param block_lower Should the upper-triangular zero constraint be enforced? Default is \code{TRUE}
#' (strongly suggested).
#' @param robust If \code{TRUE}, robust covariance matrix is used in place of the sample covariance. Default
#' is \code{FALSE}.
#' @param corr If \code{TRUE}, the analysis will employ the correlation matrix instead of the covariance matrix.
#' @param mcd If \code{TRUE}, the robust estimator used for the covariance is the same proposed in Pison et al. (2003),
#' otherwise the default value of the function \code{CovRob} of the \code{robust} library is employed. Default is
#' \code{FALSE}.
#' @param trace If \code{TRUE} then trace information is being printed every 1000 iterations of the ECM algorithm.
#' @return A list  containing the following components:
#' \item{\code{Phi},\code{Lambda_s}, \code{Psi_s}}{the estimated model matrices.}
#' \item{loglik}{the value of the log likelihood function at the final estimates.}
#' \item{\code{AIC, BIC}}{model selection criteria at the estimate.}
#' \item{\code{npar}}{number of model parameters.}
#' \item{iter}{the number of ECM iterations performed.}
#' @export
#' @import robust
#' @importFrom stats cor cov factanal prcomp
#' @references De Vito, R., Bellio, R., Trippa, L. and Parmigiani, G. (2016). Multi-study Factor Analysis,
#' arXiv:1611.06350. \url{https://arxiv.org/abs/1611.06350}
#' @references Pison, G., Rousseeuw, P.J., Filzmoser, P. and Croux, C. (2003). Robust factor analysis. Journal
#' Multivariate Analysis, 84, 145-172.
ecm_msfa <- function(X_s, start = NULL, nIt = 50000, tol = 10^-10, block_lower = TRUE, robust = FALSE,
                     corr = TRUE, mcd = FALSE, trace = TRUE)
{
  Phi <- start$Phi
  Lambda_s <- start$Lambda_s
  Psi_s <- start$Psi_s
  #######
  p <- dim(Phi)[[1]]
  k <- dim(Phi)[[2]]
  S <- length(X_s)
  j_s <- n_s <- numeric(S)
  #######defining objects
  Psi_s1 <- list()
  cov_s <- list()
  L_s <- list()

  ######1st round of cycle
  for(s in 1:S){
  	n_s[s] <-  dim(X_s[[s]])[[1]]
  	j_s[s] <-  dim(Lambda_s[[s]])[[2]]
    Psi_s1[[s]] <-  diag(1/diag(Psi_s[[s]]))
    if((!robust) & (!corr)) cov_s[[s]] <- cov(X_s[[s]])
    if((!robust) & corr) cov_s[[s]] <- cor(X_s[[s]])
    if(robust & mcd) cov_s[[s]] <- covRob(X_s[[s]], estim = "mcd", quan = .75, ntrial = 1000, corr = corr)$cov
    if(robust & (!mcd)) cov_s[[s]] <- covRob(X_s[[s]], corr = corr)$cov
    }
  ######E-step
  out <- exp_values(Phi, Lambda_s, Psi_s, Psi_s1, cov_s, getdet = TRUE)
  Sig_s1 <- out$Sig_s1; ds_s= out$ds_s;
  l_stop0 <- 0
  lm1 <- 0
  l0 <- loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)

  for (i in (1:nIt))
  {
   ###########CM1 ---------------------------------------------------------------------------------------

   ######expected values
   out <- exp_values(Phi, Lambda_s, Psi_s, Psi_s1, cov_s)
   Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs; Tfcsfcs <- out$Tfcsfcs; Tfcsfs <- out$Tfcsfs
   ######update  of Phi_s
   Psi_new <- list()
   Psi_new1 <- list()

   for(s in 1:S){
   	Psi_new[[s]]  <- diag(cov_s[[s]] + Phi %*% Tfcsfcs[[s]] %*% t(Phi) + Lambda_s[[s]] %*%
   	                 Tfsfs[[s]] %*% t(Lambda_s[[s]]) - 2*Txsfcs[[s]] %*% t(Phi) -  2*Txsfs[[s]] %*% t(Lambda_s[[s]]) +
   	                 2 * Phi %*% Tfcsfs[[s]] %*% t(Lambda_s[[s]]))
   	Psi_new[[s]] <- diag(Psi_new[[s]])
   	##########inverse
   	Psi_new1[[s]] <- diag(1/diag(Psi_new[[s]]))
   	 }

   ###########CM2 ---------------------------------------------------------------------------------------

   ######expected values
   out<- exp_values(Phi, Lambda_s, Psi_new, Psi_new1, cov_s)
   Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs;
   Tfcsfcs <- out$Tfcsfcs; Tfcsfs <- out$Tfcsfs

   ######update of Phi
   C_s <- list()
   kron_s <- list()
      for(s in 1:S){
      	C_s[[s]] <- n_s[s] * Psi_new1[[s]] %*% Txsfcs[[s]] - n_s[s] * Psi_new1[[s]] %*% Lambda_s[[s]] %*% t(Tfcsfs[[s]])
      	kron_s[[s]] <- kronecker(t(Tfcsfcs[[s]]), n_s[s] * Psi_new1[[s]])
      	}

    C <- Reduce('+', C_s)
    kron <- Reduce('+', kron_s)
    Phi_vec <- solve(kron) %*% matrix(as.vector(C))
    Phi_new <- matrix(Phi_vec, p, k)


    if (block_lower) Phi_new[upper.tri(Phi_new)] <- 0

   ########CM3 ---------------------------------------------------------------------------------------

   ######expected values
   out <- exp_values(Phi_new, Lambda_s, Psi_new, Psi_new1, cov_s)
   Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs;
   Tfcsfcs <-  out$Tfcsfcs; Tfcsfs <- out$Tfcsfs

   ######update of Phi
   Lambda_new <- list()
   for(s in 1:S){
   	Lambda_new[[s]] <- matrix(((Txsfs[[s]] - Phi_new %*% Tfcsfs[[s]]) %*% solve(Tfsfs[[s]])), p, j_s[s])
   }

  if (block_lower)  {
   	for (s in 1:S){
   		L_s[[s]] <- cbind(Phi_new, Lambda_new[[s]])
   		L_s[[s]][upper.tri(L_s[[s]])] <- 0
   		Phi_new <- matrix(L_s[[s]][,1:k], nrow=p, ncol = k)
      Lambda_new[[s]] <- L_s[[s]][,(k+1):(k+j_s[s])]
        }
    }

  ###########stopping rule
  out <- exp_values(Phi_new, Lambda_new, Psi_new, Psi_new1, cov_s, getdet = TRUE)
  Sig_s1 <- out$Sig_s1
  ds_s <- out$ds_s
  l1 <- loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
  a <- (l1 - l0)/ (l0-lm1)
  l_stop <- lm1 + (1/ (1-a)) * (l0-lm1)
  l0 <- loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
  if((trace) & (i %% 100 == 0))  cat("i=", i, "Criterion for convergence ", abs(l_stop-l_stop0),  "\n")
  if((abs(l_stop-l_stop0)<tol) & i > 1 & l_stop != Inf) break
  Psi_s <- Psi_new
  Phi <- Phi_new
  Lambda_s <- Lambda_new
  Psi_s1 <- Psi_new1
  lm1 <- l0
  l0 <- l1
  l_stop0 <- l_stop
 }
  ############return output
  # part1 <- c()
  # part2 <- c()
  # tot_s <- j_s + k
  # for(s in 1:S){
  #   part1[s] <- c(j_s[s] * p)
  #   part2[s] <- c((tot_s[s] * (tot_s[s] - 1)) / 2)
  # }
  # sum1 <- k * p + sum(part1) + p * S
  # sum2 <- sum(part2) - (k * (k - 1)) / 2
  # npar <- (sum1  - sum2)
  n_tot <- sum(n_s)
  npar <- p * S + k * (p - ( k - 1) / 2) +  sum(j_s * (p - k - (j_s - 1) / 2))
  AIC <- -2 * l1 + npar * 2
  BIC <- -2 * l1 + npar * log(n_tot)
  res <- list(Phi = Phi, Lambda_s = Lambda_s, Psi_s = Psi_s, loglik = l1,
              AIC = AIC, BIC = BIC, npar=npar,
              iter = i,  cov_s = cov_s, n_s = n_s)
 return(res)
}


#' Estimates the parameters of study-specific FA models
#'
#' Maximum likelihood estimation of study-specific FA models parameters via the ECM
#' algorithm, adopting the upper-triangular zero constraint to achieve identification
#' for each loading matrix. Note: the function can also estimate a FA model for a single
#' study, by specifiyng \code{X_s = list(data)}, where \code{data} is the data matrix.
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
#' @param tot_s Number of latent factors for each study. A vector of positive integers of length \eqn{S}{S}.
#' @param nIt Maximum number of iterations for the ECM algorithm. Default is 20000.
#' @param tol Tolerance for declaring convergence of the ECM algorithm. Default is 10^-7.
#' @param block_lower Should the upper-triangular zero constraint be enforced? Default is \code{TRUE}
#' (strongly suggested).
#' @param robust If \code{TRUE}, robust covariance matrix is used in place of the sample covariance. Default
#' is \code{FALSE}.
#' @param corr If \code{TRUE}, the analysis will employ the correlation matrix instead of the covariance matrix.
#' @param mcd If \code{TRUE}, the robust estimator used for the covariance is the same proposed in Pison et al. (2003),
#' otherwise the default value of the function \code{CovRob} of the \code{robust} library is employed. Default is
#' \code{FALSE}.
#' @param trace If \code{TRUE} then trace information is being printed every \code{traceIT} iterations of the ECM algorithm.
#' @param traceIT Frequency of tracing information.
#' @return A list  containing the following components:
#' \item{\code{Omega_s}, \code{Psi_s}}{the estimated model matrices.}
#' \item{loglik}{the value of the log likelihood function at the final estimates.}
#' \item{\code{AIC, BIC}}{model selection criteria at the estimate.}
#' \item{\code{npar}}{number of model parameters.}
#' \item{iter}{the number of ECM iterations performed.}
#' @export
#' @import robust psych
#' @references De Vito, R., Bellio, R., Trippa, L. and Parmigiani, G. (2016). Multi-study Factor Analysis,
#' arXiv:1611.06350. \url{https://arxiv.org/abs/1611.06350}
#' @references Pison, G., Rousseeuw, P.J., Filzmoser, P. and Croux, C. (2003). Robust factor analysis. Journal
#' Multivariate Analysis, 84, 145-172.
ecm_fa <- function(X_s, tot_s, nIt = 50000, tol = 10^-10, block_lower = TRUE, robust = FALSE, corr = TRUE, mcd = FALSE, trace = TRUE, traceIT = 1000)
{
  Omega_s <- list()
  Psi_s <- list()
  #######
  p <- ncol(X_s[[1]])
  S <- length(X_s)
  n_s <- numeric(S)
  #######defining objects
  Psi_s1 <- list()
  cov_s <- list()
  Phi <- matrix(0, nrow=p, ncol=1)
  ######1st round of cycle
  for(s in 1:S){
    n_s[s] <-  dim(X_s[[s]])[[1]]
    if((!robust) & (!corr)) cov_s[[s]] <- cov(X_s[[s]])
    if((!robust) & corr) cov_s[[s]] <- cor(X_s[[s]])
    if(robust & mcd) cov_s[[s]] <- covRob(X_s[[s]], estim = "mcd", quan = .75, ntrial = 1000, corr = corr)$cov
    if(robust & (!mcd)) cov_s[[s]] <- covRob(X_s[[s]], corr = corr)$cov
    FA.s <- factanal(X_s[[s]], factors = tot_s[[s]], covmat = cov_s[[s]],  n.obs=nrow(X_s[[s]]), rotation = "none")
    Omega_s[[s]] <- FA.s$loadings
    Psi_s[[s]] <- diag(FA.s$uniq)
    Psi_s1[[s]] <-  diag(1/diag(Psi_s[[s]]))
  }
  ######E-step
  out <- exp_values(Phi, Omega_s, Psi_s, Psi_s1, cov_s, getdet = TRUE)
  Sig_s1 <- out$Sig_s1; ds_s= out$ds_s;
  l_stop0 <- 0
  lm1 <- 0
  l0 <- loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
  for (i in (1:nIt))
  {
    ###########CM1 ---------------------------------------------------------------------------------------

    ######expected values
    out <- exp_values(Phi, Omega_s, Psi_s, Psi_s1, cov_s)
    Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs; Tfcsfcs <- out$Tfcsfcs; Tfcsfs <- out$Tfcsfs
    ######update  of Phi_s
    Psi_new <- list()
    Psi_new1 <- list()

    for(s in 1:S){
      Psi_new[[s]]  <- diag(cov_s[[s]] + Omega_s[[s]] %*% Tfsfs[[s]] %*% t(Omega_s[[s]]) -  2*Txsfs[[s]] %*% t(Omega_s[[s]]) )
      Psi_new[[s]] <- diag(Psi_new[[s]])
      ##########inverse
      Psi_new1[[s]] <- diag(1/diag(Psi_new[[s]]))
    }

    ###########CM2 ---------------------------------------------------------------------------------------

    ######expected values
    out<- exp_values(Phi, Omega_s, Psi_new, Psi_new1, cov_s)
    Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs;
    Tfcsfcs <- out$Tfcsfcs; Tfcsfs <- out$Tfcsfs

    ######update of Phi: not needed


    ########CM3 ---------------------------------------------------------------------------------------

    ######expected values
    out <- exp_values(Phi, Omega_s, Psi_new, Psi_new1, cov_s)
    Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs;
    Tfcsfcs <-  out$Tfcsfcs; Tfcsfs <- out$Tfcsfs

    ######update of Phi
    Omega_new <- list()
    for(s in 1:S) {
      Omega_new[[s]] <- matrix((Txsfs[[s]] %*% solve(Tfsfs[[s]])), p, tot_s[s])
      Omega_new[[s]][upper.tri(Omega_new[[s]])] <- 0
    }

    ###########stopping rule
    out <- exp_values(Phi, Omega_new, Psi_new, Psi_new1, cov_s, getdet = TRUE)

     Sig_s1 <- out$Sig_s1
    ds_s <- out$ds_s
    l1 <- loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
    a <- (l1 - l0)/ (l0-lm1)
    l_stop <- lm1 + (1/ (1-a)) * (l0-lm1)
    l0 <- loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
    if((trace) & (i %% 100 == 0))  cat("i=", i, "Criterion for convergence ", abs(l_stop-l_stop0), "\n")
    if( (abs(l_stop-l_stop0)<tol) & i > 1 & l_stop != Inf) break
    Psi_s <- Psi_new
    Omega_s <- Omega_new
    Psi_s1 <- Psi_new1
    lm1 <- l0
    l0 <- l1
    l_stop0 <- l_stop
  }
  ############return output
  npar <- p * S + sum(tot_s * (p - (tot_s - 1) / 2))
  n_tot <- sum(n_s)
  AIC <- -2 * l1 + npar * 2
  BIC <- -2 * l1 + npar * log(n_tot)
  res <- list(Omega_s = Omega_s, Psi_s = Psi_s, loglik = l1, AIC = AIC,
              BIC = BIC, npar=npar, iter = i)
  return(res)
}


#' Immune System Data
#
#'
#'
#' A data set used by De Vito et al. (2016) as an illustrative example.
#'
#' @examples
#' \dontrun{
#' #The commands below show how the dataset was obtained from libraries on the Bioconductor repository.
#' source("http://bioconductor.org/biocLite.R")
#' biocLite(c("limma", "curatedOvarianData", "RTCGAToolbox"), suppressUpdates=TRUE)
#' library(curatedOvarianData)
#' library(RTCGAToolbox)
#' data(package="curatedOvarianData")
#' data(GSE20565_eset)
#' data(GSE9891_eset)
#' im_response <- c(
#'  #"TP53",
#'  #first pathwayy
#'  "TUBA3C","PRKACG","FGF6","FGF23","FGF22","FGF20","ASB4","TUBB1","LAT","ULBP1","NCR1",
#'  "SIGLEC5","CD160","KLRD1","NCR3","TRIM9","FGF18","ICOSLG",
#'  #second pathway
#'  "MYH2","C9","MBL2",
#'  "GRIN2B","POLR2F","CSF2","IL5","CRP","C8A","SPTA1","GRIN2A","CCR6","FGA","LBP","DUSP9","FCN2","PRKCG","ADCY8","IL5RA","GRIN1","C8B",
#'  #third pathway
#'  "GH2","TNFSF18","GRIN2D","FGB","PRL","SPTBN5","CD70","FGG","RASGRF1","IFNG","SPTBN4","TRIM10","ACTN2","LTA","TNFSF11","GRIN2C","CAMK2B","SPTB","IL1A","TNFRSF13B","ITGA2B","CAMK2A","TRIM31","EREG")
#' GSE20_eset <- t(as.matrix(GSE20565_eset))
#' GSE98_eset <- t(as.matrix(GSE9891_eset))
#' leD <- length(im_response)
#' mat1 <- matrix(sapply(1:leD, function(i) which(colnames(GSE20_eset)==im_response[i])))
#' mat2 <- matrix(sapply(1:leD, function(i) which(colnames(GSE98_eset)==im_response[i])))
#' GSE20 <- matrix(c(GSE20_eset[,mat1]), 140, leD)
#' GSE98 <-matrix(c(GSE98_eset[,mat2]), 285, leD)
#' data_immune <- list(GSE98, GSE20)}
"data_immune"
