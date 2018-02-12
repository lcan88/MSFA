# MSFA
MSFA

---
 ## Bayesian Analysis of a MSFA model
author: "Roberta de Vito, Ruggero Bellio"


The following example illustrates a Bayesian analysis of a MSFA model.
Although the methodology has been developed targeting the $p>n$ case, for the sake
of simplicity we illustrate the analysis of the same data set employed for
maximum likelihood estimation. The data set is 
 available in the Bioconductor repository (www.bioconductor.org). 

### Getting the data
Some pre-processing is required to get the data into a form suitable for the
analysis. This was already done, and the resulting data frame is saved into the
`data_immune` object. The commands that were used to form it are included
in the help file for the data object.

```{r help, echo = TRUE, results = TRUE, tidy = TRUE}
library(MSFA)
data(data_immune)
help(data_immune)
```


## Sampling from the posterior distribution 
We fist estimate a model with a somewhat large dimension of the various loading matrices,
so we set a dimension 10 for both the common factor loadings and the study-specific loadings. In order to get reproducible results, we set the random seed.

```{r posterior, messages = FALSE}
set.seed(1971)
out10_1010 <- sp_msfa(data_immune,  k = 10,  j_s = c(10, 10), trace = FALSE)
```

We take as the estimated $\Sigma_\Phi$   the posterior median

```{r Phi }
p <- ncol(data_immune[[1]])
nrun <- dim(out10_1010$Phi)[3]
SigmaPhi <-  SigmaLambda1 <- SigmaLambda2 <- array(0, dim=c(p, p, nrun))
for(j in 1:nrun)
{
  SigmaPhi[,,j] <- tcrossprod(out10_1010$Phi[,,j]) 
  SigmaLambda1[,,j] <- tcrossprod(out10_1010$Lambda[[1]][,,j]) 
  SigmaLambda2[,,j] <- tcrossprod(out10_1010$Lambda[[2]][,,j]) 
}
SigmaPhi <- apply(SigmaPhi, c(1, 2), median)
SigmaLambda1 <- apply(SigmaLambda1, c(1, 2), median)
SigmaLambda2 <- apply(SigmaLambda2, c(1, 2), median)
Phi <- apply(out10_1010$Phi, c(1, 2), median)
```

## Choice of the number of factors
Then we proceed to the choice of the number of common latent factors.

```{r SigmaPhi, fig.width=4.5, fig.height=4.5}
plot(sp_eigen(SigmaPhi), pch = 16)
abline(h = 0.05, col = 2)
```

We note that 5 factors are above the $5\%$ threshold, so we choose $K=5$. We proceed in a similar way for
the two study-specific loading matrices:
```{r SigmaLambda, fig.width=6.5, fig.height=4.5}
par(mfrow=c(1, 2))
plot(sp_eigen(SigmaLambda1), pch=16)
abline(h = 0.05, col = 2)
plot(sp_eigen(SigmaLambda2), pch=16)
abline(h = 0.05, col = 2)
```

We end up with dimensions 3 and 4 for the study-specific factor loadings.

##  OP prostprocessing
We post-process the estimated loading matrix by the OP procedure
```{r OP}
Phi_OP10 <- sp_OP(out10_1010$Phi[,1:5,], itermax = 10)
```

For larger data size, we recommend to reduce the output level in the call to
```sp_msfa```.

