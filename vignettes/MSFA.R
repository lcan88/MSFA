## ----help, echo = TRUE, results = TRUE, tidy = TRUE----------------------
library(MSFA)
data(data_immune)
help(data_immune)

## ---- starting values, messages = FALSE----------------------------------
start_value <- start_msfa(X_s = data_immune, k = 4, j_s = c(2, 3))

## ---- get estimate, results = FALSE, cache = TRUE, cache.path="cache/"----
mle <-  ecm_msfa(data_immune, start_value, trace = FALSE)

## ---- heatmap, fig.show = 'hold', fig.width = 7.5, fig.height = 6.5, message = FALSE----
library(gplots)
heatmap.2(mle$Phi,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none', density.info="none", col=heat.colors(256))

