setClass("stepPFlexmix",
         representation(models="list",
                        k="integer",
                        nrep="integer",
                        logLiks="matrix",
                        call="call"))



stepPFlexmix <- function(..., nrep=3, verbose=TRUE, drop=TRUE,
                         unique=FALSE)
{
  MYCALL <- match.call()
  MYCALL1 <- MYCALL

  k = 2
  bestFlexmix <- function(...)
  {
    z = new("Pflexmix", logLik=-Inf)
    logLiks = rep(NA, length.out = nrep)
    for(m in seq_len(nrep)){
      if(verbose) cat(" *")
      x = try(Pflexmix(...))
      if (!is(x, "try-error")) {
        logLiks[m] <- logLik(x)
        if(logLik(x) > logLik(z))
          z = x
      }
    }
    return(list(z = z, logLiks = logLiks))
  }

  z = list()
  if(is.null(k)){
    RET = bestFlexmix(...)
    z[[1]] <- RET$z
    logLiks <- as.matrix(RET$logLiks)
    z[[1]]@call <- MYCALL
    z[[1]]@control@nrep <- nrep
    names(z) <- as.character(z[[1]]@k)
    if(verbose) cat("\n")
  }
  else{
    k = as.integer(k)
    logLiks <- matrix(nrow = length(k), ncol = nrep)
    for(n in seq_along(k)){
      ns <- as.character(k[n])
      if(verbose) cat(k[n], ":")
      RET <- bestFlexmix(...)
      z[[ns]] = RET$z
      logLiks[n,] <- RET$logLiks
      MYCALL1[["k"]] <- as.numeric(k[n])
      z[[ns]]@call <- MYCALL1
      z[[ns]]@control@nrep <- nrep
      if(verbose) cat("\n")
    }
  }
  logLiks <- logLiks[is.finite(sapply(z, logLik)),,drop=FALSE]
  z <- z[is.finite(sapply(z, logLik))]
  rownames(logLiks) <- names(z)
  if (!length(z)) stop("no convergence to a suitable mixture")

  if(drop & (length(z)==1)){
    return(z[[1]])
  }
  else{
    z <- return(new("stepPFlexmix",
                    models=z,
                    k=as.integer(names(z)),
                    nrep=as.integer(nrep),
                    logLiks=logLiks,
                    call=MYCALL))
    if(unique) z <- unique(z)
    return(z)
  }
}

###**********************************************************

setMethod("unique", "stepPFlexmix",
          function(x, incomparables=FALSE, ...)
          {
            z <- list()
            K <- sapply(x@models, function(x) x@k)
            logLiks <- x@logLiks
            keep <- rep(TRUE, nrow(logLiks))

            for(k in sort(unique(K))){
              n <- which(k==K)
              if(length(n)>1){
                l <- sapply(x@models[n], logLik)
                z[as.character(k)] <- x@models[n][which.max(l)]
                keep[n[-which.max(l)]] <- FALSE
              }
              else
                z[as.character(k)] <- x@models[n]
            }
            logLiks <- logLiks[keep,,drop=FALSE]
            rownames(logLiks) <- names(z)
            attr(logLiks, "na.action") <- NULL
            mycall <- x@call
            mycall["unique"] <- TRUE

            return(new("stepPFlexmix",
                       models=z,
                       k=as.integer(names(z)),
                       nrep=x@nrep,
                       logLiks=logLiks,
                       call=mycall))
          })

setMethod("logLik", "stepPFlexmix",
          function(object, ..., k = 2)
          {
            ll <- lapply(object@models, function(x) logLik(x))
            df <- sapply(ll, attr, "df")
            nobs <- sapply(ll, attr, "nobs")
            ll <- unlist(ll)
            attr(ll, "df") <- df
            attr(ll, "nobs") <- nobs
            class(ll) <- "logLik"
            ll
          })


setMethod("nobs", signature(object="Pflexmix"),
          function(object, ...) {
            if (is.null(object@weights)) nrow(object@posterior$scaled) else  sum(object@weights)
          })

setMethod("logLik", signature(object="Pflexmix"),
          function(object, newdata, ...){
            if (missing(newdata)) {
              z <- object@logLik
              attr(z, "df") <- object@df
              attr(z, "nobs") <- nobs(object)
              class(z) <- "logLik"
            } else {
              z <- sum(log(rowSums(posterior(object, newdata = newdata, unscaled = TRUE))))
              attr(z, "df") <- object@df
              attr(z, "nobs") <- nrow(newdata)
              class(z) <- "logLik"
            }
            z
          })

setMethod("ICL", signature(object="Pflexmix"),
          function(object, ...){
            -2 * clogLik(object) + object@df * log(nobs(object))
          })

setMethod("clogLik", signature(object="Pflexmix"),
          function(object, ...){
            first <- if (length(object@group)) groupFirst(object@group) else TRUE
            post <- object@posterior$unscaled[first,,drop=FALSE]
            n <- nrow(post)
            sum(log(post[seq_len(n) + (clusters(object)[first] - 1)*n]))
          })

setMethod("EIC", signature(object="Pflexmix"),
          function(object, ...) {
            first <- if (length(object@group)) groupFirst(object@group) else TRUE
            post <- object@posterior$scaled[first,,drop=FALSE]
            n <- nrow(post)
            lpost <- log(post)
            if (any(is.infinite(lpost))) lpost[is.infinite(lpost)] <- -10^3
            1 + sum(post * lpost)/(n * log(object@k))
          })
