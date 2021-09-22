# ###**********************************************************
# # add slot SIGMA compare to flexmix
setClass("Pflexmix",
         representation(posterior="ANY",
                        weights="ANY",
                        iter="numeric",
                        cluster="integer",
                        logLik="numeric",
                        df="numeric",
                        control="FLXcontrol",
                        group="factor",
                        size="integer",
                        converged="logical",
                        k0="integer",
                        SIGMA = "matrix"),
         prototype(group=(factor(integer(0))),
                   formula=.~.),
         contains="FLXdist")

# Add slot singular to mark the error and Wupdate to choose different EM update algorithm
setClass("FLXMRParticalglm",
         representation(singular = "numeric",
                        Wupdate = "numeric"),
         prototype(singular = 0,
                   Wupdate = 1),
         contains = "FLXMRglm")

# initial function to init a new class FLXMRParticalglm
FLXMRParticalglm <- function(formula = . ~ .,
                             family = c("gaussian"),Wupdate = 1,  ...) {
  family <- match.arg(family)
  new("FLXMRParticalglm", FLXMRglm(formula, family, ...),
      name = paste("FLXMRParticalglm", family, sep=":"), Wupdate = Wupdate)
}

# to re-write some function
setMethod("FLXgetModelmatrix", signature(model="FLXMRParticalglm"),
          function(model, data, formula, lhs=TRUE, ...) {
            model <- callNextMethod(model, data, formula, lhs)
            model
          })


setMethod("FLXremoveComponent", signature(model = "FLXMRParticalglm"),
          function(model, nok, ...){
            model@singular = nok
            model
          })


# m-step contain 3 method control by parameter Wupdate
setMethod("FLXmstep", signature(model = "FLXMRParticalglm"),
          function(model, weights, components, prior) {
            if(model@Wupdate == 1){
              w = (weights[, 1, drop = FALSE]) / sum(weights[, 1, drop = FALSE])
            }
            else if(model@Wupdate == 2){
              N = nrow(model@y)
              
              if(length(components[[1]]@df) == 0)  w = rep(1/N,N)
              else{
                
                logpostunscaled  = matrix(sapply(components, function(x) x@logLik(model@x, model@y)), nrow = nrow(model@y))
                
                postunscaled <- if (nrow(prior) > 1) logpostunscaled + log(prior)
                else sweep(logpostunscaled, 2, log(prior), "+")
                postunscaled <- exp(postunscaled)
                
                # update w
                w = exp(logpostunscaled[, 1, drop = FALSE])/rowSums(postunscaled)
                w = w/sum(w)
              }
            }
            else stop("Unknown W's updating way")
            defineComponent <- function(para) {
              logLik <- function(x, y, ...){
                h = bw.nrd0(y)
                KenerlW = 1/h * dnorm(as.matrix(dist(y)/h), 0, 1)
                log(as.matrix(KenerlW) %*% para)
              }
              
              new("FLXcomponent",
                  parameters=list(w = para),
                  logLik=logLik,
                  df=0)
            }
            comp.1 <- defineComponent(w)
            c(list(comp.1),
              FLXmstep(as(model, "FLXMRglm"), weights[, -1, drop=FALSE], components[2]))
          })

setGeneric("Pflexmix",
           function(formula, data=list(), k=NULL,
                    cluster=NULL, model=NULL, concomitant=NULL, control=NULL,
                    weights = NULL)
             standardGeneric("Pflexmix"))

## The following two methods only fill in and rearrange the model argument
setMethod("Pflexmix",
          signature(formula = "formula", model="missing"),
          function(formula, data=list(), k=NULL, cluster=NULL,
                   model=NULL, concomitant=NULL, control=NULL, weights=NULL)
          {
            mycall = match.call()
            z <- Pflexmix(formula=formula, data=data, k=k, cluster=cluster,
                          model=list(FLXMRglm()), concomitant=concomitant,
                          control=control, weights = weights)
            z@call <- mycall
            z
          })

setMethod("Pflexmix",
          signature(formula = "formula", model="FLXMRParticalglm"),
          function(formula, data=list(), k=NULL, cluster=NULL,
                   model=NULL, concomitant=NULL, control=NULL, weights=NULL)
          {
            mycall = match.call()
            z <- Pflexmix(formula=formula, data=data, k=k, cluster=cluster,
                          model=list(model), concomitant=concomitant,
                          control=control, weights=weights)
            z@call <- mycall
            z
          })


## This is the real thing
setMethod("Pflexmix",
          signature(formula = "formula", model="list"),
          function(formula, data=list(), k=NULL, cluster=NULL,
                   model=NULL, concomitant=NULL, control=NULL, weights=NULL)
          {
            mycall = match.call()
            control = as(control, "FLXcontrol")
            if (!is(concomitant, "FLXP")) concomitant <- FLXPconstant()
            
            groups <- .FLXgetGrouping(formula, data)
            model <- lapply(model, FLXcheckComponent, k, cluster)
            k <- unique(unlist(sapply(model, FLXgetK, k)))
            if (length(k) > 1) stop("number of clusters not specified correctly")
            if (k > 2) stop("the number of clusters must be 2")
            
            model <- lapply(model, FLXgetModelmatrix, data, formula)
            
            groups$groupfirst <-
              if (length(groups$group)) groupFirst(groups$group)
            else rep(TRUE, FLXgetObs(model[[1]]))
            
            if (is(weights, "formula")) {
              weights <- model.frame(weights, data = data, na.action = NULL)[,1]
            }
            
            ## if weights and grouping is specified the weights within each
            ## group need to be the same
            if (!is.null(weights) & length(groups$group)>0) {
              unequal <- tapply(weights, groups$group, function(x) length(unique(x)) > 1)
              if (any(unequal)) stop("identical weights within groups needed")
            }
            
            postunscaled <- initPosteriors(k, cluster, FLXgetObs(model[[1]]), groups)
            
            if (ncol(postunscaled) == 1L)
              concomitant <- FLXPconstant()
            
            concomitant <- FLXgetModelmatrix(concomitant, data = data,
                                             groups = groups)
            
            
            z <- PFLXfit(model=model, concomitant=concomitant, control=control,
                         postunscaled=postunscaled, groups=groups, weights = weights)
            
            z@formula = formula
            z@call = mycall
            z@k0 = as.integer(k)
            z
          })

###**********************************************************
setGeneric("PFLXfit",
           function(model, concomitant, control,
                    postunscaled=NULL, groups, weights)
             standardGeneric("PFLXfit"))

###**********************************************************
setMethod("PFLXfit", signature(model="list"),
          function(model, concomitant, control, postunscaled=NULL, groups, weights)
          {
            ### initialize
            k <- ncol(postunscaled)
            N <- nrow(postunscaled)
            control <- allweighted(model, control, weights)
            if(control@verbose>0)
              cat("Classification:", control@classify, "\n")
            if (control@classify %in% c("SEM", "random")) iter.rm <- 0
            group <- groups$group
            groupfirst <- groups$groupfirst
            if(length(group)>0) postunscaled <- groupPosteriors(postunscaled, group)
            
            logpostunscaled <- log(postunscaled)
            postscaled <- exp(logpostunscaled - log_row_sums(logpostunscaled))
            
            llh <- -Inf
            if (control@classify %in% c("SEM", "random")) llh.max <- -Inf
            converged <- FALSE
            components <- rep(list(rep(list(new("FLXcomponent")), k)), length(model))
            ### EM
            for(iter in seq_len(control@iter.max)) {
              ### M-Step
              postscaled = .FLXgetOK(postscaled, control, weights)
              prior <- if (is.null(weights))
                ungroupPriors(concomitant@fit(concomitant@x, postscaled[groupfirst,,drop=FALSE]),
                              group, groupfirst)
              else ungroupPriors(concomitant@fit(concomitant@x, (postscaled/weights)[groupfirst & weights > 0,,drop=FALSE], weights[groupfirst & weights > 0]),
                                 group, groupfirst)
              
              components <- lapply(seq_along(model), function(i) FLXmstep(model[[i]], postscaled, components[[i]], prior))
              postunscaled <- matrix(0, nrow = N, ncol = k)
              for (n in seq_along(model))
                postunscaled <- postunscaled + FLXdeterminePostunscaled(model[[n]], components[[n]])
              if(length(group)>0)
                postunscaled <- groupPosteriors(postunscaled, group)
              ### E-Step
              ## Code changed thanks to Nicolas Picard
              ## to avoid problems with small likelihoods
              postunscaled <- if (nrow(prior) > 1) postunscaled + log(prior)
              else sweep(postunscaled, 2, log(prior), "+")
              logpostunscaled <- postunscaled
              postunscaled <- exp(postunscaled)
              postscaled <- exp(logpostunscaled - log_row_sums(logpostunscaled))
              ##<FIXME>: wenn eine beobachtung in allen Komonenten extrem
              ## kleine postunscaled-werte hat, ist exp(-postunscaled)
              ## numerisch Null, und damit postscaled NaN
              ## log(rowSums(postunscaled)) ist -Inf
              ##</FIXME>
              if (any(is.nan(postscaled))) {
                index <- which(as.logical(rowSums(is.nan(postscaled))))
                postscaled[index,] <- if(nrow(prior)==1) rep(prior, each = length(index)) else prior[index,]
                postunscaled[index,] <- .Machine$double.xmin
              }
              
              if (any(is.na(postscaled))) {
                # browser()
                model <- lapply(model,function(x) {
                  x@singular = 3
                  x
                })
                break
              }
              ### check convergence
              # browser()
              nok <- if (nrow(postscaled) == 1) which(postscaled < control@minprior) else {
                if (is.null(weights)) which(colSums(postscaled[groupfirst,,drop=FALSE]) < control@minprior)
                else which(colSums(postscaled[groupfirst,] * weights[groupfirst]) < control@minprior)
              }
              
              if(length(nok)){
                # browser()
                model <- lapply(model,function(x) {
                  x@singular = nok
                  x
                })
                break
              }
              
              llh.old <- llh
              llh <- if (is.null(weights)) sum(log_row_sums(logpostunscaled[groupfirst,,drop=FALSE]))
              else sum(log_row_sums(logpostunscaled[groupfirst,,drop=FALSE])*weights[groupfirst])
              if(is.na(llh) | is.infinite(llh))
                # browser()
                if(is.na(llh) | is.infinite(llh))
                  stop(paste(formatC(iter, width=4),
                             "Log-likelihood:", llh))
              if (abs(llh-llh.old)/(abs(llh)+0.1) < control@tolerance){
                if(control@verbose>0){
                  printIter(iter, llh)
                  cat("converged\n")
                }
                converged <- TRUE
                break
              }
              if (control@classify=="random") {
                if (llh.max < llh) {
                  components.max <- components
                  prior.max <- prior
                  postscaled.max <- postscaled
                  postunscaled.max <- postunscaled
                  llh.max <- llh
                }
              }
              if(control@verbose && (iter%%control@verbose==0))
                printIter(iter, llh)
            }
            
            # var = T
            # calculate the variance of the paramater estimations
            # browser()
            # single = unlist(lapply(model, function(x) x@singular))
            if((!length(nok)) & (dim(prior)[1] > 1)){
              f_comp <- matrix(0, nrow = N, ncol = k)
              for (n in seq_along(model))
                f_comp <- f_comp + FLXdeterminePostunscaled(model[[n]], components[[n]])
              if(length(group)>0)
                f_comp <- groupPosteriors(f_comp, group)
              f <- if (nrow(prior) > 1) exp(f_comp + log(prior))
              beta = components[[1]][[2]]@parameters$coef
              sigma = components[[1]][[2]]@parameters$sigma
              temp_cal = (model[[1]]@y - model[[1]]@x %*% beta)
              f_diff = cbind(c(f[, 2] * temp_cal / (sigma^2)) * model[[1]]@x,
                             c(f[, 2] * (-1/sigma +temp_cal^2/(sigma^3))))
              
              h = bw.nrd0(model[[1]]@y)
              KenerlW = 1/h * dnorm(as.matrix(dist(model[[1]]@y)/h), 0, 1)
              w = components[[1]][[1]]@parameters$w
              # FI
              # browser()
              zata1 = -f_diff * c(prior[, 1])/c((rowSums(f)^2))
              zata2 =  t(c(w)*KenerlW) - exp(f_comp[, 1])
              p = length(beta)
              Fi <- matrix(0, p + 1, p+1)
              for(i in 1:N){
                theta =rep(0, p + 1)
                for(j in 1:N){
                  if(i != j)
                    theta = theta + (zata1[i, ] * zata2[i, j] + zata1[j, ] * zata2[j, i])/2
                }
                theta = theta/(N-1)
                Fi = Fi +  theta %*% t(theta)
              }
              Fi = Fi/N
              
              # V
              f_hessian = matrix(0, p + 1, p+1)
              V = matrix(0, p+1, p+1)
              S1 = f_diff / rowSums(f)
              for(i in 1:N){
                f_hessian[1:p, 1:p] <- f[i, 2]*(temp_cal[i]^2 / sigma^4 - 1/sigma^2) * (model[[1]]@x[i, ] %*% t(model[[1]]@x[i, ]))
                f_hessian[1:p, p + 1] <- f[i, 2]*(-3*temp_cal[i]/sigma^3 + temp_cal[i]^3/sigma^5)* model[[1]]@x[i, ]
                f_hessian[p+1, 1:p] <- f_hessian[1:p, p + 1]
                
                f_hessian[p+1, p+1] <- f[i, 2] * (2/sigma^2 - 5*temp_cal[i]^2/sigma^4 + temp_cal[i]^4/sigma^6)
                
                V = V + (f_hessian/sum(f[i, ]) - S1[i, ] %*% t(S1[i, ]))
              }
              V = V/N
              
              SIGMA = solve(V) %*% Fi %*% solve(V)
              SIGMA = SIGMA/N
            }else
              SIGMA = matrix(NA, 1, 1)
            
            ### Construct return object
            if (control@classify=="random") {
              components <- components.max
              prior <- prior.max
              postscaled <- postscaled.max
              postunscaled <- postunscaled.max
              llh <- llh.max
              iter <- control@iter.max - iter.rm
            }
            
            components <- lapply(seq_len(k), function(i) lapply(components, function(x) x[[i]]))
            names(components) <- paste("Comp", seq_len(k), sep=".")
            cluster <- max.col(postscaled)
            size <-  if (is.null(weights)) tabulate(cluster, nbins=k) else tabulate(rep(cluster, weights), nbins=k)
            names(size) <- seq_len(k)
            concomitant <- FLXfillConcomitant(concomitant, postscaled[groupfirst,,drop=FALSE], weights[groupfirst])
            df <- concomitant@df(concomitant@x, k) + sum(sapply(components, sapply, slot, "df"))
            control@nrep <- 1
            prior <- if (is.null(weights)) colMeans(postscaled[groupfirst,,drop=FALSE])
            else colSums(postscaled[groupfirst,,drop=FALSE] * weights[groupfirst])/sum(weights[groupfirst])
            
            retval <- new("Pflexmix", model=model, prior=prior,
                          posterior=list(scaled=postscaled,
                                         unscaled=postunscaled),
                          weights = weights,
                          iter=iter, cluster=cluster, size = size,
                          logLik=llh, components=components,
                          concomitant=concomitant,
                          control=control, df=df, group=group, k=as(k, "integer"),
                          converged=converged,
                          SIGMA = SIGMA)
            retval
          })

