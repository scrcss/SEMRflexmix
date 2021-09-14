%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: Pflexmix.Rd 5115 2017-04-07 08:18:13Z gruen $
%
\name{Pflexmix}
\alias{Pflexmix}
\alias{Pflexmix,formula,ANY,ANY,ANY,missing-method}
\alias{Pflexmix,formula,ANY,ANY,ANY,list-method}
\alias{Pflexmix,formula,ANY,ANY,ANY,FLXM-method}
\title{Flexible Mixture Modeling With Nonparametric Component}
\description{
  SEMR extend a framework for finite
  mixtures of regression models with nonparametic component from FlexMix. Parameter estimation is performed using
  the EM algorithm: the E-step is implemented by \code{flexmix}, while the M-step is . 
}
\usage{
Pflexmix(formula, data = list(), k = NULL, cluster = NULL, 
        model = NULL, concomitant = NULL, control = NULL,
        weights = NULL)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit. The
    general form is \code{y~x} where \code{y} is the response,
    \code{x} the set of predictors for repeated measurements.}
  \item{data}{An optional data frame containing the variables in the model.}
  \item{k}{Number of clusters (not needed if \code{cluster} is specified).}
  \item{cluster}{Either a matrix with \code{k} columns of initial
    cluster membership probabilities for each observation; or a factor
    or integer vector with the initial cluster
    assignments of observations at the start of the EM
    algorithm. Default is random assignment into \code{k}
    clusters.}
  \item{weights}{An optional vector of replication weights to be used in
    the fitting process.  Should be \code{NULL}, a vector or a
    formula.}
  \item{model}{Object of class \code{FLXM} or list of \code{FLXM}
    objects. Default is the object returned by calling
    \code{\link{FLXMRParticalglm}()}.}
  \item{concomitant}{Object of class \code{FLXP}. Default is the
    object returned by calling \code{\link{FLXPconstant}}.}
  \item{control}{Object of class \code{FLXcontrol} or a named list.}
  \item{object}{Object of class \code{Pflexmix}.}
  \item{\dots}{Currently not used.}
}

\details{
  FlexMix models are described by objects of class \code{FLXM},
  which in turn are created by driver functions like
  \code{\link{FLXMRParticalglm}}. Multivariate
  responses with independent components can be specified using a
  list of \code{FLXM} objects.
}
\value{
  Returns an object of class \code{Pflexmix}.
}
\author{Cyan Yu}
\references{
  Friedrich Leisch. FlexMix: A general framework for finite mixture
  models and latent class regression in R. \emph{Journal of Statistical
    Software}, \bold{11}(8), 2004. doi:10.18637/jss.v011.i08

  Bettina Gruen and Friedrich Leisch. Fitting finite mixtures of
  generalized linear regressions in R. Computational Statistics & Data
  Analysis, 51(11), 5247-5252, 2007. doi:10.1016/j.csda.2006.08.014

  Bettina Gruen and Friedrich Leisch. FlexMix Version 2: Finite
  mixtures with concomitant variables and varying and constant
  parameters Journal of Statistical Software, 28(4), 1-35, 2008.
  doi:10.18637/jss.v028.i04
  
}

\examples{
data("CCLE_Sorafenib", package = "SEMRflexmix")

y = CCLE_Sorafenib$ActArea
x = CCLE_Sorafenib$regression.covariate
z = CCLE_Sorafenib$concomitant.variable
DF <-data.frame(y = y, x, z)

Model <- FLXMRParticalglm(as.formula(paste0('~',paste(colnames(x),collapse = '+'))),Wupdate = 2)

conModel <- FLXPmultinom(as.formula(paste0('~',paste(colnames(z),collapse = '+'))))

CCLE_Pflx <- Pflexmix(y~1,k = 2,model=Model,concomitant = conModel,data=DF)
}
\keyword{regression}
\keyword{cluster}
