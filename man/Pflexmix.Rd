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
  SEMRflexmix extend a framework for finite
  mixtures of regression models with nonparametic component from FlexMix based on R package 'flexmix'. Parameter estimation is performed using
  the EM algorithm, all the step are implemented by \code{Pflexmix}. The algorithm detial can be found in article \code{Sequential Estimation for Mixture of Regression Models in Dealing with
Heterogeneous Population}
}
\usage{
Pflexmix(formula, data = list(), cluster = NULL, 
        model = NULL, concomitant = NULL, control = NULL,
        weights = NULL)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit. The
    general form is \code{y~x} where \code{y} is the response,
    \code{x} the set of predictors for repeated measurements.}
  \item{data}{An optional data frame containing the variables in the model.}
  \item{cluster}{A matrix with \code{2} columns of initial
    cluster membership probabilities for each observation; or a factor
    or integer vector with the initial cluster
    assignments of observations at the start of the EM
    algorithm. Default is random assignment into \code{2}
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
}

\value{
  Returns an object of class \code{Pflexmix}. The containing is the same as an object of class \code{flexmix} except \code{SIGMA} 
which present the variance estimation of the parameters estiomation.
}
\author{Cyan Yu}
\references{
  You, Na; Yu, Qingyun; Wang, Xueqin; Dai, Hongsheng. Sequential Estimation for Mixture of Regression Models in Dealing with
Heterogeneous Population. Submitted.
}

\examples{
data("CCLE_Sorafenib", package = "SEMRflexmix")

y = CCLE_Sorafenib$ActArea
x = CCLE_Sorafenib$regression.covariate
z = CCLE_Sorafenib$concomitant.variable
DF <-data.frame(y = y, x, z)

Model <- FLXMRParticalglm(as.formula(paste0('~',paste(colnames(x),collapse = '+'))),Wupdate = 2)

conModel <- FLXPmultinom(as.formula(paste0('~',paste(colnames(z),collapse = '+'))))

CCLE_Pflx <- Pflexmix(y~1,model=Model,concomitant = conModel,data=DF)
}
\keyword{regression}
\keyword{cluster}
