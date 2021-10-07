%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: stepFlexmix.Rd 5115 2017-04-07 08:18:13Z gruen $
%
\name{stepPFlexmix}
\alias{stepPFlexmix}
\alias{stepPFlexmix-class}
\alias{unique,stepPFlexmix-method}
\title{Run FlexMix With Nonparametric Component Repeatedly}
\description{
  Runs Pflexmix repeatedly for different initial parameters values and returns
  the maximum likelihood solution.
}
\usage{
stepPFlexmix(..., nrep = 3, verbose = TRUE, drop = TRUE,
            unique = FALSE)

\S4method{unique}{stepPFlexmix}(x, incomparables = FALSE, ...)
}
\arguments{
  \item{\dots}{Passed to \code{\link{Pflexmix}} }
  \item{nrep}{Run \code{\link{Pflexmix}}
    \code{nrep} times and keep only the solution with maximum
    likelihood. If \code{nrep} is set for the long run, it is ignored,
    because the EM algorithm is deterministic using the best solution
    discovered in the short runs for initialization.}
  \item{verbose}{If \code{TRUE}, show progress information during
    computations.}
  \item{drop}{If \code{TRUE}, then a single
    flexmix object is returned instead of a \code{"stepPFlexmix"}
    object.}
  \item{unique}{If \code{TRUE}, then \code{unique()} is called on the
    result, see below.}
  \item{x}{An object of class \code{"stepPFlexmix"}.}
  \item{incomparables}{A vector of values that cannot be
          compared. Currently, \code{FALSE} is the only possible value,
          meaning that all values can be compared.}
}
\value{
  An object of class \code{"stepPFlexmix"} containing the best models
  with respect to the log likelihood for the different intital parameters' values, else directly an object of
  class \code{"Pflexmix"}.

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

Pflx <- stepPFlexmix(y~1,model=Model,concomitant = conModel,data=DF, nrep = 3, unique = T)
}
\keyword{cluster}
\keyword{regression}
