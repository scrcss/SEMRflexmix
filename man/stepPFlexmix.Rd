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
  Runs Pflexmix repeatedly for different numbers of components and returns
  the maximum likelihood solution for each.
}
\usage{
stepPFlexmix(..., k = NULL, nrep = 3, verbose = TRUE, drop = TRUE,
            unique = FALSE)

\S4method{unique}{stepPFlexmix}(x, incomparables = FALSE, ...)
}
\arguments{
  \item{\dots}{Passed to \code{\link{Pflexmix}} }
  \item{k}{A vector of integers passed in turn to the \code{k} argument
    of \code{\link{Pflexmix}}.}
  \item{nrep}{For each value of \code{k} run \code{\link{Pflexmix}}
    \code{nrep} times and keep only the solution with maximum
    likelihood. If \code{nrep} is set for the long run, it is ignored,
    because the EM algorithm is deterministic using the best solution
    discovered in the short runs for initialization.}
  \item{verbose}{If \code{TRUE}, show progress information during
    computations.}
  \item{drop}{If \code{TRUE} and \code{k} is of length 1, then a single
    flexmix object is returned instead of a \code{"stepPFlexmix"}
    object.}
  \item{unique}{If \code{TRUE}, then \code{unique()} is called on the
    result, see below.}
  \item{x}{An object of class \code{"stepFlexmix"}.}
  \item{incomparables}{A vector of values that cannot be
          compared. Currently, \code{FALSE} is the only possible value,
          meaning that all values can be compared.}
}
\value{
  An object of class \code{"stepPFlexmix"} containing the best models
  with respect to the log likelihood for the different number of
  components in a slot if \code{length(k)>1}, else directly an object of
  class \code{"Pflexmix"}.

  If \code{unique = FALSE}, then the resulting object contains one
  model per element of \code{k} (which is the number of clusters the EM
  algorithm started with). If \code{unique = TRUE}, then the result
  is resorted according to the number of clusters contained in the
  fitted models (which may be less than the number with which the EM
  algorithm started), and only the maximum likelihood solution for each
  number of fitted clusters is kept. This operation can also be done
  manually by calling \code{unique()} on objects of class
  \code{"stepPFlexmix"}. 
}
\author{Cyan Yu}
\references{
  Friedrich Leisch. FlexMix: A general framework for finite mixture
  models and latent class regression in R. \emph{Journal of Statistical
  Software}, \bold{11}(8), 2004. doi:10.18637/jss.v011.i08

  Christophe Biernacki, Gilles Celeux and Gerard Govaert. Choosing
  starting values for the EM algorithm for getting the highest
  likelihood in multivariate Gaussian mixture models. \emph{Computational
  Statistics & Data Analysis}, \bold{41}(3--4), 561--575, 2003.

  Theresa Scharl, Bettina Gruen and Friedrch Leisch.  Mixtures of
  regression models for time-course gene expression data: Evaluation of
  initialization and random effects. \emph{Bioinformatics},
  \bold{26}(3), 370--377, 2010.
}

\examples{
data("CCLE_Sorafenib", package = "SEMRflexmix")

y = CCLE_Sorafenib$ActArea
x = CCLE_Sorafenib$regression.covariate
z = CCLE_Sorafenib$concomitant.variable
DF <-data.frame(y = y, x, z)

Model <- FLXMRParticalglm(as.formula(paste0('~',paste(colnames(x),collapse = '+'))),Wupdate = 2)

conModel <- FLXPmultinom(as.formula(paste0('~',paste(colnames(z),collapse = '+'))))

Pflx <- stepPFlexmix(y~1,k = 2,model=Model,concomitant = conModel,data=DF, nrep = 3, unique = T)
}
\keyword{cluster}
\keyword{regression}
