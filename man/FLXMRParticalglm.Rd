%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: FLXglm.Rd 5115 2017-04-07 08:18:13Z gruen $
%
\name{FLXMRParticalglm}
\alias{FLXMRParticalglm}
\title{FlexMix Interface to Generalized Linear Models With Nonparametric Component}
\description{
  This is the main driver for FlexMix interfacing the \code{\link{glm}}
  family of models.
}
\usage{
FLXMRParticalglm(formula = . ~ .,
         family = c("gaussian"),
         Wupdate = 1)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{Pflexmix}} using
    \code{\link{update.formula}}. Default is to use the original
    \code{\link{Pflexmix}} model
    formula.}
  \item{family}{A character string naming a \code{\link{glm}}
    family function. But in SEMR framework only Gussian have been extend for now.}
   \item{Wupdate}{SEMR support two way to update the parameters in the Nonparametric Comonent in M-step.}
}
\details{
  See \code{\link{Pflexmix}} for examples.
}
\value{
  Returns an object of class \code{FLXMRParticalglm}.
}
\author{Cyan Yu}
\references{
  Friedrich Leisch. FlexMix: A general framework for finite mixture
  models and latent class regression in R. \emph{Journal of Statistical
  Software}, \bold{11}(8), 2004. doi:10.18637/jss.v011.i08
}
\seealso{\code{\link{Pflexmix}}, \code{\link{glm}}}
\keyword{regression}
\keyword{models}
