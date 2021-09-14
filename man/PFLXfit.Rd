%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: FLXfit.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{PFLXfit}
\alias{PFLXfit}
\alias{PFLXfit,list-method}
\title{Fitter Function for FlexMix Models With Nonparametric Component}
\description{
  This is the basic computing engine called by \code{\link{Pflexmix}}, it
  should usually not be used directly.
}
\usage{
PFLXfit(model, concomitant, control, postunscaled = NULL, groups,
       weights)
}
\arguments{
  \item{\dots}{The same as \code{\link{FLXfit}}.}
}
\value{
  Returns an object of class \code{Pflexmix}.
}
\author{Cyan Yu}
\seealso{\code{\link{Pflexmix}}, \code{\link{Pflexmix-class}}}
\keyword{regression}
\keyword{cluster}
