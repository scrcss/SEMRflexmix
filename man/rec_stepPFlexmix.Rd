%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: stepFlexmix.Rd 5115 2017-04-07 08:18:13Z gruen $
%
\name{rec_stepPFlexmix}
\alias{rec_stepPFlexmix}
\title{Run Sequential Estimation for Mixture of Regression Models}
\description{
  Runs stepPFlexmix Sequentially for different weights of samples and returns the maximum likelihood solution for each. 
That is a sequential analysis procedure based on a mixture model with two components. 
}
\usage{
rec_stepPFlexmix(x, y, z, min.comp=10, max.coef = 10, try.times = 5, Wupdate = 1)
}
\arguments{
  \item{x}{A matrix or dataframe of regression varibales.}
  \item{y}{A vector presents the respose in the regression modle.}
  \item{z}{A matrix or dataframe of concomitant variables.}
  \item{min.comp}{The minimum of the sample size that remains in the non-parametric part to stop the sequential algorithm.}  
  \item{max.coef}{The maximum parameter estimation that we can accept otherwise we will restart the algorithm for we think the algorithm has reach at a boundary point.}
  \item{try.times}{If the algortithm reach at a boundary point we will restart the alogorithm no more than \code{try.times} times.}
  \item{Wupdate}{We provide 2 ways to update the M-step by controling the code by Wupade equal to 1 or 2.}
}
\details{
  Function rec_stepPFlexmix is a sequential analysis procedure based on a mixture model with two components.  Every parition is  a \code{\link{stepPFlexmix}}.
}
\value{
  A list containing at least the following components:
  \item{beta}{A matrix with (k -1) columns of the regression varible parameter estimation during k times sequential partitioning.}
  \item{sigma}{A vector with (k - 1) estimations of the residual standard deviation during k times sequential partitioning.}
  \item{v}{A matrix with k columns of the concomitant varible parameter estimation during k times sequential partitioning while the concomitant varible existed, otherwise, 
                 \code{v} is a vector with prior probabilities of clusters..}
  \item{U}{A matrix with one row per observation and one column per cluster..}
  \item{cluster}{Cluster assignments of observations.}
  \item{k}{Number of clusters after EM. The k-th is the non-parameteric components.}
  \item{err}{A bool varible to show whether the sequential partitioning has a case where it don't converge or not.}
  \item{singular}{A vector with k numbers to show whether the algorithm has reached at a boundary points during k times sequential partitioning. 
          \code{singular = 1} shows the entire sample was assigned to parameteric part, otherwise, \code{singular = 2} shows the entire sample was assigned to non-parametric part.}
}
\author{Cyan Yu}
\keyword{cluster}
\keyword{regression}

\examples{
data("CCLE_Sorafenib", package = "SEMRflexmix")

y = CCLE_Sorafenib$ActArea
x = CCLE_Sorafenib$regression.covariate
z = CCLE_Sorafenib$concomitant.variable

res = rec_stepPFlexmix(x, y, z)
}
\keyword{regression}
\keyword{cluster}
