%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: dmft.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{CCLE_Sorafenib}
\alias{CCLE_Sorafenib}
\docType{data}
\title{CCLE Sorafenib Data}
\description{
  The genetic profiles of 473 cell lines to explore the relationship between the molecular characterizations of patients
and their pharmacologic response, supported by the Cancer Cell Line EncyclopFedia (CCLE) project
(https://portals.broadinstitute.org/ccle)
}
\usage{data("CCLE_Sorafenib")}
\format{
  A data list with 473 observations on the following  3 componets: the response variable "ActArea", 19 regression variables and 5 concomitant variables.
  \describe{
    \item{ActArea}{The effect of drug Sorafenib on 473 cell lines.}
    \item{regression.covariate}{19 genetic profiles selected by feature selection procedure from tens of thousands of genes for drug sensitivity prediction of heterogeneous
population according to \code{Drug sensitivity prediction with high-dimensional
mixture regression.}.}
    \item{concomitant.variable}{Took the
clusters identified by flexmix without concomitant to label the samples into subgroups, and
then apply the screening and random forest method to select the concomitant variables.}
  }
}
  
\details{
   The aim of the study was to to explore the relationship between the molecular characterizations of patients
and their pharmacologic response. Because of the potential heterogeneity of patients, we selected some varible to predict the potential heterogeneity besides of the regression.
}

\source{
Li, Q., Shi, R., and Liang, F. (2019). Drug sensitivity prediction with high-dimensional
mixture regression.  \emph{PLOS one}, 14(2):e0212108.

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
\keyword{datasets}
