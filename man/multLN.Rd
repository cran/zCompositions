\name{multLN}
\alias{multLN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiplicative lognormal replacement
}
\description{
This function implements the multiplicative lognormal univariate replacement of left-censored values (e.g. values below detection limit, rounded zeros) in compositional data sets.
}
\usage{
multLN(X, label = NULL, dl = NULL, rob = FALSE, random = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
}
  \item{label}{Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote unobserved left-censored values in \code{X}.}
  \item{dl}{
Numeric vector of detection limits/thresholds (one per component/column, use e.g. 0 if no threshold for a particular one). These must be given on the same scale as \code{X}.
}
  \item{rob}{
Logical value. \code{FALSE} provides maximum-likelihood estimates of model parameters (default), \code{TRUE} provides robust estimates (see \code{NADA} package for details).
}
  \item{random}{
Logical value. \code{FALSE} imputes using the estimated geometric mean of the values < threshold (default). \code{TRUE} imputes using random values below the limit of detection.
}
}
\details{
This function depends on package \code{NADA} to produce model parameter estimates (either maximum likelihood or robust regression on order statistics) from log-transformed censored data. It produces a data set on the same scale as the input data set. If \code{X} is not closed to a constant sum, then the results are adjusted to provide a compositionally equivalent data set, expressed in the original scale, which leaves the absolute values of the observed components unaltered. Note that, following Mateu-Figueras et al. (2013), a normal distribution representation on the own positive real line is assumed for the data instead of the standard lognormal model  based on the Lebesgue measure. %Hence, estimates are obtained from the log-transformed data (space of coordinates) and transformed back by exponentiation into the positive real line.
}
\value{
A \code{\link{data.frame}} object containing the replaced compositional data set.
}
\references{
Mateu-Figueras G, Pawlowsky-Glahn V, Egozcue JJ. The normal distribution in some constrained sample spaces. SORT 2013; 37(1): 29-56.

Palarea-Albaladejo J, Martin-Fernandez JA. Values below detection limit in compositional chemical data. Analytica Chimica Acta 2013; 764: 32-43. DOI: \url{http://dx.doi.org/10.1016/j.aca.2012.12.029}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{zPatterns}}, \code{\link{multRepl}}, \code{\link{lrEM}}, \code{\link{lrDA}}, \code{\link{cmultRepl}}
}
\examples{
# Data set closed to 100 (percentages, common dl = 1%)
X <- matrix(c(26.91,8.08,12.59,31.58,6.45,14.39,
              39.73,26.20,0.00,15.22,6.80,12.05,
              10.76,31.36,7.10,12.74,31.34,6.70,
              10.85,46.40,31.89,10.86,0.00,0.00,
              7.57,11.35,30.24,6.39,13.65,30.80,
              38.09,7.62,23.68,9.70,20.91,0.00,
              27.67,7.15,13.05,32.04,6.54,13.55,
              44.41,15.04,7.95,0.00,10.82,21.78,
              11.50,30.33,6.85,13.92,30.82,6.58,
              19.04,42.59,0.00,38.37,0.00,0.00),byrow=TRUE,ncol=6)
              
X_multLN <- multLN(X,label=0,dl=rep(1,6))

# Non-closed compositional data set
data(LPdata) # data (ppm/micrograms per gram)
dl <- c(2,1,0,0,2,0,6,1,0.6,1,1,0,0,632,10) # limits of detection (0 for no limit)

# Using MV for parameter estimation
LPdata_multLN <- multLN(LPdata,label=0,dl=dl) 
# For comparison
LPdata[30:35,1:10]
round(LPdata_multLN[30:35,1:10],1)

# Using ROS for parameter estimation
LPdata_multLNrob <- multLN(LPdata,label=0,dl=dl,rob=TRUE)
round(LPdata_multLNrob[30:35,1:10],1)

# Using random values < dl
LPdata_multRLN <- multLN(LPdata,label=0,dl=dl,random=TRUE)
round(LPdata_multRLN[30:35,1:10],1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
