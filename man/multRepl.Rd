\name{multRepl}
\alias{multRepl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiplicative simple replacement
}
\description{
This function implements the multiplicative simple replacement of left-censored values (e.g. values below detection limit, rounded zeros) in compositional data sets.
}
\usage{
multRepl(X, label = NULL, dl = NULL, delta = 0.65)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Compositional vector (\code{numeric} class) or data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
}
  \item{label}{Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote unobserved left-censored values in \code{X}.}
  \item{dl}{
Numeric vector of detection limits/thresholds (one per component/column, use e.g. 0 if no threshold for a particular one). These must be given on the same scale as \code{X}.
}
  \item{delta}{
Delta parameter (fraction of the threshold) used in proportions (default = 0.65).
}
}
\details{
This function produces a replaced data set on the same scale as the input data set. If \code{X} is not closed to a constant sum, then the results are adjusted to provide a compositionally equivalent data set, expressed in the original scale, which leaves the absolute values of the observed components unaltered.
}
\value{
A \code{\link{data.frame}} object containing the replaced compositional vector or data set.
}
\references{
Martin-Fernandez JA, Barcelo-Vidal C, Pawlowsky-Glahn V. Dealing with zeros and missing values in compositional data sets using nonparametric imputation. Mathematical Geology 2003; 35: 253-78.

Palarea-Albaladejo J, Martin-Fernandez JA. Values below detection limit in compositional chemical data. Analytica Chimica Acta 2013; 764: 32-43. DOI: http://dx.doi.org/10.1016/j.aca.2012.12.029.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{zPatterns}}, \code{\link{multLN}}, \code{\link{lrEM}}, \code{\link{lrDA}}, \code{\link{cmultRepl}}
}
\examples{
# A compositional vector (NA indicates nondetect)
y <- c(0.6,NA,0.25,0.03,0.12,NA)
dl <- c(0,0.01,0,0,0,0.005)
# Using the default delta = 0.65
yr <- multRepl(y,label=NA,dl=dl)
round(yr,4)

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
              
X_multRepl <- multRepl(X,label=0,dl=rep(1,6))

# Non-closed compositional data set
data(LPdata) # data (ppm/micrograms per gram)
dl <- c(2,1,0,0,2,0,6,1,0.6,1,1,0,0,632,10) # limits of detection (0 for no limit)
LPdata_multRepl <- multRepl(LPdata,label=0,dl=dl)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
