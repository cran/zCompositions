\name{lrEM}
\alias{lrEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logratio EM algorithm
}
\description{
This function implements ordinary and robust Expectation-Maximisation algorithms to replace left-censored values (e.g. values below detection limit, rounded zeros) in compositional data sets.
}
\usage{
lrEM(X, label = NULL, dl = NULL, rob = FALSE,
        ini.cov = c("complete.obs", "multRepl"), delta = 0.65,
        tolerance = 0.0001, max.iter = 50, rlm.maxit = 150, suppress.print = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
}
  \item{label}{Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote unobserved left-censored values in \code{X}.}
  \item{dl}{
Numeric vector of detection limits/thresholds (one per component/column, use e.g. 0 if no threshold for a particular one). These must be given on the same scale as \code{X}.
}
  \item{rob}{
Logical value. \code{FALSE} provides maximum-likelihood estimates of model parameters (default), \code{TRUE} provides robust parameter estimates.
}
  \item{ini.cov}{Initial estimation of the log-ratio covariance matrix. It can be based on either complete observations (\code{"complete.obs"}, default) or multiplicative simple replacement  (\code{"multRepl"}).}
  \item{delta}{
If \code{ini.cov="multRepl"}, delta parameter for initial multiplicative replacement (\code{\link{multRepl}}) in proportions (default = 0.65).
}
  \item{tolerance}{Convergence criterion for the EM algorithm (default = 0.0001).}
  \item{max.iter}{
Maximum number of iterations for the EM algorithm (default = 50).
}
  \item{rlm.maxit}{
If \code{rob=TRUE}, maximum number of iterations for the embedded robust regression estimation (default = 150; see \code{\link{rlm}} in \code{MASS} package for details).
}
  \item{suppress.print}{
Suppress printing of number of iterations required to converge (\code{suppress.print=FALSE}, default).
}
}


\details{
\code{lrEM} produces a replaced data set on the same scale as the input data set. If \code{X} is not closed to a constant sum, then the results are adjusted to provide a compositionally equivalent data set, expressed in the original scale, which leaves the absolute values of the observed components unaltered.

Under maximum likelihood estimation (default, \code{rob=FALSE}), a correction factor based on the residual covariance obtained by censored regression is applied for the correct estimation of the conditional covariance matrix in the maximisation step of the EM algorithm. This is required in order to obtain the conditional expectation of the sum of cross-products between two components in the case that both involve imputed values. Note that the procedure is based on the oblique additive log-ratio (alr) transformation to simplify calculations and alleviates computational burden. Nonetheless, the same results would be obtained using an isometric log-ratio transformation (ilr). Note also that alr requires at least one complete column. Otherwise, a preliminary univariate imputation of the most simplest censoring pattern may be enough.

Under robust estimation (\code{rob=TRUE}), the algorithm requires ilr transformations in order to satisfy requirements for robust estimation methods (MM-estimation by default, see \code{rlm} function for more details). In the unlikely case of censoring patterns involving samples containing only one observed component, these are imputed by multiplicative simple replacement and a warning message identifying them is printed. An initial estimation of nondetects is required to get the algorithm started. This can be based on either the subset of fully observed cases (\code{ini.cov="complete.obs"}) or a multiplicative simple replacement of all nondetects in the data set (\code{ini.cov="multRepl"}). Note that the robust regression method involved includes random elements which can, occasionally, give rise to \code{\link{NaN}} values getting the routine execution halted. If this happened, we suggest to simply re-run the function once again.
}
\value{
A \code{\link{data.frame}} object containing the replaced compositional data set. The number of iterations required for convergence is also printed (this can be suppressed by setting \code{suppress.print=TRUE}).
}
\references{
Martin-Fernandez, J.A., Hron, K., Templ, M., Filzmoser, P., Palarea-Albaladejo, J. Model-based replacement of rounded zeros in compositional data: classical and robust approaches. Computational Statistics & Data Analysis 2012; 56: 2688-2704.

Palarea-Albaladejo J, Martin-Fernandez JA, Gomez-Garcia J. A parametric approach for dealing with compositional rounded zeros. Mathematical Geology 2007; 39: 625-45.

Palarea-Albaladejo J, Martin-Fernandez JA. A modified EM alr-algorithm for replacing rounded zeros in compositional data sets. Computers & Geosciences 2008; 34(8): 902-917.

Palarea-Albaladejo J, Martin-Fernandez JA. Values below detection limit in compositional chemical data. Analytica Chimica Acta 2013; 764: 32-43. DOI: \url{http://dx.doi.org/10.1016/j.aca.2012.12.029}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{zPatterns}}, \code{\link{multRepl}}, \code{\link{multLN}}, \code{\link{lrDA}}, \code{\link{cmultRepl}}
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
              
X_lrEM <- lrEM(X,label=0,dl=rep(1,6),ini.cov="multRepl")
X_roblrEM <- lrEM(X,label=0,dl=rep(1,6),ini.cov="multRepl",rob=TRUE,tolerance=0.001)

# Non-closed compositional data set
data(LPdata) # data (ppm/micrograms per gram)
dl <- c(2,1,0,0,2,0,6,1,0.6,1,1,0,0,632,10) # limits of detection (0 for no limit)
LPdata2 <- subset(LPdata,select=-c(Cu,Ni,La))  # select a subset for illustration purposes
dl2 <- dl[-c(5,7,10)]

LPdata_lrEM <- lrEM(LPdata2,label=0,dl=dl2)
LPdata_roblrEM <- lrEM(LPdata2,label=0,dl=dl2,rob=TRUE,tolerance=0.005)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line