\name{lrEM}
\alias{lrEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-ratio EM algorithm
}
\description{
This function implements model-based ordinary and robust Expectation-Maximisation algorithms to impute left-censored data (e.g. values below detection limit, rounded zeros) via coordinates representation of compositional data which incorporate the information of the relative covariance structure.

This function can be also used to impute missing data instead by setting \code{imp.missing = TRUE} (see \code{\link{lrEMplus}} to treat censored and missing data simultaneously).
}
\usage{
lrEM(X, label = NULL, dl = NULL, rob = FALSE,
        ini.cov = c("complete.obs", "multRepl"), frac = 0.65,
        tolerance = 0.0001, max.iter = 50, rlm.maxit = 150,
        imp.missing = FALSE, suppress.print = FALSE, closure = NULL,
        z.warning = 0.8, z.delete = TRUE, delta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Compositional data set (\code{\link{matrix}} or \code{\link{data.frame}} class).
}
  \item{label}{Unique label (\code{\link{numeric}} or \code{\link{character}}) used to denote zeros/unobserved values in \code{X}.}
  \item{dl}{
Numeric vector or matrix of detection limits/thresholds. These must be given on the same scale as \code{X}. If \code{NULL} the column minima are used as thresholds. 
}
  \item{rob}{
Logical value. \code{FALSE} provides maximum-likelihood estimates of model parameters (default), \code{TRUE} provides robust parameter estimates.
}
  \item{ini.cov}{Initial estimation of either the log-ratio covariance matrix (ML estimation) or unobserved data (robust estimation). It can be based on either complete observations (\code{"complete.obs"}, default) or multiplicative simple replacement (\code{"multRepl"}).}
  \item{frac}{
If \code{ini.cov="multRepl"}, parameter for initial multiplicative simple replacement of left-censored data (see \code{\link{multRepl}}) (default = 0.65).
}
  \item{tolerance}{Convergence criterion for the EM algorithm (default = 0.0001).}
  \item{max.iter}{
Maximum number of iterations for the EM algorithm (default = 50).
}
  \item{rlm.maxit}{
If \code{rob=TRUE}, maximum number of iterations for the embedded robust regression estimation (default = 150; see \code{\link[MASS]{rlm}} for details).
}
  \item{imp.missing}{If \code{TRUE} then unobserved data identified by \code{label} are treated as missing data (default = \code{FALSE}).}
  \item{suppress.print}{
Suppress printed feedback (\code{suppress.print = FALSE}, default).
}
  \item{closure}{Closure value used to add a residual part if needed when \code{ini.cov="multRepl"} is used (see \code{?multRepl}).}
 \item{z.warning}{Threshold used to identify individual rows or columns including an excess of zeros/unobserved values (to be specify in proportions, default \code{z.warning=0.8}).}
  \item{z.delete}{Logical value. If set to \code{TRUE}, rows/columns identified by \code{z.warning} are omitted in the imputed data set. Otherwise, the function stops in error when rows/columns are identified by \code{z.warning} (default \code{z.delete=TRUE}).}
  \item{delta}{This argument has been deprecated and replaced by \code{frac} (see package's NEWS for details).}
}


\details{
After convergence, this function imputes unobserved compositional data by their estimated conditional expected values through coordinates representation, given the information from the observed data. For left-censoring problems, it allows for either single (\code{vector} form) or multiple (\code{matrix} form, same size as \code{X}) limits of detection by component. Any threshold value can be set for non-censored elements (e.g. use 0 if no threshold for a particular column or element of the data matrix).

It produces an imputed data set on the same scale as the input data set. If \code{X} is not closed to a constant sum, then the results are adjusted to provide a compositionally equivalent data set, expressed in the original scale, which leaves the absolute values of the observed components unaltered.

Under maximum likelihood (ML) estimation (default, \code{rob=FALSE}), a correction factor based on the residual covariance obtained by censored regression is applied for the correct estimation of the conditional covariance matrix in the maximisation step of the EM algorithm. This is required in order to obtain the conditional expectation of the sum of cross-products between two components in the case that both involve imputed values. Note that the procedure is based on the oblique additive log-ratio (alr) transformation to simplify calculations and alleviates computational burden. Nonetheless, the same results would be obtained using an isometric log-ratio transformation (ilr). Note also that alr requires at least one complete column. Otherwise, a preliminary imputation, e.g. by \code{\link{multRepl}} or \code{\link{multLN}}, of the most simplest censoring pattern may be enough. The argument \code{ini.cov} determines how the initial estimation of the log-ratio covariance matrix required to start the EM process is worked out. Note that the estimation of the covariance matrix, and hence the lrEM routine, requires a regular data set, i.e. having more observations than variables in the data.

Under robust estimation (\code{rob=TRUE}), the algorithm requires ilr transformations in order to satisfy requirements for robust estimation methods (MM-estimation by default, see \code{\link[MASS]{rlm}} function for more details). An initial estimation of nondetects is required to get the algorithm started. This can be based on either the subset of fully observed cases (\code{ini.cov="complete.obs"}) or a multiplicative simple replacement of all nondetects in the data set (\code{ini.cov="multRepl"}). Note that the robust regression method involved includes random elements which can, occasionally, give rise to \code{\link{NaN}} values getting the routine execution halted. If this happened, we suggest to simply re-run the function once again.

Note that conditional imputation based on log-ratio coordinates cannot be conducted when there exist censoring patterns including samples with only one observed component. As a workaround, \code{lrEM} applies multiplicative simple replacement (\code{\link{multRepl}}) on those and a warning message identifying the problematic cases is printed. Alternatively, it might be sensible to simply remove those non-informative samples from the data set.

\emph{Missing data imputation}

When \code{imp.missing = TRUE}, unobserved values are treated as general missing data and imputed by their conditional expectation using the EM algorithm. Either maximum-likelihood or robust estimation can be used through the \code{rob} argument. For this case, the argument \code{label} indicates the unique label for missing values. The algorithm can be initiated using either \code{"complete.obs"} or \code{"multRepl"} (for missing data) as specified by the \code{ini.cov} argument. The argument \code{dl} is ignored.
}
\value{
A \code{\link{data.frame}} object containing the imputed compositional data set expressed in the original scale. The number of iterations required for convergence is also printed (this can be suppressed by setting \code{suppress.print=TRUE}).
}
\references{
Martin-Fernandez, J.A., Hron, K., Templ, M., Filzmoser, P., Palarea-Albaladejo, J. Model-based replacement of rounded zeros in compositional data: classical and robust approaches. Computational Statistics & Data Analysis 2012; 56: 2688-2704.

Palarea-Albaladejo J, Martin-Fernandez JA, Gomez-Garcia J. A parametric approach for dealing with compositional rounded zeros. Mathematical Geology 2007; 39: 625-45.

Palarea-Albaladejo J, Martin-Fernandez JA. A modified EM alr-algorithm for replacing rounded zeros in compositional data sets. Computers & Geosciences 2008; 34: 902-917.

Palarea-Albaladejo J, Martin-Fernandez JA. Values below detection limit in compositional chemical data. Analytica Chimica Acta 2013; 764: 32-43.

Palarea-Albaladejo J. and Martin-Fernandez JA. zCompositions -- R package for multivariate imputation of
left-censored data under a compositional approach. Chemometrics and Intelligence Laboratory Systems 2015; 143: 85-96.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{zPatterns}}, \code{\link{lrSVD}}, \code{\link{lrDA}}, \code{\link{multRepl}}, \code{\link{multLN}}, \code{\link{multKM}}, \code{\link{cmultRepl}}
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

# Multiple limits of detection by component
mdl <- matrix(0,ncol=6,nrow=10)
mdl[2,] <- rep(1,6)
mdl[4,] <- rep(0.75,6)
mdl[6,] <- rep(0.5,6)
mdl[8,] <- rep(0.5,6)
mdl[10,] <- c(0,0,1,0,0.8,0.7)

X_lrEM2 <- lrEM(X,label=0,dl=mdl,ini.cov="multRepl")

# Non-closed compositional data set
data(LPdata) # data (ppm/micrograms per gram)
dl <- c(2,1,0,0,2,0,6,1,0.6,1,1,0,0,632,10) # limits of detection (0 for no limit)
LPdata2 <- subset(LPdata,select=-c(Cu,Ni,La))  # select a subset for illustration purposes
dl2 <- dl[-c(5,7,10)]

LPdata2_lrEM <- lrEM(LPdata2,label=0,dl=dl2)
LPdata2_roblrEM <- lrEM(LPdata2,label=0,dl=dl2,rob=TRUE,tolerance=0.005)

# Two subsets of limits of detection (using e.g. robust parameter estimation)
 # Using a subset of LPdata for faster execution
data(LPdata) # data (ppm/micrograms per gram)
LPdata2 <- subset(LPdata,select=-c(Cu,Ni,La))
dl2 <- c(2,1,0,0,0,1,0.6,1,0,0,632,10)
 # DLs for first 50 samples of LPdata2
dl2a <- matrix(rep(1,50),ncol=1)\%*\%dl2
 # DLs for last 46 samples of LPdata
dl2b <- matrix(rep(1,46),ncol=1)\%*\%c(1,0.5,0,0,0,0.75,0.3,1,0,0,600,8) 

mdl <- rbind(dl2a,dl2b)
LPdata2_roblrEM <- lrEM(LPdata2,label=0,dl=mdl,rob=TRUE,tolerance=0.005)

# Treating zeros as general missing data for illustration purposes only
LPdata2_miss <- lrEM(LPdata2,label=0,imp.missing=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
