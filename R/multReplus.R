multReplus <- function(X, dl = NULL, frac = 0.65, closure = NULL, z.warning = 0.8, z.delete = TRUE, delta = NULL){
  
  if (any(X<0, na.rm=T)) stop("X contains negative values")
  if (is.character(dl)) stop("dl must be a numeric vector or matrix")
  if (is.null(dl)){ # If dl not given use min per column
    dl <- apply(X,2, function(x) min(x[x!=0]))
    warning("No dl vector or matrix provided. The minimum observed values for each column used as detection limits.")
  }
  if (is.vector(dl)) dl <- matrix(dl,nrow=1)
  dl <- as.matrix(dl) # Avoids problems when dl might be multiple classes
  if ((is.vector(X)) | (nrow(X)==1)) stop("X must be a data matrix")
  if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
  if ((nrow(dl)>1) & (nrow(dl)!=nrow(X))) stop("The number of rows in X and dl do not agree")
  if (any(is.na(X))==FALSE) stop("No missing data were found in the data set")
  if (any(X==0, na.rm=T)==FALSE) stop("No zeros were found in the data set")
  
  if (!missing("delta")){
    warning("The delta argument is deprecated, use frac instead: frac has been set equal to delta.")
    frac <- delta
  }
  
  gm <- function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x[!is.na(x)]))
  }
  
  nam <- NULL
  if (!is.null(names(X))) nam <- names(X)
  if (is.vector(X)) X <- as.data.frame(matrix(X,ncol=length(X)),stringsAsFactors=TRUE)
  
  ## Preliminaries ----

  X <- as.data.frame(X,stringsAsFactors=TRUE)
  nn <- nrow(X); D <- ncol(X)
  X <- as.data.frame(apply(X,2,as.numeric),stringsAsFactors=TRUE)
  c <- apply(X,1,sum,na.rm=TRUE)

  checkNumZerosCol <- apply(X, 2, function(x) sum(is.na(x)))
  
  if (any(checkNumZerosCol/nrow(X) > z.warning)) {    
    cases <- which(checkNumZerosCol/nrow(X) > z.warning)    
    if (z.delete == TRUE) {
      if (length(cases) > (ncol(X)-2)) {
        stop(paste("Almost all columns contain >", z.warning*100,
                   "% zeros/unobserved values (see arguments z.warning and z.delete).",
                   sep=""))
      }      
      X <- X[,-cases]      
      action <- "deleted"
      
      warning(paste("Column no. ",cases," containing >", z.warning*100,
                    "% zeros/unobserved values ", action, " (see arguments z.warning and z.delete).\n",
                    sep=""))
    } else {      
      action <- "found"      
      warning(paste("Column no. ",cases," containing >", z.warning*100,
                    "% zeros/unobserved values ", action, " (see arguments z.warning and z.delete. Check out with zPatterns()).\n",
                    sep=""))      
    }
  }
  
  checkNumZerosRow <- apply(X, 1, function(x) sum(is.na(x)))  
  if (any(checkNumZerosRow/ncol(X) > z.warning)) {    
    cases <- which(checkNumZerosRow/ncol(X) > z.warning)    
    if (z.delete == TRUE) {
      if (length(cases) > (nrow(X)-2)) {
        stop(paste("Almost all rows contain >", z.warning*100,
                   "% zeros/unobserved values (see arguments z.warning and z.delete).",
                   sep=""))
      }
      X <- X[-cases,]      
      action <- "deleted"
      
      warning(paste("Row no. ",cases," containing >", z.warning*100,
                    "% zeros/unobserved values ", action, " (see arguments z.warning and z.delete).\n",
                    sep=""))
    } else {      
      action <- "found"      
      warning(paste("Row no. ", cases," containing >", z.warning*100,
                    "% zeros/unobserved values ", action,
                    " (see arguments z.warning and z.delete. Check out with zPatterns()).\n",
                    sep=""))      
    }
  }
  
  if (nrow(dl)==1) dl <- matrix(rep(1,nn),ncol=1)%*%dl

  # Check for closure
  closed <- 0
  if (all( abs(c - mean(c)) < .Machine$double.eps^0.3 )) closed <- 1

  Y <- X
  
  if (!is.null(closure)){
    if (closed == 1) {stop("closure: The data are already closed to ",c[1])}
    resid <- apply(X,1, function(x) closure-sum(x, na.rm = TRUE))
    Xresid <- cbind(X,resid,stringsAsFactors=TRUE)
    c <- rep(closure,nn)
    Y <- Xresid
  }
  
  ## Imputation of missing (ignoring 0s in column if any) ----
  
  gms <- apply(X,2,function(x) gm(x[x!=0]))
  for (i in 1:nn){
    if (any(is.na(X[i,]))){
      z <- which(is.na(X[i,]))
      Y[i,z] <- gms[z]
      if (!is.null(closure)){
        Y[i,-z] <- ((c[i]-(sum(Y[i,z])))/sum(Xresid[i,-z]))*Xresid[i,-z]
        tmp <- Y[i,-(D+1)]
        nz_idx <- which(tmp[-z]!=0)[1] # Use a non-zero part for adjustment
        X[i,z] <- as.numeric((X[i,-z][nz_idx]/tmp[-z][nz_idx]))*Y[i,z]
      }
      else{
        Y[i,-z] <- ((c[i]-(sum(Y[i,z])))/sum(X[i,-z]))*X[i,-z]
        nz_idx <- which(Y[i,-z]!=0)[1] # Use a non-zero part for adjustment
        X[i,z] <- as.numeric((X[i,-z][nz_idx]/Y[i,-z][nz_idx]))*Y[i,z]
      }
    }
  } 
  
  if (closed==1){
    X <- t(apply(X,1,function(x) x/sum(x)*c[1]))
  } 
  
  if (any(X < 0, na.rm = T)) stop("multRepl: negative imputed values were generated (please check out help for advice)")
  
  ## Imputation of zeros ----
  
  X <- multRepl(X,label=0,dl=dl,frac=frac,closure=closure,z.warning=z.warning,z.delete=z.delete)
  
  ## Final section ----
  
  if (!is.null(nam)) names(X) <- nam
  
  return(as.data.frame(X,stringsAsFactors=TRUE))
  
}