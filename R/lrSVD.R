lrSVD <- function(X, label = NULL, dl = NULL, frac = 0.65, ncp = 2,
                  imp.missing = FALSE, beta = 0.5, method = c("ridge", "EM"),
                  row.w = NULL, coeff.ridge = 1, threshold = 1e-4, seed = NULL, nb.init = 1,
                  max.iter = 1000, z.warning = 0.8, z.delete=TRUE, ...) {
  
  if (any(X < 0, na.rm = T)) stop("X contains negative values")
  
  if ((is.vector(X)) | (nrow(X) == 1)) stop("X must be a data matrix")
  if (is.null(label)) stop("A value for label must be given")
  if (!is.na(label)) {
    if (!any(X == label, na.rm = T))
      stop(paste("Label", label, "was not found in the data set"))
    if (label != 0 &
        any(X == 0, na.rm = T))
      stop("Zero values not labelled as censored values were found in the data set")
    if (any(is.na(X)))
      stop(paste("NA values not labelled as censored values were found in the data set"))
  }
  if (is.na(label)) {
    if (any(X == 0, na.rm = T))
      stop("Zero values not labelled as censored values were found in the data set")
    if (!any(is.na(X), na.rm = T))
      stop(paste("Label", label, "was not found in the data set"))
  }
  
  if (imp.missing==FALSE){
    if (is.character(dl)) stop("dl must be a numeric vector or matrix")
    if (is.null(dl)){ # If dl not given use min per column
      dl <- apply(X,2, function(x) min(x[x!=label]))
      warning("No dl vector or matrix provided. The minimum observed values for each column used as detection limits.")
    }
    if (is.vector(dl)) dl <- matrix(dl,nrow=1)
    dl <- as.matrix(dl) # Avoids problems when dl might be multiple classes
  }
  
  if (imp.missing == FALSE) {
    if (ncol(dl) != ncol(X)) stop("The number of columns in X and dl do not agree")
    if ((nrow(dl) > 1) & (nrow(dl) != nrow(X))) stop("The number of rows in X and dl do not agree")
  }
  if (is.numeric(ncp)){
    if (ncp > min(nrow(X) - 2, ncol(X) - 2)) stop("ncp is too large for the size of the data matrix")
  }
  if (is.null(row.w)) row.w = rep(1, nrow(X)) / nrow(X) # Equal weight for all rows
  
  svd.triplet <- function(X, row.w = NULL, col.w = NULL, ncp = Inf) # From FactoMineR package
  {
    tryCatch.W.E <- function(expr) {
      W <- NULL
      w.handler <- function(w) {
        W <<- w
        invokeRestart("muffleWarning")
      }
      list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), 
                                       warning = w.handler), warning = W)
    }
    if (is.null(row.w)) 
      row.w <- rep(1/nrow(X), nrow(X))
    if (is.null(col.w)) 
      col.w <- rep(1, ncol(X))
    ncp <- min(ncp, nrow(X) - 1, ncol(X))
    row.w <- row.w/sum(row.w)
    X <- t(t(X) * sqrt(col.w)) * sqrt(row.w)
    if (ncol(X) < nrow(X)) {
      svd.usuelle <- tryCatch.W.E(svd(X, nu = ncp, nv = ncp))$val
      if (names(svd.usuelle)[[1]] == "message") {
        svd.usuelle <- tryCatch.W.E(svd(t(X), nu = ncp, nv = ncp))$val
        if (names(svd.usuelle)[[1]] == "d") {
          aux <- svd.usuelle$u
          svd.usuelle$u <- svd.usuelle$v
          svd.usuelle$v <- aux
        }
        else {
          bb <- eigen(crossprod(X, X), symmetric = TRUE)
          svd.usuelle <- vector(mode = "list", length = 3)
          svd.usuelle$d[svd.usuelle$d < 0] = 0
          svd.usuelle$d <- sqrt(svd.usuelle$d)
          svd.usuelle$v <- bb$vec[, 1:ncp]
          svd.usuelle$u <- t(t(crossprod(t(X), svd.usuelle$v))/svd.usuelle$d[1:ncp])
        }
      }
      U <- svd.usuelle$u
      V <- svd.usuelle$v
      if (ncp > 1) {
        mult <- sign(as.vector(crossprod(rep(1, nrow(V)), 
                                         as.matrix(V))))
        mult[mult == 0] <- 1
        U <- t(t(U) * mult)
        V <- t(t(V) * mult)
      }
      U <- U/sqrt(row.w)
      V <- V/sqrt(col.w)
    }
    else {
      svd.usuelle <- tryCatch.W.E(svd(t(X), nu = ncp, nv = ncp))$val
      if (names(svd.usuelle)[[1]] == "message") {
        svd.usuelle <- tryCatch.W.E(svd(X, nu = ncp, nv = ncp))$val
        if (names(svd.usuelle)[[1]] == "d") {
          aux <- svd.usuelle$u
          svd.usuelle$u <- svd.usuelle$v
          svd.usuelle$v <- aux
        }
        else {
          bb <- eigen(crossprod(t(X), t(X)), symmetric = TRUE)
          svd.usuelle <- vector(mode = "list", length = 3)
          svd.usuelle$d[svd.usuelle$d < 0] = 0
          svd.usuelle$d <- sqrt(svd.usuelle$d)
          svd.usuelle$v <- bb$vec[, 1:ncp]
          svd.usuelle$u <- t(t(crossprod(X, svd.usuelle$v))/svd.usuelle$d[1:ncp])
        }
      }
      U <- svd.usuelle$v
      V <- svd.usuelle$u
      mult <- sign(as.vector(crossprod(rep(1, nrow(V)), as.matrix(V))))
      mult[mult == 0] <- 1
      V <- t(t(V) * mult)/sqrt(col.w)
      U <- t(t(U) * mult)/sqrt(row.w)
    }
    vs <- svd.usuelle$d[1:min(ncol(X), nrow(X) - 1)]
    num <- which(vs[1:ncp] < 1e-15)
    if (length(num) == 1) {
      U[, num] <- U[, num, drop = FALSE] * vs[num]
      V[, num] <- V[, num, drop = FALSE] * vs[num]
    }
    if (length(num) > 1) {
      U[, num] <- t(t(U[, num]) * vs[num])
      V[, num] <- t(t(V[, num]) * vs[num])
    }
    res <- list(vs = vs, U = U, V = V)
    return(res)
  }
  
  impute <- function(X = NULL, dl = NULL, bal = NULL, frac = 0.65, ncp = 2, beta = 0.5, method=c("ridge","EM"),
                     row.w = NULL, coeff.ridge = 1, threshold = 1e-4, seed = NULL, max.iter = 1000, init = 1, ...) {
    
    # (scale argument removed: no scaling of olr columns by (weighted) variance allowed)
    
    # Weighted AVERAGE = moyenne poids
    moy.p <- function(V, poids) {
      res <- sum(V * poids,na.rm = TRUE)/sum(poids[!is.na(V)])
    }
    
    # Geometric mean by columns (for missing data case)
    gm <- function(x, na.rm = TRUE){
      exp(sum(log(x), na.rm = na.rm)/length(x))
    }
    
    nb.iter <- 1
    old <- Inf
    objective <- 0
    if (!is.null(seed)) {set.seed(seed)} # fix seed to have same results
    
    # OLR of initial data matrix
    
    # Missing pattern
    missRaw <- which(is.na(X))
    obsRaw <- which(!is.na(X))
    X <- as.matrix(X)
    Xaux <- X # copy original raw data with NA
    caux <- apply(Xaux, 1, sum, na.rm = TRUE)
    
    # Initial imputation
    if (imp.missing == FALSE) {
      if (init == 1) {X[missRaw] <- frac*dl[missRaw]} # mult repl
      else {X[missRaw] <- runif(1,0.50,0.8)*dl[missRaw]} # random initialisations if nb.init/init > 1
    }
    else{
      # Initial geo mean imputation of missing
      gmeans <- apply(Xaux,2,function(x) gm(x))
      nn <- nrow(X)
      gmeans <- matrix(rep(1, nn), ncol = 1) %*% gmeans
      X[missRaw] <- gmeans[missRaw]
    }
    
    # Xhat: OLR-coordinates
    Xhat <- t(bal %*% t(log(X)))
    
    # Number of components
    ncp <- min(ncp, ncol(Xhat), nrow(Xhat) - 1)
    # Weighted column mean  
    mean.p.or=mean.p <- apply(Xhat, 2, moy.p,row.w)
    # Matrix centring
    Xhat <- t(t(Xhat) - mean.p)
    # Update X: olr.inv(Xhat)
    X <- exp(t(t(bal) %*% t(Xhat)))
    X <- X / apply(X, 1, sum)
    # Aux data matrix for observed and non-observed data
    fittedX <- fittedXus <- Xhat
    fittedXRaw <- fittedXusRaw <- X
    if (ncp == 0) {nb.iter <- 0}
    
    while (nb.iter > 0) {
      # Update data matrix  
      X[missRaw] <- fittedXRaw[missRaw]
      # Xhat: OLR-coordinates
      Xhat <- t(bal %*% t(log(X)))
      # Recover the centre
      Xhat <- t(t(Xhat) + mean.p)
      mean.p <- apply(Xhat, 2, moy.p, row.w)
      # violations check
      X <- exp(t(t(bal) %*% t(Xhat)))
      X <- X/apply(X, 1, sum)
      Xaux2 <- X * caux
      viol <- which(Xaux2 > dl)
      Xaux2[viol] <- dl[viol]
      Xhat <- t(bal %*% t(log(Xaux2)))
      mean.p <- apply(Xhat, 2, moy.p, row.w)
      # Update X: olr.inv(Xhat)
      X <- exp(t(t(bal) %*% t(Xhat)))
      X <- X / apply(X, 1, sum)
      # Impute observed values
      fittedXusRC <- t(t(fittedXus) + mean.p)# recover centre
      # RAW values
      # INV-OLR #
      fittedXusRCRaw <- exp(t(t(bal) %*% t(fittedXusRC)))
      fittedXusRCRaw <- fittedXusRCRaw / apply(fittedXusRCRaw, 1, sum)
      
      X[obsRaw] <- ((fittedXusRCRaw[obsRaw]) ^ (1 - beta)) * ((X[obsRaw]) ^ beta)
      
      # check DL
      X <- X/apply(X,1,sum)
      Xaux2 <- X*caux # re-scaled to original
      viol <- which(Xaux2 > dl)
      Xaux2[viol] <- dl[viol]
      
      # Xhat: OLR-coordinates
      Xhat <- t(bal %*% t(log(Xaux2)))
      # Update mean
      mean.p <- apply(Xhat, 2, moy.p, row.w)
      # Centring
      Xhat <- t(t(Xhat) - mean.p)
      # Update X
      # INV-OLR
      X <- exp(t(t(bal) %*% t(Xhat)))
      X <- X / apply(X, 1, sum)
      # SVD calculation WEIGHTED by row.w and rank ncp
      svd.res <- svd.triplet(Xhat, row.w = row.w, ncp = ncp)
      sigma2 <-
        nrow(Xhat) * ncol(Xhat) / min(ncol(Xhat), nrow(Xhat) - 1) * sum((svd.res$vs[-c(1:ncp)] ^
                                                                           2) / ((nrow(Xhat) - 1) * ncol(Xhat) - (nrow(Xhat) - 1) * ncp - ncol(Xhat) * ncp + ncp ^
                                                                                   2))
      sigma2 <- min(sigma2 * coeff.ridge, svd.res$vs[ncp + 1] ^ 2)
      if (method == "EM")
        sigma2 <- 0
      # usual lambda 
      lambda.us <- svd.res$vs[1:ncp]
      # calculate the usual new matrix
      fittedXus <-
        tcrossprod(t(t(svd.res$U[, 1:ncp, drop = FALSE] * row.w) * lambda.us), svd.res$V[, 1:ncp, drop =
                                                                                           FALSE])
      fittedXus <- fittedXus / row.w
      # lambda for regularisation
      lambda.shrinked <- (svd.res$vs[1:ncp] ^ 2 - sigma2) / svd.res$vs[1:ncp]
      # calculate the new matrix for regularisation
      fittedX <- tcrossprod(t(t(svd.res$U[, 1:ncp, drop = FALSE] * row.w) * lambda.shrinked), svd.res$V[, 1:ncp, drop = FALSE])
      fittedX <- fittedX / row.w
      # calculate the Frobenius norm of the difference between iterations (convergence)
      # INV-OLR
      fittedXRaw <- exp(t(t(bal) %*% t(fittedX)))
      fittedXRaw <- fittedXRaw / apply(fittedXRaw, 1, sum)
      
      diffRaw <- X / fittedXRaw
      diffRaw[missRaw] <- 1
      # OLR-coordinates
      diff <- t(bal %*% t(log(diffRaw)))
      
      objective <- sum(diff^2*row.w)
      # objective <- mean((Xhat[-missing]-fittedX[-missing])^2)
      
      # Convergence
      criterion <- abs(1 - objective/old)
      old <- objective
      nb.iter <- nb.iter + 1
      if (!is.nan(criterion)) {
        if ((criterion < threshold) && (nb.iter > 5))  nb.iter <- 0
        if ((objective < threshold) && (nb.iter > 5))  nb.iter <- 0
      }
      if (nb.iter > max.iter) {
        nb.iter <- 0
        warning(paste("Stopped after ",max.iter," iterations"))
      }
    }
    # END LOOP WHILE
    
    # Preparing the results
    Xhat <- t(t(Xhat) + mean.p.or)
    
    # Update X
    # INV-OLR
    X <- exp(t(t(bal) %*% t(Xhat)))
    X <- X/apply(X,1,sum)
    
    # completeObs
    completeObs <- Xaux / apply(Xaux, 1, sum, na.rm = TRUE)
    completeObs[missRaw] <- X[missRaw]
    completeObs <- completeObs / apply(completeObs, 1, sum)
    # violation dl check
    completeObs <- completeObs * caux
    viol <- which(completeObs > dl)
    completeObs[viol] <- dl[viol]
    completeObs[obsRaw] <- Xaux[obsRaw]
    completeObs <- completeObs/apply(completeObs, 1, sum)
    #
    fittedX <- t(t(fittedX) + mean.p)
    
    # INV-OLR
    fittedXRaw <- exp(t(t(bal) %*% t(fittedX)))
    fittedXRaw <- fittedXRaw/apply(fittedXRaw,1,sum)
    
    # Return complete matrix and imputed matrix
    result <- list()
    result$completeObs <- completeObs
    result$fittedX <- fittedXRaw
    
    return(result) 
  }
  
  method <- match.arg(method)
  
  ## Preliminaries ----  
  
  X <- as.data.frame(X, stringsAsFactors = TRUE)
  nn <- nrow(X)
  D <- ncol(X)
  X[X == label] <- NA
  X <- as.data.frame(apply(X, 2, as.numeric), stringsAsFactors = TRUE)
  c <- apply(X, 1, sum, na.rm = TRUE)
  
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
  
  # Check for closure
  closed <- 0
  if (all(abs(c - mean(c)) < .Machine$double.eps^0.3)) closed <- 1
  
  # Sort columns decreasingly according to observed cells
  Xaux <- as.matrix(X)
  XauxClosed <- Xaux / apply(Xaux, 1, sum, na.rm = TRUE) #as.matrix(X) # copy original data set
  # Replace zeros by NA
  Xna <- X # copy of original data
  #colnames(Xna)
  pz <- apply(apply(Xna, 2, is.na), 2, sum) / nn
  # Ordered decreasingly by number of zeros (NA)
  X <- Xna[, order(-pz)]
  
  if (imp.missing == FALSE) {
    if (nrow(dl) == 1) {dl <- matrix(dl[, order(-pz)], nrow = 1)}
    else{dl <- dl[, order(-pz)]}
  }
  
  # Balance matrix for olr
  Smat <- diag(rep(1,D))
  Smat[upper.tri(Smat)] <- -1
  Smat <- Smat[-D,]
  bal <- Smat
  numsbp <- dim(Smat)[1]
  for (f in 1:numsbp) {
    den <- sum(bal[f, ] == -1)
    num <- sum(bal[f, ] == 1)
    bal[f, bal[f, ] == 1] <- sqrt(den / ((den + num) * num))
    bal[f, bal[f, ] == -1] <- -sqrt(num / ((den + num) * den))
  }
  
  # Build dl matrix for SVD imputation
  if (imp.missing == FALSE) {
    if (nrow(dl) == 1)
      dl <- matrix(rep(1, nn), ncol = 1) %*% dl
  }
  else {dl <- matrix(0, nrow = nn, ncol = D)} # fake dl matrix for the case of missing
  
  # Set observed data as upper bound for estimates of observed values
  observedRaw <- which(!is.na(X))
  missingRaw <- which(is.na(X))
  Xaux2 <- as.matrix(X)
  # DL observed and missing values = maximum value = observed = infinity upper bound
  Xmax <- apply(X, 2, max, na.rm = TRUE)
  Xmax <- matrix(rep(1, nn), ncol = 1) %*% Xmax
  dl[observedRaw] <- Xmax[observedRaw]
  if (imp.missing == TRUE) {
    dl[missingRaw] <- Xmax[missingRaw]
  }
  colnames(dl) <- colnames(X)
  
  ## Imputation ---
  
  for (i in 1:nb.init) {
    
    if (!any(is.na(X))) return(X)
    
    res.impute <- impute(X = X, dl = dl, bal = bal, frac = frac, ncp = ncp, beta = beta, method = method,
                         row.w = row.w, coeff.ridge = coeff.ridge, threshold = threshold, seed = if (!is.null(seed)) {
                           (seed * (i - 1))} else {NULL},max.iter = max.iter,init = i)
  }
  
  ## Final section ---
  
  # Re-scale to original units
  Y <- res.impute$completeObs
  XauxClosed <- XauxClosed[, order(-pz)]
  XauxClosed[missingRaw] <- Y[missingRaw]
  X <- XauxClosed * c
  
  if (closed == 1) {
    X <- (X / apply(X, 1, sum)) * c[1]
  }
  
  # Original order
  X <- X[, colnames(Xna)]
  
  return(as.data.frame(X, stringsAsFactors = TRUE))
}