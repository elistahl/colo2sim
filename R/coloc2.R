# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}


var_mle_from_z = function(z,n,maf){
  var_mle= 1/(2*maf*(1-maf) * ( n + z^2))
  return(var_mle)
}

approx.bf.estimates.pvalue  <- function(z, n, maf, suffix=NULL, sdY=1){
  V = var_mle_from_z(z,n,maf)
  return(approx.bf.estimates.ave(z, V, type, suffix=suffix))
}

approx.bf.estimates.ave <- function (z, V, type, suffix=NULL, sdY=1) {
  listVec <- list(lABF_sd1 = lABF.fn(z, V, sd.prior=sqrt(0.01)*sdY), lABF_sd2 = lABF.fn(z, V, sd.prior=sqrt(0.1)*sdY), lABF_sd3 = lABF.fn(z, V, sd.prior=sqrt(0.5)*sdY))
  m <- do.call(cbind, listVec)
  lABF <- apply(m, 1, function(x) logsum(x) -log(3))
  ret <- data.frame(V, z, m, lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}


lABF.fn <- function (z, V, sd.prior=0.15) {
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  return(lABF)
  #return(list("lABF" = lABF, "r" = r))
}














symmetricize <-
  function(matrix, method=c("max", "min","avg", "ld", "ud"), adjacencyList=FALSE){
    method <- match.arg(method)

    if (missing(matrix)){
      stop("You must provide a matrix to symmetricize.")
    }

    x <- matrix
    if (method=="ld"){
      #solution from Michael Conklin, http://www.biostat.wustl.edu/archives/html/s-news/2000-03/msg00127.html
      x[matrix(c(col(x)[lower.tri(x)], row(x)[lower.tri(x)]), ncol = 2)] <-
        x[matrix(c(row(x)[lower.tri(x)], col(x)[lower.tri(x)]), ncol = 2)]
    }
    if (method=="ud"){
      x[matrix(c(col(x)[upper.tri(x)], row(x)[upper.tri(x)]), ncol = 2)] <-
        x[matrix(c(row(x)[upper.tri(x)], col(x)[upper.tri(x)]), ncol = 2)]
    }
    if (method=="max" || method=="min" || method=="avg"){
      #retrieve the (row,col) indices for the lower-diagonal entries in the matrix
      ldi <- matrix(c(row(x)[lower.tri(x)], col(x)[lower.tri(x)]), ncol = 2)
      #already column-major, so we can leave this

      #retrieve the (row,col) indices for the upper-diagonal entries in the matrix
      udi <- matrix(c(row(x)[upper.tri(x)], col(x)[upper.tri(x)]), ncol = 2)
      #in order for these to be in the symmetrical order as ldi, we need to sort.
      udi <- udi[order(udi[,1], udi[,2]),]

      #extract the upper and lower diagonal elements in a way that's symmetrical in their indexces
      ud <- x[udi]
      ld <- x[ldi]

      #replace with either the min, max, or mean, depending on the selection
      if (method=="max"){
        x[ldi] <- apply(rbind(ld,ud),2,max);
        x[udi] <- apply(rbind(ld,ud),2,max);
      }
      if (method=="min"){
        x[ldi] <- apply(rbind(ld,ud),2,min);
        x[udi] <- apply(rbind(ld,ud),2,min);
      }
      if (method=="avg"){
        x[ldi] <- apply(rbind(ld,ud),2,mean);
        x[udi] <- apply(rbind(ld,ud),2,mean);
      }
    }

    if (!adjacencyList){
      #return the adjacency matrix
      return(x)
    }
    else{
      #convert to adjacency list and add addressing
      if (is.null(rownames(matrix))){
        stop("You requested adjacency list format, but the matrix provided has no (row) names.")
      }
      names <- rownames(matrix)
      add <- getTableAddressing(names)
      toReturn <- (cbind(add, x[upper.tri(x)]))
      colnames(toReturn) <- c("Source", "Dest", "Symmetric")
      return(toReturn)
    }
  }
