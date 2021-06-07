#' Parametric Bootstrap Mean Squared Error of EBLUP's for domain means using Geoadditive Small Area Model
#'
#' This function calculates MSE of EBLUP's based on unit level using Geoadditive Small Area Model
#'
#' @param formula the model that to be fitted
#' @param zspline matrix that used in model for random effect of spline-2
#' @param dom domain codes
#' @param xmean matrix of auxiliary variables means for each domains
#' @param zmean matrix of spline-2 means for each domains
#' @param data data unit level that used as data frame that containing the variables named in formula and dom
#' @param B iteration bootstraping
#'
#' @export pbmsegeo

pbmsegeo<-function(formula, zspline, dom, xmean, zmean, data, B=100)
{
  result <- list(est = NA, mse=NA)
  namedom<-deparse(substitute(dom))
  if (!missing(data)) {
    formuladata <- model.frame(formula, na.action = na.omit,data)
    X <- model.matrix(formula, data)
    dom <- data[,"area"]
  }
  else {
    formuladata <- model.frame(formula, na.action = na.omit)
    X <- model.matrix(formula)
  }

  if (is.factor(dom))
    dom <- as.vector(dom)
  if (sum(c(nrow(formuladata) != length(dom)) , (nrow(formuladata) != nrow(zspline)))!=0)
    stop("   Arguments formula [rows=", nrow(formuladata),
         "] , dom [rows=", length(dom), "] and zmatrix [rows=", nrow(zspline), "] must be the same length.\n")

  y <- formuladata[, 1]
  z <- zspline

  result$est<-eblupgeo(y~x1+x2+x3, zspline, dom, xmean, zmean, data)

  udom<-unique(dom)
  m<-length(udom)
  k<-ncol(z)
  n<-length(y)
  D<-matrix(0,nrow=n,ncol=m)
  for(i in 1:n){
    domain<-	dom[i]
    D[i,]<-rep(c(0,1,0),c(domain-1,1,m-domain))
  }

  beta<-result$est$fit$coefficients$fixed
  sigma2.gamma<-result$est$sigma2[1]
  sigma2.u<-result$est$sigma2[2]
  sigma2.e<-result$est$sigma2[3]

  cat("\nBootstrap procedure with B =", B, "iterations starts.\n")
  summse.pb<-NULL
  error<-NULL
  b<-1
  seed<-1000
  while(b<=B){
    set.seed(seed)
    u.boot <- rnorm(m,0,sqrt(sigma2.u))
    gamma.boot <- rnorm(k, 0, sqrt(sigma2.gamma))
    e.boot <- rnorm(n, 0, sqrt(sigma2.e))
    ebar.boot<-aggregate(e.boot,list(dom),mean)[,2]
    theta.boot <- as.matrix(xmean) %*% as.vector(beta) + as.matrix(zmean) %*% as.vector(gamma.boot) + u.boot+ebar.boot
    y_ij.boot <- as.matrix(X) %*% as.vector(beta) + as.matrix(z) %*% as.vector(gamma.boot) + D %*% u.boot + e.boot
    model.boot <- eblupgeo(y_ij.boot~ x1+x2+x3 , z, dom, xmean, zmean, data)
    selisih<-(model.boot$eblup - theta.boot)^2
    error<-rbind(error,t(selisih))
    b=b+1
    seed<-seed+10*b
  }

  result$mse<-colMeans(error, na.rm=T)
  return(result)
}
