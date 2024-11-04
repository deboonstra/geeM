#' @title Correlation information criteria
#'
#' @description Calculating the correlation information criteria (CIC) based on
#' Hardin and Hilbe (2012).
#'
#' @param object An object of class \code{"geem"} representing the fit.
#' @param tol A numeric value specifying the tolerance when calculating the
#' inverse of the model-based covariance matrix.
#' Default is \code{.Machine$double.eps}.
#' @param envir An environment object to which find the \code{"geem"} object.
#'
#' @details The calculation of CIC proposed by Hin and Wang (2009) is a
#' modification of the quasi-likehood under the working independence model
#' (QIC) criteria, which was proposed by Pan (2001a), by only using the penalty
#' term of QIC. This version of CIC uses the mean structure parameter estimates
#' and the scale parameter estimate based on the working correlation structure
#' to calculate the model-based covariance matrix assuming the correlation
#' structure is independent. However, in 2012, Hardin and Hilbe proposed to
#' calculate CIC using the model-based covariance matrix given the mean
#' structure parameter estimates and the scale parameter estimate produced by
#' the independence model. Hin and Wang noted that the difference between the
#' two calculations of CIC is \eqn{O_{p}(n^{-1/2})}; thus, they do not expect
#' there to be much of a difference between the two variants of CIC.
#'
#' If the MASS package is loaded then the \code{\link{ginv}} function is used
#' for matrix inversion. Otherwise the standard \code{\link{solve}} function is
#' used.
#'
#' @return A numeric vector containing the CIC value for the \code{"geem"}
#' object under consideration.
#'
#' @references
#' Pan, W. (2001a). \emph{Akaike's information criterion in generalized
#' estimating equations}. Biometrics, 57, 120-125.\cr
#' Hardin, J.W.  and Hilbe, J.M. (2012). \emph{Generalized
#' Estimating Equations, 2nd Edition}, Chapman and Hall/CRC: New York. \cr
#' Hin, L.-Y. and Wang, Y-G. (2009). \emph{Working-correlation-structure
#' identification in generalized estimating equations},
#' Statistics in Medicine 28: 642-658.
#'
#' @examples
#' data(ohio, package = "geepack")
#' fit <- geeM::geem(
#'   resp ~ age + smoke + age:smoke,
#'   id = id, data = ohio, family = stats::binomial(link = "logit"),
#'   corstr = "exchangeable"
#' )
#' geeM::cic(fit)
#'
#' @export
cic <- function(object, tol = .Machine$double.eps, envir = parent.frame()) {
  # Checking function parameter types ####
  if (!methods::is(object, "geem")) {
    stop("object must be an object produced by geeM::geem")
  }
  if (!(class(tol) %in% c("numeric", "integer"))) {
    stop("tol must be a numeric or an integer value")
  }
  if (!is.environment(envir)) {
    stop("envir must be an environment object")
  }

  # Obtaining covariance matrices ####

  ## Pulling modeling call to re-fit models ####
  model_call <- object$call

  ## Assigning corstr to independence ####
  model_call$corstr <- "independence"

  ## Re-fitting model
  model_indep <- eval(model_call, envir = envir)

  ## Getting covariance matrices ####

  ### Inverse of model-based covariance under independence ####
  omega <- invert(as.matrix(model_indep$naiv.var), tol = tol)

  ### Robust covariance given working correlation structure ####
  vr <- as.matrix(object$var)

  # Calculating CIC ####
  trace <- sum(diag(x = crossprod(x = omega, y = vr)))

  # Returning CIC
  return(trace)
}