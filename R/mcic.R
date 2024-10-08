#' @title Modified correlation information criteria
#' @description Calculating the modified correlation information criteria (mCIC)
#' based on Boonstra and Cavanaugh (in works).
#'
#' @param object An object of class \code{"geem"} representing the fit.
#' @param type How mCIC is calculated by using either the model-based covariance
#' matrix produced by the "unstructured" working correlation structure or the
#' "model-based" covariance matrix produced by the working correlation
#' structure.
#' @param tol A numeric value specifying the tolerance when calculating the
#' inverse of the model-based covariance matrix.
#' Default is \code{.Machine$double.eps}.
#' @param envir An environment object to which find the \code{"geem"} object.
#'
#' @details (Out-dated) mCIC is calculated based on the trace of the product of
#' model-based covariance matrices, where the referent variances are produced
#' by the model using the unstructured (i.e., general) correlation structure.
#'
#' @return A numeric vector containing the mCIC value for the \code{"geem"}
#' object under consideration.
#'
#' @examples
#' data(ohio, package = "geepack")
#' fit <- geeM::geem(
#'   resp ~ age + smoke + age:smoke,
#'   id = id, data = ohio, family = stats::binomial(link = "logit"),
#'   corstr = "exchangeable"
#' )
#' geeM::mcic(fit)
#'
#' @export
mcic <- function(
  object,
  type = "unstructured",
  tol = .Machine$double.eps,
  envir = parent.frame()
) {
  # Checking function parameter types ####
  if (!methods::is(object, "geem")) {
    stop("object must be an object produced by geeM::geem")
  }
  if (!(class(tol) %in% c("numeric", "integer"))) {
    stop("tol must be a numeric or an integer value")
  }
  if (!(type %in% c("unstructured", "model-based"))) {
    stop("type must be unstructured or type")
  }
  if (!is.environment(envir)) {
    stop("envir must be an environment object")
  }

  # Obtaining the referent covariance matrix ####
  if (type == "unstructured") {
    # Pulling modeling call to re-fit models ####
    model_call <- object$call

    # Assigning corstr to unstructured ####
    model_call$corstr <- "unstructured"

    # Re-fitting model ####
    # returning NULL if error or warning occurrs with unstructured model fit
    model_unstructured <- tryCatch(
      expr = {
        model_unstructured <- eval(model_call, envir = envir)
      },
      warning = function(w) {
        NULL
      },
      error = function(e) {
        NULL
      }
    )

    if (is.null(model_unstructured)) {
      stop(
        paste0(
          "Unstructured model fit resulted in a warning or error. To examine",
          " warning or error, fit model outside of function call and specify",
          " corstr = 'unstructured'."
        )
      )
    }

    # Getting inverse of covariance matrices ####
    omega <- invert(as.matrix(model_unstructured$naiv.var), tol = tol)
  } else if (type == "model-based") {
    # Inverse of model-based covariance under working correlation ####
    omega <- invert(object$naiv.var, tol = tol)
  } else {
    stop("type specified does not allow calculation of mCIC")
  }

  # Obtaining the robust covariance matrix ####
  vr <- object$var

  # Calculating the measure ####
  measure <- sum(diag(x = crossprod(x = omega, y = vr)))

  # Returning mCIC
  return(measure)
}