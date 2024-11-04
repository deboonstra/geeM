#' @title Modified correlation information criteria
#' @description Calculating the modified correlation information criteria (mCIC)
#' based on Boonstra and Cavanaugh (in works).
#'
#' @param object An object of class \code{"geem"} representing the fit.
#' @param tol A numeric value specifying the tolerance when calculating the
#' inverse of the model-based covariance matrix.
#' Default is \code{.Machine$double.eps}.
#' @param envir An environment object to which find the \code{"geem"} object.
#'
#' @details mCIC is calculated based on the trace of the product of
#' robust covariance matrices produced from the fitted model and a model-based
#' referent covariance matrix, where the referent covariances are produced
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
  if (!is.environment(envir)) {
    stop("envir must be an environment object")
  }

  # Obtaining the referent covariance matrix ####
  ## Pulling modeling call to re-fit models ####
  model_call <- object$call

  ## Assigning corstr to unstructured ####
  model_call$corstr <- "unstructured"

  ## Re-fitting model ####
  ## returning NULL if error or warning occurrs with unstructured model fit
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

  ## Getting inverse of covariance matrices ####
  omega <- invert(as.matrix(model_unstructured$naiv.var), tol = tol)

  # Obtaining the robust covariance matrix ####
  vr <- object$var

  # Calculating the measure ####
  measure <- sum(diag(x = crossprod(x = omega, y = vr)))

  # Returning mCIC
  return(measure)
}