#' @title Variance information criteria
#' @description Calculating the variance information criteria (VIC) based on
#' Boonstra and Cavanaugh (in works).
#'
#' @param object An object of class \code{"geem"} representing the fit.
#' @param robust A logical denoting whether the robust covariance matrices
#' should be used. Default is `FALSE`.
#' @param envir An environment object to which find the `"geem"` object.
#' Default is `parent.frame()`.
#'
#' @details VIC is calculated based on the sum of the maximum of variance
#' ratios, where the referent variances are produced by the model using
#' the unstructured (i.e., general) correlation structure.
#'
#' When the robust covariance matrices are used through VIC may be affected
#' by the variability and efficiency, while VIC based on the model-based
#' covariance matrices may be affected by bias. *In the future, other robust
#' covariance estimators may be used.*
#'
#' @return A numeric vector containing the VIC value for the \code{"geem"}
#' object under consideration.
#'
#' @examples
#' data(ohio, package = "geepack")
#' fit <- geeM::geem(
#'   resp ~ age + smoke + age:smoke,
#'   id = id, data = ohio, family = stats::binomial(link = "logit"),
#'   corstr = "exchangeable"
#' )
#' geeM::vic(fit)
#'
#' @export
vic <- function(object, robust = FALSE, envir = parent.frame()) {
  # Checking function parameter types ####
  if (!methods::is(object, "geem")) {
    stop("object must be an object produced by geeM::geem")
  }
  if (!is.logical(robust)) {
    stop("robust must a be a logical.")
  }
  if (!is.environment(envir)) {
    stop("envir must be an environment object")
  }

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

  # Calculating VIC ####
  ## Getting covariance matrices ####

  ### Covariance under unstructured ####
  if (robust) {
    cov_unstructured <- model_unstructured$var
  } else {
    cov_unstructured <- model_unstructured$naiv.var
  }

  ### Covariance under fitted model ####
  if (robust) {
    cov_object <- object$var
  } else {
    cov_object <- object$naiv.var
  }

  ## Calculation ratio of variances ####
  ratio1 <- diag(cov_unstructured) / diag(cov_object)
  ratio2 <- diag(cov_object) / diag(cov_unstructured)

  ## Checking ratios ####
  ## Error if the covariance matrices are non-p.d.
  if (any(c(ratio1 < 0, ratio2 < 0))) {
    stop("Non-positive definite covariance matrices present.")
  }

  ## Calculating VIC
  measure <- sum(pmax(ratio1, ratio2))

  # Returning VIC
  return(measure)
}