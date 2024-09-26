#' @title Variance information criteria
#' @description Calculating the variance information criteria (VIC) based on
#' Boonstra and Cavanaugh (in works).
#'
#' @param object An object of class \code{"geem"} representing the fit.
#' @param envir An environment object to which find the \code{"geem"} object.
#'
#' @details VIC is calculated based on the sum of the maximum of model-based
#' variance ratios, where the referent variances are produced by the model using
#' the unstructured (i.e., general) correlation structure.
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
vic <- function(object, envir = parent.frame()) {
  # Checking function parameter types ####
  if (!methods::is(object, "geem")) {
    stop("object must be an object produced by geeM::geem")
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

  ### Model-based covariance under unstructured ####
  mb_unstructured <- model_unstructured$naiv.var

  ### Model-based covariance under fitted model ####
  mb_object <- object$naiv.var

  ## Calculation ratio of variances ####
  ratio1 <- diag(mb_unstructured) / diag(mb_object)
  ratio2 <- diag(mb_object) / diag(mb_unstructured)

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