#-------------------------------------------------------------------------------
#
#           summary.estCovParSSNbd
#
#-------------------------------------------------------------------------------

#' Summary of S3 Method for Class "estCovParSSNbd"
#'
#' Summary of S3 Method for Class "estCovParSSNbd"
#'
#' @param object of class \code{estCovParSSNbd}
#'
#' @return covariance parameter estimates for each covariance model.
#'
#' @author Erin Peterson
#' @seealso \code{\link{estCovParSSNbd}}
#'
#' @method summary estCovParSSNbd
#' @export


summary.estCovParSSNbd <-
    function(object, ...)

{
  catCall =  object$mfcall


  covmodels <- object[[1]]
  covmodels <- data.frame(Covariance.Model=attributes(
        	  covmodels)$terms,Parameter=attributes(covmodels)$type,
        		Estimate=covmodels)

 outpt = list(catCall = catCall,
    covariance.parameter.estimates = covmodels)
  class(outpt) <- "summary.estCovParSSNbd"
  outpt
}

#' Print S3 method for estCovParSSNbd
#'
#' @param x
#' @param ...
#'
#' @export
print.estCovParSSNbd <- function(x,...) {
    print(summary(x,...))
}

#' Print S3 method for estCovParSSNbd
#'
#' @param x
#'
#' @export
#'
print.summary.estCovParSSNbd <- function(x,
  digits = max(3L, getOption("digits") - 3L))
{
  cat("\nCall:\n", paste(deparse(x$catCall), sep = "\n", collapse = "\n"),
        "\n", sep = "")
  cat("\nCovariance Parameters:\n")
  cpe = x$covariance.parameter.estimates
  print(cpe, row.names = FALSE)
  rse = sqrt(sum(cpe[cpe[,"Parameter"] == "parsill","Estimate"]))
  cat("\nResidual standard error:", rse)

  cat("\n")
  invisible(x)
}
