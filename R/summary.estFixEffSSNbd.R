#-------------------------------------------------------------------------------
#
#           summary.estFixEffSSNbd
#
#-------------------------------------------------------------------------------

#' Summary of S3 Method for Class "estFixEffSSNbd"
#'
#' @param object of class \code{estFixEffSSNbd}
#'
#' @return Table of estimates for fixed effects.
#'
#' @author Erin Peterson
#' @seealso \code{\link{estFixEffSSNbd}}
#'
#' @method summary estFixEffSSNbd
#' @export

summary.estFixEffSSNbd <-
function(object, ...)
{
    effnames <- rownames(object$betaHat)
    b.hat <- as.vector(object$betaHat)
    bhat.se <- sqrt(diag(object[[2]]))
## n.allxy <- object$sampinfo$obs.sample.size
   p <- length(effnames)

  if(any(rownames(b.hat) %in% effnames == FALSE)) {
            ## dataXY issue
  stop(cat("glmssn has computed estimates for",rownames(b.hat),"but the summary command expects estimates for",effnames,collapse=" "))
        }

	## bhat.0.NA <- rep(NA, times = length(effnames))
	## bhat.0.NA[setzero] <- 0
	## bhat.0.NA[!setzero & !setNA2] <- b.hat
	## NAvec <- rep(NA, times = length(effnames))
	## bhatse.0.NA <- NAvec
	## bhatse.0.NA[!setzero & !setNA2] <- bhat.se
	## tvec <- NAvec
        ## tvec[!setzero & !setNA2] <- b.hat/bhat.se
        tvec <- as.vector(b.hat/bhat.se)
	##pvec <- NAvec
	## pvec[!setzero & !setNA2] <-
	##	round(100000*(1 - pt(abs(b.hat/bhat.se), df = n.allxy - p))*2)/100000
        pvec <- as.vector(round(100000*(1 - pt(abs(b.hat/bhat.se), df = 1000 - p))*2)/100000)
        fixed.eff.est <- data.frame(FactorLevel = effnames, Estimate = b.hat,
			std.err = bhat.se, t.value = tvec, prob.t = pvec)

  outpt = list(fixed.effects.estimates =fixed.eff.est)
  class(outpt) <- "summary.estFixEffSSNbd"
  outpt
}


#' Print method for estFixEffSSNbd
#'
#' @param x
#' @param ...
#'
#' @export

print.estFixEffSSNbd <- function(x,...) {
    print(summary(x,...))
}

#' Print method for estFixEffSSNbd
#'
#' @param x
#' @param ...
#' @export
#'
print.summary.estFixEffSSNbd <- function(x,
  signif.stars = getOption("show.signif.stars"), ...)
{

  cat("\nCoefficients:\n")
  coef1 = x[[1]]
  coefs = coef1[,2:5]
  ##coefs <- x
  ##row.names(coefs) = coef1[,1]
  colnames(coefs) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  printCoefmat(coefs, signif.stars = signif.stars,
            na.print = "NA", ...)

  cat("\n")
  invisible(x)

}
