###############################################
#### AUTHOR:    Arnost Komarek             ####
####            (2003)                     ####
####                                       ####
#### FILE:      residuals.smoothSurvReg.R  ####
####                                       ####
#### FUNCTIONS: residuals.smoothSurvReg    ####
###############################################

### =================================================================================
### residuals.smoothSurvReg: Compute residuals for objects of class 'smoothSurvReg'
### =================================================================================
## object .......... object of class 'smoothSurvReg'
## ... ........ other arguments passed to 'residuals' function 
##              (it'shere only for compatibility with a generic function)
residuals.smoothSurvReg <- function(object, ...){
   ny <- ncol(object$y)
   nx <- ncol(object$x)

   regres <- object$regres[, "Value"]
   beta <- regres[1:nx]
   if (length(regres) > nx) scale <- regres["Scale"]
   else                     scale <- object$init.regres["Scale", "Value"]   ## scale was fixed

   eta <- object$x %*% beta

   y1 <- (object$y[,1] - eta)/scale
   if (ny > 2){
      y2 <- (object$y[,2] - eta)/scale
      y2[object$y[,3] == 2] <- 0
      y2[object$y[,3] == 0] <- 0
   }

   if (ny <= 2){
      out <- cbind(y1, object$y[,2])
      colnames(out) <- c("res", "censor")
   }
   else{
      out <- cbind(y1, y2, object$y[,3])
      colnames(out) <- c("res", "res2", "censor")
   }

   return(out)
}

