###########################################
#### AUTHOR:    Arnost Komarek         ####
####            25/02/2004             ####
###             03/05/2004             ####
####                                   ####
#### FILE:      print.estimTdiff.R     ####
####                                   ####
#### FUNCTIONS: print.estimTdiff       ####
###########################################

### ====================================================================
### print.estimTdiff: Print objects of class 'estimTdiff'
### ====================================================================
## x .......... object of class 'estimTdiff'
## digits ..... # of printed digits
## ... ........ other arguments passed to 'print' function
print.estimTdiff <- function(x, digits = min(options()$digits, 4), ...)
{

    if(is.null(digits))
        digits <- min(options()$digits, 4)

    cat("\nCovariate Values Compared:\n")
    if (is.null(attr(x, "cov1"))){
      cat("   Only intercept was in the model.\n")
    }
    else{
      cat("   Covariate values for T1:\n")
      print(attr(x, "cov1"), digits = digits, ...)
      cat("\n")
      cat("   Covariate values for T2:\n")
      print(attr(x, "cov2"), digits = digits, ...)            
    }
    cat("\n")

    if (!is.null(attr(x, "logscale.cov1"))){
      cat("\nLog-Scale Covariate Values Compared:\n")
      cat("   Log-scale covariate values for T1:\n")
      print(attr(x, "logscale.cov1"), digits = digits, ...)
      cat("\n")
      cat("   Log-scale covariate values for T2:\n")
      print(attr(x, "logscale.cov2"), digits = digits, ...)
      cat("\n")
    }         
    
    cat("\nEstimates of Expectations:\n")
    
    Z1 <- x$ET1 / x$sd.ET1
    Z2 <- x$ET2 / x$sd.ET2
    Zdiff <- x$diffT / x$sd.diffT

    p1 <- 2 * pnorm(-abs(Z1))
    p2 <- 2 * pnorm(-abs(Z2))
    pdiff <- 2 * pnorm(-abs(Zdiff))

    show1 <- data.frame(x$ET1, x$sd.ET1, Z1, p1)
    show2 <- data.frame(x$ET2, x$sd.ET2, Z2, p2)
    show3 <- data.frame(x$diffT, x$sd.diffT, Zdiff, pdiff)        

    colnames(show1) <- c("T1", "Std.Error", "Z", "p")
    colnames(show2) <- c("T2", "Std.Error", "Z", "p")
    colnames(show3) <- c("T1 - T2", "Std.Error", "Z", "p")    
    rownames(show1) <- paste("Value ", 1:length(x$ET1), sep = "")
    rownames(show2) <- paste("Value ", 1:length(x$ET1), sep = "")
    rownames(show3) <- paste("Value ", 1:length(x$ET1), sep = "")    
    print(show1, digits = digits, ...); cat("\n")
    print(show2, digits = digits, ...); cat("\n")
    print(show3, digits = digits, ...); cat("\n")    
}

