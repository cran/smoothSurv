###############################################
#### AUTHOR:    Arnost Komarek             ####
####            (2003)                     ####
####                                       ####
#### FILE:      survfit.smoothSurvReg.R    ####
####                                       ####
#### FUNCTIONS: survfit.smoothSurvReg      ####
###############################################

### ===================================================================================
### survfit.smoothSurvReg: Compute survivor curves for objects of class 'smoothSurvReg'
### ===================================================================================
## formula ... object of class 'smoothSurvReg' (name of the parameter is a little bit ambiguous
##             but I have to call in this way due to compatibility with a generic function)
## cov
## plot
## cdf
## by
## xlim
## ylim
## xlab
## ylab
## type
## lty
## main
## sub
## legend
## bty
## ... ....... other parameters passed to plot function
survfit.smoothSurvReg <- function(formula,
                                  cov = NULL,
                                  plot = TRUE,
                                  cdf = FALSE,
                                  by = NULL,
                                  xlim = NULL,
                                  ylim = c(0, 1),
                                  xlab = "t",
                                  ylab = expression(paste("S(","t",")", sep = "")),
                                  type = "l",
                                  lty = NULL,
                                  main = NULL,
                                  sub = NULL,
                                  legend = NULL,
                                  bty = "n",
                                  ...
                         )
{
   x <- formula
   if (x$fail >= 99){
        cat("No survivor curve, smoothSurvReg failed.\n")
        return(invisible(x))
   }

   ## Some fitted values
   ccoef <- x$spline[["c coef."]]
   knots <- x$spline$Knot
   sigma0 <- x$spline[["SD basis"]][1]
   alpha <- x$adjust["(Intercept)","Value"]
   sigma <- x$adjust["Scale","Value"]
   shift <- x$error.dist$Mean[1]
   scale <- x$error.dist$SD[1]

   ## Standardized survivor function of the fitted error distribution
   sfitted <- function(u){
      normals <- pnorm(scale*u + shift, mean = knots, sd = sigma0)
      value <- 1 - (t(ccoef) %*% normals)[1]
      return(value)
   }

   ## Unstandardized survivor function of the fitted error distribution
   ## (survivor function of alpha + sigma epsilon)
   sfitted.un <- function(u){
      normals <- pnorm((u - alpha)/sigma, mean = knots, sd = sigma0)
      value <- 1 - (t(ccoef) %*% normals)[1]
      return(value)
   }

   ## Regression parameters
   beta <- x$regres[,"Value"]
   if (x$estimated["Scale"]){
      ncov <- length(beta) - 3   # minus Intercept, Log(scale) and Scale
   }
   else{
      ncov <- length(beta) - 1   # minus Intercept
   }

   if (is.null(cov) && ncov > 0) cov <- matrix(rep(0, ncov), nrow = 1)
   if (ncov == 1)                cov <- matrix(cov, ncol = 1)
   if (ncov == 0)                cov <- NULL

   ## Different covariates combinations
   row.cov <- ifelse(is.null(dim(cov)), 1, dim(cov)[1])
   col.cov <- ifelse(is.null(dim(cov)),
                     ifelse(is.null(cov), 0, length(cov)),
                     dim(cov)[2])

   ## Linear predictor
   if (col.cov != ncov) stop("Incorrect cov parameter ")
   if (ncov > 0){
      beta <- matrix(beta[2:(ncov + 1)], nrow = ncov, ncol = 1)    # remove intercept and Log(scale) with scale
      cov <- matrix(cov, nrow = row.cov, ncol = col.cov)
      eta <- as.numeric(cov %*% beta)
   }
   else{
      eta <- 0
   }

   ## Grid
   if (is.null(xlim)){
      xmin <- 0
      xmax <- exp(max(x$y[,1]))
      xlim <- c(xmin, xmax)
   }
   if (is.null(by)){
      by <- (xlim[2] - xlim[1])/100
   }
   if (xlim[1] < 0) xlim[1] <- 0
   if (xlim[2] < 0) xlim[2] <- xlim[1] + 0.01

   grid <- seq(xlim[1], xlim[2], by) + 0.01

   ## Values
   etas <- matrix(rep(eta, rep(length(grid), row.cov)), ncol = row.cov)
   grid2 <- matrix(rep(grid, row.cov), ncol = row.cov)
   grid2 <- log(grid2) - etas
   Sfun <- list()
   for (i in 1:row.cov){
      grid3 <- matrix(grid2[,i], ncol = 1)
      Sfun[[i]] <- apply(grid3, 1, "sfitted.un")
      if (cdf) Sfun[[i]] <- 1 - Sfun[[i]]
   }

   ## lty
   if (is.null(lty)){
      lty <- 1:row.cov
   }

   ## main and sub
   if (is.null(main)){
      main <- ifelse(cdf, "Fitted Cum. Distribution Function", "Fitted Survivor Function")
   }
   if (is.null(sub)){
      aic <- round(x$aic, digits = 3)
      df <- round(x$degree.smooth$df, digits = 2)
      nparam <- x$degree.smooth[["Number of parameters"]]
      sub <- paste("AIC = ", aic, ",   df = ", df, ",   nParam = ", nparam, sep="")
   }

   ## Plot it
   if (plot){
      plot(grid, Sfun[[1]],
           type = type, lty = lty[1], ylim = ylim, xlab = xlab, ylab = ylab, bty = bty, ...)
      title(main = main, sub = sub)
      if (row.cov > 1){
         for (i in 2:row.cov){
            lines(grid, Sfun[[i]], lty = lty[i])
         }
      }
      if (is.null(legend)){
         leg <- numeric(2)
         leg[1] <- ifelse(cdf, xlim[1], xlim[2])
         leg[2] <- ylim[2]
         legjust <- numeric(2)
         legjust[1] <- ifelse(cdf, 0, 1)
         legjust[2] <- 1
         legend(leg[1], leg[2], legend = paste("cov", 1:row.cov, sep = ""), lty = lty, bty = "n",
            xjust = legjust[1], yjust = legjust[2])
     }
   }

   to.return <- data.frame(grid, Sfun[[1]])
   if (row.cov > 1)
   for (i in 2:row.cov){
      to.return <- cbind(to.return, Sfun[[i]])
   }
   names(to.return) <- c("x", paste("y", 1:row.cov, sep = ""))

   if (plot) return(invisible(to.return))
   else      return(to.return)
}

