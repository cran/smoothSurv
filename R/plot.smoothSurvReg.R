###########################################
#### AUTHOR:    Arnost Komarek         ####
####            (2003)                 ####
####                                   ####
#### FILE:      plot.smoothSurvReg.R   ####
####                                   ####
#### FUNCTIONS: plot.smoothSurvReg     ####
###########################################

### ====================================================================
### plot.smoothSurvReg: Plot objects of class 'smoothSurvReg'
### ====================================================================
## x ......... object of class 'smoothSurvReg'
## plot ...... T/F, do I want to plot it?
## resid ..... T/F, do I want to plot residuals on the x axe?
##             (midpoints are plotted for interval censored observ.)
## knots ..... T/F do I want to plot knots?
## compare ... T/F, do I want to draw standardized normal, logistic and extreme value densities?
## components. T/F, do I want to plot components of the mixture?
##             (if both compare and components are true than compare is set to FALSE)
## standard .. T/F
##               T ... I want to plot standardized fitted distrib. (with zero mean and unit variance)
##               F ... I want to plot distribution of alpha + sigma epsilon
## by ........ distance between two points of the grid for plotting
## toler.c ... tolerance to determine which G-spline coeff. are zero
## xlim
## ylim
## xlab
## ylab
## type
## lty
## main ...... standard arguments for 'plot' function
## sub
## bty
## ... ....... other arguments passed to 'plot' function
##
## RETURN: data.frame(x,y) to be used to produce the plot later on
plot.smoothSurvReg <- function(x,
                          plot = TRUE,
                          resid = TRUE,
                          knots = TRUE,
                          compare = TRUE,
                          components = FALSE,
                          standard = TRUE,
                          by = NULL,
                          toler.c = 1e-5,
                          xlim = NULL,
                          ylim = NULL,
                          xlab = expression(epsilon),
                          ylab = expression(paste("f(",epsilon,")", sep = "")),
                          type = "l",
                          lty = 1,
                          main = NULL,
                          sub = NULL,
                          bty = "n",
                          ...)
{
   if (x$fail >= 99){
        cat("No summary, smoothSurvReg failed.\n")
        return(invisible(x))
   }

   if (compare && components) compare <- FALSE
   if (!standard){
      compare <- FALSE
      components <- FALSE
   }

   ## Density function of extreme value distribution
   dextreme <- function(u){
      value <- exp(u-exp(u))
      return(value)
   }

   ## Some fitted values
   ccoef <- x$spline[["c coef."]]
   knots <- x$spline$Knot
   sigma0 <- x$spline[["SD basis"]][1]
   alpha <- x$adjust["(Intercept)","Value"]
   sigma <- x$adjust["Scale","Value"]
   shift <- x$error.dist$Mean[1]
   scale <- x$error.dist$SD[1]

   ## Standardized density of the fitted distribution
   dfitted <- function(u){
      normals <- scale*dnorm(scale*u + shift, mean = knots, sd = sigma0)
      value <- t(ccoef) %*% normals
      return(value)
   }

   ## Unstandardized density of the fitted distribution
   dfitted.un <- function(u){
      normals <- (1/sigma)*dnorm((u - alpha)/sigma, mean = knots, sd = sigma0)
      value <- t(ccoef) %*% normals
      return(value)
   }

   ## xlim
   small.c <- ccoef < toler.c
   if (is.null(xlim)){
      knot.min <- min(knots[!small.c])
      knot.max <- max(knots[!small.c])
      if (standard) xlim <- c(knot.min - 3*sigma0, knot.max + 3*sigma0)
      else          xlim <- c(sigma*knot.min + alpha - 3*sigma0, sigma*knot.max + alpha + 3*sigma0)
   }

   ## y values
   if (is.null(by)) by <- (xlim[2] - xlim[1])/100
   rooster <- seq(xlim[1], xlim[2], by = by)
   rooster2 <- matrix(rooster, ncol = 1)
   mean.extr <- -0.5772
   sig.extr <- pi/sqrt(6)
   sig.logis <- pi/sqrt(3)
   y.extreme <- sig.extr*dextreme(sig.extr*rooster + mean.extr)
   y.logis <- sig.logis*dlogis(sig.logis*rooster)
   y.normal <- dnorm(rooster)
   dens.use <- ifelse(standard, "dfitted", "dfitted.un")
   y.fitted <- apply(rooster2, 1, dens.use)

#   mm <- cumsum(y.fitted*rooster*by)
#   vv <- cumsum(y.fitted*(rooster^2)*by)
#   cat("Mean: "); print(mm[length(mm)])
#   cat("Variance: "); print(vv[length(vv)])

   ## ylim
   if (is.null(ylim)){
      if (standard){
         ymax <- max(y.extreme, y.logis, y.normal, y.fitted) + 0.05
      }
      else{
         ymax <- max(y.fitted) + 0.05
      }
      ylim <- c(-0.02, ymax)
   }


   ## main
   if (is.null(main)){
      ll <- round(x$degree.smooth$Lambda, digits = 3)
      logll <- round(x$degree.smooth[, "Log(Lambda)"], digits = 3)
      main <- paste("Error distribution,   ", "Log(Lambda) = ", logll, sep="")
   }


   ## sub
   if (is.null(sub)){
      aic <- round(x$aic, digits = 3)
      df <- round(x$degree.smooth$df, digits = 2)
      nparam <- x$degree.smooth[["Number of parameters"]]
      sub <- paste("AIC = ", aic, ",   df = ", df, ",   nParam = ", nparam, sep="")
   }

   ## Plot it
   if (plot){
     ltyplot <- c(lty,2,3,4)
     plot(rooster, y.fitted, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, 
            type=type, lty=lty, bty = bty, ...)
     title(sub = sub)
     if (compare && standard){
       lines(rooster, y.normal, lty=ltyplot[2])
       lines(rooster, y.extreme, lty=ltyplot[3])
       lines(rooster, y.logis, lty=ltyplot[4])
     }
   }

   ## Legend
   if (plot && compare){
      legend <- c("Fitted", "Normal", "Extreme val.", "Logistic")
      legend(xlim[1], ylim[2], legend, bty="n", lty=ltyplot)
   }

   ## Plot components
   if (plot && components){
      rooster2 <- matrix(rep(rooster, length(knots)), nrow = length(knots), byrow = TRUE)
      knots2 <- matrix(rep(knots, length(rooster)), ncol = length(rooster))
      y.comp <- dnorm(rooster2, mean = knots2, sd = sigma0)
      for (i in 1:length(knots)){
         lines(rooster, ccoef[i]*y.comp[i,], lty = 2)
      }
   }

   ## Plot residuals (if wanted)
   if (plot && resid){
      res <- resid(x)
      if (!standard){
         res[,1] <- alpha + sigma * res[,1]
         if (ncol(res) == 3) res[,2] <- alpha + sigma * res[,2]
      }

      if (ncol(res) == 3){   ## compute mid-points for interval censored residuals
         midp <- 0.5*(res[,1] + res[,2])
         res.use <- rep(NA, nrow(res))
         res.use[res[,3] == 3] <- midp[res[,3] == 3]
         res.use[res[,3] != 3] <- res[res[,3] != 3, 1]
         res <- cbind(res.use, res[,3])
      }
      n0 <- sum(res[,2] == 0); n1 <- sum(res[,2] == 1); n2 <- sum(res[,2] == 2); n3 <- sum(res[,2] == 3)
      ref0 <- 0.05; ref1 <- 0.01; ref2 <- 0.07; ref3 <- 0.03
      if (n0 > 0) points(res[res[,2] == 0, 1], rep(ref0, n0), pch = 4)
      if (n1 > 0) points(res[res[,2] == 1, 1], rep(ref1, n1), pch = 3)
      if (n2 > 0) points(res[res[,2] == 2, 1], rep(ref2, n2), pch = 2)
      if (n3 > 0) points(res[res[,2] == 3, 1], rep(ref3, n3), pch = 5)
   }

   if (plot && knots){
      if (!standard) knots <- sigma * knots + alpha
      knots.plot <- knots[!small.c]
      knots.plot <- knots.plot[knots.plot >= xlim[1] & knots.plot <= xlim[2]]
      zero <- rep(ylim[1], length(knots.plot))
      points(knots.plot, zero, pch = 19)
   }

   to.return <- data.frame(x = rooster, y = y.fitted)

   if (plot) return(invisible(to.return))
   else      return(to.return)

}

