pkgname <- "smoothSurv"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('smoothSurv')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("a2c")
### * a2c

flush(stderr()); flush(stdout())

### Name: a2c
### Title: Work Function for 'smoothSurvReg'
### Aliases: a2c
### Keywords: internal utilities

### ** Examples

ccoef <- c(0.1, 0.2, 0.15, 0.3, 0.25)

### Compute 'a' counterparts
acoef <- c2a(ccoef, 1)
print(acoef)

### And back 'c', ccoef2 should be same as ccoef
ccoef2 <- a2c(acoef)
print(ccoef2)



cleanEx()
nameEx("c2a")
### * c2a

flush(stderr()); flush(stdout())

### Name: c2a
### Title: Work Function for 'smoothSurvReg'
### Aliases: c2a
### Keywords: internal utilities

### ** Examples

ccoef <- c(0.1, 0.2, 0.15, 0.3, 0.25)

### Compute 'a' counterparts
acoef <- c2a(ccoef, 1)
print(acoef)

### And back 'c', ccoef2 should be same as ccoef
ccoef2 <- a2c(acoef)
print(ccoef2)



cleanEx()
nameEx("dextreme")
### * dextreme

flush(stderr()); flush(stdout())

### Name: extreme value
### Title: Density of the Extreme Value Distribution of a Minimum.
### Aliases: dextreme dstextreme
### Keywords: distribution

### ** Examples

dextreme(1, (sqrt(6)/pi)*0.5772, sqrt(6)/pi)
dstextreme(1)        ## approximately same result as on the previous row



cleanEx()
nameEx("dstlogis")
### * dstlogis

flush(stderr()); flush(stdout())

### Name: standardized logistic
### Title: Density of Standardized Logistic Distribution.
### Aliases: dstlogis
### Keywords: distribution

### ** Examples

dstlogis(0)
dstlogis(seq(-3, 3, 0.2))



cleanEx()
nameEx("eval.Gspline")
### * eval.Gspline

flush(stderr()); flush(stdout())

### Name: eval.Gspline
### Title: Evaluate a G-spline in a grid of values
### Aliases: eval.Gspline
### Keywords: dplot

### ** Examples

  spline <- minPenalty(knots=seq(-4.2, 4.2, by=0.3), sdspline=0.2, difforder=3)$spline
  values <- eval.Gspline(spline, seq(-4.5, 4.5, by=0.05))
  plot(values, type="l", bty="n", lwd=3)



cleanEx()
nameEx("find.c")
### * find.c

flush(stderr()); flush(stdout())

### Name: find.c
### Title: Work Function for 'smoothSurvReg'
### Aliases: find.c
### Keywords: internal utilities

### ** Examples

knots <- seq(-4, 4, 0.5)
sd0 <- 0.3
ccoef <- find.c(knots, sd0, dist = "dstlogis")

### We plot the approximation together with the truth
###
grid <- seq(-4, 4, 0.05)
truth <- dstlogis(grid)

### Following lines compute the values of the approximation
grid.big <- matrix(grid, nrow = length(grid), ncol = length(knots))
knots.big <- matrix(knots, nrow = length(grid), ncol = length(knots), byrow = TRUE)
normals <- dnorm(grid.big, mean = knots.big, sd = sd0)
approx <- normals %*% ccoef

### Plot it
plot(grid, approx, type = "l", xlab = "y", ylab = "f(y)", bty = "n")
lines(grid, truth, lty = 2)
legend(-4, 0.35, c("approx", "truth"), lty = 1:2, bty = "n")



cleanEx()
nameEx("give.c")
### * give.c

flush(stderr()); flush(stdout())

### Name: give.c
### Title: Work Function for 'smoothSurvReg'
### Aliases: give.c
### Keywords: internal utilities

### ** Examples

knots <- seq(-4, 4, 0.5)
sd0 <- 0.3
ccoef <- find.c(knots, sd0, dist = "dstlogis")

last.three <- c(3, 7, 10)
c.rest <- ccoef[-last.three]
ccoef2 <- give.c(knots, sd0, last.three, c.rest)

print(ccoef)
print(ccoef2)    ## Almost no change



cleanEx()
nameEx("minPenalty")
### * minPenalty

flush(stderr()); flush(stdout())

### Name: minPenalty
### Title: Minimize the penalty term under the two (mean and variance)
###   constraints
### Aliases: minPenalty
### Keywords: optimize

### ** Examples

optimum <- minPenalty(knots=seq(-4.2, 4.2, by = 0.3), sdspline=0.2, difforder=3)
where <- optimum$spline
print(where)
show <- eval.Gspline(where, seq(-4.2, 4.2, by=0.05))
plot(show, type="l", bty="n", lwd=2)



cleanEx()
nameEx("piece")
### * piece

flush(stderr()); flush(stdout())

### Name: piece
### Title: Left Continuous Piecewise Constant Function with a Finite
###   Support.
### Aliases: piece
### Keywords: utilities

### ** Examples

my.breaks <- c(-2, 1.5, 4, 7)
my.values <- c(0.5, 0.9, -2)
grid <- seq(-3, 8, by = 0.25)
piece(grid, my.breaks, my.values)



cleanEx()
nameEx("smoothSurvReg")
### * smoothSurvReg

flush(stderr()); flush(stdout())

### Name: smoothSurvReg
### Title: Regression for a Survival Model with Smoothed Error Distribution
### Aliases: smoothSurvReg
### Keywords: survival smooth

### ** Examples

##### EXAMPLE 1:  Common scale
##### ========================
### We generate interval censored data and fit a model with few artificial covariates
set.seed(221913282)
x1 <- rbinom(50, 1, 0.4)                                         ## binary covariate
x2 <- rnorm(50, 180, 10)                                         ## continuous covariate
y1 <- 0.5*x1 - 0.01*x2 + 0.005 *x1*x2 + 1.5*rnorm(50, 0, 1)      ## generate log(T), left limit
t1 <- exp(y1)                                                    ## left limit of the survival time
t2 <- t1 + rgamma(50, 1, 1)                                      ## right limit of the survival time
surv <- Surv(t1, t2, type = "interval2")                         ## survival object

## Fit the model with an interaction
fit1 <- smoothSurvReg(surv ~ x1 * x2, logscale = ~1, info = FALSE, lambda = exp(2:(-1)))

## Print the summary information
summary(fit1, spline = TRUE)

## Plot the fitted error distribution
plot(fit1)

## Plot the fitted error distribution with its components
plot(fit1, components = TRUE)

## Plot the cumulative distribution function corresponding to the error density
survfit(fit1, cdf = TRUE)

## Plot survivor curves for persons with (x1, x2) = (0, 180) and (1, 180)
cov <- matrix(c(0, 180, 0,   1, 180, 180), ncol = 3, byrow = TRUE)
survfit(fit1, cov = cov)

## Plot hazard curves for persons with (x1, x2) = (0, 180) and (1, 180)
cov <- matrix(c(0, 180, 0,   1, 180, 180), ncol = 3, byrow = TRUE)
hazard(fit1, cov = cov)

## Plot densities for persons with (x1, x2) = (0, 180) and (1, 180)
cov <- matrix(c(0, 180, 0,   1, 180, 180), ncol = 3, byrow = TRUE)
fdensity(fit1, cov = cov)

## Compute estimates expectations of survival times for persons with
## (x1, x2) = (0, 180), (1, 180), (0, 190), (1, 190), (0, 200), (1, 200)
## and estimates of a difference of these expectations:
## T(0, 180) - T(1, 180), T(0, 190) - T(1, 190), T(0, 200) - T(1, 200),
cov1 <- matrix(c(0, 180, 0,   0, 190, 0,   0, 200, 0), ncol = 3, byrow = TRUE)
cov2 <- matrix(c(1, 180, 180,   1, 190, 190,   1, 200, 200), ncol = 3, byrow = TRUE)
print(estimTdiff(fit1, cov1 = cov1, cov2 = cov2))


##### EXAMPLE 2:  Scale depends on covariates
##### =======================================
### We generate interval censored data and fit a model with few artificial covariates
set.seed(221913282)
x1 <- rbinom(50, 1, 0.4)                                         ## binary covariate
x2 <- rnorm(50, 180, 10)                                         ## continuous covariate
x3 <- runif(50, 0, 1)                                            ## covariate for the scale parameter
logscale <- 1 + x3
scale <- exp(logscale)
y1 <- 0.5*x1 - 0.01*x2 + 0.005 *x1*x2 + scale*rnorm(50, 0, 1)    ## generate log(T), left limit
t1 <- exp(y1)                                                    ## left limit of the survival time
t2 <- t1 + rgamma(50, 1, 1)                                      ## right limit of the survival time
surv <- Surv(t1, t2, type = "interval2")                         ## survival object

## Fit the model with an interaction
fit2 <- smoothSurvReg(surv ~ x1 * x2, logscale = ~x3, info = FALSE, lambda = exp(2:(-1)))

## Print the summary information
summary(fit2, spline = TRUE)

## Plot the fitted error distribution
plot(fit2)

## Plot the fitted error distribution with its components
plot(fit2, components = TRUE)

## Plot survivor curves for persons with (x1, x2) = (0, 180) and (1, 180)
## x3 = 0.8 and 0.9
cov <- matrix(c(0, 180, 0,   1, 180, 180), ncol = 3, byrow = TRUE)
logscale.cov <- c(0.8, 0.9)
survfit(fit2, cov = cov, logscale.cov = logscale.cov)

## Plot hazard curves for persons with (x1, x2) = (0, 180) and (1, 180)
## x3 = 0.8 and 0.9
cov <- matrix(c(0, 180, 0,   1, 180, 180), ncol = 3, byrow = TRUE)
logscale.cov <- c(0.8, 0.9)
hazard(fit2, cov = cov, logscale.cov=c(0.8, 0.9))

## Plot densities for persons with (x1, x2) = (0, 180) and (1, 180)
## x3 = 0.8 and 0.9
cov <- matrix(c(0, 180, 0,   1, 180, 180), ncol = 3, byrow = TRUE)
logscale.cov <- c(0.8, 0.9)
fdensity(fit2, cov = cov, logscale.cov = logscale.cov)


## More involved examples can be found in script files
## used to perform analyses  and draw pictures 
## presented in above mentioned references.
## These scripts and some additional files can be found as *.tar.gz files
## in the /inst/doc directory of this package.
##



cleanEx()
nameEx("std.data")
### * std.data

flush(stderr()); flush(stdout())

### Name: std.data
### Title: Standardization of the Data
### Aliases: std.data
### Keywords: manip

### ** Examples

variable1 <- rnorm(30)
variable2 <- rbinom(30, 1, 0.4)
variable3 <- runif(30)
data.example <- data.frame(variable1, variable2, variable3)
## We standardize only the first and the third column.
data.std <- std.data(data.example, c("variable1", "variable3"))
print(data.std)
print(c(mean(data.std$variable1), sd(data.std$variable1)))
print(c(mean(data.std$variable3), sd(data.std$variable3)))



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
