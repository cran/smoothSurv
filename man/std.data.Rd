\name{std.data}
\alias{std.data}
\title{
   Standardization of the Data
}
\description{
   Chosen columns of a \code{data.frame} are standardized by
   using a sample mean and sample standard deviation.
}
\usage{
std.data(datain, cols)
}
\arguments{
   \item{datain}{
      Input data frame.
   }
   \item{cols}{
      A character vector with names of the columns to be standardized.
   }
}
\details{
   For each chosen column the sample mean and the sample standard deviation
   is computed and then used to standardize the column. Missing values are 
   ignored when computing the mean and the standard deviation.
}
\value{
   A \code{data.frame} with same variables as the input \code{data.frame}.
   Chosen columns are standardized.
 }
 \seealso{
   \code{\link{scale}}
 }  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek@med.kuleuven.be}
}
\examples{
variable1 <- rnorm(30)
variable2 <- rbinom(30, 1, 0.4)
variable3 <- runif(30)
data.example <- data.frame(variable1, variable2, variable3)
## We standardize only the first and the third column.
data.std <- std.data(data.example, c("variable1", "variable3"))
print(data.std)
print(c(mean(data.std$variable1), sd(data.std$variable1)))
print(c(mean(data.std$variable3), sd(data.std$variable3)))
}
\keyword{manip}
