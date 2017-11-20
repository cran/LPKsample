\name{LP.comean}
\alias{LP.comean}
\alias{LP.smooth}
\title{Function to find LP-comeans}
\description{
The function computes the LP comeans between \code{x} and \code{y}.
}
\usage{
LP.comean(x, y, type = "smooth", smooth.method = "BIC")
}

\arguments{
  \item{x}{vector, observations of an univariate random variable}
  \item{y}{vector, observations of another univariate random variable}
  \item{type}{set to "smooth" to enable smoothing of LP comeans}
  \item{smooth.method}{\code{"AIC"} or \code{"BIC"}; Method used for smoothing the LP comeans }
}

\value{
  A list containing:
  \item{LPINFOR }{The test statistics based on LP comeans}
  \item{p.val}{Test p-value}
  \item{LP.matrix}{LP comean matrix}
}
\references{Mukhopadhyay, S. and Wang, K. (2017), "A Nonparametric Approach to High-dimensional K-sample Comparison Problem".

Parzen, E. and Mukhopadhyay, S. (2012)  "Modeling, Dependence, Classification, United Statistical Science, Many Cultures".
}
\author{ Mukhopadhyay, S. and Wang, K.}

\examples{
#example: LP-comean for two simple vectors:
 y<-c(1,2,3,4,5)
 z<-c(0,-1,-1,3,4)
 comeanYZ=LP.comean(y,z)
#sum square statistics of LP comean:
 comeanYZ$LPINFOR
#p-value:
 comeanYZ$p.val
#comean matrix:
 comeanYZ$LP.matrix
}

\keyword{ GLP other functions }