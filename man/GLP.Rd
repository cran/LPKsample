\name{GLP}
\alias{GLP}
\alias{graph.clust.test}

\title{
A function to perform K-sample test using GLP algorithm
}
\description{
This function performs the GLP multivariate K-sample learning.
}
\usage{
GLP(X,y,m.max=4,alpha=0.05,c=0.5,clust.alg='kmeans',LPtype='smooth',
smooth.method='BIC',return.LPT=F,return.clust=F)
}

\arguments{
  \item{X}{ A \eqn{n}-by-\eqn{d} matrix of the observations, the observations should be grouped by their respective classes}
  \item{y}{ A length \eqn{n} vector indicating the sample class}
  \item{m.max}{ An integer, maximum order of LP component to investigate, default: 4 }
  \item{alpha}{ Numeric, confidence level \eqn{\alpha} , default: 0.05}
  \item{c}{ Numeric, parameter for polynomial kernel, default: 0.5 }
  \item{clust.alg}{ \code{"mclust"} or \code{"kmeans"}; algorithm used for clustering in graph community detection }
  \item{LPtype}{ set to \code{"smooth"} to smooth the LP comeans for final testing phase, default: "smooth"}
  \item{smooth.method}{ \code{"AIC"} or \code{"BIC"}; Method used for smoothing the LP comeans, default: "BIC"}
  \item{return.LPT}{ logical, whether or not to return the data driven covariate matrix, default: FALSE}
  \item{return.clust}{ logical, whether or not to return the class labels assigned by graph community detection, default: FALSE}
}

\value{
 A list containing the following items:
   \item{GLP}{Overall GLP statistics}
   \item{pval}{Overall P-value}
   \item{table}{The GLP component table indicating the significance of each component}
   \item{components}{significant eLP components for the data set}
   \item{LPT}{(optional) matrix of data driven covariates}
   \item{clust}{(optional) class labels assigned by graph community detection}
}

\references{
Mukhopadhyay, S. and Wang, K. (2017), "A Nonparametric Approach to High-dimensional K-sample Comparison Problem".
}
\author{ Mukhopadhyay, S. and Wang, K.}

\examples{


  ##1.muiltivariate normal distribution with only mean difference:
  ##generate data, n1=n2=10, dimension 25
   X1<-matrix(rnorm(250,mean=0,sd=1),10,25)
   X2<-matrix(rnorm(250,mean=0.5,sd=1),10,25)
   y<-c(rep(1,10),rep(2,10))
   X<-rbind(X1,X2)
  ##GLP test:
   locdiff.test<-GLP(X,y,m.max=4)

  \dontrun{
  ##2.Leukemia data example
   data(leukemia)
   attach(leukemia)
   leukemia.test<-GLP(X,class,clust.alg="mclust")
  ##confirmatory results:
   leukemia.test$GLP  # overall statistic
   #[1] 0.2723404
   leukemia.test$pval # overall p-value
   #[1] 9.504488e-06
  ##exploratory outputs:
   leukemia.test$table  # rows as shown in Table 3 of reference
   #     component  comp.GLP       pvalue
   #[1,]         1 0.2723404 9.504488e-06
   #[2,]         2 0.0000000 1.000000e+00
   #[3,]         3 0.0000000 1.000000e+00
   #[4,]         4 0.0000000 1.000000e+00
  }
}


\keyword{ GLP main function }