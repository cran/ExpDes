\name{samiuddin}
\alias{samiuddin}

\title{Test for homogeneity of variances of Samiuddin}
\description{Performs the test for homogeneity of variances of Samiuddin (1976).}

\usage{
samiuddin(trat, resp, t, r)
}

\arguments{
\item{trat}{
   Numeric or complex vector containing treatments.
  }
  \item{resp}{
    Numeric or complex vector containing the response variable.
  }
  \item{t}{
Number of treatments.
  }
\item{r}{
Numeric or complex vector containing the number of replications of each treatment.
}
}

%\details{ }

\value{
Returns the p-value of Samiuddin's test of homogeneity of variances and its practical interpretation for significance level of 5\%.
}

\references{
SAMIUDDIN, M. Bayesian test of homogeneity of variance. \emph{Journal of the American Statistical Association}, 71(354):515-517, Jun. 1976.

NOGUEIRA, D, P.; PEREIRA, G, M. Desempenho de testes para homogeneidade de variâncias em delineamentos inteiramente casualizados. \emph{Sigmae}, Alfenas, v.2, n.1, p. 7-22. 2013.
}

\author{
Denismar Alves Nogueira
    
Eric Batista Ferreira
    
Marcos Costa de Paula
    
Mateus Pimenta Siqueira Lima
}

%\note{ }
%%  ~~further notes~~
  
  
  %% ~Make other sections like Warning with \section{Warning }{....} ~
  
\seealso{
See also: \code{\link{bartlett}}, \code{\link{layard}}, \code{\link{levene}}, \code{\link{oneillmathews}}
}
    
\examples{
data(ex1)
attach(ex1)
crd(trat, ig, quali = FALSE, hvar='samiuddin', sigF = 0.05)
}
