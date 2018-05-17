\name{layard}
\alias{layard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test for homogeneity of variances of Layard}
\description{Performs the test for homogeneity of variances of Layard for Jackknife (1973).}
\usage{
layard(trat, resp, t, r)
}
%- maybe also 'usage' for other objects documented here.
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
%%  ~~ If necessary, more details than the description above ~~
  
\value{
Returns the p-value of the Layard test of homogeneity of variances and its practical interpretation for the significance level of 5\%.
}

\references{
LAYARD, M. N. J. Robust large-sample tests for homogeneity of variances. \emph{Journal of the American Statistical Association}, v.68, n.341, p.195-198, 1973.

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
See also: \code{\link{bartlett}}, \code{\link{samiuddin}}, \code{\link{levene}}, \code{\link{oneillmathews}}
}
    
\examples{
data(ex1)
attach(ex1)
crd(trat, ig, quali = FALSE, hvar='layard')
}