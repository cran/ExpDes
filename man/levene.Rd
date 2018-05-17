\name{levene}
\alias{levene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test for homogeneity of variances of Levene}
\description{Performs the test for homogeneity of variances of Levene (1960).}
\usage{
levene(trat, resp, t, r)
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
  
\value{
Returns the p-value of Levene's test of homogeneity of variances and its practical interpretation for significance level of 5\%.
}

\references{
LEVENE, H. Robust tests for equality of variances. in: Olkin, i.; ghurye, s.g.; hoeffding, w.; madow, w.g.; mann, h.b. (eds.). \emph{Contribution to Probability and Statistics. Stanford}, CA: Stanford University Press, pages 278-292, 1960.

NOGUEIRA, D, P.; PEREIRA, G, M. Desempenho de testes para homogeneidade de variâncias em delineamentos inteiramente casualizados. \emph{Sigmae}, Alfenas, v.2, n.1, p. 7-22. 2013.
}
  
\author{
Denismar Alves Nogueira
    
Eric Batista Ferreira
    
Marcos Costa de Paula
    
Mateus Pimenta Siqueira Lima
}

%\note{ }
  
  
\seealso{
See also: \code{\link{bartlett}}, \code{\link{layard}}, \code{\link{oneillmathews}}, \code{\link{samiuddin}}
}
    
\examples{
data(ex1)
attach(ex1)
crd(trat, ig, quali = FALSE, hvar='levene')
}
