\name{oneillmathews}
\alias{oneillmathews}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test for homogeneity of variances of ONeill and Mathews (CRD)}
\description{Performs the test for homogeneity of variances of ONeill and Mathews (2000).}
\usage{
oneillmathews(trat, resp, t, r)
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
Returns the p-value of ONeill and Mathews' test of homogeneity of variances and its practical interpretation for significance level of 5\%.
  }

\references{
O´NEILL, M. E.; MATHEWS, K. L. A weighted least squares approach to levene test of homogeneity of variance. \emph{Australian e New Zealand Journal Statistical}, 42(1):81-100, 2000.

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
See also: \code{\link{bartlett}}, \code{\link{layard}}, \code{\link{levene}}, \code{\link{samiuddin}}
}
    
\examples{
data(ex1)
attach(ex1)
crd(trat, ig, quali = FALSE, hvar='oneillmathews', sigF = 0.05)
}
