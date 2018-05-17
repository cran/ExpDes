\name{bartlett}
\alias{bartlett}

\title{Test for Homogeneity of Variances: Bartlett}
\description{Performs the test for homogeneity of variances of Bartlett (1937).}
\usage{
bartlett(trat, resp, t, r)
}

\arguments{
\item{trat}{
Numeric or complex vector containing the treatments.
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
Returns the p-value of Bartlett's test of homogeneity of variances and its practical interpretation for 5\% of significance.
}

\references{
BARTLETT, M. S. Properties of sufficiency and statistical tests. \emph{Proceedings of the Royal Statistical Society - Serie A}, 60:268-282, 1937.

NOGUEIRA, D, P.; PEREIRA, G, M. Desempenho de testes para homogeneidade de variâncias em delineamentos inteiramente casualizados. \emph{Sigmae}, Alfenas, v.2, n.1, p. 7-22. 2013.
}

%% ~put references to the literature/web site here ~

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
See also: \code{\link{levene}}, \code{\link{layard}}, \code{\link{oneillmathews}}, \code{\link{samiuddin}}
}

\examples{
data(ex1)
attach(ex1)
crd(trat, ig, quali = FALSE, hvar='bartlett')
}