% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oneillmathews.R
\name{oneillmathews}
\alias{oneillmathews}
\title{Test for homogeneity of variances of ONeill and Mathews (CRD)}
\usage{
oneillmathews(trat, resp, t, r)
}
\arguments{
\item{trat}{Numeric or complex vector containing treatments.}

\item{resp}{Numeric or complex vector containing the response
variable.}

\item{t}{Number of treatments.}

\item{r}{Numeric or complex vector containing the number of
replications of each treatment.}
}
\value{
Returns the p-value of ONeill and Mathews' test of
homogeneity of variances and its practical interpretation
for significance level of 5\%.
}
\description{
\code{oneillmathews} Performs the test for homogeneity of
variances of ONeill and Mathews (2000).
}
\examples{
data(ex1)
attach(ex1)
crd(trat, ig, quali = FALSE, hvar = "oneillmathews",
sigF = 0.05)
}
\references{
O'NEILL, M. E.; MATHEWS, K. L. A weighted least
squares approach to levene test of homogeneity of variance.
\emph{Australian e New Zealand Journal Statistical},
42(1):81-100, 2000.

NOGUEIRA, D, P.; PEREIRA, G, M. Desempenho de testes para
homogeneidade de variancias em delineamentos inteiramente
casualizados. \emph{Sigmae}, Alfenas, v.2, n.1, p. 7-22.
2013.
}
\seealso{
\code{\link{bartlett}}, \code{\link{layard}},
\code{\link{levene}}, \code{\link{samiuddin}}.
}
\author{
Eric B Ferreira,
 \email{eric.ferreira@unifal-mg.edu.br}
 @author Denismar Alves Nogueira
 @author Marcos Costa de Paula
 @author Mateus Pimenta Siqueira Lima
}
