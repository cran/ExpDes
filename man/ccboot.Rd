% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccboot.r
\name{ccboot}
\alias{ccboot}
\title{Multiple comparison: Bootstrap}
\usage{
ccboot(
  y,
  trt,
  DFerror,
  SSerror,
  alpha = 0.05,
  group = TRUE,
  main = NULL,
  B = 1000
)
}
\arguments{
\item{y}{Numeric or complex vector containing the response
varible.}

\item{trt}{Numeric or complex vector containing the
treatments.}

\item{DFerror}{Error degrees of freedom.}

\item{SSerror}{Error sum of squares.}

\item{alpha}{Significance of the test.}

\item{group}{TRUE or FALSE}

\item{main}{Title}

\item{B}{Number of bootstrap resamples.}
}
\value{
Multiple means comparison for the bootstrap test.
}
\description{
\code{ccboot} Performs the Ramos and Ferreira (2009)
multiple comparison bootstrap test.
}
\examples{
data(ex1)
attach(ex1)
crd(trat, ig, quali = TRUE, mcomp='ccboot', sigF = 0.05)
}
\references{
RAMOS, P. S., FERREIRA, D. F. Agrupamento de
medias via bootstrap de populacoes normais e nao-normais,
Revista Ceres, v.56, p.140-149, 2009.
}
\author{
Eric B Ferreira,
 \email{eric.ferreira@unifal-mg.edu.br}

Patricia de Siqueira Ramos

Daniel Furtado Ferreira
}
