% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccF.r
\name{ccF}
\alias{ccF}
\title{Multiple comparison: Calinski and Corsten}
\usage{
ccF(y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE, main = NULL)
}
\arguments{
\item{y}{Numeric or complex vector containing the response
varible.}

\item{trt}{Numeric or complex vector containing the
treatments.}

\item{DFerror}{Error degrees of freedom.}

\item{SSerror}{Error sum of squares.}

\item{alpha}{Significance of the test.}

\item{group}{TRUE or FALSE.}

\item{main}{Title.}
}
\value{
Multiple means comparison for the Calinski and
Corsten test.
}
\description{
\code{ccF} Performs the Calinski and Corsten test based on
the F distribution.
}
\examples{
data(ex2)
attach(ex2)
rbd(trat, provador, aparencia, quali = TRUE, mcomp='ccf',
sigT = 0.05, sigF = 0.05)
}
\references{
CALI\'{N}SKI, T.; CORSTEN, L. C. A. Clustering
means in ANOVA by Simultaneous Testing. Biometrics. v. 41,
p. 39-48, 1985.
}
\author{
Eric B Ferreira,
 \email{eric.ferreira@unifal-mg.edu.br}

Patricia de Siqueira Ramos

Daniel Furtado Ferreira
}
