% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snk.R
\name{snk}
\alias{snk}
\title{Multiple comparison: Student-Newman-Keuls test}
\usage{
snk(y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE, main = NULL)
}
\arguments{
\item{y}{Numeric or complex vector containing the response
variable.}

\item{trt}{Numeric or complex vector containing the
treatments.}

\item{DFerror}{Error degrees of freedom.}

\item{SSerror}{Error sum of squares.}

\item{alpha}{Significance level.}

\item{group}{TRUE or FALSE.}

\item{main}{Title.}
}
\value{
Returns the multiple comparison of means according
to the test of SNK.
}
\description{
\code{snk} Performs the test of SNK, for multiple
comparison of means.
}
\examples{
data(ex1)
attach(ex1)
crd(trat, ig, quali = TRUE, mcomp = "snk", sigT = 0.05)
}
\seealso{
\code{\link{scottknott}}, \code{\link{duncan}},
\code{\link{lsd}}, \code{\link{lsdb}}, \code{\link{ccboot}},
\code{\link{tukey}}, \code{\link{ccF}}.
}
\author{
Eric B Ferreira,
 \email{eric.ferreira@unifal-mg.edu.br}

Denismar Alves Nogueira

Portya Piscitelli Cavalcanti
}
