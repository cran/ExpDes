\name{ccf}
\alias{ccf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiple comparison: Calinski and Corsten}
\description{Performs the Calinski and Corsten test based on the F distribution.}
\usage{
ccf(y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE, main = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Numeric or complex vector containing the response varible.}
  \item{trt}{Numeric or complex vector containing the treatments.}
  \item{DFerror}{Error degrees of freedom.}
  \item{SSerror}{Error sum of squares.}
  \item{alpha}{Significance of the test.}
  \item{group}{TRUE or FALSE}
  \item{main}{Title}
}

\value{
Multiple means comparison for the Calinski and Corsten test.
}

\references{
CALI\'{N}SKI, T.; CORSTEN, L. C. A. Clustering means in ANOVA by Simultaneous Testing. Biometrics. v. 41, p. 39-48, 1985.
}
\author{
Patricia de Siqueira Ramos

Daniel Furtado Ferreira

Eric Batista Ferreira
}
%\note{
%%  ~~further notes~~

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{snk}}, \code{\link{duncan}}, \code{\link{lsd}}, \code{\link{lsdb}}, \code{\link{tukey}}, \code{\link{scottknott}}, \code{\link{ccboot}} 
}

\examples{
data(ex2)
attach(ex2)
rbd(trat, provador, aparencia, quali = TRUE, mcomp='ccf', sigT = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
