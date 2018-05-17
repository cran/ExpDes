\name{scottknott}
\Rdversion{1.1}
\alias{scottknott}

\title{Multiple comparison: Scott-Knott test}
\description{Performs the test of Scott-Knott, for multiple comparison of means.}
\usage{
scottknott(y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE, main = NULL)
}

\arguments{
  \item{y}{Numeric or complex vector containing the response variable.}
  \item{trt}{Numeric or complex vector containing the treatments.}
  \item{DFerror}{Error degrees of freedom.}
  \item{SSerror}{Error sum of squares.}
  \item{alpha}{Significance level.}
  \item{group}{TRUE or FALSE}
  \item{main}{Title}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~

\value{
Returns the multiple comparison of means according to the test of Scott-Knott.
}
\references{
RAMALHO, M. A. P.; FERREIRA, D. F.; OLIVEIRA, A. C. de. Experimentacao em Genetica e Melhoramento de Plantas. 2a ed. Lavras: UFLA. 2005. 300p.
}
\author{
Denismar Alves Nogueira

Eric Batista Ferreira

Portya Piscitelli Cavalcanti

(Adapted from Laercio Junio da Silva - GPL(>=2))
}
%\note{
%%  ~~further notes~~

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{snk}}, \code{\link{duncan}}, \code{\link{lsd}}, \code{\link{lsdb}}, \code{\link{ccboot}}, \code{\link{tukey}}, \code{\link{ccf}} 
}

\examples{
data(ex1)
attach(ex1)
crd(trat, ig, quali = TRUE, mcomp='sk', sigT = 0.05)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
