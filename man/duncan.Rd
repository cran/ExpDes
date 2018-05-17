\name{duncan}
\Rdversion{1.1}
\alias{duncan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiple comparison: Duncan test}
\description{Performs the test of Duncan for multiple comparison of means.}
\usage{
duncan(y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE, main = NULL)
}
%- maybe also 'usage' for other objects documented here.
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
Returns the multiple comparison of means according to the test of Duncan.
}
%\references{
%% ~put references to the literature/web site here ~

\author{
Denismar Alves Nogueira

Eric Batista Ferreira

Portya Piscitelli Cavalcanti
}
%\note{
%%  ~~further notes~~

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{snk}}, \code{\link{tukey}}, \code{\link{lsd}}, \code{\link{lsdb}}, \code{\link{ccboot}}, \code{\link{scottknott}}, \code{\link{ccf}} 
}

\examples{
data(ex1)
attach(ex1)
crd(trat, ig, quali = TRUE, mcomp='duncan', sigT = 0.05)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
