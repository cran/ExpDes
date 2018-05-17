\name{tukey}
\Rdversion{1.1}
\alias{tukey}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiple comparison: Tukey's test}
\description{Performs the test of Tukey, for multiple comparison of means.}
\usage{
tukey(y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE, main = NULL)
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
\details{
It is necessary first makes a analysis of variance.
}
\value{
y  Numeric 
trt  factor 
DFerror  Numeric 
MSerror  Numeric 
alpha  Numeric 
group Logic 
main  Text 
}
\references{
Principles and procedures of statistics a biometrical approach Steel and Torry and Dickey. Third Edition 1997 
}
\author{
Denismar Alves Nogueira

Eric Batista Ferreira

Portya Piscitelli Cavalcanti

(Adapted from Felipe de Mendiburu - GPL)
}
%\note{
%%  ~~further notes~~

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{snk}}, \code{\link{duncan}}, \code{\link{lsd}}, \code{\link{lsdb}}, \code{\link{ccboot}}, \code{\link{scottknott}}, \code{\link{ccf}} 
}

\examples{
data(ex1)
attach(ex1)
crd(trat, ig, quali = TRUE, mcomp='tukey', sigT = 0.05)
}
%\examples{

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
