\name{reg.poly}
\Rdversion{1.1}
\alias{reg.poly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Polinomial Regression
}
\description{
Fits sequencial regression models until the third power.
}
\usage{
reg.poly(resp, treat, DFerror, SSerror, DFtreat, SStreat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resp}{Numeric or complex vector containing the response variable.}
  \item{treat}{Numeric or complex vector containing the treatments.}
  \item{DFerror}{Error degrees of freedom}
  \item{SSerror}{Error sum of squares}
  \item{DFtreat}{Treatments' dregrees of freedom}
  \item{SStreat}{Treatments' sum of squares}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~

\value{
Returns coefficients, significance and ANOVA of the fitted regression models.
}
\references{
GOMES, F. P. Curso de Estatistica Experimental. 10a ed. Piracicaba: ESALQ/USP. 1982. 430.
}
\author{
Denismar Alves Nogueira

Eric Batista Ferreira

Portya Piscitelli Cavalcanti
}
%\note{
%%  ~~further notes~~

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~

%\examples{

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
