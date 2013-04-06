\name{ccboot}
\alias{ccboot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiple comparison: Bootstrap
}
\description{
Performs the Ramos and Ferreira (2009) multiple comparison bootstrap test.
}
\usage{
ccboot(y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE, main = NULL, B = 1000)
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
  \item{B}{Number of bootstrap resamples.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~

\value{
Multiple means comparison for the bootstrap test.
}
\references{
RAMOS, P. S., FERREIRA, D. F. Agrupamento de medias via bootstrap de populacoes normais e nao-normais, Revista Ceres, v.56, p.140-149, 2009.
}
\author{
Patricia de Siqueira Ramos

Daniel Furtado Ferreira

Eric Batista Ferreira
}
%\note{
%%  ~~further notes~~

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~

\examples{
data(ex1)
attach(ex1)
crd(trat, ig, quali = TRUE, mcomp='ccboot', sigF = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
