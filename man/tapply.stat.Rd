\name{tapply.stat}
\Rdversion{1.1}
\alias{tapply.stat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistics of data grouped by factors
}
\description{
This process lies in finding statistics which consist of more than one variable, grouped or crossed by factors. The table must be organized by columns between variables and factors. 
}
\usage{
tapply.stat(y, x, stat = "mean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{data.frame variables}
  \item{x}{data.frame factors}
  \item{stat}{Method}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~

\value{
y Numeric  

x Numeric  

stat method = "mean", ...  
}
%\references{
%% ~put references to the literature/web site here ~

\author{
Denismar Alves Nogueira

Eric Batista Ferreira

Portya Piscitelli Cavalcanti

(Adapted from Felipe de Mendiburu - GPL)
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
