\name{lastC}
\Rdversion{1.1}
\alias{lastC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Setting the last character of a chain}
\description{A special function for the group of treatments in the multiple comparison tests. Use order.group.}
\usage{
lastC(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{letters}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~

\value{
x  character 
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

\seealso{
order.group 
}
\examples{
x<-c("a","ab","b","c","cd")
lastC(x)
# "a" "b" "b" "c" "d"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
