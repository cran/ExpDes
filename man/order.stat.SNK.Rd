\name{order.stat.SNK}
\Rdversion{1.1}
\alias{order.stat.SNK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Grouping the treatments averages in a comparison with a minimum value}
\description{Orders the groups of means according to the test of SNK.}
\usage{
order.stat.SNK(treatment, means, minimum)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{treatment}{treatment}
  \item{means}{means of treatment}
  \item{minimum}{minimum value for the comparison}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~

\value{
trt Factor  
means Numeric  
minimum Numeric 
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
%\examples{

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
