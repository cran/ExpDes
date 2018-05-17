\name{order.group}
\Rdversion{1.1}
\alias{order.group}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Ordering the treatments according to the multiple comparison}
\description{Orders the groups os means.}
\usage{
order.group(trt, means, N, MSerror, Tprob, std.err, parameter = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trt}{Treatments}
  \item{means}{Means of treatment}
  \item{N}{Replications}
  \item{MSerror}{Mean square error}
  \item{Tprob}{minimum value for the comparison}
  \item{std.err}{standard error}
  \item{parameter}{Constante 1 (Sd), 0.5 (Sx)}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~

\value{
trt Factor  
means Numeric  
N Numeric  
MSerror Numeric  
Tprob value between 0 and 1  
std.err Numeric  
parameter Constant  
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
order.stat 
}
%\examples{

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
