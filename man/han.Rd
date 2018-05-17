\name{han}
\alias{han}

\title{Test for homogeneity of variances of Han}
\description{Performs the test for homogeneity of variances of Han (1969).}
\usage{
han(resp, trat, block)
}

\arguments{
\item{resp}{
Numeric or complex vector containing the response variable.
}

\item{trat}{
Numeric or complex vector containing the treatments.
}
  
\item{block}{
Numeric or complex vector containing the blocks.
}
}

%\details{ }
%%  ~~ If necessary, more details than the description above ~~
  
\value{
Returns the p-value of Han's test of homogeneity of variances and its practical interpretation for 5\% of significance.
}

\references{
HAN, C. P. Testing the homogeneity of variances in a two-way classification. \emph{Biometrics}, 25:153-158, Mar. 1969.

RIBEIRO, R. \emph{Proposta e comparacao do desempenho de testes para homogeneidade de variancia de modelos de classicacao one-way e two-way}. Iniciacao Cientifica. (Iniciacao Cientifica) - Universidade Federal de Alfenas. 2012.
}

  
\author{
Denismar Alves Nogueira
    
Eric Batista Ferreira
    
Marcos Costa de Paula
    
Mateus Pimenta Siqueira Lima
}

%\note{ }
%%  ~~further notes~~
  
  
%% ~Make other sections like Warning with \section{Warning }{....} ~
  
\seealso{
See also: \code{\link{anscombetukey}}, \code{\link{oneillmathews}}
}
    
\examples{
data(ex2)
attach(ex2)
rbd(trat, provador, aparencia, hvar='han')
    }