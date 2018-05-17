\name{oneilldbc}
\alias{oneilldbc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test for homogeneity of variances of ONeill and Mathews (RBD)}
\description{Performs the test for homogeneity of variances of ONeill and Mathews (2002).}
\usage{
oneilldbc(resp, trat, block)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resp}{
   Numeric or complex vector containing the response variable.
  }
  \item{trat}{
   Numeric or complex vector containing treatments.
  }
  \item{block}{
   Numeric or complex vector containing blocks.
  }
}

%\details{ Return p-value }
%%  ~~ If necessary, more details than the description above ~~
  
\value{
Returns the p-value of ONeill and Mathews' test of homogeneity of variances and its practical interpretation for significance level of 5\%.
}

\references{
O´NEILL, M. E.; MATHEWS, K. L. Levene tests of homogeneity of variance for general block and treatment designs. \emph{Biometrics}, 58:216-224, Mar. 2002.

RIBEIRO, R. \emph{Proposta e comparacao do desempenho de testes para homogeneidade de variancia de modelos de classi cação one-way e two-way}. Iniciacao Cientifica. (Iniciacao Cientifica) - Universidade Federal de Alfenas. 2012.
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
    see: \code{\link{anscombetukey}}, \code{\link{han}}
  }

\examples{
data(ex2)
attach(ex2)
rbd(trat, provador, aparencia, hvar='oneillmathews')
}