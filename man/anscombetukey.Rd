\name{anscombetukey}
\alias{anscombetukey}

\title{Test for homogeneity of variances of Anscombe and Tukey}

\description{Performs the test for homogeneity of variances of Anscombe and Tukey (1963).}

\usage{
anscombetukey(resp, trat, block, glres, msres, sstrat, ssblock, residuals, fitted.values)
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

\item{glres}{
Residual degrees of freedom.
}

\item{msres}{
Residual Mean Square.
}

\item{sstrat}{
Residual Sum of Squares.
}

\item{ssblock}{
Sum of Squares for blocks.
}

\item{residuals}{
Numeric or complex vector containing the residuals.
}

\item{fitted.values}{
Numeric or complex vector containing the fitted values.
}
}

%\details{ }
  
\value{
Returns the p-value of Anscombe and Tukey's test of homogeneity of variances and its practical interpretation for 5\% of significance.
}

\references{
ANSCOMBE, F. J.; TUKEY, J. W. \emph{The examination and analysis of residuals.} Technometrics, 5:141-160, 1963.

RIBEIRO, R. \emph{Proposta e comparação do desempenho de testes para homogeneidade de variância de modelos de classi cação one-way e two-way}. Iniciação Científica. (Iniciação Científica) - Universidade Federal de Alfenas. 2012.
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
See also: \code{\link{han}}, \code{\link{oneillmathews}}
}

\examples{
data(ex2)
attach(ex2)
rbd(trat, provador, aparencia, hvar='anscombetukey')
}
