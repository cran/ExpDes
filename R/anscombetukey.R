#' Test for homogeneity of variances of Anscombe and Tukey
#'
#' \code{anscombetukey} Performs the test for homogeneity of
#'  variances of Anscombe and Tukey (1963).
#' @param resp Numeric or complex vector containing the response
#' variable.
#' @param Trat Numeric or complex vector containing the
#' treatments.
#' @param Bloco Numeric or complex vector containing the blocks.
#' @param glres Residual degrees of freedom.
#' @param msres Residual Mean Square.
#' @param sstrat Residual Sum of Squares.
#' @param ssbloco Sum of Squares for blocks.
#' @param residuals Numeric or complex vector containing the
#' residuals.
#' @param fitted.values Numeric or complex vector containing the
#' fitted values.
#' @return Returns the p-value of Anscombe and Tukey's test of
#' homogeneity of variances and its practical interpretation for
#' 5\% of significance.
#' @author Eric B Ferreira, \email{eric.ferreira@@unifal-mg.edu.br}
#' @author Denismar Alves Nogueira
#' @author Marcos Costa de Paula
#' @author Mateus Pimenta Siqueira Lima
#' @seealso \code{\link{han}}, \code{\link{oneillmathews}}.
#' @references ANSCOMBE, F. J.; TUKEY, J. W. \emph{The
#' examination and analysis of residuals.} Technometrics,
#' 5:141-160, 1963.
#'
#' RIBEIRO, R. \emph{Proposta e comparacao do desempenho de
#' testes para homogeneidade de variancia de modelos de
#' classificacao one-way e two-way}. Iniciacao Cientifica.
#' (Iniciacao Cientifica) - Universidade Federal de Alfenas.
#' 2012.
#' @examples
#' data(ex2)
#' attach(ex2)
#' rbd(trat, provador, aparencia, quali = TRUE, mcomp = "tukey",
#' hvar='anscombetukey', sigT = 0.05, sigF = 0.05)
#' @export

anscombetukey<-function(resp, Trat, Bloco, glres, msres, sstrat,
                        ssbloco, residuals, fitted.values)
{
Trat<-length(Trat)
Bloco<-length(Bloco)

div1<-(2*glres*(msres)^2)/glres+2
div2<-((((Trat-2)*(Bloco-1))/Trat*Bloco)*sstrat)+((((Trat-1)*(Bloco-2))/Trat*Bloco)*ssbloco)
Fc17<-((sum((residuals^2)*(fitted.values-mean(resp))))^2)/(div1*div2)
pvalor.hvar<-1-pf(Fc17,1,((Trat-1)*(Bloco-1)))
output <- pvalor.hvar
return(output)
}
