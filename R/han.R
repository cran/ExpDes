#' Test for homogeneity of variances of Han
#'
#' \code{han} Performs the test for homogeneity of variances of
#' Han (1969).
#' @param resp Numeric or complex vector containing the response
#' variable.
#' @param trat Numeric or complex vector containing the
#' treatments.
#' @param block Numeric or complex vector containing the blocks.
#' @return Returns the p-value of Han's test of homogeneity of
#' variances and its practical interpretation for 5\% of
#' significance.
#' @references HAN, C. P. Testing the homogeneity of variances
#' in a two-way classification. \emph{Biometrics}, 25:153-158,
#' Mar. 1969.
#'
#' RIBEIRO, R. \emph{Proposta e comparacao do desempenho de
#' testes para homogeneidade de variancia de modelos de
#' classicacao one-way e two-way}. Iniciacao Cientifica.
#' (Iniciacao Cientifica) - Universidade Federal de Alfenas.
#' 2012.
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#'  @author Denismar Alves Nogueira
#'  @author Marcos Costa de Paula
#'  @author Mateus Pimenta Siqueira Lima
#' @seealso \code{\link{anscombetukey}},
#' \code{\link{oneillmathews}}.
#' @examples
#' data(ex2)
#' attach(ex2)
#' rbd(trat, provador, aparencia, hvar = "han")
#' @export

han <-
function(resp, trat, block)
{
  Trat<-length(levels(trat)) # numero de tratamentos (Trat deve ser factor)
  Block<-length(levels(block)) # numero de blocos    (Bloco deve ser factor)
  if (Block>Trat){   # se numero de blocos for menor ou igual ao numero de tratamentos n?o rodar.
   dife<-matrix(0,Block,Trat)
   ymedia<-matrix(0,Block,1)
   rp<-0
   for(j in 1:Block) {
    for(i in 1:Trat) {
      ymedia[j]<-mean(resp[(rp+1):(rp+Trat)])
      dife[j,i]<-(resp[rp+i]-ymedia[j])
      }
    rp<-Trat*j
    }
   modelohan<-lm(ymedia ~ dife[,2:Trat])
   pvalor.hvar<-1-pf(summary(modelohan)[[10]][1],summary(modelohan)[[10]][2],summary(modelohan)[[10]][3])
   output <- pvalor.hvar
   return(output)
   }
}
