#' Test for homogeneity of variances of ONeill and Mathews (RBD)
#'
#' \code{oneilldbc} Performs the test for homogeneity of
#' variances of ONeill and Mathews (2002).
#' @param trat Numeric or complex vector containing treatments.
#' @param resp Numeric or complex vector containing the response
#' variable.
#' @param block Numeric or complex vector containing blocks.
#' @return Returns the p-value of ONeill and Mathews' test of
#' homogeneity of variances and its practical interpretation
#' for significance level of 5\%.
#' @references O'NEILL, M. E.; MATHEWS, K. L. Levene tests of
#' homogeneity of variance for general block and treatment
#' designs. \emph{Biometrics}, 58:216-224, Mar. 2002.
#'
#' RIBEIRO, R. \emph{Proposta e comparacao do desempenho de
#' testes para homogeneidade de variancia de modelos de
#' classificacao one-way e two-way}. Iniciacao Cientifica.
#' (Iniciacao Cientifica) - Universidade Federal de Alfenas.
#' 2012.
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#'  @author Denismar Alves Nogueira
#'  @author Marcos Costa de Paula
#'  @author Mateus Pimenta Siqueira Lima
#' @seealso \code{\link{anscombetukey}}, \code{\link{han}}.
#' @examples
#' data(ex2)
#' attach(ex2)
#' rbd(trat, provador, aparencia, hvar = "oneillmathews")
#' @export

oneilldbc <-
function(resp, trat, block){
  #
  ntrat<-length(levels(factor(trat)))
  nbloc<-length(levels(factor(block)))
  data<-data.frame(trat,block,resp)
  data<-data[order(trat),]
  zdados<-y<-matrix(0,ntrat,nbloc)
    for(i in 1:ntrat) {
      y[i,]<-data$resp[((i-1)*nbloc + 1) : (i*nbloc)]
      }
  trat.mean<-apply(y,1,mean)
  bloc.mean<-apply(y,2,mean)
  g.mean<-mean(resp)
    for(i in 1:ntrat){
      for(j in 1:nbloc){
        zdados[i,j]=abs(y[i,j]- trat.mean[i] - bloc.mean[j] + g.mean)
      }
    }
  zdados<-as.vector(zdados)
  dadosz<-data.frame('z'=zdados,'blocagem'=rep(1:nbloc, each=ntrat),'tratamento'=rep(seq(1:ntrat),nbloc))
  Fc6<-summary(aov(dadosz$z ~ factor(dadosz$tratamento) + factor(dadosz$blocagem)))[[1]][1,4] # pvalor = posicao [1,5]
  rho<-c(-1/(ntrat-1), -1/(nbloc-1), 1/((nbloc-1)*(ntrat-1)))
  w0<-1-(2/pi)
  w1<-(2/pi)*(sqrt(1-rho[1]^2)+rho[1]*asin(rho[1])-1)
  w2<-(2/pi)*(sqrt(1-rho[2]^2)+rho[2]*asin(rho[2])-1)
  w3<-(2/pi)*(sqrt(1-rho[3]^2)+rho[3]*asin(rho[3])-1)
  m<-(w0-w1-w2+w3)/(w0-w1+(nbloc-1)*(w2-w3))
  Fc18<-m*Fc6
  pvalor.hvar<-(1-pf(Fc18, (ntrat-1), (nbloc-1)*(ntrat-1)))
  output <- pvalor.hvar
  return(output)
}
