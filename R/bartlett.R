#' Test for Homogeneity of Variances: Bartlett
#'
#' \code{bartlett} Performs the test for homogeneity of
#' variances of Bartlett (1937).
#' @param trat Numeric or complex vector containing the
#' treatments.
#' @param resp Numeric or complex vector containing the
#' response variable.
#' @param t Number of treatments.
#' @param r Numeric or complex vector containing the number of
#' replications of each treatment.
#' @return Returns the p-value of Bartlett's test of
#' homogeneity of variances and its practical interpretation
#' for 5\% of significance.
#' @references BARTLETT, M. S. Properties of sufficiency and
#' statistical tests. \emph{Proceedings of the Royal
#' Statistical Society - Serie A}, 60:268-282, 1937.
#'
#' NOGUEIRA, D, P.; PEREIRA, G, M. Desempenho de testes para
#' homogeneidade de vari?ncias em delineamentos inteiramente
#' casualizados. \emph{Sigmae}, Alfenas, v.2, n.1, p. 7-22.
#' 2013.
#' @author Eric B Ferreira, \email{eric.ferreira@@unifal-mg.edu.br}
#' @author Denismar Alves Nogueira
#' @author Marcos Costa de Paula
#' @author Mateus Pimenta Siqueira Lima
#' @seealso \code{\link{levene}},
#' \code{\link{oneillmathews}}, \code{\link{samiuddin}}
#' @examples
#' data(ex1)
#' attach(ex1)
#' crd(trat, ig, quali = FALSE, hvar='bartlett', sigF = 0.05)
#' @export

bartlett <-function(trat, resp, t, r)
{
  vari<-tapply(resp,trat,var)
  S2p<-sum((r-1)*vari)/(length(resp)-t)
  A<-(length(resp)-t)*log(S2p)-sum((r-1)*log(vari))
  B<-(1/(3*(t-1)))*(sum(1/(r-1))-(1/(length(resp)-t)))
  Xc1<-A/(1+B)
  pvalue<-1-pchisq(Xc1, t-1)
  output <- pvalue
  return(output)
}
