#' Test for homogeneity of variances of ONeill and Mathews (CRD)
#'
#' \code{oneillmathews} Performs the test for homogeneity of
#' variances of ONeill and Mathews (2000).
#' @param trat Numeric or complex vector containing treatments.
#' @param resp Numeric or complex vector containing the response
#' variable.
#' @param t Number of treatments.
#' @param r Numeric or complex vector containing the number of
#' replications of each treatment.
#' @return Returns the p-value of ONeill and Mathews' test of
#' homogeneity of variances and its practical interpretation
#' for significance level of 5\%.
#' @references O'NEILL, M. E.; MATHEWS, K. L. A weighted least
#' squares approach to levene test of homogeneity of variance.
#' \emph{Australian e New Zealand Journal Statistical},
#' 42(1):81-100, 2000.
#'
#' NOGUEIRA, D, P.; PEREIRA, G, M. Desempenho de testes para
#' homogeneidade de variancias em delineamentos inteiramente
#' casualizados. \emph{Sigmae}, Alfenas, v.2, n.1, p. 7-22.
#' 2013.
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#'  @author Denismar Alves Nogueira
#'  @author Marcos Costa de Paula
#'  @author Mateus Pimenta Siqueira Lima
#' @seealso \code{\link{bartlett}}, \code{\link{layard}},
#' \code{\link{levene}}, \code{\link{samiuddin}}.
#' @examples
#' data(ex1)
#' attach(ex1)
#' crd(trat, ig, quali = FALSE, hvar = "oneillmathews",
#' sigF = 0.05)
#' @export

oneillmathews <-
function(trat, resp, t, r)
{
  Trat<-factor(trat)
  zdados1.1<-matrix(0,length(resp),1)
  rr<-t/sum(1/r)
  rp<-0
   for(k in 1:length(resp)) {
     zdados1.1[k]<-abs(resp[k]-mean(resp[(rp+1):(rp+r[Trat[k]])]))/sqrt(1-(1/rr))
     if(k<length(resp)){if(trat[k]<trat[k+1]){rp<-sum(r[1:Trat[k]])}}
    }
  Fc5.1<-summary(aov(zdados1.1 ~ trat))[[1]][1,4] # pvalor = posicao [1,5]
  b<-(1-2/pi)
  c<-(2/pi)*(1/(rr-1))*(sqrt(rr*(rr-2))+asin(1/(rr-1))-(rr-1))
  m<-(b-c)/(b+(rr-1)*c)
  Fc13<-m*Fc5.1
  pvalor<-(1-pf(Fc13, (t-1), summary(aov(zdados1.1 ~ trat))[[1]][2,1]))
  output <- pvalor
  return(output)
}
