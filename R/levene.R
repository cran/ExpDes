#' Test for homogeneity of variances of Levene
#'
#' \code{levene} Performs the test for homogeneity of variances
#' of Levene (1960).
#' @param trat Numeric or complex vector containing treatments.
#' @param resp Numeric or complex vector containing the response
#' variable.
#' @param t Number of treatments.
#' @param r Numeric or complex vector containing the number of
#' replications of each treatment.
#' @return Returns the p-value of Levene's test of homogeneity
#' of variances and its practical interpretation for
#' significance level of 5\%.
#' @references LEVENE, H. Robust tests for equality of
#' variances. In: Olkin, I.; Ghurye, S.G.; Hoeffding, W.;
#' Madow, W.G.; Mann, H.B. (eds.). \emph{Contribution to
#' Probability and Statistics. Stanford}, CA: Stanford
#' University Press, pages 278-292, 1960.
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
#' @seealso \code{\link{bartlett}}, \code{\link{samiuddin}},
#' \code{\link{layard}}, \code{\link{oneillmathews}}.
#' @examples
#' data(ex1)
#' attach(ex1)
#' crd(trat, ig, quali = FALSE, hvar = "levene")
#' @export

levene <-
function(trat, resp, t, r)
{
  Trat<-factor(trat)
   zdados1<-matrix(0,length(resp),1)
   rp<-0
   for(k in 1:length(resp)) {
    zdados1[k]<-abs(resp[k]-mean(resp[(rp+1):(rp+r[Trat[k]])]))
    if(k<length(resp)){if(trat[k]!=trat[k+1]){rp<-sum(r[1:Trat[k]])}}
     }
   pvalor<-summary(aov(zdados1 ~ Trat))[[1]][1,5]
   output <- pvalor
   return(output)
}
