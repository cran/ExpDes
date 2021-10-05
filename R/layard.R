#' Test for homogeneity of variances of Layard
#'
#' \code{layard} Performs the test for homogeneity of variances
#' of Layard for Jackknife (1973).
#' @param trat Numeric or complex vector containing treatments.
#' @param resp Numeric or complex vector containing the response
#' variable.
#' @param t Number of treatments.
#' @param r Numeric or complex vector containing the number of
#' replications of each treatment.
#' @return Returns the p-value of the Layard test of homogeneity
#' of variances and its practical interpretation for the
#' significance level of 5\%.
#' @references LAYARD, M. N. J. Robust large-sample tests for
#' homogeneity of variances. \emph{Journal of the American
#' Statistical Association}, v.68, n.341, p.195-198, 1973.
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
#' \code{\link{levene}}, \code{\link{oneillmathews}}.
#' @examples
#' data(ex1)
#' attach(ex1)
#' crd(trat, ig, quali = FALSE, hvar = "layard")
#' @export

layard <-
function(trat, resp, t, r)
{
 vari<-matrix(0,t,1)
 varij<-matrix(0,max(r),t)
 U<-matrix(0,max(r),t)
 rp<-0
 for(i in 1:t) {
   vv<-resp[(rp+1):(rp+r[i])]
   vari[i]<-var(vv)
   for(j in 1:r[i]) {
     varij[j,i]<-var(vv[-j])
     U[j,i]<-(r[i]*log(vari[i]))-((r[i]-1)*log(varij[j,i]))
    }
   rp<-sum(r[1:i])
 }
 Uij<-as.vector(U)
 dadosUij<-cbind(trat,Uij)
 dadosUij<-as.data.frame(dadosUij)
 pvalor<-summary(aov(dadosUij$Uij ~ trat))[[1]][1,5]
 output <- pvalor
 return(output)
}
