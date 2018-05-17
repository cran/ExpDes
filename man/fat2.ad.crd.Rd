\name{fat2.ad.crd}
\Rdversion{1.1}
\alias{fat2.ad.crd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Double factorial scheme plus one additional treatment in CRD}
\description{Analyses experiments in balanced Completely Randomized Design in double factorial scheme with an additional treatment, considering a fixed model.}
\usage{
fat2.ad.crd(factor1, factor2, repet, resp, respAd, quali = c(TRUE, TRUE), 
mcomp = "tukey", fac.names = c("F1", "F2"), sigT = 0.05, sigF = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{factor1}{Numeric or complex vector containing the factor 1 levels.}
  \item{factor2}{Numeric or complex vector containing the factor 2 levels.}
  \item{repet}{Numeric or complex vector containing the replications.}
  \item{resp}{Numeric or complex vector containing the response variable.}
  \item{respAd}{Numeric or complex vector containing the additional treatment.}
  \item{quali}{Logic. If TRUE (default), the treatments are assumed qualitative, if FALSE, quantitatives.}
  \item{mcomp}{Allows choosing the multiple comparison test; the \emph{default} is the test of Tukey, however, the options are: the LSD test ('lsd'), the LSD test with Bonferroni protection ('lsdb'), the test of Duncan ('duncan'), the test of Student-Newman-Keuls ('snk'), the test of Scott-Knott ('sk'), the Calinski and Corsten test ('ccf') and bootstrap multiple comparison's test ('ccboot').}
  \item{fac.names}{Allows labeling the factors 1 and 2.}
  \item{sigT}{The signficance to be used for the multiple comparison test; the default is 5\%.}
  \item{sigF}{The signficance to be used for the F test of ANOVA; the default is 5\%.}
}

\details{
The arguments sigT and mcomp will be used only when the treatment are qualitative.
}
\value{
The output contains the ANOVA of the referred CRD, the Shapiro-Wilk normality test for the residuals of the model, the fitted regression models
(when the treatments are quantitative) and/or the multiple comparison tests (when the treatments are qualitative).
}
\references{
HEALY, M. J. R. The analysis of a factorial experiment with additional treatments. Journal of Agricultural Science, Cambridge, v. 47, p. 205-206. 1956.

FERREIRA, E. B.; CAVALCANTI, P. P.; NOGUEIRA D. A. Funcao para analisar experimentos em fatorial duplo com um tratamento adicional, em uma so rodada. 
In: CONGRESSO DE POS-GRADUACAO DA UNIVERSIDADE FEDERAL DE LAVRAS, 19., 2010, Lavras. Resumos... Lavras: UFLA, 2010.
}
\author{
Denismar Alves Nogueira

Eric Batista Ferreira
   
Portya Piscitelli Cavalcanti
}

\note{
The \code{\link{graphics}} can be used to construct regression plots and \code{\link{plotres}} for residuals plots
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For more examples, see: \code{\link{fat2.crd}}, \code{\link{fat2.rbd}}, \code{\link{fat3.crd}}, \code{\link{fat3.rbd}}, \code{\link{fat2.ad.crd}}, \code{\link{fat2.ad.rbd}}, \code{\link{fat3.ad.crd}} and \code{\link{fat3.ad.rbd}}.
}

\examples{
data(ex8)
attach(ex8)
data(secaAd)
fat2.ad.crd(inoculante, biodiesel, vaso, seca, secaAd, quali = c(TRUE,FALSE), mcomp = "tukey", 
 fac.names = c("Inoculant", "Biodiesel"), sigT = 0.05, sigF = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
