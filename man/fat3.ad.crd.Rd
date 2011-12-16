\name{fat3.ad.crd}
\Rdversion{1.1}
\alias{fat3.ad.crd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Triple factorial scheme plus an additional treatment in CRD
}
\description{
Analyses experiments in balanced Completely Randomized Design in triple factorial scheme with an additional treatment, considering a fixed model.
}
\usage{
fat3.ad.crd(factor1, factor2, factor3, repet, resp, respAd, quali = c(TRUE, TRUE, TRUE), mcomp = "tukey", fac.names = c("F1", "F2", "F3"), sigT = 0.05, sigF = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{factor1}{Numeric or complex vector containing the factor 1 levels.}
  \item{factor2}{Numeric or complex vector containing the factor 2 levels.}
  \item{factor3}{Numeric or complex vector containing the factor 3 levels.}
  \item{repet}{Numeric or complex vector containing the replications.}
  \item{resp}{Numeric or complex vector containing the response variable.}
  \item{respAd}{Numeric or complex vector containing the additional treatment.}
  \item{quali}{Logic. If TRUE (default), the treatments are assumed qualitative, if FALSE, quantitatives.}
  \item{mcomp}{Allows choosing the multiple comparison test; the \emph{default} is the test of Tukey, however, the options are: the LSD test ('lsd'),
  the LSD test with Bonferroni protection ('lsdb'), the test of Duncan ('duncan'), the test of Student-Newman-Keuls ('snk'), the test of Scott-Knott
  ('sk') and bootstrap multiple comparison's test ('ccboot').}
  \item{fac.names}{Allows labeling the factors 1, 2 and 3.}
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
}
\author{
Denismar Alves Nogueira

Eric Batista Ferreira
   
Portya Piscitelli Cavalcanti
}
%\note{
%%  ~~further notes~~

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For more examples, see: fat2.crd, fat2.rbd, fat3.crd, fat3.rbd, fat2.ad.crd, fat2.ad.rbd and fat3.ad.rbd.
}
\examples{
data(ex6)
attach(ex6)
data(respAd)
fat3.ad.crd(fatorA, fatorB, fatorC, rep, resp, respAd, quali = c(TRUE, TRUE, TRUE), mcomp = "duncan", fac.names = c("Factor A", "Factor B", "Factor C"), sigT = 0.05, sigF = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
