\name{crd}
\alias{crd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{One factor Completely Randomized Design}
\description{Analyses balanced experiments in Completely Randomized Design under one single factor, considering a fixed model.}
\usage{
crd(treat, resp, quali = TRUE, mcomp = "tukey", nl = FALSE, hvar='bartlett', 
sigT = 0.05, sigF = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{treat}{Numeric or complex vector containing the treatments.}

\item{resp}{Numeric or complex vector containing the response variable.}

\item{quali}{Logic. If TRUE (default), the treatments are assumed qualitative, if FALSE, quantitatives.}

\item{mcomp}{Allows choosing the multiple comparison test; the \emph{default} is the test of Tukey, however, the options are: the LSD test ('lsd'), the LSD test with Bonferroni protection ('lsdb'), the test of Duncan ('duncan'), the test of Student-Newman-Keuls ('snk'), the test of Scott-Knot ('sk'), the Calinski and Corsten test ('ccf') and bootstrap multiple comparison's test ('ccboot').}

\item{nl}{Logic. If FALSE (\emph{default}) linear regression models are adjusted. IF TRUE, non-linear regression models are adjusted.}

\item{hvar}{Allows choosing the test for homogeneity of variances; the \emph{default} is the test of Bartlett, however there are other options: test of Levene ('levene'), test of Samiuddin ('samiuddin'), test of ONeill and Mathews ('oneillmathews') and the Layard test ('layard').}

\item{sigT}{The signficance to be used for the multiple comparison test; the default is 5\%.}

\item{sigF}{The signficance to be used for the F test of ANOVA; the default is 5\%.}
}

\details{
The arguments sigT and mcomp will be used only when the treatment are qualitative.
}
\value{
The output contains the ANOVA of the CRD, the Shapiro-Wilk normality test for the residuals of the model, the fitted regression models
(when the treatments are quantitative) and/or the multiple comparison tests (when the treatments are qualitative).
}
\references{
BANZATTO, D. A.; KRONKA, S. N. Experimentacao Agricola. 4 ed. Jaboticabal: Funep. 2006. 237 p.

FERREIRA, E. B.; CAVALCANTI, P. P. Funcao em codigo R para analisar experimentos em DIC simples, em uma so rodada. In: REUNIAO ANUAL
DA REGIAO BRASILEIRA DA SOCIEDADE INTERNACIONAL DE BIOMETRIA, 54./SIMPOSIO DE ESTATISTICA APLICADA A EXPERIMENTACAO AGRONOMICA, 13.,
2009, Sao Carlos. Programas e resumos... Sao Carlos, SP: UFSCar, 2009. p. 1-5.
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
For more examples, see: \code{\link{fat2.crd}}, \code{\link{fat3.crd}}, \code{\link{split2.crd}}, \code{\link{fat2.ad.crd}} and \code{\link{fat3.ad.crd}}.
}

\examples{
data(ex1)
attach(ex1)
crd(trat, ig, quali = FALSE, sigF = 0.05)
}
%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
