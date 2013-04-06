\name{rbd}
\Rdversion{1.1}
\alias{rbd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Randomized Blocks Design
}
\description{
Analyses experiments in balanced Randomized Blocks Designs under one single factor, considering a fixed model.
}
\usage{
rbd(treat, block, resp, quali = TRUE, mcomp = "tukey", sigT = 0.05, sigF = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{treat}{Numeric or complex vector containing the treatments.}
  \item{block}{Numeric or complex vector containing the blocks.}
  \item{resp}{Numeric or complex vector containing the response variable.}
  \item{quali}{Logic. If TRUE (default), the treatments are assumed qualitative, if FALSE, quantitatives.}
  \item{mcomp}{Allows choosing the multiple comparison test; the \emph{default} is the test of Tukey, however, the options are: the LSD test ('lsd'),
  the LSD test with Bonferroni protection ('lsdb'), the test of Duncan ('duncan'), the test of Student-Newman-Keuls ('snk'), the test of Scott-Knott
  ('sk') and bootstrap multiple comparison's test ('ccboot').}
  \item{sigT}{The signficance to be used for the multiple comparison test; the default is 5\%.}
  \item{sigF}{The signficance to be used for the F test of ANOVA; the default is 5\%.}
}
\details{
The arguments sigT and mcomp will be used only when the treatment are qualitative.
}
\value{                                                                                                
The output contains the ANOVA of the RBD, the Shapiro-Wilk normality test for the residuals of the model, the fitted regression models
(when the treatments are quantitative) and/or the multiple comparison tests (when the treatments are qualitative).
}

\references{
BANZATTO, D. A.; KRONKA, S. N. Experimentacao Agricola. 4 ed. Jaboticabal: Funep. 2006. 237 p.

FERREIRA, E. B.; CAVALCANTI, P. P.; NOGUEIRA D. A. Funcao em codigo R para analisar experimentos em DBC simples, em uma so rodada. 
In: JORNADA CIENTIFICA DA UNIVERSIDADE FEDERAL DE ALFENAS-MG, 2., 2009, Alfenas. Annals... ALfenas: Unifal-MG, 2009.
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
For more examples, see: fat2.rbd, fat3.rbd, split2.rbd, split2.ad.rbd and fat3.ad.rbd.
}
\examples{
data(ex2)
attach(ex2)
rbd(trat, provador, aparencia, quali = TRUE, mcomp='lsd', sigT = 0.05, sigF = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
