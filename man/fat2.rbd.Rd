\name{fat2.rbd}
\Rdversion{1.1}
\alias{fat2.rbd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Double factorial scheme in RBD
}
\description{
Analyses experiments in balanced Randomized Blocks Designs in double factorial scheme, considering a fixed model.
}
\usage{
fat2.rbd(factor1, factor2, block, resp, quali = c(TRUE, TRUE), mcomp = "tukey", fac.names = c("F1", "F2"), sigT = 0.05, sigF = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{factor1}{Numeric or complex vector containing the factor 1 levels.}
  \item{factor2}{Numeric or complex vector containing the factor 2 levels.}
  \item{block}{Numeric or complex vector containing the blocks.}
  \item{resp}{Numeric or complex vector containing the response variable.}
  \item{quali}{Logic. If TRUE (default), the treatments are assumed qualitative, if FALSE, quantitatives.}
  \item{mcomp}{Allows choosing the multiple comparison test; the \emph{default} is the test of Tukey, however, the options are: the LSD test ('lsd'),
  the LSD test with Bonferroni protection ('lsdb'), the test of Duncan ('duncan'), the test of Student-Newman-Keuls ('snk'), the test of Scott-Knott]
  ('sk') and bootstrap multiple comparison's test ('ccboot').}
  \item{fac.names}{Allows labeling the factors 1 and 2.}
  \item{sigT}{The signficance to be used for the multiple comparison test; the default is 5\%.}
  \item{sigF}{The signficance to be used for the F test of ANOVA; the default is 5\%.}
}
\details{
The arguments sigT and mcomp will be used only when the treatment are qualitative.
}
\value{
The output contains the ANOVA of the referred RBD, the Shapiro-Wilk normality test for the residuals of the model, the fitted regression models
(when the treatments are quantitative) and/or the multiple comparison tests (when the treatments are qualitative).
}
\references{
BANZATTO, D. A.; KRONKA, S. N. Experimentacao Agricola. 4 ed. Jaboticabal: Funep. 2006. 237 p.
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
For more examples, see: fat2.crd, fat3.crd, fat3.rbd, fat2.ad.crd, fat2.ad.rbd, fat3.ad.crd and fat3.ad.rbd.
}
\examples{
data(ex5)
attach(ex5)
fat2.rbd(trat, genero, bloco, sabor ,quali=c(TRUE,TRUE), mcomp="lsd", fac.names=c("Samples","Gender"), sigT = 0.05, sigF = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
