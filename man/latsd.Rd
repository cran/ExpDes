\name{latsd}
\Rdversion{1.1}
\alias{latsd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Latin Square Design}                        
\description{Analyses experiments in balanced Latin Square Design, considering a fixed model.}
\usage{
latsd(treat, row, column, resp, quali = TRUE, mcomp = "tukey", sigT = 0.05, sigF = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{treat}{Numeric or complex vector containing the treatments.}
  \item{row}{Numeric or complex vector containing the rows.}
  \item{column}{Numeric or complex vector containing the columns.}
  \item{resp}{Numeric or complex vector containing the response variable.}
  \item{quali}{Logic. If TRUE (default), the treatments are assumed qualitative, if FALSE, quantitatives.}
  \item{mcomp}{Allows choosing the multiple comparison test; the \emph{default} is the test of Tukey, however, the options are: the LSD test ('lsd'), the LSD test with Bonferroni protection ('lsdb'), the test of Duncan ('duncan'), the test of Student-Newman-Keuls ('snk'), the test of Scott-Knott ('sk'), the Calinski and Corsten test ('ccf') and bootstrap multiple comparison's test ('ccboot').}
  \item{sigT}{The signficance to be used for the multiple comparison test; the default is 5\%.}
  \item{sigF}{The signficance to be used for the F test of ANOVA; the default is 5\%.}
}
\details{
The arguments sigT and mcomp will be used only when the treatment are qualitative.
}
\value{
The output contains the ANOVA of the LSD, the Shapiro-Wilk normality test for the residuals of the model, the fitted regression models
(when the treatments are quantitative) and/or the multiple comparison tests (when the treatments are qualitative).
}

\references{
GOMES, F. P. Curso de Estatistica Experimental. 10a ed. Piracicaba: ESALQ/USP. 1982. 430.

FERREIRA, E. B.; CAVALCANTI, P. P.; NOGUEIRA D. A. Funcao em codigo R para analisar experimentos em DQL simples, em uma so rodada. 
In: CONGRESSO DE POS-GRADUACAO DA UNIVERSIDADE FEDERAL DE LAVRAS, 18., 2009, Lavras. Annals... Lavras: UFLA, 2009.
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
For more examples, see: \code{\link{crd}}, \code{\link{rbd}}
}

\examples{
data(ex3)
attach(ex3)
latsd(trat, linha, coluna, resp, quali=TRUE, mcomp="snk", sigT=0.05, sigF=0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
