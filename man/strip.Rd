\name{strip}
\alias{strip}

\title{Strip-plot experiments}
\description{Analysis Strip-plot experiments}

\usage{
strip(factor1, factor2, block, resp, quali = c(TRUE, TRUE), mcomp = "tukey", 
fac.names = c("F1", "F2"), sigT = 0.05, sigF = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{factor1}{Numeric or complex vector containing the factor 1 levels.}
  \item{factor2}{Numeric or complex vector containing the factor 2 levels.}
  \item{block}{Numeric or complex vector containing the blocks.}
  \item{resp}{Numeric or complex vector containing the response variable.}
  \item{quali}{Logic. If TRUE (default), the treatments are assumed qualitative, if FALSE, quantitatives.}
  \item{mcomp}{Allows choosing the multiple comparison test; the \emph{default} is the test of Tukey, however, the options are: the LSD test ('lsd'),
  the LSD test with Bonferroni protection ('lsdb'), the test of Duncan ('duncan'), the test of Student-Newman-Keuls ('snk'), the test of Scott-Knott
  ('sk'), the Calinski and Corsten test ('ccf') and bootstrap multiple comparison's test ('ccboot').}
  \item{fac.names}{Allows labeling the factors 1 and 2.}
  \item{sigT}{The signficance to be used for the multiple comparison test; the default is 5\%.}
  \item{sigF}{The signficance to be used for the F test of ANOVA; the default is 5\%.}
}

%\details{}

%\value{}

%\references{}

\author{
Eric B Ferreira

Lais B S Ferreira

}

\note{
The \code{\link{graphics}} can be used to construct regression plots and \code{\link{plotres}} for residuals plots
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rbd}} and \code{\link{split2.rbd}}
}

\examples{
data(ex5)
attach(ex5)
strip(trat, genero, bloco, sabor, quali=c(TRUE,TRUE), mcomp="tukey", 
 fac.names=c("Amostras","Genero"), sigT = 0.05, sigF = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
