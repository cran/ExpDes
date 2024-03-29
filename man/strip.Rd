% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strip.R
\name{strip}
\alias{strip}
\title{Strip-plot experiments}
\usage{
strip(
  factor1,
  factor2,
  block,
  resp,
  quali = c(TRUE, TRUE),
  mcomp = "tukey",
  fac.names = c("F1", "F2"),
  sigT = 0.05,
  sigF = 0.05,
  unfold = NULL
)
}
\arguments{
\item{factor1}{Numeric or complex vector containing the
factor 1 levels.}

\item{factor2}{Numeric or complex vector containing the
factor 2 levels.}

\item{block}{Numeric or complex vector containing the blocks.}

\item{resp}{Numeric or complex vector containing the
response variable.}

\item{quali}{Logic. If TRUE (default), the treatments
are assumed qualitative, if FALSE, quantitatives.}

\item{mcomp}{Allows choosing the multiple comparison
test; the \emph{default} is the test of Tukey, however,
the options are: the LSD test ('lsd'), the LSD test
with Bonferroni protection ('lsdb'), the test of Duncan
('duncan'), the test of Student-Newman-Keuls ('snk'),
the test of Scott-Knott ('sk'), the Calinski and
Corsten test ('ccF') and bootstrap multiple comparison's
test ('ccboot').}

\item{fac.names}{Allows labeling the factors 1 and 2.}

\item{sigT}{The signficance to be used for the multiple
comparison test; the default is 5\%.}

\item{sigF}{The signficance to be used for the F test
of ANOVA; the default is 5\%.}

\item{unfold}{Says what must be done after the ANOVA.
If NULL (\emph{default}), recommended tests are performed;
if '0', just ANOVA is performed; if '1', the simple effects
are tested; if '2', the double interaction is unfolded.}
}
\value{
The output contains the ANOVA of the referred
RBD, the Shapiro-Wilk normality test for the residuals
of the model, the fitted regression models (when the
treatments are quantitative) and/or the multiple
comparison tests (when the treatments are qualitative).
}
\description{
\code{strip} Analysis Strip-plot experiments.
}
\details{
The arguments sigT and mcomp will be used only
when the treatment are qualitative.
}
\note{
The \code{\link{graphics}} can be used to
construct regression plots and \code{\link{plotres}}
for residuals plots.
}
\examples{
data(ex5)
attach(ex5)
strip(trat, genero, bloco, sabor, quali = c(TRUE,TRUE),
mcomp = "tukey", fac.names = c("Amostras","Genero"),
sigT = 0.05, sigF = 0.05, unfold=NULL)
}
\seealso{
\code{\link{split2.rbd}} and \code{\link{rbd}}.
}
\author{
Eric B Ferreira,
 \email{eric.ferreira@unifal-mg.edu.br}

Laís Brambilla Storti Ferreira
}
