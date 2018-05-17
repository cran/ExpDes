\name{reg.nl}
\alias{reg.nl}

\title{Non-linear Regression}

\description{Adjusts non-linear regression models in Anova (Models: Power, Exponential, Logistic, Gompertz).}

\usage{reg.nl(resp, treat)}

\arguments{
\item{resp}{Numeric or complex vector containing the response variable.}

\item{treat}{Numeric or complex vector containing the treatments.}
}

%\details{

\value{Returns coefficients, significance and ANOVA of the fitted regression models.}

\references{DRAPER, N.R.; SMITH, H. \emph{Apllied regression analysis}. 3ed. New York : John Wiley, 1998. 706p.}

\author{
Luiz Alberto Beijo

Eric Batista Ferreira
}

%\note{

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ See too: \code{\link{graphics}}}

\examples{
data(ex_nl)
attach(ex_nl)
x<-crd(trat, resp, quali=FALSE, nl=TRUE)
graphics(x,degree = "log")
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
