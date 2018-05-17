\name{plotres}
\alias{plotres}

\title{Set of plots residuals}
\description{Residual plots for a output model. Four sets of plots are produced: (1) Histogram, (2) normal probability plot for the residual, (3) Standardized Residuals versus Fitted Values, and (4) box-plot (Standardized Residuals).}

\usage{
plotres(x)
}

\arguments{
\item{x}{Output from anova (performed in ExpDes)}
}

%\details{

%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

\references{
STEEL, R. G. D.; TORRIE, J. H. \emph{Principles and procedures in Statistics: a biometrical approach}. McGraw-Hill, New York, NY. 1980.
}

\author{Denismar Alves Nogueira}

\note{
The default "printing" of the result will produce four plots of any model of the
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See also:\code{\link{graphics}}}

\examples{
data(ex1)
attach(ex1)
a<-crd(trat, ig)
plotres(a)
}