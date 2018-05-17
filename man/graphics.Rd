\name{graphics}
\alias{graphics}

\title{Regression model plots}

\description{Plots from regression models fitted in ANOVA}

\usage{graphics(a, degree = 1, mod = TRUE, main = " ", sub = " ",
         xlab = "Levels (X)", ylab = "Response var (Y)", pch = 19,
         xlim = NULL, ylim = NULL, bty = "o") }

\arguments{
\item{a}{Output from anova (performed in ExpDes)}
\item{degree}{For polynomial models, 1 (linear model) is the \emph{default}, 2 (quadratic model), 3 (cubic model), "pot" (Power model), "log" (Logistic model), "gom" (Gompertz model) and "exp" (Exponential model).}
\item{mod}{Logic. Print the model expression and its R2 on the top of the graphic. The \emph{default} is TRUE.}
\item{main}{Title of the plot. Empty is the \emph{default}.}
\item{sub}{Subtitle of the plot. Empty is the \emph{default}.}
\item{xlab}{Name for axis X.}
\item{ylab}{Name for axis Y.}
\item{pch}{Caracter type to be used on the observed values.}
\item{xlim}{Limits for axis X.}
\item{ylim}{Limits for axis Y.}
\item{bty}{Type of box the plot is fitted in.}
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

\author{Eric Batista Ferreira}

%\note{
%%  ~~further notes~~


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{reg.poly}}, \code{\link{plotres}}}

\examples{
data(ex1)
attach(ex1)
a<-crd(trat, ig, quali=FALSE, nl=FALSE)
graphics(a, degree=1)
graphics(a, degree=2)
graphics(a, degree=3)
}