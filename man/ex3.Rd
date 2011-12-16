\name{ex3}
\alias{ex3}
\docType{data}
\title{
Forage: LSD
}
\description{
Data from an experiment aiming to select forage for minimizing the intake problem of feeding cattle in the sub-region of Paiaguas.
}
\usage{data(ex3)}
\format{
  A data frame with 49 observations on the following 4 variables.
  \describe{
    \item{\code{trat}}{a factor with levels \code{A} \code{B} \code{C} \code{D} \code{E} \code{F} \code{G}}
    \item{\code{linha}}{a numeric vector}
    \item{\code{coluna}}{a numeric vector}
    \item{\code{resp}}{a numeric vector}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~

%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~

\references{
COMASTRI FILHO, J. A. Avaliacao de especies de forrageiras nativas e exoticas na sub-regiao dos paiaguas no pantanal mato-grossense. 
Pesq. Agropec. Bras., Brasilia, v.29, n.6, p. 971-978, jun. 1994.
}
\examples{
data(ex3)
## maybe str(ex3) ; plot(ex3) ...
}
%\keyword{datasets}
