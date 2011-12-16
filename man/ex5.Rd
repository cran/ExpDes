\name{ex5}
\alias{ex5}
\docType{data}
\title{
Food bars: Double Factorial scheme in RBD
}
\description{
Data adapted from a sensorial experiment where panelists of different genders evaluated the taste of food bars.
}
\usage{data(ex5)}
\format{
  A data frame with 160 observations on the following 4 variables.
  \describe{
    \item{\code{trat}}{a factor with levels \code{10g} \code{15g} \code{15t} \code{20t}}
    \item{\code{genero}}{a factor with levels \code{F} \code{M}}
    \item{\code{bloco}}{a numeric vector}
    \item{\code{sabor}}{a numeric vector}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~

%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~

\references{
MOREIRA, D. K. T. Extrudados Expandidos de Arroz, Soja e Gergelim para Uso em Barras Alimenticias. 2010. 166p. Dissertation
(Master in Food Science) - Universidade Federal de Lavras, UFLA, Lavras, 2010.
}
\examples{
data(ex5)
## maybe str(ex5) ; plot(ex5) ...
}
%\keyword{datasets}
