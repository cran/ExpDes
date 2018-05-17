\name{ExpDes-package}
\alias{ExpDes-package}
\alias{ExpDes}
\docType{package}

\title{Experimental Designs}

\description{
Package for analysis of simple experimental designs (CRD, RBD and LSD), experiments in double factorial schemes (in CRD and RBD), experiments in a split plot in time schemes (in CRD and RBD), experiments in double factorial schemes with an additional treatment (in CRD and RBD), experiments in triple factorial scheme (in CRD and RBD) and experiments in triple factorial schemes with an additional treatment (in CRD and RBD), performing the analysis of variance and means comparison by fitting regression models until the third power (quantitative treatments) or by a multiple comparison test, Tukey test, test of Student-Newman-Keuls (SNK), Scott-Knott, Duncan test, t test (LSD), Bonferroni t test (protected LSD) and bootstrap multiple comparison's test (qualitative treaments); residual analysis (Ferreira, Cavalcanti and Nogueira, 2014) <doi:10.4236/am.2014.519280>.}
\details{
\tabular{ll}{
Package: \tab ExpDes\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2018-05-17\cr
License:  \tab GPL 2\cr
LazyLoad: \tab yes\cr
}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}
\author{
Eric Batista Ferreira

Portya Piscitelli Cavalcanti

Denismar Alves Nogueira

Maintainer: Eric Batista Ferreira <eric@unifal-mg.edu.br>
}

\references{
BANZATTO, D. A.; KRONKA, S. N. Experimentacao Agricola. 4 ed. Jaboticabal: Funep. 2006. 237 p.

COMASTRI FILHO, J. A. Avaliacao de especies de forrageiras nativas e exoticas na sub-regiao dos paiaguas no pantanal mato-grossense. Pesq. Agropec.
Bras., Brasilia, v.29, n.6, p. 971-978, jun. 1994.

FERREIRA, E. B.; CAVALCANTI, P. P. Funcao em codigo R para analisar experimentos em DIC simples, em uma so rodada.
In: REUNIAO ANUAL DA REGIAO BRASILEIRA DA SOCIEDADE INTERNACIONAL DE BIOMETRIA, 54o SIMPOSIO DE ESTATISTICA APLICADA A EXPERIMENTACAO AGRONOMICA,
13, 2009, Sao Carlos. Programas e resumos... Sao Carlos, SP: UFSCar, 2009. p. 1-5.

FERREIRA, E. B.; CAVALCANTI, P. P.; NOGUEIRA D. A. Funcao em codigo R para analisar experimentos em DBC simples, em uma so rodada.
In: JORNADA CIENTIFICA DA UNIVERSIDADE FEDERAL DE ALFENAS-MG, 2, 2009, Alfenas. Resumos... ALfenas: Unifal-MG, 2009.

FERREIRA, E. B.; CAVALCANTI, P. P.; NOGUEIRA D. A. Funcao em codigo R para analisar experimentos em DQL simples, em uma so rodada.
In: CONGRESSO DE POS-GRADUACAO DA UNIVERSIDADE FEDERAL DE LAVRAS, 18., 2009, Lavras. Resumos... Lavras: UFLA, 2009.

FERREIRA, E. B.; CAVALCANTI, P. P.; NOGUEIRA D. A. Funcao para analisar experimentos em fatorial duplo com um tratamento adicional, em uma so rodada.
In: CONGRESSO DE POS-GRADUACAO DA UNIVERSIDADE FEDERAL DE LAVRAS, 19., 2010, Lavras. Resumos... Lavras: UFLA, 2010.

GOMES, F. P. Curso de Estatistica Experimental. 10a ed. Piracicaba: ESALQ-USP. 1982. 430 p.

GUERRA, A. R. Atributos de solo sob coberturas vegetais em sistema silvipastoril em Lavras - MG. 2010. 141 p. Dissertacao (Mestrado em Engenharia
 Florestal) - Universidade Federal de Lavras, UFLA, Lavras, 2010.

HEALY, M. J. R. The analysis of a factorial experiment with additional treatments. Journal of Agricultural Science, Cambridge, v. 47, p. 205-206. 1956.

MOREIRA, D. K. T. Extrudados expandidos de arroz, soja e gergelim para uso em barras alimenticias. 2010. 166p. Dissertacao (Mestrado em Ciencias dos
 Alimentos) - Universidade Federal de Lavras, UFLA, Lavras, 2010.

PAIVA, A. P. de. Estudos tecnologicos, quimico, fisico-quimico e sensorial de barras alimenticias elaboradas com subprodutos e residuos agoindustriais.
 2008. 131p. Dissertacao (Mestrado em Ciencias dos Alimentos) - Universidade Federal de Lavras, UFLA, Lavras, 2008.

RAMALHO, M. A. P.; FERREIRA, D. F.; OLIVEIRA, A. C. de. Experimentacao em Genetica e Melhoramento de Plantas. 2a ed. Lavras: UFLA. 2005. 300p.

RAMOS, P. S., FERREIRA, D. F. Agrupamento de medias via bootstrap de populacoes normais e nao-normais, Revista Ceres, v.56, p.140-149, 2009.

REZENDE, F. A. de. Aproveitamento  da casca de cafe e borra da purificacao de gorduras e oleos residuarios em compostagem. 2010. 74p. Tese (Doutorado
 em Agronomia - Fitotecnia) - Universidade Federal de Lavras, UFLA, Lavras, 2010.

RIBEIRO, J. de A. Estudos quimicos e bioquimicos do Yacon (Samallanthus sonchifolius) in natura e processado e influencia do seu consumo sobre niveis
glicemicos e lipideos fecais de ratos. 2008. 166p. Dissertacao (Mestrado em Ciencias dos Alimentos) - Universidade Federal de Lavras, UFLA, Lavras,
 2008.

RODRIGUES, R. B. Danos do percevejo-barriga-verde Dichelops melacanthus (Dallas, 1851) (Hemiptera: Pentatomidae) na cultura do milho. 2011. 105f.
Dissertacao (Mestrado em Agronomia - Universidade Federal de Santa Maria, Santa Maria, 2011.

STEEL, R. G. D.; TORRIE, J. H.; DICKEY, D. A. Principles and procedures of statistics: a biometrical approach. 3rd Edition. 1997. 666 p.

VENABLES, W. N.; RIPLEY, B. D. Modern Applied Statistics with S-PLUS. Third Edition. Springer. 1999. p. 100.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }

\seealso{
see: \code{\link{crd}}, \code{\link{rbd}}, \code{\link{fat2.crd}}, \code{\link{fat3.crd}}, \code{\link{split2.crd}}, \code{\link{fat2.ad.crd}}, \code{\link{fat3.ad.crd}}, \code{\link{fat2.rbd}}, \code{\link{fat3.rbd}}, \code{\link{split2.rbd}}, \code{\link{strip}}, \code{\link{fat2.ad.rbd}} e \code{\link{fat3.ad.rbd}} 
}
%\examples{
%~~ simple examples of the most important functions ~~
