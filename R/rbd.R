#' Randomized Blocks Design
#'
#' \code{rbd} Analyses experiments in balanced Randomized
#' Blocks Designs under one single factor, considering a fixed
#' model.
#' @param treat Numeric or complex vector containing the
#' treatments.
#' @param block Numeric or complex vector containing the blocks.
#' @param resp Numeric or complex vector containing the
#' response variable.
#' @param quali Logic. If TRUE (default), the treatments are
#' assumed qualitative, if FALSE, quantitatives.
#' @param mcomp Allows choosing the multiple comparison test;
#' the \emph{default} is the test of Tukey, however, the
#' options are: the LSD test ('lsd'), the LSD test with
#' Bonferroni protection ('lsdb'), the test of Duncan
#' ('duncan'), the test of Student-Newman-Keuls ('snk'),
#' the test of Scott-Knot ('sk'), the Calinski and Corsten
#' test ('ccF') and bootstrap multiple comparison's test
#' ('ccboot').
#' @param nl Logic. If FALSE (\emph{default}) linear regression
#' models are adjusted. IF TRUE, non-linear regression models
#' are adjusted.
#' @param hvar Allows choosing the test for homogeneity of
#' variances; the \emph{default} is the test of Bartlett,
#' however there are other options: test of Levene ('levene'),
#' test of Samiuddin ('samiuddin'), test of ONeill and Mathews
#' ('oneillmathews') and the Layard test ('layard').
#' @param sigT The signficance to be used for the multiple
#' comparison test; the default is 5\%.
#' @param sigF The signficance to be used for the F test of
#' ANOVA; the default is 5\%.
#' @details The arguments sigT and mcomp will be used only when
#' the treatment are qualitative.
#' @return The output contains the ANOVA of the RBD, the
#' Shapiro-Wilk normality test for the residuals of the model,
#' the fitted regression models (when the treatments are
#' quantitative) and/or the multiple comparison tests (when
#' the treatments are qualitative).
#' @references BANZATTO, D. A.; KRONKA, S. N. Experimentacao
#' Agricola. 4 ed. Jaboticabal: Funep. 2006. 237 p.
#'
#' FERREIRA, E. B.; CAVALCANTI, P. P.; NOGUEIRA D. A. Funcao
#' em codigo R para analisar experimentos em DBC simples, em
#' uma so rodada. In: JORNADA CIENTIFICA DA UNIVERSIDADE
#' FEDERAL DE ALFENAS-MG, 2., 2009, Alfenas. Annals...
#' ALfenas: Unifal-MG, 2009.
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#' @author Denismar Alves Nogueira
#' @author Portya Piscitelli Cavalcanti
#' @note The \code{\link{graphics}} can be used to construct
#' regression plots and \code{\link{plotres}} for residuals
#' plots.
#' @seealso \code{\link{fat2.rbd}}, \code{\link{fat3.rbd}},
#' \code{\link{split2.rbd}}, \code{\link{strip}},
#' \code{\link{fat2.ad.rbd}} and \code{\link{fat3.ad.rbd}}.
#' @examples
#' data(ex2)
#' attach(ex2)
#' rbd(trat, provador, aparencia, quali = TRUE, mcomp = "lsd",
#' hvar = "oneillmathews", sigT = 0.05, sigF = 0.05)
#' @export

rbd <-function(treat, block, resp, quali=TRUE, mcomp='tukey', nl=FALSE,
              hvar='oneillmathews', sigT=0.05, sigF=0.05) {

Trat<-factor(treat)
Bloco<-factor(block)
anava<-aov(resp~Trat+Bloco)
tab<-summary(anava)
colnames(tab[[1]])<-c('DF','SS','MS','Fc','Pr>Fc')
tab[[1]]<-rbind(tab[[1]],c(apply(tab[[1]],2,sum)))
rownames(tab[[1]])<-c('Treatament','Block','Residuals','Total')
cv<-round(sqrt(tab[[1]][3,3])/mean(resp)*100, 2)
tab[[1]][4,3]=NA
cat('------------------------------------------------------------------------
Analysis of Variance Table\n------------------------------------------------------------------------\n')
print(tab[[1]])
cat('------------------------------------------------------------------------\nCV =',cv,'%\n')

#Normality test
pvalor.shapiro<-shapiro.test(anava$residuals)$p.value
cat('\n------------------------------------------------------------------------\nShapiro-Wilk normality test\n')
cat('p-value: ',pvalor.shapiro, '\n')
if(pvalor.shapiro<0.05){cat('WARNING: at 5% of significance, residuals can not be considered normal!
------------------------------------------------------------------------\n')}
else{cat('According to Shapiro-Wilk normality test at 5% of significance, residuals can be considered normal.
------------------------------------------------------------------------\n')}

#Homogeneity of variances test
if(hvar=='oneillmathews') pvalor.hvar<-oneilldbc(resp, Trat, Bloco)
if(hvar=='han') pvalor.hvar<-han(resp, Trat, Bloco)
if(hvar=='anscombetukey') pvalor.hvar<-anscombetukey(resp, Trat, Bloco, tab[[1]][3,1], as.numeric(tab[[1]][3,3]), tab[[1]][1,2], tab[[1]][2,2], anava$residuals, anava$fitted.values)
cat('\n------------------------------------------------------------------------\nHomogeneity of variances test\n')
cat('p-value: ',pvalor.hvar, '\n')
if(pvalor.hvar<0.05){cat('WARNING: at 5% of significance, residuals can not be considered homocedastic!
------------------------------------------------------------------------\n')}
else{cat('According to the test of',hvar,'at 5% of significance, the variances can be considered homocedastic.
------------------------------------------------------------------------\n')}

if(tab[[1]][1,5]<sigF){

if(quali==TRUE) {

  if(mcomp=='tukey') tukey(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)
  if(mcomp=='duncan')duncan(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)
  if(mcomp=='lsd')   lsd(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)
  if(mcomp=='lsdb')  lsdb(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)
  if(mcomp=='sk')    scottknott(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)
  if(mcomp=='snk')   snk(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)
  if(mcomp=="ccboot")ccboot(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)
  if(mcomp=="ccF")   ccF(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)

                }
else if(nl==FALSE) reg<-reg.poly(resp, treat, tab[[1]][3,1], tab[[1]][3,2], tab[[1]][1,1], tab[[1]][1,2])
else if(nl==TRUE)  reg<-reg.nl(resp, treat)
                       }
else {
    cat('\nAccording to the F test, the means can not be considered distinct.\n')
mean.table<-tapply.stat(resp,treat,mean)
colnames(mean.table)<-c('Levels','Means')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}

#Saida
out<-list()
out$residuals<-anava$residuals
out$df.residual<-anava$df.residual
out$coefficients<-anava$coefficients
out$effects<-anava$effects
out$fitted.values<-anava$fitted.values
out$means<-tapply.stat(resp,treat,mean)
if(quali==FALSE && tab[[1]][1,5]<sigF) {out$reg<-reg}
invisible(out)
}
