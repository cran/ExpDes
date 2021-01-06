#' One factor Completely Randomized Design
#'
#' \code{crd} Analyses balanced experiments in Completely
#' Randomized Design under one single factor, considering a
#' fixed model.
#' @param treat Numeric or complex vector containing the
#' treatments.
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
#' @return The output contains the ANOVA of the CRD, the
#' Shapiro-Wilk normality test for the residuals of the model,
#' the fitted regression models (when the treatments are
#' quantitative) and/or the multiple comparison tests (when
#' the treatments are qualitative).
#' @references BANZATTO, D. A.; KRONKA, S. N. Experimentacao
#' Agricola. 4 ed. Jaboticabal: Funep. 2006. 237 p.
#'
#' FERREIRA, E. B.; CAVALCANTI, P. P. Funcao em codigo R para
#' analisar experimentos em DIC simples, em uma so rodada. In:
#' REUNIAO ANUAL DA REGIAO BRASILEIRA DA SOCIEDADE
#' INTERNACIONAL DE BIOMETRIA, 54./SIMPOSIO DE ESTATISTICA
#' APLICADA A EXPERIMENTACAO AGRONOMICA, 13., 2009, Sao Carlos.
#' Programas e resumos... Sao Carlos, SP: UFSCar, 2009. p. 1-5.
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#' @author Denismar Alves Nogueira
#' @author Portya Piscitelli Cavalcanti
#' @seealso \code{\link{fat2.crd}}, \code{\link{fat3.crd}},
#' \code{\link{split2.crd}}, \code{\link{fat2.ad.crd}} and
#' \code{\link{fat3.ad.crd}}.
#' @examples
#' data(ex1)
#' attach(ex1)
#' crd(trat, ig, quali = FALSE, sigF = 0.05)
#' @export

crd <-function(treat, resp, quali=TRUE, mcomp='tukey', nl=FALSE,
              hvar='bartlett', sigT=0.05, sigF=0.05) {

Trat<-factor(treat)
anova<-aov(resp~Trat)
tab<-summary(anova)

#### Computing the number of replications for each treatment ####
#i<-0
#ii<-1
#rr<-1
t<-length(levels(Trat))
#r<-matrix(0,t,1)
#for(i in 1:(length(treat)-1)) {
#  if (Trat[i]==Trat[i+1]) {rr<-rr+1} else {r[ii]<-rr}
#  if (Trat[i]!=Trat[i+1]) rr<-1
#  if (Trat[i]!=Trat[i+1]) ii<-ii+1
#  if ((i+1)==length(treat)){r[ii]<-rr}
#}
r<-as.numeric(table(Trat))
#########################################

colnames(tab[[1]])<-c('DF','SS','MS','Fc','Pr>Fc')
tab[[1]]<-rbind(tab[[1]],c(apply(tab[[1]],2,sum)))
rownames(tab[[1]])<-c('Treatament','Residuals','Total')
cv<-round(sqrt(tab[[1]][2,3])/mean(resp)*100, 2)
#######print
#stargazer(tab[[1]], type = "text", title="Analysis of Variance Table", digits=4, flip=FALSE, summary=F, digits.extra=4)
#cat('CV =',cv,'%\n')
#######
tab[[1]][3,3]=NA
cat('------------------------------------------------------------------------
Analysis of Variance Table\n------------------------------------------------------------------------\n')
print(tab[[1]])
cat('------------------------------------------------------------------------\nCV =',cv,'%\n')


#Normality test
pvalor.shapiro<-shapiro.test(anova$residuals)$p.value
cat('\n------------------------------------------------------------------------\nShapiro-Wilk normality test\n')
cat('p-value: ',pvalor.shapiro, '\n')
if(pvalor.shapiro<0.05){cat('WARNING: at 5% of significance, residuals can not be considered normal!
------------------------------------------------------------------------\n')}
else{cat('According to Shapiro-Wilk normality test at 5% of significance, residuals can be considered normal.
------------------------------------------------------------------------\n')}

#Homogeneity of variances test
if(hvar=='bartlett') pvalor.hvar<-bartlett(treat,resp,t,r)
if(hvar=='levene') pvalor.hvar<-levene(treat,resp,t,r)
if(hvar=='samiuddin') pvalor.hvar<-samiuddin(treat,resp,t,r)
if(hvar=='oneillmathews') pvalor.hvar<-oneillmathews(treat,resp,t,r)
if(hvar=='layard') pvalor.hvar<-layard(treat,resp,t,r)

cat('\n------------------------------------------------------------------------\nHomogeneity of variances test\n')
cat('p-value: ',pvalor.hvar, '\n')
if(pvalor.hvar<0.05){cat('WARNING: at 5% of significance, residuals can not be considered homocedastic!
------------------------------------------------------------------------\n')}
else{cat('According to the test of',hvar,'at 5% of significance, residuals can be considered homocedastic.
------------------------------------------------------------------------\n')}


if(tab[[1]][1,5]<sigF) {


if(quali==TRUE) {

  if(mcomp=='tukey') tukey(resp,Trat,tab[[1]][2,1],tab[[1]][2,2],sigT)
  if(mcomp=='duncan')duncan(resp,Trat,tab[[1]][2,1],tab[[1]][2,2],sigT)
  if(mcomp=='lsd')   lsd(resp,Trat,tab[[1]][2,1],tab[[1]][2,2],sigT)
  if(mcomp=='lsdb')  lsdb(resp,Trat,tab[[1]][2,1],tab[[1]][2,2],sigT)
  if(mcomp=='sk')    scottknott(resp,Trat,tab[[1]][2,1],tab[[1]][2,2],sigT)
  if(mcomp=='snk')   snk(resp,Trat,tab[[1]][2,1],tab[[1]][2,2],sigT)
  if(mcomp=="ccboot")ccboot(resp,Trat,tab[[1]][2,1],tab[[1]][2,2],sigT)
  if(mcomp=="ccF")   ccF(resp,Trat,tab[[1]][2,1],tab[[1]][2,2],sigT)

               }

else if(nl==FALSE) reg<-reg.poly(resp, treat, tab[[1]][2,1],
                                   tab[[1]][2,2], tab[[1]][1,1], tab[[1]][1,2])
else if(nl==TRUE)  reg<-reg.nl(resp, treat)

                       }

else {
    cat('\nAccording to the F test, the means can not be considered distinct.\n')
    cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp,treat,mean)
colnames(mean.table)<-c('Levels','Means')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}

#Saida
out<-list()
out$residuals<-anova$residuals
out$df.residual<-anova$df.residual
out$coefficients<-anova$coefficients
out$effects<-anova$effects
out$fitted.values<-anova$fitted.values
out$means<-tapply.stat(resp,treat,mean)
if(quali==FALSE && tab[[1]][1,5]<sigF) {out$reg<-reg}
invisible(out)
}
