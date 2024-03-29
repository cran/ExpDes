#' Double factorial scheme in RBD
#'
#' \code{fat2.rbd} Analyses experiments in balanced
#' Randomized Blocks Designs in double factorial scheme,
#' considering a fixed model.
#' @param factor1 Numeric or complex vector containing the
#' factor 1 levels.
#' @param factor2 Numeric or complex vector containing the
#' factor 2 levels.
#' @param block Numeric or complex vector containing the
#' blocks.
#' @param resp Numeric or complex vector containing the
#' response variable.
#' @param quali Logic. If TRUE (default), the treatments
#' are assumed qualitative, if FALSE, quantitatives.
#' @param mcomp Allows choosing the multiple comparison
#' test; the \emph{default} is the test of Tukey, however,
#' the options are: the LSD test ('lsd'), the LSD test with
#' Bonferroni protection ('lsdb'), the test of Duncan
#' ('duncan'), the test of Student-Newman-Keuls ('snk'),
#' the test of Scott-Knott ('sk'), the Calinski and Corsten
#' test ('ccF') and bootstrap multiple comparison's test
#' ('ccboot').
#' @param fac.names Allows labeling the factors 1 and 2.
#' @param sigT The signficance to be used for the multiple
#' comparison test; the default is 5\%.
#' @param sigF The signficance to be used for the F test of
#' ANOVA; the default is 5\%.
#' @param unfold Says what must be done after the ANOVA.
#' If NULL (\emph{default}), recommended tests are performed;
#' if '0', just ANOVA is performed; if '1', the simple effects
#' are tested; if '2', the double interaction is unfolded.
#' @details The arguments sigT and mcomp will be used only
#' when the treatment are qualitative.
#' @return The output contains the ANOVA of the referred
#' RBD, the Shapiro-Wilk normality test for the residuals
#' of the model, the fitted regression models (when the
#' treatments are quantitative) and/or the multiple
#' comparison tests (when the treatments are qualitative).
#' @references BANZATTO, D. A.; KRONKA, S. N.
#' Experimentacao Agricola. 4 ed. Jaboticabal: Funep.
#' 2006. 237 p.
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#' @author Denismar Alves Nogueira
#' @author Portya Piscitelli Cavalcanti
#' @note The \code{\link{graphics}} can be used to
#' construct regression plots and \code{\link{plotres}}
#' for residuals plots.
#' @seealso \code{\link{fat3.rbd}},
#' \code{\link{split2.rbd}}, \code{\link{strip}},
#' \code{\link{fat2.ad.rbd}} and \code{\link{fat3.ad.rbd}}.
#' @examples
#' data(ex5)
#' attach(ex5)
#' fat2.rbd(trat, genero, bloco, sabor ,quali =
#' c(TRUE,TRUE), mcomp = "lsd", fac.names = c("Samples",
#' "Gender"), sigT = 0.05, sigF = 0.05, unfold=NULL)
#' @export

fat2.rbd <- function(factor1,
                     factor2,
                     block,
                     resp,
                     quali=c(TRUE,TRUE),
                     mcomp='tukey',
                     fac.names=c('F1','F2'),
                     sigT=0.05,
                     sigF=0.05,
                     unfold=NULL) {

cat('------------------------------------------------------------------------\nLegend:\n')
cat('FACTOR 1: ',fac.names[1],'\n')
cat('FACTOR 2: ',fac.names[2],'\n------------------------------------------------------------------------\n\n')

fatores<-cbind(factor1,factor2)
Fator1<-factor(factor1)
Fator2<-factor(factor2)
Bloco<-factor(block)
nv1<-length(summary(Fator1))   #Diz quantos niveis tem o fator 1.
nv2<-length(summary(Fator2))   #Diz quantos niveis tem o fator 2.
J<-length(summary(Bloco))
lf1<-levels(Fator1)
lf2<-levels(Fator2)
lfB<-levels(Bloco)
anava<-aov(resp~ Bloco + Fator1*Fator2)
tab<-summary(anava)
colnames(tab[[1]])<-c('DF','SS','MS','Fc','Pr>Fc')
tab[[1]]<-rbind(tab[[1]],c(apply(tab[[1]],2,sum)))
rownames(tab[[1]])<-c('Block',fac.names[1],fac.names[2],paste(fac.names[1],'*',fac.names[2],sep=''),'Residuals','Total')
cv<-round(sqrt(tab[[1]][5,3])/mean(resp)*100, 2)
tab[[1]][6,3]=' '
cat('\nAnalysis of Variance Table\n------------------------------------------------------------------------\n')
print(tab[[1]])
cat('------------------------------------------------------------------------\nCV =',cv,'%\n')

#Teste de normalidade
pvalor.shapiro<-shapiro.test(anava$residuals)$p.value
cat('\n------------------------------------------------------------------------\nShapiro-Wilk normality test\n')
cat('p-value: ',pvalor.shapiro, '\n')
if(pvalor.shapiro<0.05){cat('WARNING: at 5% of significance, residuals can not be considered normal!
------------------------------------------------------------------------\n')}
else{cat('According to Shapiro-Wilk normality test at 5% of significance, residuals can be considered normal.
------------------------------------------------------------------------\n')}

# Creating unfold #########################################
if(is.null(unfold)){
  if(tab[[1]][4,5]>sigF)  {unfold<-c(unfold,1)}
  if(tab[[1]][4,5]<=sigF) {unfold<-c(unfold,2)}
}

#Para interacao nao significativa, fazer...
if(any(unfold==1)) {
cat('\nNo significant interaction: analyzing the simple effect
------------------------------------------------------------------------\n')
fatores<-data.frame('fator 1'=factor1,'fator 2' = factor2)

for(i in 1:2){

#Para os fatores QUALITATIVOS, teste de Tukey
if(quali[i]==TRUE && tab[[1]][i+1,5]<=sigF) {
    cat(fac.names[i])
      if(mcomp=='tukey'){
    tukey(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                    }
  if(mcomp=='duncan'){
    duncan(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                    }
  if(mcomp=='lsd'){
    lsd(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                    }
  if(mcomp=='lsdb'){
    lsdb(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                    }
  if(mcomp=='sk'){
    scottknott(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                    }
  if(mcomp=='snk'){
    snk(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                    }
  if(mcomp=='ccboot'){
    ccboot(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                    }
  if(mcomp=="ccF"){
    ccF(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                  }
                   }
if(quali[i]==TRUE && tab[[1]][i+1,5]>sigF) {
    cat(fac.names[i])
    cat('\nAccording to the F test, the means of this factor are statistical equal.\n')
    cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Levels','Means')
print(mean.table)
cat('------------------------------------------------------------------------')
}

#Para os fatores QUANTITATIVOS, regressao
if(quali[i]==FALSE && tab[[1]][i+1,5]<=sigF){
    cat(fac.names[i])
    reg.poly(resp, fatores[,i], tab[[1]][5,1], tab[[1]][5,2], tab[[1]][i+1,1], tab[[1]][i+1,2])
}

if(quali[i]==FALSE && tab[[1]][i+1,5]>sigF) {
    cat(fac.names[i])
    cat('\nAccording to the F test, the means of this factor are statistical equal.\n')
    cat('------------------------------------------------------------------------\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Levels','Means')
print(mean.table)
cat('------------------------------------------------------------------------')
}
cat('\n')
}
}

#Se a interacao for significativa, desdobrar a interacao
if(any(unfold==2)) {
cat("\n\n\nSignificant interaction: analyzing the interaction
------------------------------------------------------------------------\n")

#Desdobramento de FATOR 1 dentro do niveis de FATOR 2
cat("\nAnalyzing ", fac.names[1], ' inside of each level of ', fac.names[2], '
------------------------------------------------------------------------\n')

des1<-aov(resp~ Bloco + Fator2/Fator1)

l1<-vector('list',nv2)
names(l1)<-names(summary(Fator2))
v<-numeric(0)
for(j in 1:nv2) {
        for(i in 0:(nv1-2)) v<-cbind(v,i*nv2+j)
        l1[[j]]<-v
        v<-numeric(0)
                }
des1.tab<-summary(des1,split=list('Fator2:Fator1'=l1))[[1]]
#Montando a tabela de ANAVA do des1
glB=tab[[1]][1,1]
glb=nv2-1
glf1=c(as.numeric(des1.tab[4:(nv2+3),1]))
glE=tab[[1]][5,1]
glT=tab[[1]][6,1]

SQB=tab[[1]][1,2]
SQb=tab[[1]][3,2]
SQf1=c(as.numeric(des1.tab[4:(nv2+3),2]))
SQE=tab[[1]][5,2]
SQT=tab[[1]][6,2]

QMB=SQB/glB
QMb=SQb/glb
QMf1=SQf1/glf1
QME=SQE/glE
QMT=SQT/glT

FcB=QMB/QME
Fcb=QMb/QME
Fcf1=QMf1/QME

rn<-numeric(0)
for(j in 1:nv2){ rn<-c(rn, paste(paste(fac.names[1],':',fac.names[2],sep=''),lf2[j]))}

anavad1<-data.frame("DF"=c(round(c(glB, glb, glf1, glE, glT))),
"SS"=c(round(c(SQB,SQb,SQf1,SQE,SQT),5)),
"MS"=c(round(c(QMB,QMb,QMf1,QME),5),''),
"Fc"=c(round(c(FcB,Fcb,Fcf1),4),'',''),
"Pr>Fc"=c(round(c(1-pf(FcB,glB,glE),1-pf(Fcb,glb,glE),1-pf(Fcf1,glf1,glE)),4),'', ''))
rownames(anavad1)=c("Block",fac.names[2],rn,"Residuals","Total")
cat('------------------------------------------------------------------------
Analysis of Variance Table\n------------------------------------------------------------------------\n')
print(anavad1)
cat('------------------------------------------------------------------------\n\n')


for(i in 1:nv2) {
  if(des1.tab[(i+3),5]<=sigF){
    if(quali[1]==TRUE){
                      cat('\n\n',fac.names[1],' inside of the level ',lf2[i],' of ',fac.names[2],'
------------------------------------------------------------------------')
                        if(mcomp=='tukey'){
                          tukey(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                          }
                        if(mcomp=='duncan'){
                          duncan(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                           }
                        if(mcomp=='lsd'){
                          lsd(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                        }
                        if(mcomp=='lsdb'){
                          lsdb(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                         }
                        if(mcomp=='sk'){
                          scottknott(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                       }
                        if(mcomp=='snk'){
                          snk(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                        }
                        if(mcomp=='ccboot'){
                          ccboot(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                        }
                        if(mcomp=="ccF"){
                          ccF(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                        }
                      }
    else{  #regressao
    cat('\n\n',fac.names[1],' inside of the level ',lf2[i],' of ',fac.names[2],'
------------------------------------------------------------------------')
    reg.poly(resp[Fator2==lf2[i]], factor1[Fator2==lf2[i]], tab[[1]][5,1], tab[[1]][5,2], des1.tab[i+3,1], des1.tab[i+3,2][[1]])
        }
                              }
    else{cat('\n\n',fac.names[1],' inside of the level ',lf2[i],' of ',fac.names[2],'\n')
    cat('\nAccording to the F test, the means of this factor are statistical equal.\n')
    cat('------------------------------------------------------------------------\n')
        mean.table<-tapply.stat(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],mean)
        colnames(mean.table)<-c('  Levels','    Means')
        print(mean.table)
        cat('------------------------------------------------------------------------\n')
        }
                 }
cat('\n\n')

#Desdobramento de FATOR 2 dentro do niveis de FATOR 1
cat("\nAnalyzing ", fac.names[2], ' inside of each level of ', fac.names[1], '
------------------------------------------------------------------------\n')

des2<-aov(resp~ Bloco + Fator1/Fator2)

l2<-vector('list',nv1)
names(l2)<-names(summary(Fator1))
v<-numeric(0)
for(j in 1:nv1) {
        for(i in 0:(nv2-2)) v<-cbind(v,i*nv1+j)
        l2[[j]]<-v
        v<-numeric(0)
                }
des2.tab<-summary(des2,split=list('Fator1:Fator2'=l2))[[1]]
#Montando a tabela de ANAVA do des2
gla=nv1-1
glf2=c(as.numeric(des2.tab[4:(nv1+3),1]))

SQa=tab[[1]][2,2]
SQf2=c(as.numeric(des2.tab[4:(nv1+3),2]))

QMa=SQa/gla
QMf2=SQf2/glf2

Fca=QMa/QME
Fcf2=QMf2/QME

rn<-numeric(0)
for(i in 1:nv1){ rn<-c(rn, paste(paste(fac.names[2],':',fac.names[1],sep=''),lf1[i]))}

anavad2<-data.frame("DF"=c(round(c(glB, gla, glf2, glE, glT))),
"SS"=c(round(c(SQB,SQa,SQf2,SQE,SQT),5)),
"MS"=c(round(c(QMB,QMa,QMf2,QME),5),''),
"Fc"=c(round(c(FcB,Fca,Fcf2),4),'',''),
"Pr>Fc"=c(round(c(1-pf(FcB,glB,glE),1-pf(Fca,gla,glE),1-pf(Fcf2,glf2,glE)),4),'', ''))
rownames(anavad2)=c("Block",fac.names[1],rn,"Residuals","Total")
cat('------------------------------------------------------------------------
Analysis of Variance Table\n------------------------------------------------------------------------\n')
print(anavad2)
cat('------------------------------------------------------------------------\n\n')


for(i in 1:nv1) {
  if(des2.tab[(i+3),5]<=sigF){
    if(quali[2]==TRUE){
                      cat('\n\n',fac.names[2],' inside of the level ',lf1[i],' of ',fac.names[1],'
------------------------------------------------------------------------')
                      if(mcomp=='tukey'){
                          tukey(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                          }
                        if(mcomp=='duncan'){
                          duncan(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                           }
                        if(mcomp=='lsd'){
                          lsd(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                        }
                        if(mcomp=='lsdb'){
                          lsdb(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                         }
                        if(mcomp=='sk'){
                          scottknott(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                       }
                        if(mcomp=='snk'){
                          snk(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                        }
                        if(mcomp=='ccboot'){
                          ccboot(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                        }
                        if(mcomp=="ccF"){
                          ccF(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                        }
                        }
    else{  #regressao
        cat('\n\n',fac.names[2],' inside of the level ',lf1[i],' of ',fac.names[1],'
------------------------------------------------------------------------')
        reg.poly(resp[Fator1==lf1[i]], factor2[Fator1==lf1[i]], tab[[1]][5,1], tab[[1]][5,2], des2.tab[i+3,1], des2.tab[i+3,2][[1]])
        }
                             }
    else{cat('\n\n',fac.names[2],' inside of the level ',lf1[i],' of ',fac.names[1],'\n')
    cat('\nAccording to the F test, the means of this factor are statistical equal.\n')
    cat('------------------------------------------------------------------------\n')
        mean.table<-tapply.stat(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],mean)
        colnames(mean.table)<-c('  Levels','    Means')
        print(mean.table)
        cat('------------------------------------------------------------------------\n')
        }

                }
}
#Saida
out<-list()
out$residuals<-anava$residuals
out$df.residual<-anava$df.residual
out$coefficients<-anava$coefficients
out$effects<-anava$effects
out$fitted.values<-anava$fitted.values
out$means.factor1<-tapply.stat(resp,fatores[,1],mean)
out$means.factor2<-tapply.stat(resp,fatores[,2],mean)
tabmedia<-model.tables(anava, "means")
out$means.inside<-tabmedia$tables$`Fator1:Fator2`
#if(quali==FALSE && tab[[1]][1,5]<sigF) {out$reg<-reg}
invisible(out)
}
