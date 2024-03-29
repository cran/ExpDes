#' Strip-plot experiments
#'
#' \code{strip} Analysis Strip-plot experiments.
#' @param factor1 Numeric or complex vector containing the
#' factor 1 levels.
#' @param factor2 Numeric or complex vector containing the
#' factor 2 levels.
#' @param block Numeric or complex vector containing the blocks.
#' @param resp Numeric or complex vector containing the
#' response variable.
#' @param quali Logic. If TRUE (default), the treatments
#' are assumed qualitative, if FALSE, quantitatives.
#' @param mcomp Allows choosing the multiple comparison
#' test; the \emph{default} is the test of Tukey, however,
#' the options are: the LSD test ('lsd'), the LSD test
#' with Bonferroni protection ('lsdb'), the test of Duncan
#' ('duncan'), the test of Student-Newman-Keuls ('snk'),
#' the test of Scott-Knott ('sk'), the Calinski and
#' Corsten test ('ccF') and bootstrap multiple comparison's
#' test ('ccboot').
#' @param fac.names Allows labeling the factors 1 and 2.
#' @param sigT The signficance to be used for the multiple
#' comparison test; the default is 5\%.
#' @param sigF The signficance to be used for the F test
#' of ANOVA; the default is 5\%.
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
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#' @author Laís Brambilla Storti Ferreira
#' @note The \code{\link{graphics}} can be used to
#' construct regression plots and \code{\link{plotres}}
#' for residuals plots.
#' @seealso \code{\link{split2.rbd}} and \code{\link{rbd}}.
#' @examples
#' data(ex5)
#' attach(ex5)
#' strip(trat, genero, bloco, sabor, quali = c(TRUE,TRUE),
#' mcomp = "tukey", fac.names = c("Amostras","Genero"),
#' sigT = 0.05, sigF = 0.05, unfold=NULL)
#' @export

strip <- function(factor1,
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
cat('FACTOR 1 (Whole plot): ',fac.names[1],'\n')
cat('FACTOR 2 (strip-plot): ',fac.names[2],'\n------------------------------------------------------------------------\n\n')

cont<-c(2,4)                   #endereco das linhas dos fatores F1 e F2
Fator1<-factor(factor1)
Fator2<-factor(factor2)
bloco<-factor(block)
nv1<-length(summary(Fator1))   #Diz quantos niveis tem o fator 1.
nv2<-length(summary(Fator2))   #Diz quantos niveis tem o fator 2.

anava<-aov(resp ~ Fator1*Fator2 + Fator1*bloco + Fator2:bloco)
tab1<-summary(anava)
colnames(tab1[[1]])<-c('DF', 'SS', 'MS', 'Fc', 'Pr(>Fc)')
tab<-rbind(tab1[[1]],apply(tab1[[1]],2,sum))
tab<-rbind(tab[3,],tab[1,],tab[5,],tab[2,],tab[6,],tab[4,],tab[7,],tab[8,])
rownames(tab)<-c('Block',fac.names[1],'Error a',fac.names[2],'Error b',paste(fac.names[1],'*',fac.names[2],sep=''),'Error c','Total')
QMerrobloco<-tab[3,3] + tab[5,3] - tab[7,3]
tab[,4]<-c(tab[1,3]/QMerrobloco,tab[2,3]/tab[3,3],NA,tab[4,3]/tab[5,3],NA,tab[6,3]/tab[7,3],NA,NA)
tab[,5]<-c(1-pf(tab[1,4],tab[1,1],tab[3,1]),1-pf(tab[2,4],tab[2,1],tab[3,1]),NA,
           1-pf(tab[4,4],tab[4,1],tab[5,1]),NA,1-pf(tab[6,4],tab[6,1],tab[7,1]),NA,NA)
tab[8,3]<-NA

cv1=sqrt(as.numeric(tab[3,3]))/mean(resp)*100
cv2=sqrt(as.numeric(tab[5,3]))/mean(resp)*100
cv3=sqrt(as.numeric(tab[7,3]))/mean(resp)*100
tab<-round(tab,6)

cat('------------------------------------------------------------------------
Analysis of variance table\n------------------------------------------------------------------------\n')
print(tab)
cat('------------------------------------------------------------------------
CV 1 =',cv1,'%\nCV 2 =', cv2,'%\nCV 3 =', cv3,'%\n')

fatores<-data.frame('fator 1' = factor1,'fator 2' = factor2)

###############################################################################################################
#Teste de normalidade
pvalor.shapiro<-shapiro.test(anava$residuals)$p.value
cat('\n------------------------------------------------------------------------\nShapiro-Wilk normality test (Error b)\n')
cat('p-value: ',pvalor.shapiro, '\n')
if(pvalor.shapiro<0.05){cat('WARNING: at 5% of significance, residuals can not be considered normal!
                            ------------------------------------------------------------------------\n')}
else{cat('According to Shapiro-Wilk normality test at 5% of significance, residuals can be considered normal.
------------------------------------------------------------------------\n')}

# Creating unfold #########################################
if(is.null(unfold)){
  if(tab[6,5]>sigF)  {unfold<-c(unfold,1)}
  if(tab[6,5]<=sigF) {unfold<-c(unfold,2)}
}

#Para interacao nao significativa, fazer...
if(any(unfold==1)) {
cat('\nNo significant interaction: analyzing the main effects
------------------------------------------------------------------------\n')

for(i in 1:2){

#Para os fatores QUALITATIVOS, teste de medias
if(quali[i]==TRUE && as.numeric(tab[cont[i],5])<=sigF) {
    cat(fac.names[i])

  if(mcomp=='tukey'){
    tukey(resp,fatores[,i],as.numeric(tab[cont[i]+1,1]), as.numeric(tab[cont[i]+1,2]),sigT)
                    }
  if(mcomp=='duncan'){
    duncan(resp,fatores[,i],as.numeric(tab[cont[i]+1,1]), as.numeric(tab[cont[i]+1,2]),sigT)
                    }
  if(mcomp=='lsd'){
    lsd(resp,fatores[,i],as.numeric(tab[cont[i]+1,1]), as.numeric(tab[cont[i]+1,2]),sigT)
                    }
  if(mcomp=='lsdb'){
    lsdb(resp,fatores[,i],as.numeric(tab[cont[i]+1,1]), as.numeric(tab[cont[i]+1,2]),sigT)
                    }
  if(mcomp=='sk'){
    scottknott(resp,fatores[,i],as.numeric(tab[cont[i]+1,1]), as.numeric(tab[cont[i]+1,2]),sigT)
                    }
  if(mcomp=='snk'){
    snk(resp,fatores[,i],as.numeric(tab[cont[i]+1,1]), as.numeric(tab[cont[i]+1,2]),sigT)
                  }
  if(mcomp=="ccboot"){
  ccboot(resp,fatores[,i],as.numeric(tab[cont[i]+1,1]), as.numeric(tab[cont[i]+1,2]),sigT)
                     }
  if(mcomp=="ccF"){
  ccF(resp,fatores[,i],as.numeric(tab[cont[i]+1,1]), as.numeric(tab[cont[i]+1,2]),sigT)
                     }
          }

if(quali[i]==TRUE && as.numeric(tab[cont[i],5]>sigF)) {
    cat(fac.names[i])
    cat('\nAccording to F test, the means of this factor are not different.\n')
    cat('------------------------------------------------------------------------\n')
    mean.table<-tapply.stat(resp,fatores[,i],mean)
    colnames(mean.table)<-c('Levels','Means')
    print(mean.table)
    cat('------------------------------------------------------------------------')
}

#Para os fatores QUANTITATIVOS, regressao
if(quali[i]==FALSE && as.numeric(tab[cont[i],5])<=sigF){
    cat(fac.names[i])
    reg.poly(resp, fatores[,i], tab[cont[i]+1,1], as.numeric(tab[cont[i]+1,2]), as.numeric(tab[cont[i],1]),
             as.numeric(tab[cont[i],2]))
}

if(quali[i]==FALSE && as.numeric(tab[cont[i],5])>sigF) {
    cat(fac.names[i])
    cat('\nAccording to the F test, the means of this factor are not different.\n')
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

#Desdobramento de FATOR 1 dentro dos niveis de FATOR 2
cat("\nAnalyzing ", fac.names[1], ' inside of each level of ', fac.names[2], '
------------------------------------------------------------------------\n')

#Somas de quadrados do fator 1 dentro dos niveis de fator 2
l2<-names(summary(Fator2))

sq<-numeric(0)

for(k in 1:nv2) {
soma<-numeric(0)
for(j in 1:nv1) {
sub<-resp[Fator1==levels(Fator1)[j] & Fator2==levels(Fator2)[k]]
q.som<-length(sub)
soma<-c(soma, sum(sub))
                 }
sq<-c(sq, sum(soma^2)/q.som - sum(soma)^2/(q.som*length(soma)))
                 }
gl.sattert<-(as.numeric(tab[3,3])+(nv2-1)*as.numeric(tab[7,3]))^2/((as.numeric(tab[3,3])^2/as.numeric(tab[3,1]))
             + (((nv2-1)*as.numeric(tab[7,3]))^2/as.numeric(tab[7,1])))
gl.f1f2<-c(rep(nv1-1,nv2),gl.sattert)
sq<-c(sq, NA)
qm.f1f2<-sq[1:nv2]/gl.f1f2[1:nv2]
qm.ecomb<-(as.numeric(tab[3,3])+(nv2-1)*as.numeric(tab[7,3]))/nv2
qm.f1f2<-c(qm.f1f2,qm.ecomb)
fc.f1f2<-c(qm.f1f2[1:nv2]/qm.f1f2[nv2+1],NA)
p.f1f2<-c(1-pf(fc.f1f2,gl.f1f2,gl.sattert))
tab.f1f2<-data.frame('DF'=gl.f1f2,'SS'=sq,'MS'=qm.f1f2,'Fc'=fc.f1f2, 'p-value'=p.f1f2)
nome.f1f2<-numeric(0)
for(j in 1:nv2){
nome.f1f2<-c(nome.f1f2, paste(fac.names[1], ' : ', fac.names[2],' ',l2[j],' ',sep=''))
                }
nome.f1f2<-c(nome.f1f2,'Pulled error')
rownames(tab.f1f2)<-nome.f1f2
tab.f1f2<-round(tab.f1f2,6)
tab.f1f2[nv2+1,2]<-tab.f1f2[nv2+1,3]*tab.f1f2[nv2+1,1]
tab.f1f2[nv2+1,5]<-tab.f1f2[nv2+1,4]<-''
print(tab.f1f2)
    cat('------------------------------------------------------------------------\n\n')

for(i in 1:nv2) {

    cat('\n',fac.names[1], 'inside', fac.names[2], l2[i] )
    cat('\n------------------------------------------------------------------------')

  if(quali[1]==TRUE & as.numeric(tab.f1f2[i,5])<=sigF) {

    if(mcomp=='tukey'){
    tukey(resp[fatores[,2]==l2[i]], fatores[,1][fatores[,2]==l2[i]], as.numeric(tab.f1f2[nv2+1,1]),as.numeric(tab.f1f2[nv2+1,2]), sigT)
                      }

  if(mcomp=='duncan'){
    duncan(resp[fatores[,2]==l2[i]],fatores[,1][fatores[,2]==l2[i]],as.numeric(tab.f1f2[nv2+1,1]),as.numeric(tab.f1f2[nv2+1,2]),sigT)
                    }

  if(mcomp=='lsd'){
    lsd(resp[fatores[,2]==l2[i]],fatores[,1][fatores[,2]==l2[i]],as.numeric(tab.f1f2[nv2+1,1]),as.numeric(tab.f1f2[nv2+1,2]),sigT)
                    }

  if(mcomp=='lsdb'){
    lsdb(resp[fatores[,2]==l2[i]],fatores[,1][fatores[,2]==l2[i]],as.numeric(tab.f1f2[nv2+1,1]),as.numeric(tab.f1f2[nv2+1,2]),sigT)
                    }

  if(mcomp=='sk'){
    scottknott(resp[fatores[,2]==l2[i]],fatores[,1][fatores[,2]==l2[i]],as.numeric(tab.f1f2[nv2+1,1]),as.numeric(tab.f1f2[nv2+1,2]),sigT)
                    }

  if(mcomp=='snk'){
    snk(resp[fatores[,2]==l2[i]],fatores[,1][fatores[,2]==l2[i]],as.numeric(tab.f1f2[nv2+1,1]),as.numeric(tab.f1f2[nv2+1,2]),sigT)
                   }
  if(mcomp=="ccboot"){
  ccboot(resp[fatores[,2]==l2[i]],fatores[,1][fatores[,2]==l2[i]],as.numeric(tab.f1f2[nv2+1,1]),as.numeric(tab.f1f2[nv2+1,2]),sigT)
                     }
  if(mcomp=="ccF"){
  ccF(resp[fatores[,2]==l2[i]],fatores[,1][fatores[,2]==l2[i]],as.numeric(tab.f1f2[nv2+1,1]),as.numeric(tab.f1f2[nv2+1,2]),sigT)
                     }
                                                   }

if(quali[1]==FALSE & as.numeric(tab.f1f2[i,5])<sigF) {             #Fazer regressao
    reg.poly(resp[fatores[,2]==l2[i]], fatores[,1][fatores[,2]==l2[i]], as.numeric(tab.f1f2[nv2+1,1]),
    as.numeric(tab.f1f2[nv2+1,2]), as.numeric(tab.f1f2[i,1]), as.numeric(tab.f1f2[i,2]))
                                                   }

if(as.numeric(tab.f1f2[i,5])>sigF) {
    cat('\nAccording to F test, the means of this factor are not distinct.\n')
    cat('------------------------------------------------------------------------\n')
    mean.table<-tapply.stat(resp[fatores[,2]==l2[i]],fatores[,1][fatores[,2]==l2[i]],mean)
    colnames(mean.table)<-c('Levels','Means')
    print(mean.table)
    cat('------------------------------------------------------------------------\n')
                                   }
                        }


#Desdobramento de FATOR 2 dentro dos niveis de FATOR 1
cat("\n\nAnalysing ", fac.names[2], ' inside each level of ', fac.names[1], '
------------------------------------------------------------------------\n')

#Somas de quadrados do fator 2 dentro dos niveis de fator 1
l1<-names(summary(Fator1))

sq<-numeric(0)

for(k in 1:nv1) {
soma<-numeric(0)
for(j in 1:nv2) {
parc<-resp[Fator1==levels(Fator1)[k] & Fator2==levels(Fator2)[j]]
q.som<-length(parc)
soma<-c(soma, sum(parc))
                 }
sq<-c(sq, sum(soma^2)/q.som - sum(soma)^2/(q.som*length(soma)))
                 }
gl.sattert<-(as.numeric(tab[5,3])+(nv2-1)*as.numeric(tab[7,3]))^2/((as.numeric(tab[5,3])^2/as.numeric(tab[5,1]))
            + (((nv2-1)*as.numeric(tab[7,3]))^2/as.numeric(tab[7,1])))
gl.f2f1<-c(rep(nv2-1,nv1),gl.sattert)
sq<-c(sq, NA)
qm.f2f1<-sq[1:nv1]/gl.f2f1[1:nv1]
qm.ecomb<-(as.numeric(tab[5,3])+(nv1-1)*as.numeric(tab[7,3]))/nv1
qm.f2f1<-c(qm.f2f1,qm.ecomb)
fc.f2f1<-c(qm.f2f1[1:nv1]/qm.f2f1[nv1+1],NA)
p.f2f1<-c(1-pf(fc.f2f1,gl.f2f1,gl.sattert))
tab.f2f1<-data.frame('DF'=gl.f2f1,'SS'=sq,'MS'=qm.f2f1,'Fc'=fc.f2f1, 'p-value'=p.f2f1)
nome.f2f1<-numeric(0)
for(j in 1:nv1){
nome.f2f1<-c(nome.f2f1, paste(fac.names[2], ' : ', fac.names[1],' ',l1[j],' ',sep=''))
                }
nome.f2f1<-c(nome.f2f1,'Pulled error')
rownames(tab.f2f1)<-nome.f2f1
tab.f2f1<-round(tab.f2f1,6)
tab.f2f1[nv1+1,2]<-tab.f2f1[nv1+1,3]*tab.f2f1[nv1+1,1]
tab.f2f1[nv1+1,5]<-tab.f2f1[nv1+1,4]<-''
print(tab.f2f1)
    cat('------------------------------------------------------------------------\n\n')


for(i in 1:nv1) {

    cat('\n',fac.names[2], 'inside', fac.names[1], l1[i] )
    cat('\n------------------------------------------------------------------------')


  if(quali[2]==TRUE & as.numeric(tab.f2f1[i,5])<sigF) {             #Fazer teste de comparacao multipla

    if(mcomp=='tukey'){
    tukey(resp[fatores[,1]==l1[i]], fatores[,2][fatores[,1]==l1[i]], as.numeric(tab.f2f1[nv1+1,1]),as.numeric(tab.f2f1[nv1+1,2]),sigT)
                    }

  if(mcomp=='duncan'){
    duncan(resp[fatores[,1]==l1[i]],fatores[,2][fatores[,1]==l1[i]],as.numeric(tab.f2f1[nv1+1,1]),as.numeric(tab.f2f1[nv1+1,2]),sigT)
                    }

  if(mcomp=='lsd'){
    lsd(resp[fatores[,1]==l1[i]],fatores[,2][fatores[,1]==l1[i]],as.numeric(tab.f2f1[nv1+1,1]),as.numeric(tab.f2f1[nv1+1,2]),sigT)
                    }

  if(mcomp=='lsdb'){
    lsdb(resp[fatores[,1]==l1[i]],fatores[,2][fatores[,1]==l1[i]],as.numeric(tab.f2f1[nv1+1,1]),as.numeric(tab.f2f1[nv1+1,2]),sigT)
                    }

  if(mcomp=='sk'){
    scottknott(resp[fatores[,1]==l1[i]],fatores[,2][fatores[,1]==l1[i]],as.numeric(tab.f2f1[nv1+1,1]),as.numeric(tab.f2f1[nv1+1,2]),sigT)
                    }

  if(mcomp=='snk'){
    snk(resp[fatores[,1]==l1[i]],fatores[,2][fatores[,1]==l1[i]],as.numeric(tab.f2f1[nv1+1,1]),as.numeric(tab.f2f1[nv1+1,2]),sigT)
                  }
  if(mcomp=="ccboot"){
  ccboot(resp[fatores[,1]==l1[i]],fatores[,2][fatores[,1]==l1[i]],as.numeric(tab.f2f1[nv1+1,1]),as.numeric(tab.f2f1[nv1+1,2]),sigT)
                     }
  if(mcomp=="ccF"){
  ccF(resp[fatores[,1]==l1[i]],fatores[,2][fatores[,1]==l1[i]],as.numeric(tab.f2f1[nv1+1,1]),as.numeric(tab.f2f1[nv1+1,2]),sigT)
                     }
    cat('------------------------------------------------------------------------\n\n')
                                                      }


  if(quali[2]==FALSE & as.numeric(tab.f2f1[i,5])<sigF){            #Fazer regressao
    reg.poly(resp[fatores[,1]==l1[i]], fatores[,2][fatores[,1]==l1[i]], as.numeric(tab.f2f1[nv1+1,1]),
    as.numeric(tab.f2f1[nv1+1,2]), as.numeric(tab.f2f1[i,1]), as.numeric(tab.f2f1[i,2]))
                                                   }


if(as.numeric(tab.f2f1[i,5])>sigF) {
    cat('\nAccording to F test, the means of this factor are not distinct.\n')
    cat('------------------------------------------------------------------------\n')
    mean.table<-tapply.stat(resp[fatores[,1]==l1[i]],fatores[,2][fatores[,1]==l1[i]],mean)
    colnames(mean.table)<-c('Levels','Means')
    print(mean.table)
    cat('------------------------------------------------------------------------\n')
                                                }

}

}
## error a ##
tabmedia<-model.tables(anava, "means")
#error.plot<-as.vector(t(as.matrix(tabmedia$tables$`Fator1:bloco`)-as.vector(tabmedia$tables$Fator1))-as.vector(tabmedia$tables$bloco))
#Saida
out<-list()
out$residuals<-anava$residuals
#out$residuals.a<-error.plot
out$df.residual<-anava$df.residual
#out$df.residual.a<-as.numeric(tab[2,1])
out$coefficients<-anava$coefficients
out$effects<-anava$effects
out$fitted.values<-anava$fitted.values
out$means.factor1<-tapply.stat(resp,fatores[,1],mean)
out$means.factor2<-tapply.stat(resp,fatores[,2],mean)
out$means.inside<-tabmedia$tables$`Fator1:Fator2`
#if(quali==FALSE && tab[[1]][1,5]<sigF) {out$reg<-reg}
invisible(out)
}
