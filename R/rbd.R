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
tab[[1]][4,3]=' '
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
  if(mcomp=="ccf")   ccf(resp,Trat,tab[[1]][3,1],tab[[1]][3,2],sigT)
  
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
