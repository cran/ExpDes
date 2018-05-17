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
tab[[1]][3,3]=' '
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
  if(mcomp=="ccf")   ccf(resp,Trat,tab[[1]][2,1],tab[[1]][2,2],sigT)
  
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
