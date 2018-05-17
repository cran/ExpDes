samiuddin <-
function(trat, resp, t, r)
{
  Trat<-factor(trat)
  t=length(levels(Trat))
  m<-matrix(0,t,1)
  somma<-matrix(0,t,1)
  a2<-matrix(0,t,1)
  rp<-0
  for(i in 1:t) {
    dife<-0  
    soma<-0
    for(j in 1:r[i]) {
      dife<-(resp[rp+j]-mean(resp[(rp+1):(rp+r[i])]))^2
      soma<-soma+dife
    }
    somma[i]<-soma
    rp<-sum(r[1:i])
    m[i]<-(((r[i]-1)/somma[i])^(1/3))*(1-(2/(9*(r[i]-1))))
    a2[i]<-2/(9*(somma[i]^(2/3))*(r[i]-1)^(1/3))
  }
  mm<-sum(m/a2)/sum(1/a2)
  pvalor<-pchisq(sum(((m-mm)^2)/a2), (t-1))
  output <- pvalor
  return(output)
}
