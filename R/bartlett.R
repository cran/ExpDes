bartlett <-
function(trat, resp, t, r)
{
  vari<-matrix(0,t,1)
  rp<-0
  for(i in 1:t) {
    vari[i]<-var(resp[(rp+1):(rp+r[i])])
    rp<-sum(r[1:i])
  }
  S2p<-sum((r-1)*vari)/(length(resp)-t)
  A<-(length(resp)-t)*log(S2p)-sum((r-1)*log(vari))
  B<-(1/(3*(t-1)))*(sum(1/(r-1))-(1/(length(resp)-t)))
  Xc1<-A/(1+B)
  pvalor<-1-pchisq(Xc1, t-1)
  output <- pvalor
  return(output)
}
