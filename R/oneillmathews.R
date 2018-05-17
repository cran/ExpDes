oneillmathews <-
function(trat, resp, t, r)
{
  Trat<-factor(trat)
  zdados1.1<-matrix(0,length(resp),1)
  rr<-t/sum(1/r)
  rp<-0
   for(k in 1:length(resp)) {
     zdados1.1[k]<-abs(resp[k]-mean(resp[(rp+1):(rp+r[Trat[k]])]))/sqrt(1-(1/rr))
     if(k<length(resp)){if(trat[k]<trat[k+1]){rp<-sum(r[1:Trat[k]])}}
    }
  Fc5.1<-summary(aov(zdados1.1 ~ trat))[[1]][1,4] # pvalor = posicao [1,5]
  b<-(1-2/pi)
  c<-(2/pi)*(1/(rr-1))*(sqrt(rr*(rr-2))+asin(1/(rr-1))-(rr-1))
  m<-(b-c)/(b+(rr-1)*c)
  Fc13<-m*Fc5.1
  pvalor<-(1-pf(Fc13, (t-1), summary(aov(zdados1.1 ~ trat))[[1]][2,1]))   
  output <- pvalor
  return(output)
}
