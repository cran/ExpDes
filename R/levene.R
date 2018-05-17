levene <-
function(trat, resp, t, r)
{
  Trat<-factor(trat)
   zdados1<-matrix(0,length(resp),1)
   rp<-0
   for(k in 1:length(resp)) {
    zdados1[k]<-abs(resp[k]-mean(resp[(rp+1):(rp+r[Trat[k]])]))
    if(k<length(resp)){if(trat[k]!=trat[k+1]){rp<-sum(r[1:Trat[k]])}}
     }
   pvalor<-summary(aov(zdados1 ~ Trat))[[1]][1,5]
   output <- pvalor
   return(output)
}
