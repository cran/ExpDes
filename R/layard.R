layard <-
function(trat, resp, t, r)
{
 vari<-matrix(0,t,1)
 varij<-matrix(0,max(r),t)
 U<-matrix(0,max(r),t)
 rp<-0
 for(i in 1:t) {
   vv<-resp[(rp+1):(rp+r[i])]
   vari[i]<-var(vv)
   for(j in 1:r[i]) {
     varij[j,i]<-var(vv[-j])
     U[j,i]<-(r[i]*log(vari[i]))-((r[i]-1)*log(varij[j,i]))
    }
   rp<-sum(r[1:i])
 }
 Uij<-as.vector(U)
 dadosUij<-cbind(trat,Uij)
 dadosUij<-as.data.frame(dadosUij)
 pvalor<-summary(aov(dadosUij$Uij ~ trat))[[1]][1,5]
 output <- pvalor
 return(output)
}
