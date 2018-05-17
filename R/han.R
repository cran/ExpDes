han <-
function(resp, trat, block)
{
  Trat<-length(levels(trat)) # numero de tratamentos (Trat deve ser factor)
  Block<-length(levels(block)) # numero de blocos    (Bloco deve ser factor)
  if (Block>Trat){   # se numero de blocos for menor ou igual ao numero de tratamentos não rodar.
   dife<-matrix(0,Block,Trat)
   ymedia<-matrix(0,Block,1)
   rp<-0
   for(j in 1:Block) {
    for(i in 1:Trat) {
      ymedia[j]<-mean(resp[(rp+1):(rp+Trat)])
      dife[j,i]<-(resp[rp+i]-ymedia[j])
      }
    rp<-Trat*j
    }
   modelohan<-lm(ymedia ~ dife[,2:Trat])
   pvalor.hvar<-1-pf(summary(modelohan)[[10]][1],summary(modelohan)[[10]][2],summary(modelohan)[[10]][3])
   output <- pvalor.hvar
   return(output)
   }
}
