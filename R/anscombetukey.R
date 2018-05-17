anscombetukey <-
function(resp, trat, block, glres, msres, sstrat, ssblock, residuals, fitted.values)
{
 Trat<-length(trat)
 Bloco<-length(block)
 div1<-(2*glres*(msres)^2)/glres+2
 div2<-((((Trat-2)*(Bloco-1))/Trat*Bloco)*sstrat)+((((Trat-1)*(Bloco-2))/Trat*Bloco)*ssblock)
 Fc17<-((sum((residuals^2)*(fitted.values-mean(resp))))^2)/(div1*div2)
 pvalor.hvar<-1-pf(Fc17,1,((Trat-1)*(Bloco-1)))
 output <- pvalor.hvar
 return(output)
}
