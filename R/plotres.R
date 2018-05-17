plotres <-function(x){
  resid<-x$residuals
  df.resid<-x$df.residual
  fitted.val<-x$fitted.values
  var.res<-sum(resid^2)/df.resid
  respad<-resid/sqrt(var.res)
  par(mfrow=c(2,2))
  # Grafico1
  hist(respad, xlab="Standardized Residuals", main="Histogram", freq=FALSE)
  x<-c()
  curve(dnorm(x,mean=0,sd=1),col=2,lty=1,lwd=1,add=TRUE)
  #plot(fitted.val, resid, xlab="Fitted Values", ylab="Residuals")
  #abline(h=0, col = "lightgray", lty = 3)
  #title("Residuals vs Fitted Values")
  #
  # Grafico2 P-Pplot com as bandas de confiança https://stats.stackexchange.com/questions/111288/confidence-bands-for-qq-line
  good<-!is.na(resid)
  ord<-order(resid[good])
  ord.x<-resid[good][ord]
  n<-length(ord.x)
  P<-ppoints(n)
  z<-qnorm(P)
  plot(z,ord.x, xlab="z", ylab="Residuals")
  coef<-coef(lm(ord.x~z)) #rlm
  b0<-coef[1]
  b1<-coef[2]
  abline(b0,b1,col="red",lwd=2)
  conf<-0.95
  zz<-qnorm(1-(1-conf)/2)
  SE<-(b1/dnorm(z))*sqrt(P*(1-P)/n)     
  fit.value<-b0+b1*z
  upper<-fit.value+zz*SE
  lower<-fit.value-zz*SE
  lines(z,upper,lty=2,lwd=2,col="red")
  lines(z,lower,lty=2,lwd=2,col="red")
  title("Normal Q-Q (95%)")
  #
  # Grafico3 Resíduos standardizados versus Valores Ajustados da Variável Resposta
  plot(fitted.val, respad, xlab="Fitted Values", ylab="Standardized Residuals")
  abline(h=0,  col = "red", lty = 3)
  title("Standardized Residuals vs Fitted Values")
  #
  # Grafico4 
  boxplot(respad)
  title("Standardized Residuals")
}
