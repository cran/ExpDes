#' Residual plots
#'
#' \code{plotres} Residual plots for a output model. Four sets
#' of plots are produced: (1) Histogram, (2) normal probability
#' plot for the residual, (3) Standardized Residuals versus
#' Fitted Values, and (4) box-plot (Standardized Residuals).
#' @param x Output from anova (performed in ExpDes).
#' @references STEEL, R. G. D.; TORRIE, J. H. \emph{Principles
#' and procedures in Statistics: a biometrical approach}.
#' McGraw-Hill, New York, NY. 1980.
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#'  @author Denismar Alves Nogueira
#'  @note The default produces four plots regarding the ANOVA
#'  assumptions.
#' @seealso \code{\link{graphics}}.
#' @examples
#' data(ex1)
#' attach(ex1)
#' a<-crd(trat, ig)
#' plotres(a)
#' @export

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
  # Grafic3 
  plot(fitted.val, respad, xlab="Fitted Values", ylab="Standardized Residuals")
  abline(h=0,col="red",lty = 3)
  title("Standardized Residuals vs Fitted Values")
  # Grafic4
  boxplot(respad)
  title("Standardized Residuals")
  par(mfrow=c(1,1))
}
