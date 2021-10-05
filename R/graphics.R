#' Regression model plots
#'
#' \code{graphics} Plots from regression models fitted in ANOVA.
#' @param a Output from anova (performed in ExpDes).
#' @param degree For polynomial models, 1 (linear model) is the
#' \emph{default}, 2 (quadratic model), 3 (cubic model), "pot"
#' (Power model), "log" (Logistic model), "gom" (Gompertz model)
#' and "exp" (Exponential model).
#' @param mod Logic. Print the model expression and its R2 on
#' the top of the graphic. The \emph{default} is TRUE.
#' @param main Title of the plot. Empty is the \emph{default}.
#' @param sub Subtitle of the plot. Empty is the \emph{default}.
#' @param xlab Name for axis X.
#' @param ylab Name for axis Y.
#' @param pch Caracter type to be used on the observed values.
#' @param xlim Limits for axis X.
#' @param ylim Limits for axis Y.
#' @param bty Type of box the plot is fitted in.
#' @references STEEL, R. G. D.; TORRIE, J. H. \emph{Principles
#' and procedures in Statistics: a biometrical approach}.
#' McGraw-Hill, New York, NY. 1980.
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#' @seealso \code{\link{reg.poly}}, \code{\link{plotres}}.
#' @examples
#' data(ex1)
#' attach(ex1)
#' a<-crd(trat, ig, quali=FALSE, nl=FALSE)
#' graphics(a, degree=1)
#' graphics(a, degree=2)
#' graphics(a, degree=3)
#' @importFrom "grDevices" "dev.new"
#' @importFrom "graphics" "abline" "boxplot" "curve" "hist"
#' "lines" "mtext" "par" "points" "title"
#' @export

graphics <-
function(a, degree = 1, mod = TRUE, main = " ", sub = " ",
         xlab = "Levels (X)", ylab = "Response var (Y)", pch = 19,
         xlim = NULL, ylim = NULL, bty = "o"){

a<-a$reg
xob<-as.numeric(as.vector(a$'Table of means'[,1]))
x<-seq(min(xob), max(xob), by=0.1)

if(degree==1) {
dev.new()
b0<-a$'Coefficients linear reg'[1]; b1<-a$'Coefficients linear reg'[2]
y<-b0 + b1*x
yob<-as.numeric(as.vector(a$'Table of means'[,2]))
if(is.null(ylim)==TRUE) ylim=c(min(y,yob), max(y,yob))
plot(x,y,'l', main=main, sub=sub, bty=bty, xlab=xlab, ylab=ylab,
     xlim=xlim, ylim=ylim)
if(mod==TRUE) mtext(paste('y =',round(b0,3),'+',round(b1,3),
     'x  ', ' R^2 = ', round(a$'R2 linear reg'*100,2),'%'),side=3)
points(xob, yob, pch=pch)
             }

if(degree==2) {
dev.new()
b0<-a$'Coefficients quadratic reg'[1]; b1<-a$'Coefficients quadratic reg'[2]
b2<-a$'Coefficients quadratic reg'[3]
y<-b0 + b1*x + b2*x^2
yob<-as.numeric(as.vector(a$'Table of means'[,2]))
if(is.null(ylim)==TRUE) ylim=c(min(y,yob), max(y,yob))
plot(x,y,'l', main=main, sub=sub, bty=bty, xlab=xlab, ylab=ylab,
     xlim=xlim, ylim=ylim)
if(mod==TRUE) mtext(paste('y = ',round(b0,3),'+',round(b1,3),'x+',
     round(b2,3),'x^2  ',' R^2 = ', round(a$'R2 quadratic reg'*100,2),'%'),
     side=3)
points(xob, yob, pch=pch)
             }

if(degree==3) {
dev.new()
b0<-a$'Coefficients cubic reg'[1]; b1<-a$'Coefficients cubic reg'[2]
b2<-a$'Coefficients cubic reg'[3]; b3<-a$'Coefficients cubic reg'[4]
y<-b0 + b1*x + b2*x^2 + b3*x^3
yob<-as.numeric(as.vector(a$'Table of means'[,2]))
if(is.null(ylim)==TRUE) ylim=c(min(y,yob), max(y,yob))
plot(x,y,'l', main=main, sub=sub, bty=bty, xlab=xlab, ylab=ylab,
     xlim=xlim, ylim=ylim)
if(mod==TRUE) mtext(paste('y = ',round(b0,3),'+',round(b1,3),'x+',
     round(b2,3),'x^2+',round(b3,3),'x^3  ',' R^2 = ',
     round(a$'R2 cubic reg'*100,2),'%'),side=3)
points(xob, yob, pch=pch)

}

#reg nlinear

if(degree=='pow' && is.numeric(a$'Power model coefficients')) {
  if(is.numeric(a$'Power model coefficients')==FALSE) { print(a$'power model AIC') } # estava por fora do if
  dev.new()
  b0<-a$'Power model coefficients'[1]; b1<-a$'Power model coefficients'[2]
  y<-b0*x^b1
  yob<-as.numeric(as.vector(a$'Table of means'[,2]))
  if(is.null(ylim)==TRUE) ylim=c(min(y,yob), max(y,yob))
  plot(x,y,'l', main=main, sub=sub, bty=bty, xlab=xlab, ylab=ylab,
       xlim=xlim, ylim=ylim)
  if(mod==TRUE) mtext(paste('y = ',round(b0,3),'* x^',round(b1,3),' Approx R2 = ',
                            round(a$'Approximate R2 of power Model'*100,2),'%'),side=3)
  points(xob, yob, pch=pch)
}

if(degree=='exp' && is.numeric(a$'Exponential model coefficients')) {
  if(is.numeric(a$'Exponential model coefficients')==FALSE) { print(a$'exponential model AIC') } # estava por fora do if
  dev.new()
  b0<-a$'Exponential model coefficients'[1]; b1<-a$'Exponential model coefficients'[2]
  y<-b0*exp(b1*x)
  yob<-as.numeric(as.vector(a$'Table of means'[,2]))
  if(is.null(ylim)==TRUE) ylim=c(min(y,yob), max(y,yob))
  plot(x,y,'l', main=main, sub=sub, bty=bty, xlab=xlab, ylab=ylab,
       xlim=xlim, ylim=ylim)
  if(mod==TRUE) mtext(paste('y = ',round(b0,3),'* exp(',round(b1,3),'x)   Approx R2 = ',
                            round(a$'Approximate R2 of exponential Model'*100,2),'%'),side=3)
  points(xob, yob, pch=pch)
}

if(degree=='log' && is.numeric(a$'Logistic model coefficients')) {
  if(is.numeric(a$'Logistic model coefficients')==FALSE) { print(a$'logistic model AIC') } # estava por fora do if
  dev.new()
  b0<-a$'Logistic model coefficients'[1]
  b1<-a$'Logistic model coefficients'[2]
  b2<-a$'Logistic model coefficients'[3]
  y<-b0/(1+exp(b1-(b2*x)))
  yob<-as.numeric(as.vector(a$'Table of means'[,2]))
  if(is.null(ylim)==TRUE) ylim=c(min(y,yob), max(y,yob))
  plot(x,y,'l', main=main, sub=sub, bty=bty, xlab=xlab, ylab=ylab,
       xlim=xlim, ylim=ylim)
  if(mod==TRUE) mtext(paste('y = ',
                            round(b0,3),'/(1+exp(',round(b1,3),'-(',round(b2,3),'*x)))   Approx R2 = ',
                            round(a$'Approximate R2 of logistic Model'*100,2),'%'),side=3)
  points(xob, yob, pch=pch)
}

if(degree=='gom' && is.numeric(a$'Gompertz model coefficients')) {
  if(is.numeric(a$'Gompertz model coefficients')==FALSE) { print(a$'Gompertz model AIC') } # estava por fora do if
  dev.new()
  b0<-a$'Gompertz model coefficients'[1]
  b1<-a$'Gompertz model coefficients'[2]
  b2<-a$'Gompertz model coefficients'[3]
  y<-b0*exp(-exp(b1-(b2*x)))
  yob<-as.numeric(as.vector(a$'Table of means'[,2]))
  if(is.null(ylim)==TRUE) ylim=c(min(y,yob), max(y,yob))
  plot(x,y,'l', main=main, sub=sub, bty=bty, xlab=xlab, ylab=ylab,
       xlim=xlim, ylim=ylim)
  if(mod==TRUE) mtext(paste('y = ',
                            round(b0,3),'*exp(-exp(',round(b1,3),'-(',round(b2,3),'*x)))   Approx R2 = ',
                            round(a$'Approximate R2 of Gompertz Model'*100,2),'%'),side=3)
  points(xob, yob, pch=pch)
}

}
