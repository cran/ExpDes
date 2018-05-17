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
