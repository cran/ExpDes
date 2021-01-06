#' Non-linear Regression
#'
#' \code{reg.nl} Adjusts non-linear regression models in Anova
#' (Models: Power, Exponential, Logistic, Gompertz).
#' @param treat Numeric or complex vector containing the
#' treatments.
#' @param resp Numeric or complex vector containing the
#' response variable.
#' @return Returns coefficients, significance and ANOVA of the
#' fitted regression models.
#' @references DRAPER, N.R.; SMITH, H. \emph{Apllied regression
#' analysis}. 3ed. New York : John Wiley, 1998. 706p.
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#' @author Luiz Alberto Beijo
#' @seealso \code{\link{graphics}}.
#' @examples
#' data(exnl)
#' attach(exnl)
#' x<-crd(trat, resp, quali = FALSE, nl = TRUE)
#' graphics(x, degree = "log")
#' @importFrom "stargazer" "stargazer"
#' @export

reg.nl<-function(resp, treat) {

  ###========================================================###
  #        -  Mod Potencia
  ###========================================================###
  #Obtencao de estimativas iniciais
  #Protecao para o dominio da funcao log
  if(any(resp==0)) {resp1<-resp+1e-10; yli=log(resp1)}
  if(any(resp==0)==FALSE) yli=log(resp)
  if(any(treat==0)) {treat1<-treat+1e-10; xli=log(treat1)}
  if(any(treat==0)==FALSE) xli=log(treat)

  modli=lm(yli~xli)

  apinit=exp(coef(modli)[1]) # Estimativa inicial do parametro beta
  bpinit=(coef(modli)[2])    # Estimativa inicial do parametro gama

  api<-(apinit)[[1]][1:1]
  bpi<-(bpinit)[[1]][1:1]

  mod_pot<-numeric(0)
  try( mod_pot <- nls(resp~alfa*treat^beta, start=c(alfa=api, beta=bpi)),
       silent=TRUE)
  cat('------------------------------------------------------------------------\n')

  if(length(mod_pot)==0) {cat('
------------------------------------------------------------------------
Power Model
------------------------------------------------------------------------
Error in model fit! (singularity, convergence, etc)
')

tm1<-data.frame('Error in model fit!')
aic1 <- rs1 <- c('Error in model fit!')

  }

  if(length(mod_pot)!=0) {

    b1<-summary(mod_pot)
    tm1<-data.frame('Estimate' = round(b1[[10]][,1],8),
                    'Standard Error' = round(b1[[10]][,2],5),
                    'tc'=round(b1[[10]][,3],5),
                    'p-value' = round(b1[[10]][,4],5))
    rownames(tm1)<-c('Alpha','Beta')

    #AIC
    aic1<-AIC(mod_pot)

    # R2 Aproximado
    f=fitted.values(mod_pot)
    r=residuals(mod_pot)
    dr=sum((r/resp)^2)
    sr=sum(r^2)
    qe=mean(resp)
    d=(f-qe)
    dq=d^2
    rs1=sum(dq)/(sum(dq)+sr)

    output1<-list('
------------------------------------------------------------------------\n
Power Model
------------------------------------------------------------------------\n' = tm1,
                  '
------------------------------------------------------------------------\n
AIC
------------------------------------------------------------------------\n' = aic1,
'
------------------------------------------------------------------------\n
Approximate R2 of Power Model
------------------------------------------------------------------------\n' = rs1)
    #print(output1,right=TRUE)
    stargazer(tm1, type = "text", title="Power Model", summary=F, digits=4)
    stargazer(aic1, type = "text", title="AIC", digits=4, style="apsr")
    stargazer(rs1, type = "text", title="Approximate R2 of Power Model", digits=4, flip=FALSE, summary=F, style="apsr")
    cat('------------------------------------------------------------------------\n')

  }

  ###========================================================###
  #        Modelo Exponencial
  ###========================================================###

  #Obtencao de estimativas iniciais
  yeli=log(resp)
  xeli=(treat)
  modlie=lm(yeli~xeli)

  aeinit=exp(coef(modlie)[1]) # Estimativa inicial do parametro beta
  beinit=(coef(modlie)[2])    # Estimativa inicial do parametro gama

  aei<-(aeinit)[[1]][1:1]
  bei<-(beinit)[[1]][1:1]

  mod_exp<-numeric(0)
  try( mod_exp<-nls(resp~alfa*exp(beta*treat), start=c(alfa=aei, beta=bei)),
       silent=TRUE)

  if(length(mod_exp)==0) {cat('
------------------------------------------------------------------------\n
Exponential Model
------------------------------------------------------------------------\n
Error in model fit! (singularity, convergence, etc)
                              ')

                          tm2<-data.frame('Error in model fit!')
                          aic2 <- rs2 <- c('Error in model fit!')

  }

  if(length(mod_exp)!=0) {

    b2<-summary(mod_exp)

    tm2<-data.frame('Estimate' = round(b2[[10]][,1],8),
                    'Standard Error' = round(b2[[10]][,2],5),
                    'tc'=round(b2[[10]][,3],5),
                    'p-value' = round(b2[[10]][,4],5))
    rownames(tm2)<-c('Alpha','Beta')

    aic2<-AIC(mod_exp)

    # R2 Aproximado
    f=fitted.values(mod_exp)
    r=residuals(mod_exp)
    dr=sum((r/resp)^2)
    sr=sum(r^2)
    qe=mean(resp)
    d=(f-qe)
    dq=d^2
    rs2=sum(dq)/(sum(dq)+sr)

    output2<-list('
------------------------------------------------------------------------\n
Exponential Model
------------------------------------------------------------------------\n' = tm2,
'
------------------------------------------------------------------------\n
AIC
------------------------------------------------------------------------\n' = aic2,
'
------------------------------------------------------------------------\n
Approximate R2 of Exponential Model
------------------------------------------------------------------------\n' = rs2)
    #print(output2,right=TRUE)
    stargazer(tm2, type = "text", title="Exponential Model", summary=F, digits=4)
    stargazer(aic2, type = "text", title="AIC", digits=4, style="apsr")
    stargazer(rs2, type = "text", title="Approximate R2 of Exponential Model", digits=4, flip=FALSE, summary=F, style="apsr")
    cat('------------------------------------------------------------------------\n')
  }

  ###========================================================###
  #                      Modelo logistico
  ###========================================================###

  #Obtencao de estimativas iniciais
  ylinit=log(((max(resp)+2)/resp)-1)
  modlil=lm(ylinit~treat)
  summary(modlil)

  b=coef(modlil)[1]
  gm=coef(modlil)[2]

  glinit=-1*gm             # Estimativa inicial do parametro gama
  alinit=max(resp)+2       # Estimativa inicial do parametro alfa
  blinit=(coef(modlil)[1]) # Estimativa inicial do parametro beta

  ali<-(alinit)[[1]][1:1]
  bli<-(blinit)[[1]][1:1]
  gli<-(glinit)[[1]][1:1]

  mod_logi<-numeric(0)
  try( mod_logi<-nls(resp~alfa/(1+exp(beta-(gama*treat))),
                     start=c(alfa=ali,beta=bli,gama=gli)),
       silent=TRUE)

  if(length(mod_logi)==0) {cat('
------------------------------------------------------------------------\n
Logistic Model
------------------------------------------------------------------------\n
Error in model fit! (singularity, convergence, etc)
                               ')

                           tm3<-data.frame('Error in model fit!')
                           aic3 <- rs3 <- c('Error in model fit!')

  }

  if(length(mod_logi)!=0) {

    b3<-summary(mod_logi)

    tm3<-data.frame('Estimate' = round(b3[[10]][,1],8),
                    'Standard Error' = round(b3[[10]][,2],5),
                    'tc'=round(b3[[10]][,3],5),
                    'p-value' = round(b3[[10]][,4],5))
    rownames(tm3)<-c('Alpha','Beta','Gamma')

    yest=fitted(mod_logi)
    aic3<-AIC(mod_logi)

    # Coeficiente de Determinacao Aproximado
    fl=fitted.values(mod_logi)
    rl=residuals(mod_logi)
    drl=sum((rl/resp)^2)
    srl=sum(rl^2)
    qel=mean(resp)
    dl=(fl-qel)
    dql=dl^2
    rs3=sum(dql)/(sum(dql)+srl)

    output3<-list('
------------------------------------------------------------------------\n
Logistic Model
------------------------------------------------------------------------\n' = tm3,
'
------------------------------------------------------------------------\n
AIC
------------------------------------------------------------------------\n' = aic3,
'
------------------------------------------------------------------------\n
Approximate R2 of Logistic Model
------------------------------------------------------------------------\n' = rs3)
    #print(output3,right=TRUE)
    stargazer(tm3, type = "text", title="Logistic Model", summary=F, digits=4)
    stargazer(aic3, type = "text", title="AIC", digits=4, style="apsr")
    stargazer(rs3, type = "text", title="Approximate R2 of Logistic Model", digits=4, flip=FALSE, summary=F, style="apsr")
    cat('------------------------------------------------------------------------\n')
  }


  ###========================================================###
  #         Modelo Gompertz
  ###========================================================###

  #Obtencao de estimativas iniciais
  ylig=log(-1*log(resp/(max(resp)+2)))
  modlig=lm(ylig~treat)

  bg=coef(modlig)[1]
  gmg=coef(modlig)[2]

  gginit=-1*gmg            # Estimativa inicial do parametro gama
  aginit=max(resp)+2       # Estimativa inicial do parametro alfa
  bginit=(coef(modlig)[1]) # Estimativa inicial do parametro beta

  agi<-(aginit)[[1]][1:1]
  bgi<-(bginit)[[1]][1:1]
  ggi<-(gginit)[[1]][1:1]

  #modelo gompertz
  mod_gomp<-numeric(0)
  try( mod_gomp<-nls(resp~(alfa*exp(-exp(beta-(gama*treat)))),
                     start=c(alfa=22,beta=1.5,gama=0.4)),
       silent=TRUE)

  if(length(mod_gomp)==0) {cat('
------------------------------------------------------------------------\n
Gompertz Model
------------------------------------------------------------------------\n
Error in model fit! (singularity, convergence, etc)
                               ')

                           tm4<-data.frame('Error in model fit!')
                           aic4 <- rs4 <- c('Error in model fit!')

  }

  if(length(mod_gomp)!=0) {

    b4<-summary(mod_gomp)
    tm4<-data.frame('Estimate' = round(b4[[10]][,1],8),
                    'Standard Error' = round(b4[[10]][,2],5),
                    'tc'=round(b4[[10]][,3],5),
                    'p-value' = round(b4[[10]][,4],5))
    rownames(tm4)<-c('Alpha','Beta','Gamma')

    yest_gomp=fitted(mod_gomp)
    aic4<-AIC(mod_gomp)

    # Coeficiente de Determinacao Aproximado
    fgo=fitted.values(mod_gomp)
    rg=residuals(mod_gomp)
    drg=sum((rg/resp)^2)
    srg=sum(rg^2)
    qeg=mean(resp)
    dg=(fgo-qeg)
    dqg=dg^2
    rs4=sum(dqg)/(sum(dqg)+srg)

    output4<-list('
------------------------------------------------------------------------\n
Gompertz Model
------------------------------------------------------------------------\n' = tm4,
'
------------------------------------------------------------------------\n
AIC
------------------------------------------------------------------------\n' = aic4,
'
------------------------------------------------------------------------\n
Approximate R2 of Gompertz Model
------------------------------------------------------------------------\n' = rs4)
    #print(output4,right=TRUE)
    stargazer(tm4, type = "text", title="Gompertz Model", summary=F, digits=4)
    stargazer(aic4, type = "text", title="AIC", digits=4, style="apsr")
    stargazer(rs4, type = "text", title="Approximate R2 of Gompertz Model", digits=4, flip=FALSE, summary=F, style="apsr")
    cat('------------------------------------------------------------------------\n')
  }



  ###################### Output ####################################
  mean.table<-tapply.stat(resp,treat,mean)
  colnames(mean.table)<-c('  Levels','   Observed Means')


  regout<-list("Table of means" = mean.table,
               "Power model coefficients" = tm1[,1], "power model AIC" = aic1,
               "Approximate R2 of power Model"= rs1,
               "Exponential model coefficients" = tm2[,1], "exponential model AIC" = aic2,
               "Approximate R2 of exponential Model"= rs2,
               "Logistic model coefficients" = tm3[,1], "logistic model AIC" = aic3,
               "Approximate R2 of logistic Model"= rs3,
               "Gompertz model coefficients" = tm4[,1], "Gompertz model AIC" = aic4,
               "Approximate R2 of Gompertz Model"= rs4)
  invisible(regout)

}
