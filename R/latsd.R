#' Latin Square Design
#'
#' \code{lastd} Analyses experiments in balanced Latin Square
#' Design, considering a fixed model.
#' @param treat Numeric or complex vector containing the
#' treatments.
#' @param row Numeric or complex vector containing the rows.
#' @param column Numeric or complex vector containing the
#' columns.
#' @param resp Numeric or complex vector containing the
#' response variable.
#' @param quali Logic. If TRUE (default), the treatments are
#' assumed qualitative, if FALSE, quantitatives.
#' @param  mcomp Allows choosing the multiple comparison test;
#' the \emph{default} is the test of Tukey, however, the
#' options are: the LSD test ('lsd'), the LSD test with
#' Bonferroni protection ('lsdb'), the test of Duncan
#' ('duncan'), the test of Student-Newman-Keuls ('snk'), the
#' test of Scott-Knott ('sk'), the Calinski and Corsten test
#' ('ccF') and bootstrap multiple comparison's test ('ccboot').
#' @param sigT The signficance to be used for the multiple
#' comparison test; the default is 5\%.
#' @param sigF The signficance to be used for the F test of
#' ANOVA; the default is 5\%.
#' @details The arguments sigT and mcomp will be used only
#' when the treatment are qualitative.
#' @return The output contains the ANOVA of the LSD, the
#' Shapiro-Wilk normality test for the residuals of the model,
#' the fitted regression models (when the treatments are
#' quantitative) and/or the multiple comparison tests (when the
#' treatments are qualitative).
#' @references GOMES, F. P. Curso de Estatistica Experimental.
#' 10a ed. Piracicaba: ESALQ/USP. 1982. 430.
#'
#' FERREIRA, E. B.; CAVALCANTI, P. P.; NOGUEIRA D. A. Funcao
#' em codigo R para analisar experimentos em DQL simples, em
#' uma so rodada. In: CONGRESSO DE POS-GRADUACAO DA
#' UNIVERSIDADE FEDERAL DE LAVRAS, 18., 2009, Lavras.
#' Annals... Lavras: UFLA, 2009.
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#'  @author Denismar Alves Nogueira
#'  @author Portya Piscitelli Cavalcanti
#'  @note The \code{\link{graphics}} can be used to construct
#'  regression plots and \code{\link{plotres}} for residuals
#'  plots.
#' @seealso \code{\link{crd}}, \code{\link{rbd}}.
#' @examples
#' data(ex3)
#' attach(ex3)
#' latsd(trat, linha, coluna, resp, quali = TRUE, mcomp = "snk",
#' sigT = 0.05, sigF = 0.05)
#' @export

latsd <-
function(treat, row, column, resp, quali=TRUE, mcomp='tukey', sigT=0.05, sigF=0.05) {

Trat<-factor(treat)
Linha<-factor(row)
Coluna<-factor(column)
anava<-aov(resp~Trat+Linha+Coluna)
tab<-summary(anava)

colnames(tab[[1]])<-c('DF','SS','MS','Fc','Pr>Fc')
tab[[1]]<-rbind(tab[[1]],c(apply(tab[[1]],2,sum)))
rownames(tab[[1]])<-c('Treatament','Row','Column','Residuals','Total')
cv<-round(sqrt(tab[[1]][4,3])/mean(resp)*100, 2)
tab[[1]][5,3]=NA
cat('------------------------------------------------------------------------\nAnalysis of Variance Table
------------------------------------------------------------------------\n')
print(tab[[1]])
cat('------------------------------------------------------------------------\nCV =',cv,'%\n')


#Teste de normalidade
pvalor.shapiro<-shapiro.test(anava$residuals)$p.value
cat('\n------------------------------------------------------------------------\nShapiro-Wilk normality test\n')
cat('p-value: ',pvalor.shapiro, '\n')
if(pvalor.shapiro<0.05){cat('WARNING: at 5% of significance, residuals can not be considered normal!
------------------------------------------------------------------------\n')}
else{cat('According to Shapiro-Wilk normality test at 5% of significance, residuals can be considered normal.
------------------------------------------------------------------------\n')}

if(tab[[1]][1,5]<sigF){

if(quali==TRUE) {

  if(mcomp=='tukey'){
    tukey(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
                    }
  if(mcomp=='duncan'){
    duncan(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
                    }
  if(mcomp=='lsd'){
    lsd(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
                    }
  if(mcomp=='lsdb'){
    lsdb(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
                    }
  if(mcomp=='sk'){
    scottknott(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
                    }
  if(mcomp=='snk'){
    snk(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
                    }
  if(mcomp=='ccboot'){
    ccboot(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
                    }
    if(mcomp=='ccF'){
    ccF(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
                    }
                }
else{
    reg<-reg.poly(resp, treat, tab[[1]][4,1], tab[[1]][4,2], tab[[1]][1,1], tab[[1]][1,2])
}
                       }
else {
    cat('\nAccording to the F test, the means can not be considered distinct.\n')
mean.table<-tapply.stat(resp,treat,mean)
colnames(mean.table)<-c('Levels','Means')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}

#Saida
out<-list()
out$residuals<-anava$residuals
out$df.residual<-anava$df.residual
out$coefficients<-anava$coefficients
out$effects<-anava$effects
out$fitted.values<-anava$fitted.values
out$means<-tapply.stat(resp,treat,mean)
if(quali==FALSE && tab[[1]][1,5]<sigF) {out$reg<-reg}
invisible(out)
}
