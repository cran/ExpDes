#' Multiple comparison: Bonferroni's Least Significant
#' Difference test
#'
#' \code{lsdb} Performs the t test (LSD) with Bonferroni's
#' protection, for multiple comparison of means.
#' @param y Numeric or complex vector containing the response
#' variable.
#' @param trt Numeric or complex vector containing the
#' treatments.
#' @param DFerror Error degrees of freedom.
#' @param SSerror Error sum of squares.
#' @param alpha Significance level.
#' @param group TRUE or FALSE.
#' @param main Title.
#' @return Returns the multiple comparison of means according
#' to the LSDB test.
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#' @author Denismar Alves Nogueira
#' @author Portya Piscitelli Cavalcanti
#' @seealso \code{\link{snk}}, \code{\link{duncan}},
#' \code{\link{ccboot}}, \code{\link{lsd}},
#' \code{\link{scottknott}}, \code{\link{tukey}},
#' \code{\link{ccF}}.
#' @examples
#' data(ex1)
#' attach(ex1)
#' crd(trat, ig, quali = TRUE, mcomp = "lsdb", sigT = 0.05)
#' @importFrom "utils" "combn" "read.table"
#' @export

lsdb <-
function (y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE,    main = NULL)
{
    MSerror <- SSerror/DFerror
    name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    means <- tapply.stat(junto[, 1], junto[, 2], stat = "mean")
    sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
    nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
    means <- data.frame(means, std.err = sds[, 2]/sqrt(nn[, 2]),
        replication = nn[, 2])
    names(means)[1:2] <- c(name.t, name.y)
    ntr <- nrow(means)
    #Tprob <- qtukey(1 - alpha, ntr, DFerror)
    alphap=(2*alpha)/(ntr*(ntr-1))
    Tprob <- qt(1-(alphap/2),DFerror)*sqrt(2)
    nr <- unique(nn[, 2])
    nfila <- c("Alpha", "Error Degrees of Freedom", "Error Mean Square",
        "Critical Value of Studentized Range")
    nvalor <- c(alpha, DFerror, MSerror, Tprob)
    #cat("\nStudy:", main)
    #cat("\n\nHSD Test for", name.y, "\n")
    xtabla <- data.frame(...... = nvalor)
    row.names(xtabla) <- nfila
    #print(xtabla)
    #cat("\nTreatment Means\n")
    #print(data.frame(row.names = NULL, means))
    if (group) {
        if (length(nr) == 1) {
            HSD <- Tprob * sqrt(MSerror/nr)
            #cat("\nHonestly Significant Difference", HSD)
        }
        else {
            nr1 <- 1/mean(1/nn[, 2])
            HSD <- Tprob * sqrt(MSerror/nr1)
            #cat("\nHonestly Significant Difference", HSD)
            #cat("\nHarmonic Mean of Cell Sizes ", nr1)
            #cat("\n\nDifferent HSD for each comparison")
        }
        cat("\nBonferroni's T test (protected LSD)\n------------------------------------------------------------------------")
        cat("\nGroups  Treatments  Means\n")
        output <- order.group(means[, 1], means[, 2], means[,4], MSerror, Tprob, means[, 3], parameter = 0.5)
        cat('------------------------------------------------------------------------\n')
    }
    if (!group) {
        comb <- combn(ntr, 2)
        nn <- ncol(comb)
        dif <- rep(0, nn)
        pvalue <- rep(0, nn)
        for (k in 1:nn) {
            i <- comb[1, k]
            j <- comb[2, k]
            dif[k] <- abs(means[i, 2] - means[j, 2])
            sdtdif <- sqrt(MSerror * (1/means[i, 4] + 1/means[j,
                4]))
            pvalue[k] <- round(1 - ptukey(dif[k] * sqrt(2)/sdtdif,
                ntr, DFerror), 4)
        }
        tr.i <- comb[1, ]
        tr.j <- comb[2, ]
        #cat("\nComparison between treatments means\n\n")
        #print(data.frame(row.names = NULL, tr.i, tr.j, diff = dif, pvalue = pvalue))
        output <- data.frame(trt = means[, 1], means = means[,2], M = "", N = means[, 4], std.err = means[, 3])
    }
#    return(output)
}
