#' Ordering the treatments according to the multiple comparison
#'
#' \code{order.group} It orders the groups of means.
#' @param trt Treatments.
#' @param means Means of treatment.
#' @param N Replications.
#' @param MSerror Mean square error.
#' @param Tprob Minimum value for the comparison.
#' @param std.err Standard error.
#' @param parameter Constante 1 (Sd), 0.5 (Sx).
#' @return
#' \itemize{
#' \item trt Factor
#' \item means Numeric
#' \item N Numeric
#' \item MSerror Numeric
#' \item Tprob value between 0 and 1
#' \item std.err Numeric
#' \item parameter Constant
#' }
#' @author Eric B Ferreira,
#'  \email{eric.ferreira@@unifal-mg.edu.br}
#'  @author Denismar Alves Nogueira
#'  @author Portya Piscitelli Cavalcanti
#'  (Adapted from Felipe de Mendiburu - GPL)
#' @seealso \code{\link{order.stat.SNK}}.
#' @export

order.group <-
function(trt, means, N, MSerror, Tprob, std.err, parameter = 1)
{
    N <- rep(1/mean(1/N), length(N))
    n <- length(means)
    z <- data.frame(trt, means, N, std.err)
    w <- z[order(z[, 2], decreasing = TRUE), ]
    M <- rep("", n)
    k <- 1
    j <- 1
    k <- 1
    cambio <- n
    cambio1 <- 0
    chequeo = 0
    M[1] <- letters[k]
    while (j < n) {
        chequeo <- chequeo + 1
        if (chequeo > n)
            break
        for (i in j:n) {
            minimo <- Tprob * sqrt(parameter * MSerror * (1/N[i] +
                1/N[j]))
            s <- abs(w[i, 2] - w[j, 2]) <= minimo
            if (s) {
                if (lastC(M[i]) != letters[k])
                  M[i] <- paste(M[i], letters[k], sep = "")
            }
            else {
                k <- k + 1
                cambio <- i
                cambio1 <- 0
                ja <- j
                for (jj in cambio:n) M[jj] <- paste(M[jj], " ",
                  sep = "")
                M[cambio] <- paste(M[cambio], letters[k], sep = "")
                for (v in ja:cambio) {
                  if (abs(w[v, 2] - w[cambio, 2]) > minimo) {
                    j <- j + 1
                    cambio1 <- 1
                  }
                  else break
                }
                break
            }
        }
        if (cambio1 == 0)
            j <- j + 1
    }
    w <- data.frame(w, stat = M)
    trt <- as.character(w$trt)
    means <- as.numeric(w$means)
    N <- as.numeric(w$N)
    std.err <- as.numeric(w$std.err)
    for (i in 1:n) {
        cat(M[i], "\t", trt[i], "\t", means[i], "\n")
    }
    output <- data.frame(trt, means, M, N, std.err)
    return(output)
}
