#' Generalized inverse
#'
#' \code{ginv} Computes the Moore-Penrose generalized inverse
#' of a matrix X.
#' @param X Matrix for which the Moore-Penrose inverse is
#' required.
#' @param tol A relative tolerance to detect zero singular
#' values.
#' @return A MP generalized inverse matrix for X.
#' @references Venables, W. N. and Ripley, B. D. (1999) Modern
#' Applied Statistics with S-PLUS. Third Edition. Springer.
#' p.100.
#' @seealso \code{\link{solve}}, \code{\link{svd}},
#' \code{\link{eigen}}
#' @export

ginv <-
function(X, tol = sqrt(.Machine$double.eps))
{
    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
        stop("'X' must be a numeric or complex matrix")
    if (!is.matrix(X))
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X))
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive))
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive))
        array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
        t(Xsvd$u[, Positive, drop = FALSE]))
}
