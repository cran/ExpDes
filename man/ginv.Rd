\name{ginv}
\alias{ginv}

\title{
Generalized inverse
}
\description{
Calculates the Moore-Penrose generalized inverse of a matrix X. 
}
\usage{
ginv(X, tol = sqrt(.Machine$double.eps))
}

\arguments{
  \item{X}{Matrix for which the Moore-Penrose inverse is required.}
  \item{tol}{A relative tolerance to detect zero singular values.}
}
\value{
A MP generalized inverse matrix for X.
}
\references{
Venables, W. N. and Ripley, B. D. (1999) Modern Applied Statistics with S-PLUS. Third Edition. Springer. p.100. 
}

\seealso{
See also: solve, svd, eigen.
}
\examples{
## Not run: 
# The function is currently defined as
function(X, tol = sqrt(.Machine$double.eps))
{
## Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] \%*\% (t(s$u[, nz])/s$d[nz]) else X, dimnames = dnx[2:1])
}

## End(Not run)
}
