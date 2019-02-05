#' Variance of Bn
#'
#'
#' Estimates the variance of the Bn statistic using the resampling
#'  procedure described in Cybis, Gabriela B., Marcio Valk, and Sílvia RC Lopes. "Clustering and classification problems in genetics through U-statistics."
#' Journal of Statistical Computation and Simulation 88.10 (2018)
#' and Valk, Marcio, and Gabriela Bettella Cybis. "U-statistical inference for hierarchical clustering." arXiv preprint arXiv:1805.12179 (2018).
#'
#' Either \code{data} or \code{md} should be provided.
#' If data are entered directly, Bn will be computed considering the squared Euclidean
#'  distance, which is compatible with \code{\link{is_homo}}, \code{\link{uclust}} and
#'   \code{\link{uhclust}}.
#'
#' @param group_sizes A vector with two entries: size of group 1 and size of group 2.
#' @param md Matrix of distances between all data points.
#' @param data Data matrix. Each row represents an observation.
#' @param numB Number of resampling iterations. Only used if no groups are of size 1.
#'
#' @return Variance of Bn
#'
#' @seealso \code{\link{bn}}
#'
#' @examples
#'
#' n=5
#' x=matrix(rnorm(n*20),ncol=20)
#' # option (a) entering the data matrix directly and considering a group of size 1
#' var_bn(c(1,4),data=x)
#'
#' # option (b) entering the distance matrix and considering a groups of size 2 and 3
#' md=as.matrix(dist(x))^2
#' var_bn(c(2,3),md)
#'
#'
#'
# NAMESPACE INFO
#' @import robcor
#' @export
#'
#' @importFrom grDevices rainbow
#' @importFrom graphics plot text
#' @importFrom stats as.dendrogram as.hclust dist pnorm runif var
#' @importFrom utils combn
#'
#'
var_bn <- function(group_sizes, md = NULL, data = NULL, numB = 2000) {
    if (is.null(md)) {
        # Computing data matrix if one is not provided
        if (is.null(data)) {
            stop("No data provided")
        }
        md <- as.matrix(dist(data) ^ 2)
    }

    if (class(md) != "matrix") {
        stop("md is not of class matrix")
    }

    ngv <- group_sizes

    if (min(ngv) == 1) {
        var <- boot_sigma1(ngv, md)
    } else {
        var <- boot_sigma(ngv, numB, md)
    }

    return(var)
}

###########################################################
# Bootstrap to find sigma (variance of Bn) for n1>2       #
###########################################################
boot_sigma <- function(ngv, numB = 2000, md)
{
    n <- sum(ngv)
    if (n > 6) {
        B <- rep(0, numB)
        for (i in 1:numB) {
            vaux1 <- sample(n, n, replace = FALSE)
            mataux <- md[vaux1, vaux1]
            B[i] <- Bn(ngv, mataux)#Bn
        }

        varBboot <- var(B)
    }
    else {
        B <- vector()
        cb <- combn(n, 2)
        for (i in 1:(dim(cb)[2])) {
            vaux1 <- cb[, i]
            vaux2 <- c(1:n)[-vaux1]
            vaux <- c(vaux1, vaux2)
            mataux <- md[vaux, vaux]
            B[i] <- Bn(ngv, mataux)#Bn
        }

        a <- robcor::robacf(B, lag.max = 0, type = "covariance", plot = FALSE)
        varBboot <- as.numeric(a$acf)
    }

    return(varBboot)
}

###########################################################
# Bootstrap to find sigma (variance of Bn) for n1=1       #
###########################################################
boot_sigma1 <- function(ngv, md)
{
    n <- sum(ngv)
    B <- rep(0, n)
    for (i in 1:n) {
        vaux1 <- i
        vaux2 <- c(1:n)[-i]
        vaux <- c(vaux1, vaux2)
        mataux <- md[vaux, vaux]
        B[i] <- Bn(c(1, n - 1), mataux)#Bn
    }

    a <- robcor::robacf(B, lag.max = 0, type = "covariance", plot = FALSE)

    varBboot <- as.numeric(a$acf)

    return(varBboot)
}
