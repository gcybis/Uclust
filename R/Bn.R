#' Computes Bn Statistic.
#'
#' Returns the value for the Bn statistic that measures the degree of separation between two groups.
#' The statistic is computed through the difference of average within group distances to average between
#' group distances. Large values of Bn indicate large group separation. Under overall sample homogeneity
#' we have E(Bn)=0.
#'
#' Either \code{data} OR \code{md} should be provided.
#' If data are entered directly, Bn will be computed considering the squared Euclidean distance, which is compatible with
#' \code{\link{is_homo}}, \code{\link{uclust}} and \code{\link{uhclust}}.
#'
#' For more detail see Cybis, Gabriela B., Marcio Valk, and Sílvia RC Lopes. "Clustering and classification problems in genetics through U-statistics."
#' Journal of Statistical Computation and Simulation 88.10 (2018)
#' and Valk, Marcio, and Gabriela Bettella Cybis. "U-statistical inference for hierarchical clustering." arXiv preprint arXiv:1805.12179 (2018).
#' @param group_id A vector of 0s and 1s indicating to which group the samples belong. Must be in the same order as data or md.
#' @param md Matrix of distances between all data points.
#' @param data Data matrix. Each row represents an observation.
#' @return Value of the Bn statistic.
#'
#' @examples
#' n=5
#' x=matrix(rnorm(n*10),ncol=10)
#' bn(c(1,0,0,0,0),data=x)     # option (a) entering the data matrix directly
#' md=as.matrix(dist(x))^2
#' bn(c(0,1,1,1,1),md)         # option (b) entering the distance matrix
#'
#' @export
bn <- function(group_id, md = NULL, data = NULL) {
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

    group1 <- (group_id == 0)
    group2 <- group_id == 1
    ngv <- c(sum(group1), sum(group2))

    ng <- dim(md)[1]
    if (ng != sum(ngv)) {
        stop("Incorrect dimension or group_id")
    }

    if (min(ngv) > 1) {
        #for groups with n1>2 (oiriginal definition of Bn)

        s11 <- sum(md[group1, group1]) / 2
        s22 <- sum(md[group2, group2]) / 2
        s12 <- sum(md[group1, group2])

        a1 <- (1 / (ngv[1] * ngv[2])) * s12
        a2 <- (2 / (ngv[1] * (ngv[1] - 1))) * s11
        a3 <- (2 / (ngv[2] * (ngv[2] - 1))) * s22
        sBn <- (ngv[1] * ngv[2] / (ng * (ng - 1))) * (2 * a1 - a2 - a3)
    } else{
        #if n1=1 (extended definition of Bn)
        if (ngv[1] == 1) {
            s22 <- sum(md[group2, group2]) / 2
            s12 <- sum(md[group1, group2])

            a1 <- (1 / (ngv[1] * ngv[2])) * s12
            a2 <- 0
            a3 <- (2 / (ngv[2] * (ngv[2] - 1))) * s22
            sBn <- (ngv[1] * ngv[2] / (ng * (ng - 1))) * (a1 - a2 - a3)
        } else{
            s11 <- sum(md[group1, group1]) / 2
            s12 <- sum(md[group1, group2])

            a1 <- (1 / (ngv[1] * ngv[2])) * s12
            a2 <- (2 / (ngv[1] * (ngv[1] - 1))) * s11
            a3 <- 0
            sBn <- (ngv[1] * ngv[2] / (ng * (ng - 1))) * (a1 - a2 - a3)
        }
    }
    sBn
}

############################################
# Internal function: computes Bn
############################################
#Computes Bn
## ngv is a vector with 2 entries: size of group 1 (n1) and size of group 2 (n2)
## md is the distance matrix, ordered so that the first n1 elements are from group 1t

Bn <- function(ngv, md) {
    ng <- sum(ngv)
    maux1 <- matrix(0, nrow = ng, ncol = ng)

    if (min(ngv) > 1) {
        #for groups with n1>2 (oiriginal definition of Bn)

        s11 <- sum(md[(1:ngv[1]), 1:ngv[1]]) / 2
        s22 <- sum(md[(ngv[1] + 1):ng, (ngv[1] + 1):ng]) / 2
        s12 <- sum(md[1:ngv[1], (ngv[1] + 1):ng])

        a1 <- (1 / (ngv[1] * ngv[2])) * s12
        a2 <- (2 / (ngv[1] * (ngv[1] - 1))) * s11
        a3 <- (2 / (ngv[2] * (ngv[2] - 1))) * s22
        sBn <- (ngv[1] * ngv[2] / (ng * (ng - 1))) * (2 * a1 - a2 - a3)
    }
    else {
        #if n1=1 (extended definition of Bn)
        if (ngv[1] == 1) {
            s22 <- sum(md[(ngv[1] + 1):ng, (ngv[1] + 1):ng]) / 2
            s12 <- sum(md[1:ngv[1], (ngv[1] + 1):ng])

            a1 <- (1 / (ngv[1] * ngv[2])) * s12
            a2 <- 0
            a3 <- (2 / (ngv[2] * (ngv[2] - 1))) * s22
            sBn <- (ngv[1] * ngv[2] / (ng * (ng - 1))) * (a1 - a2 - a3)
        }
        else {
            s11 <- sum(md[(1:ngv[1]), 1:ngv[1]]) / 2
            s22 <- sum(md[(ngv[1] + 1):ng, (ngv[1] + 1):ng]) / 2
            s12 <- sum(md[1:ngv[1], (ngv[1] + 1):ng])


            a1 <- (1 / (ngv[1] * ngv[2])) * s12
            a2 <- (2 / (ngv[1] * (ngv[1] - 1))) * s11
            a3 <- 0
            sBn <- (ngv[1] * ngv[2] / (ng * (ng - 1))) * (a1 - a2 - a3)
        }
    }

    sBn
}
