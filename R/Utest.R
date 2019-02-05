


#' U test
#'
#'
#' Test for the separation of two groups.
#' The null hypothesis states that the groups are homogeneous and the alternative hypothesis states that they are separate.
#'
#'
#'
#' Either \code{data} or \code{md} should be provided.
#' If data are entered directly, Bn will be computed considering the squared Euclidean
#'  distance, which is compatible with \code{\link{is_homo}}, \code{\link{uclust}} and
#'   \code{\link{uhclust}}.
#'
#' For more details see Cybis, Gabriela B., Marcio Valk, and Sílvia RC Lopes. "Clustering and classification problems in genetics through U-statistics."
#' Journal of Statistical Computation and Simulation 88.10 (2018)
#'
#' @param group_id A vector of 0s and 1s indicating to which group the samples belong. Must be in the same order as data or md.
#' @param md Matrix of distances between all data points.
#' @param data Data matrix. Each row represents an observation.
#' @param numB Number of resampling iterations.
#'
#' @return Returns a list with the following elements:\describe{
#'   \item{Bn}{Test Statistic}
#'   \item{Pvalue}{Replication based p-value}
#'   \item{Replication}{Number of replications used to compute p-value}
#' }
#'
#' @seealso \code{\link{bn}},\code{\link{is_homo}}
#'
#' @examples
#'
#'# Simulate a dataset with two separate groups, the first 5 rows have mean 0 and
#'# the last 5 rows have mean 5.
#' data <- matrix(c(rnorm(75, 0), rnorm(75, 5)), nrow = 10, byrow=TRUE)
#'
#' # U test for mixed up groups
#' utest(group_id=c(1,0,1,0,1,0,1,0,1,0), data=data, numB=3000)
#' # U test for correct group definitions
#' utest(group_id=c(1,1,1,1,1,0,0,0,0,0), data=data, numB=3000)
#'
#'
#'@export
## Public version of the U test

utest <- function(group_id, md = NULL, data = NULL, numB=1000)
{
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

    group1 <- which(group_id == 0)
    group2 <- which(group_id == 1)
    ngv <- c(length(group1), length(group2))

    ng <- dim(md)[1]
    if (ng != sum(ngv)) {
        stop("Incorrect dimension or group_id")
    }

    md = md[c(group1,group2),c(group1,group2)]

   pvalue= Utest(ngv,numB,md)
   ans=list(Bn(ngv, md),pvalue, numB)
   names(ans)=c("Bn","Pvalue", "Replications")


   message(
       paste(
           "\t U-test for group separation  \n\nTest Statistic Bn =",
           round(ans$Bn, digits = 4),
           "\t p-value = ",
           round(ans$Pvalue, digits = 4),
           "\nAlternative hypothesis: The groups are not homogeneous, \nthere exists some separation between groups. "
       )
   )



   invisible(ans)

   }


############################################
# Internal function: U test
############################################
# test for homogeneity of two groups

# Riquire  Bn.R

# imput
   # ngv=c(n1,n2), n1=number of elements in the group 1 and n2=number of elements in the group 2
   # md is a nxn matrix of distances. n=n1+n2.
   # numB = number of bootstraps (we recomend 1000<numB<10000)
# output
   # p.value based on bootstraps
Utest <- function(ngv, numB, md)
{
    ng <- sum(ngv)
    comp <- Bn(ngv, md)
    vaux <- rep(0, ng)
    mataux <- matrix(0, ng, ng)
    absoluto <- rep(0, numB)

    for (i in 1:numB) {
        vaux <- floor(runif(ng, 1, (ng + 1)))
        mataux <- md[vaux, vaux]
        aux <- Bn(ngv, mataux)
        if (aux <= comp) {
            absoluto[i] <- 1
        }
        else {
            absoluto[i] <- 0
        }
    }

    p.value <- 1 - (sum(absoluto)) / numB
    p.value
}
