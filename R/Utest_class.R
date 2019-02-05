#' Test for classification of a sample in one of two groups.
#'
#' The null hypothesis is that the new data is not well classified into the first group when compared to the second group. The
#' alternative hypothesis is that the data is well classified into the first group.
#'
#' The test is performed considering the squared Euclidean distance.
#'
#' For more detail see Cybis, Gabriela B., Marcio Valk, and Sílvia RC Lopes. "Clustering and classification problems in genetics through U-statistics."
#' Journal of Statistical Computation and Simulation 88.10 (2018)
#' and Valk, Marcio, and Gabriela Bettella Cybis. "U-statistical inference for hierarchical clustering." arXiv preprint arXiv:1805.12179 (2018).
#' @param x A numeric vector to be classified.
#' @param data Data matrix. Each row represents an observation.
#' @param group_id A vector of 0s (first group) and 1s indicating to which group the samples belong. Must be in the same order as data.
#' @param bootstrap_iter Numeric scalar. The number of bootstraps. It's recommended
#'   \eqn{1000 < bootstrap_iter < 10000}.
#' @return A list with class "utest_classify" containing the following components:
#'     \item{statistic}{the value of the test statistic.}
#'     \item{p_value}{The p-value for the test.}
#'     \item{bootstrap_iter}{the number of bootstrap iterations.}
#' @examples
#' # Example 1
#' # Five observations from each group, G1 and G2. Each observation has 60 dimensions.
#' data <- matrix(c(rnorm(300, 0), rnorm(300, 10)), ncol = 60, byrow=TRUE)
#' # Test data comes from G1.
#' x <- rnorm(60, 0)
#' # The test correctly indicates that the test data should be classified into G1 (p < 0.05).
#' utest_classify(x, data, group_id = c(rep(0,times=5),rep(1,times=5)))
#'
#' # Example 2
#' # Five observations from each group, G1 and G2. Each observation has 60 dimensions.
#' data <- matrix(c(rnorm(300, 0), rnorm(300, 10)), ncol = 60, byrow=TRUE)
#' # Test data comes from G2.
#' x <- rnorm(60, 10)
#' # The test correctly indicates that the test data should be classified into G2 (p > 0.05).
#' utest_classify(x, data, group_id = c(rep(1,times=5),rep(0,times=5)))
#' @export
#'
utest_classify <- function(x, data, group_id, bootstrap_iter = 1000)
{


    ng1 <- sum(group_id==0)
    ng2 <- sum(group_id==1)

    if((ng1+ng2)!=nrow(data)){stop("Incorrect group_id element.")}


    data <- rbind(data[group_id==0,],x,data[group_id==1,])


    distances <- as.matrix(dist(data))
    B <- rep(0, bootstrap_iter)
    B1_0 <- Bn(c(ng1 + 1, ng2), distances)
    B2_0 <- Bn(c(ng1, 1 + ng2), distances)

    D = B1_0 - B2_0

    for (i in 1:bootstrap_iter)
    {
        vaux1 <- floor(runif(ng1, 1, ng1 + 1))
        vaux2 <- floor(runif(ng2, ng1 + 1, ng2 + ng1 + 2))
        vaux3 <- floor(runif(1, ng1 + 1, ng2 + ng1 + 2))
        mataux <- distances[c(vaux1, vaux3, vaux2), c(vaux1, vaux3, vaux2)]
        ngv1 <- c(ng1 + 1, ng2)
        B1 <- Bn(ngv1, mataux)
        ngv2 <- c(ng1, ng2 + 1)
        B2 <- Bn(ngv2, mataux)
        B[i] <- B1 - B2
    }

    utest_classify_obj <- list(test_statistic = D,
                               p_value = mean(B > D),
                               #groups = groups,
                               bootstrap_iter = bootstrap_iter)

    class(utest_classify_obj) <- "utest_classify"

    return(utest_classify_obj)
}

#' Simple print method for utest_classify objects.
#'
#'
#'
#'
#' @param x  utest_classify object
#' @param ...  additional parameters passed to the function
#'
#' @export
print.utest_classify <- function(x, ...) {
    cat("\n")
    cat("\tU test for classification\n\n")
    cat("Test statistic:", x$test_statistic, "\n")
    cat("p-value:", x$p_value, (.level_symbol(x$p_value)), "\n")
    cat("Alternative hypothesis: it should be classified into the first group\n")
    #cat("First group:", obj$groups[1], ", second group:", obj$groups[2], "\n")
    cat("Bootstrap iterations:", x$bootstrap_iter, "\n")
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
    cat("\n")

    invisible(x)
}

.level_symbol <- function(p_value) {
    symbols <- c("***", "**", "*", ".", "")
    sig_levels <- c(0.001, 0.01, 0.05, 0.1, 1)

    return(symbols[which.max(p_value < sig_levels)])
}
