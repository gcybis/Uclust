#############################################################
#### Homogeneity test main function
#############################################################
#' U-statistic based homogeneity test
#'
#' Homogeneity test based on the statistic \code{bn}. The test assesses whether there exists a data partition
#' for which group separation is statistically significant according to the U-test. The null hypothesis
#' is overall sample homogeneity, and a sample is considered homogeneous if it cannot be divided into
#' two statistically significant subgroups.
#'
#' This is the homogeneity test of Cybis et al. (2017) extended to account for groups of size 1.
#' The test is performed through two steps: an optimization procedure that finds the data partition that
#' maximizes the standardized Bn and a test for the resulting maximal partition. Should be used in high dimension small sample size settings.
#'
#'
#'
#' Either \code{data} or \code{md} should be provided.
#' If data are entered directly, Bn will be computed considering the squared Euclidean distance.
#' It is important that if a distance matrix is entered, it consists of squared Euclidean distances, otherwise test results are
#' invalid.
#'
#' Variance of \code{bn} is estimated through resampling, and thus, p-values may vary a bit in different runs.
#'
#' For more detail see Cybis, Gabriela B., Marcio Valk, and Sílvia RC Lopes. "Clustering and classification problems in genetics through U-statistics."
#' Journal of Statistical Computation and Simulation 88.10 (2018)
#' and Valk, Marcio, and Gabriela Bettella Cybis. "U-statistical inference for hierarchical clustering." arXiv preprint arXiv:1805.12179 (2018).
#'
#'
#' @param md Matrix of squared Euclidean distances between all data points.
#' @param data Data matrix. Each row represents an observation.
#' @param rep Number of times to repeat optimization procedure. Important for problems with
#' multiple optima.
#' @return Returns a list with the following elements:\describe{
#'   \item{minFobj}{Test statistic. Minimum of the objective function for optimization (-stdBn).}
#'   \item{group1}{Elements in group 1 in the maximal partition. (obs: this is not the best
#'   partition for the data, see \code{uclust})}
#'   \item{group2}{Elements in group 2 in the maximal partition.}
#'   \item{p.MaxTest}{P-value for the homogeneity test.}
#'   \item{Rep.Fobj}{Values for the minimum objective function on all \code{rep} optimization runs.}
#'   \item{bootB}{Resampling variance estimate for partitions with groups of size n/2 (or (n-1)/2 and (n+1)/2 if n is odd).}
#'   \item{bootB1}{Resampling variance estimate for partitions with one group of size 1.}
#'
#' }
#'
#'
#' @examples
#'x = matrix(rnorm(500000),nrow=50)  #creating homogeneous Gaussian dataset
#'res = is_homo(data=x)
#'
#'x[1:30,] = x[1:30,]+0.15   #Heterogeneous dataset (first 30 samples have different mean)
#'res = is_homo(data=x)
#'
#'md = as.matrix(dist(x)^2)   #squared Euclidean distances for the same data
#'res = is_homo(md)
#'
#'# Multidimensional sacling plot of distance matrix
#'fit <- cmdscale(md, eig = TRUE, k = 2)
#'x <- fit$points[, 1]
#'y <- fit$points[, 2]
#'plot(x,y, main=paste("Homogeneity test: p-value =",res$p.MaxTest))
#'
#'@export
is_homo <- function(md = NULL, data = NULL, rep = 10) {
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

    n <- dim(md)[1]
    if (n<=3){
        stop("samples size n is too small for homogeneity test")
    }


    Fobj <- vector()
    grupo1 <- list()

    ans <- optimstdBn(md)
    Fobj[1] <- ans$Fobj[length(ans$Fobj)]
    grupo1[[1]] <- ans$grupo1
    bootB <- ans$bootB

    bootB1 <- ans$bootB1

    for (i in 2:rep) {
        ans <- optimstdBn(md, bootB = bootB, bootB1 = bootB1)
        Fobj[i] <- ans$Fobj[length(ans$Fobj)]
        grupo1[[i]] <- ans$grupo1
    }

    minFobj <- min(Fobj)
    g1 <- grupo1[[which(Fobj == min(Fobj))[1]]]
    g2 <- (1:n)[-g1]
    p.homo <- test_max_stdBn(minFobj, n)

   # if (length(g1) == 1) {
#        varBn <- bootB1
#    } else{
#        varBn <- bootB * smile(n)[length(g1)] / smile(n)[floor(n / 2)]    ############### Vamos remover o Var Bn - pois nao é usado e tem erro - CONFIRMAR
#    }

    ans <- list(minFobj, g1, g2, p.homo, Fobj, bootB, bootB1)#, varBn)
    names(ans) <- list("minFobj", "group1", "group2", "p.MaxTest", "Rep.Fobj", "bootB", "bootB1")#, "varBn")

    message(
        paste(
            "\t Homogeneity test - Max distribution \n\nmax(stdBn) =",
            round(-ans$minFobj, digits = 4),
            "\t p-value = ",
            round(ans$p.MaxTest, digits = 4),
            "\nAlternative Hypothesis: The data can be divided into two significant subgroups."
        )
    )

    invisible(ans)
}

#############################################################
####function in old smile.R
#############################################################
## Computes smile function - used for relation between Bn variance for different group sizes
smile <- function(n) {
    Cn <- vector()

    for (n11 in 2:(n - 2)) {
        n22 <- n - n11
        Cn[n11] <- (((n11 * n22) / (n * (n - 1)) ^ 2) * (2 * n ^ 2 - 6 * n + 4) / ((n11 - 1) * (n22 - 1)))
    }

    return(Cn)
}

#############################################################
##### Max test functions
#############################################################
# Gumbel Correction
gumbel_correction <- function(nt, Fobj) {
    bn <- sqrt(2 * log(nt)) - (log(log(nt)) + log(4 * pi * log(2) ^ 2)) / (2 * sqrt(2 * log(nt)))
    an <- (log((4 * log(2) ^ 2) / log(4 / 3) ^ 2)) / (2 * sqrt(2 * log(nt)))
    bgumbel <- (1 / an) * (Fobj - bn)

    return(1 - exp(-exp(-bgumbel))) # Gumbel cdf
}

#Returns the p-value for the maximum standardized Bn
test_max_stdBn <- function(Fobj, n, nt = 0) {
    if (nt == 0) {
        #computes number of tests if not provided
        nt <- 2 ^ (n - 1) - 1
    }
    if (nt < 2 ^ (28)) {
        # For small n (n<30) uses regular max test
        p <- 1 - exp((nt) * pnorm(-Fobj, log.p = TRUE)) #p-value
    }
    else {
        p <- gumbel_correction(nt, -Fobj)
    } # For larger n uses Gumbel correction

    return(p)
}





#############################################################
##### Computes the standardized Bn objective function
#############################################################
# for use in standardized Bn optimization only
# returns stdBn for particular assignment
objstdBn <- function(assign, varBn, mdm) {
    n <- length(assign)
    n1 <- sum(assign == 0)
    n2 <- n - n1

    if (n1 == 0  | n2 == 0) {
        return(Inf)
    }
    else {
        vaux1 <- which(assign == 0)
        vaux2 <- c(1:n)[-vaux1]
        m <- mdm[c(vaux1, vaux2), c(vaux1, vaux2)]
        Bns <- Bn(c(n1, n2), m)

        return(-Bns / sqrt(varBn[n1]))
    }
}

#############################################################
##### Standardized Bn Optimization
##### Returns the grouping for which objective function (-stdBn)
##### converges to minumim
#############################################################
optimstdBn <- function(mdm, itmax = 200, centers = -1, bootB = -1, bootB1 = -1) {
    md <- mdm # distance matrix

    n <- dim(md)[1] # sample size

    #initiates data structures for keeping tabs on  assignments
    it <- 1
    ass <- vector()
    ass_old <- rep(2, n)
    ASS <- matrix(ncol = n, nrow = itmax)   # Matriz keeps optimization history
    Fobj <- vector()

    #smile function
    Cn <- vector()
    varBn <- vector()
    numB <- 2000

    if (bootB1 == -1) {
        bootB1 <- boot_sigma1(c(1, (n - 1)), md)
    }

    if (bootB == -1) {
        bootB <- boot_sigma(c(floor(n / 2), (n - floor(n / 2))), numB, md) # returns variance of Bn with group size  c(floor(n/2),(n-floor(n/2))
    }

    for (n1 in 2:(n - 2)) {
        n2 <- n - n1
        Cn[n1] <- (((n1 * n2) / (n * (n - 1)) ^ 2) * (2 * n ^ 2 - 6 * n + 4) / ((n1 - 1) * (n2 - 1)))
    }

    for (n1 in 2:(n - 2)) {
        n2 <- n - n1
        varBn[n1] <- Cn[n1] * bootB / Cn[floor(n / 2)]
    }

    varBn[1] <- bootB1
    varBn[n - 1] <- bootB1


    ### Initializing Parameters for optimization

    # random cluster centers
    if (centers == -1) {
        centers <- sample(n, 2)
    }
    # print(centers)
    # initializa assignments - closest centers

    for (i in 1:n) {
        ass[i] <- (md[i, centers[1]]) > (md[i, centers[2]])
    }

    ass

    ASS[1, ] <- ass

    # TRUE bellongs to group 2
    #### Start iterations

    while (it < itmax && !prod(ass == ass_old)) {
        ass_old <- ass

        ord <- sample(1:n)
        for (i in ord) {
            ass[i] <- 0
            f0 <- objstdBn(ass, varBn, md)
            ass[i] <- 1
            f1 <- objstdBn(ass, varBn, md)

            if (f0 < f1) {
                ass[i] <- 0
            }
        }

        Fobj[it] <- objstdBn(ass, varBn, md)

        it <- it + 1
        ASS[it, ] <- ass
    }

    # Assmble answer
    ans <- list(which(ass == ass[1]), Fobj, it - 1, ASS[1:(it + 1), ], varBn, bootB, bootB1)
    names(ans) <- c("grupo1", "Fobj", "numIt", "history", "varBn", "bootB", "bootB1")

    ans
}
