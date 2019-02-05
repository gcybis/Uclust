#########################################################
###   Main Uclust function
#########################################################
#' U-statistic based significance clustering
#'
#' Partitions the sample into the two significant subgroups with the largest Bn statistic. If no significant
#' partition exists, the test will return "homogeneous".
#'
#' This is the significance clustering procedure of Valk and Cybis (2018).
#' The method first performs a homogeneity test to verify whether the data can be significantly
#' partitioned. If the hypothesis of homogeneity is rejected, then the method will search, among all
#' the significant partitions, for the partition that better separates the data, as measured by larger
#' \code{bn} statistic. This function should be used in high dimension small sample size settings.
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
#'See also \code{is_homo}, \code{uhclust}, \code{Utest_class}.
#'
#' @param md Matrix of squared Euclidean distances between all data points.
#' @param data Data matrix. Each row represents an observation.
#' @param alpha Significance level.
#' @param rep Number of times to repeat optimization procedures. Important for problems with
#' multiple optima.
#' @return  Returns a list with the following elements:\describe{
#'   \item{cluster1}{Elements in group 1 in the final partition. This is the significant partition with
#'   maximal Bn, if sample is heterogeneous.}
#'   \item{cluster2}{Elements in group 2 in the final partition.}
#'   \item{p.value}{P-value for the test that renders the final partition, if heterogeneous.
#'   Homogeneity test p-value, if homogeneous.}
#'   \item{alpha_corrected}{Bonferroni corrected significance level for the test that renders the final
#'   partition, if heterogeneous. Homogeneity test significance level, if homogeneous.}
#'   \item{n1}{Size of the smallest cluster}
#'   \item{ishomo}{Logical, returns \code{TRUE} when the sample is homogeneous.}
#'   \item{Bn}{Value of Bn statistic for the final partition, if heterogeneous.
#'   Value of Bn statistic for the maximal homogeneity test partition, if homogeneous.}
#'   \item{varBn}{Variance estimate for final partition, if heterogeneous.
#'   Variance estimate for the maximal homogeneity test partition, if homogeneous.}
#'   \item{ishomoResult}{Result of homogeneity test (see \code{is_homo}).}
#' }
#'
#' @examples
#' set.seed(17161)
#' x = matrix(rnorm(100000),nrow=50)  #creating homogeneous Gaussian dataset
#' res = uclust(data=x)
#'
#' x[1:30,] = x[1:30,]+0.25   #Heterogeneous dataset (first 30 samples have different mean)
#' res = uclust(data=x)
#'
#' md = as.matrix(dist(x)^2)   #squared Euclidean distances for the same data
#' res = uclust(md)
#'
#' # Multidimensional scaling plot of distance matrix
#' fit <- cmdscale(md, eig = TRUE, k = 2)
#' x <- fit$points[, 1]
#' y <- fit$points[, 2]
#' col=rep(3,dim(md)[1])
#' col[res$cluster2]=2
#' plot(x,y, main=paste("Multidimensional scaling plot of data:
#'                     homogeneity p-value =",res$ishomoResult$p.MaxTest),col=col)
#'
#'
#'@export
uclust <- function(md = NULL, data = NULL, alpha = 0.05, rep = 15) {
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
    is.h <- suppressMessages(is_homo(rep = rep, md)) # finds the config with max std Bn and applies max test.
    ResultadoTesteIsHomo <- is.h   #tirar depois

    if (is.h$p.MaxTest < alpha) {
        # if non homogeneous
        bootB <- is.h$bootB
        bootB1 <- is.h$bootB1

        oBn <- rep.optimBn(md, bootB = bootB, rep = rep) # finds max Bn
        maxBn <- -oBn$minFobj

        n1 <- length(oBn$grupo1)
        is.group <- t.Bnbonf(maxBn, n, n1, alpha, bootB = bootB, bootB1 = bootB1)

        if (is.group$significant == TRUE) {
            clust <- oBn$grupo1
            p <- is.group$p.value
        }
        else {
            #1-perform otimization with n1=1
            #2-perform restricted otimization with central group sizes
            #3-Test for largest Bn feom 1 or 2 - if signif - stop!
            # else  repeat 2 e 3.
            n2 <- n - n1
            minsize <- (n1 == floor(n / 2) || n2 == floor(n / 2)) # the group sizes with smaller variance, for which we want to do restricted optimization
            oBn1 <- maxBnsize1(md)
            maxBn1 <- oBn1$maxBn

            if (minsize == TRUE) {
                clust <- oBn1$grupo1
                is.group <- t.Bnbonf(maxBn1, n, 1, alpha, bootB = bootB, bootB1 = bootB1)
                p <- is.group$p.value
                n1 <- 1
            }

            while (minsize == FALSE && is.group$significant == FALSE) {
                # stop when significant or reach min group size
                # center of simle function
                n1m <- min(n1, n2)
                n1_min <- n1m + 1
                n1_max <- n - n1_min

                oBn <- rep.optimBn_restrict(md, n1_max, n1_min, rep = rep)

                maxBn <- -oBn$minFobj # compare to groups of size 1
                maxBnmax <- max(maxBn1, maxBn)
                if (maxBnmax == maxBn1) {
                    n1 <- 1
                    is.group <- t.Bnbonf(maxBnmax, n, n1, alpha, bootB = bootB, bootB1 = bootB1)
                    if (is.group$significant == FALSE) {
                        n1 <- n1_min
                    }
                }
                else {
                    n1 <- length(oBn$grupo1) # dfine n1 for group with smaller Bn
                    is.group <- t.Bnbonf(maxBnmax, n, n1, alpha, bootB = bootB, bootB1 = bootB1)
                }

                n2 <- n - n1
                minsize <- (n1 == floor(n / 2) || n2 == floor(n / 2))
                if (is.group$significant == TRUE) {
                    clust <- oBn$grupo1
                    p <- is.group$p.value
                }
                if (minsize == TRUE && is.group$significant == FALSE) {
                    clust <- oBn1$grupo1
                    is.group <- t.Bnbonf(maxBn1, n, 1, alpha, bootB = bootB, bootB1 = bootB1)
                    p <- is.group$p.value
                    n1 <- 1
                }
            }
        }

        alpha_correct <- alpha / (2 ^ (n - 1) - 1)
    }
    else {
        # if homogeneous
        clust <- 1:n
        p <- is.h$p.MaxTest
        n1 <- n
        alpha_correct <- alpha
    }
    if (n1 != n) {
        if (is.group$significant == FALSE) {
            # rare case: not homogeneous, but didn`t find any signif partition
            n1 <- length(is.h$grupo1)
            ishomo <- FALSE
            p <- is.h$p.MaxTest
            alpha_correct <- alpha
            bn <- Bn(c(length(is.h$grupo1), length(is.h$grupo2)), md[c(is.h$grupo1, is.h$grupo2), c(is.h$grupo1, is.h$grupo2)])
            varBn <- (bn / is.h$minFobj) ^ 2
            clust <- is.h$grupo1
        }
        else {
            #regular case: not homogeneous and finds partition
            ord <- c(clust, (1:n)[-clust])
            bn <- Bn(c(n1, (n - n1)), md[ord, ord])
            varBn <- is.group$varBn
            ishomo <- FALSE
        }

        if (n < 25) {
            message(
                paste(
                    "\t U-Statstics Significance Clustering  \n\nmax(stdBn) =",
                    round(-ResultadoTesteIsHomo$minFobj, digits = 4),
                    "\t homogeneity p-value = ",
                    round(ResultadoTesteIsHomo$p.MaxTest, digits = 4),
                    "\t alpha = ",
                    alpha,
                    "\nReject H0: Significant partition found \nGroup 1:",
                    paste(clust, collapse = " "),
                    "\nGroup 2:",
                    paste(c(1:n)[-clust], collapse = " ")
                )
            )
        }
        else {
            message(
                paste(
                    "\t U-Statstics Significance Clustering  \n\nmax(stdBn) =",
                    round(-ResultadoTesteIsHomo$minFobj, digits = 4),
                    "\t homogeneity p-value = ",
                    round(ResultadoTesteIsHomo$p.MaxTest, digits = 4),
                    "\t alpha = ",
                    alpha,
                    "\nReject H0: Significant partition found: \nGroup 1 of size",
                    n1,
                    "\nGroup 2 of size",
                    n - n1
                )
            )
        }
    }
    else {
        n1 <- length(is.h$grupo1)
        ishomo <- TRUE
        bn <- Bn(c(length(is.h$group1), length(is.h$group2)), md[c(is.h$group1, is.h$group2), c(is.h$group1, is.h$group2)])
        varBn <- (bn / is.h$minFobj) ^ 2

        message(
            paste(
                "\t U-Statstics Significance Clustering  \n\nmax(stdBn) =",
                round(-ResultadoTesteIsHomo$minFobj, digits = 4),
                "\t homogeneity p-value = ",
                round(ResultadoTesteIsHomo$p.MaxTest, digits = 4),
                "\t alpha = ",
                alpha,
                "\nAccept H0: the sample is homogeneous."
            )
        )
    }
    clust2 <- c(1:n)[-clust]
    ans <- list(clust, clust2, p, alpha_correct, n1, ishomo, bn, varBn, ResultadoTesteIsHomo)

    names(ans) <- c("cluster1", "cluster2", "p.value", "alpha_corrected", "n1", "ishomo", "Bn", "varBn", "ishomoResult")

    invisible(return(ans))
}

#############################################################################################
###   Computes objective function for obtimization based on Bn statistic (non standardized)
#############################################################################################
# internal function: only used in Bn optimzation
objBn <- function(assign, mdm) {
    n <- length(assign)
    n1 <- sum(assign == 0)
    n2 <- n - n1

    if (n1 == 0  | n2 == 0) {
        return(Inf)
    }
    else {
        vaux1 <- which(assign == assign[1])
        vaux2 <- c(1:n)[-vaux1]
        n1 <- length(vaux1)
        n2 <- n - n1
        m <- mdm[c(vaux1, vaux2), c(vaux1, vaux2)]
        Bns <- Bn(c(n1, n2), m)
        return(-Bns)
    }
}

#############################################################################################
### Finds the configuration with max Bn among all configurations
#############################################################################################
optimBn <- function(md, itmax = 200, centers = -1, bootB = -1) {
    n <- dim(md)[1] # number of samples

    #creating structures for registring assignments
    it <- 1
    ass <- vector()
    ass_old <- rep(2, n)
    ASS <- matrix(ncol = n, nrow = itmax)   # matrix for optimization history
    Fobj <- vector()

    ### Initializing optimization paramters

    # choose random centers, if not defined in options
    if (centers == -1) {
        centers <- sample(n, 2)
    }

    # inicialize assignments to closest center
    for (i in 1:n) {
        ass[i] <- (md[i, centers[1]]) > (md[i, centers[2]])
    }

    ass

    ASS[1, ] <- ass

    # TRUE belongs to group 2
    #### Start Iteration
    while (it < itmax && !prod(ass == ass_old)) {
        ass_old <- ass
        ord <- sample(n, n)

        for (i in ord) {
            ass[i] <- 0
            f0 <- objBn(ass, md)
            ass[i] <- 1
            f1 <- objBn(ass, md)

            if (f0 < f1) {
                ass[i] <- 0
            }
        }

        Fobj[it] <- objBn(ass, md)

        it <- it + 1
        ASS[it, ] <- ass
    }

    #assemble answer
    ans <- list(which(ass == ass[1]), Fobj, it - 1, ASS[1:(it + 1), ], bootB)
    names(ans) <- c("grupo1", "Fobj", "numIt", "history", "bootB")

    ans
}

#############################################################################################
### Optimization function with multiple staring ponts for local optima
### Finds the configuration with max Bn among all configurations: rep.optimBn
#############################################################################################
rep.optimBn <- function(mdm, rep = 15, bootB = -1) {
    ans <- optimBn(mdm, bootB = bootB)
    Fobj <- vector()
    n <- dim(mdm)[1]
    grupo1 <- list()
    Fobj[1] <- ans$Fobj[length(ans$Fobj)]
    grupo1[[1]] <- ans$grupo1
    bootB <- ans$bootB

    for (i in 2:rep) {
        ans <- optimBn(mdm, bootB = bootB)
        Fobj[i] <- ans$Fobj[length(ans$Fobj)]
        grupo1[[i]] <- ans$grupo1
    }

    minFobj <- min(Fobj)
    g1 <- grupo1[[which(Fobj == min(Fobj))[1]]]
    g2 <- (1:n)[-g1]

    ans <- list(minFobj, g1, g2, Fobj, bootB)
    names(ans) <- list("minFobj", "grupo1", "grupo2", "Fobj", "bootB")

    return(ans)
}



#############################################################################################
### Finds the configuration with max Bn among all configurations with one group of size one
#############################################################################################
maxBnsize1 <- function(mdm) {
    n <- dim(mdm)[1]
    vecBn <- vector()

    for (i in 1:n) {
        vecBn[i] <- Bn(c(1, n - 1), mdm[c(i, c(1:n)[-i]), c(i, c(1:n)[-i])])
    }
    maxBn <- max(vecBn)
    ans <- list(maxBn, which(vecBn == maxBn))
    names(ans) <- c("maxBn", "grupo1")

    return(ans)
}

######################################################################
# OptimBn_resctrict - restricted optimization function used in uclust
######################################################################
optimBn_restrict <- function(md, n1_max, n1_min, itmax = 200) {
    if (n1_max < n1_min) {
        print("ERROR: n1_max  must be larger than n1_min")
    }
    else {
        n <- dim(md)[1] # sample size n

        #set data structures to keep tabs on assignments
        it <- 1
        ass_old <- rep(2, n)
        ASS <- matrix(ncol = n, nrow = (itmax + 1))   # matrix to register optimzation history

        Fobj <- vector()
        ### initialize optimization parameters

        ass <- rep(0, n)
        ass[sample(n, floor(n / 2))] <- 1
        ASS[1, ] <- ass
        # true belongs to group  2

        #### Start iteration
        if (n1_max == n1_min) {
            # case  n1_max==n1_min
            count <- 0
            while (it < itmax && count < max(n / 2, 10)) {
                ass_old <- ass
                ord <- sample(n, n)
                for (i in ord) {
                    v <- (1:n)[-i]
                    j <- sample(v, 1)
                    f0 <- objBn(ass, md)
                    temp_ass <- ass
                    ass[j] <- ass[i]
                    ass[i] <- temp_ass[j]

                    if (sum(ass) <= n1_max && sum(ass) >= n1_min) {
                        f1 <- objBn(ass, md)
                    }
                    else {
                        f1 <- Inf
                    }

                    if (f0 < f1) {
                        ass <- temp_ass
                    }
                }

                Fobj[it] <- objBn(ass, md)

                it <- it + 1
                ASS[it, ] <- ass
                if (prod(ass == ass_old)) {
                    count <- count + 1
                }
                else{
                    count <- 0
                }
            }
        }
        else {
            ################################    case: n1_max different to n1_min
            while (it < itmax) {
                ass_old <- ass
                ord <- sample(n, n)
                for (i in ord) {
                    ass[i] <- 0
                    if (sum(ass) <= n1_max && sum(ass) >= n1_min) {
                        f0 <- objBn(ass, md)
                    }
                    else {
                        f0 <- Inf
                    }

                    ass[i] <- 1
                    if (sum(ass) <= n1_max && sum(ass) >= n1_min) {
                        f1 <- objBn(ass, md)
                    }
                    else {
                        f1 <- Inf
                    }

                    if (f0 < f1) {
                        ass[i] <- 0
                    }
                }

                Fobj[it] <- objBn(ass, md)

                it <- it + 1
                ASS[it, ] <- ass
            }
        }
        #assemble response

        ans <- list(which(ass == ass[1]), Fobj, it - 1, ASS[1:(it + 1), ])
        names(ans) <- c("grupo1", "Fobj", "numIt", "history")

        ans
    }
}

#############################################################################################
### Optimization function with multiple staring ponts for local optima
### Finds the configuration with max Bn among  restricted set: rep.optimBn_restrict
#############################################################################################
rep.optimBn_restrict <- function(mdm, n1_max, n1_min, rep = 15) {
    ans <- optimBn_restrict(mdm, n1_max, n1_min)
    Fobj <- vector()
    n <- dim(mdm)[1]
    grupo1 <- list()
    Fobj[1] <- ans$Fobj[length(ans$Fobj)]
    grupo1[[1]] <- ans$grupo1

    for (i in 2:rep) {
        ans <- optimBn_restrict(mdm, n1_max, n1_min)
        Fobj[i] <- ans$Fobj[length(ans$Fobj)]
        grupo1[[i]] <- ans$grupo1
    }

    minFobj <- min(Fobj)
    g1 <- grupo1[[which(Fobj == min(Fobj))[1]]]
    g2 <- (1:n)[-g1]

    ans <- list(minFobj, g1, g2, Fobj)
    names(ans) <- list("minFobj", "grupo1", "grupo2", "Fobj")

    return(ans)
}

#############################################################################################
### Bonferroni corected U test used for restricted optimization problems
#############################################################################################
t.Bnbonf <- function(Bn, n, n1, alpha, bootB = -1, bootB1 = -1, md = 1, std = FALSE) {
    if (std == FALSE) {
        # if Bn is not standardized (std=FALSE)
        Cn <- vector()
        numB <- 2000
        if (n1 == 1 || n1 == (n - 1)) {
            if (bootB1 == -1) {
                bootB1 <- boot_sigma1(c(1, (n - 1)), md)
            }
            varBn <- bootB1
        }
        else {
            if (bootB == -1) {
                bootB <- boot_sigma(c(floor(n / 2), (n - floor(n / 2))), numB, md) # computes var of Bn with group sizes c(floor(n/2),(n-floor(n/2))
            }
            n1b <- floor(n / 2)
            n2 <- n - n1b
            Cnd <- (((n1b * n2) / (n * (n - 1)) ^ 2) * (2 * n ^ 2 - 6 * n + 4) / ((n1b - 1) * (n2 - 1)))

            n1b <- n1
            n2 <- n - n1b
            Cnn <- (((n1b * n2) / (n * (n - 1)) ^ 2) * (2 * n ^ 2 - 6 * n + 4) / ((n1b - 1) * (n2 - 1)))

            varBn <- Cnn * bootB / Cnd #variance of Bn for test
        }
        p <- pnorm(Bn / sqrt(varBn), lower.tail = FALSE)
    }
    else {
        p <- pnorm(Bn, lower.tail = FALSE) # standardized Bn is used here. Nesse caso std=TRUE
    }
    significant <- (p < (alpha / (2 ^ (n - 1) - 1)))
    ans <- list(p, significant, varBn)
    names(ans) <- c("p.value", "significant", "varBn")

    return(ans)
}
