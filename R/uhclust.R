#########################################################
###   Main Uhclust function
#########################################################
#' U-statistic based significance hierarchical clustering
#'
#' Hierarchical clustering method that partitions the data only when these partitions are statistically significant.
#'
#'
#' This is the significance hierarchical clustering procedure of Valk and Cybis (2018). The data are
#' repeatedly partitioned into two subgroups, through function \code{uclust}, according to a hierarchical scheme.
#' The procedure stops when resulting subgroups are homogeneous or have fewer than 3 elements.
#' This function should be used in high dimension small sample size settings.
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
#'  See also \code{is_homo}, \code{uclust} and \code{Utest_class}.
#'
#' @param md Matrix of squared Euclidean distances between all data points.
#' @param data Data matrix. Each row represents an observation.
#' @param alpha Significance level.
#' @param rep Number of times to repeat optimization procedures. Important for problems with
#' multiple optima.
#' @param plot Logical, \code{TRUE} if p-value annotated dendrogram should be plotted.
#' @return  Returns an object of class \code{hclust} with three additional attribute arrays:\describe{
#' \item{Pvalues}{ P-values from uclust for the final data partition at each node of the dendrogram. This
#' array is in the same order of \code{height}, and only contains values for tests that were performed.}
#' \item{alpha}{ Bonferroni corrected significance levels for uclust for the data partitions at each node
#' of the dendrogram. This array is in the same order of \code{height}, and only contains values for tests that were performed.}
#' \item{groups}{ Final group assignments.}
#' }
#'
#'
#' @examples
#'
#' x = matrix(rnorm(100000),nrow=50)  #creating homogeneous Gaussian dataset
#' res = uhclust(data=x)
#'
#'
#' x[1:30,] = x[1:30,]+0.7   #Heterogeneous dataset
#' x[1:10,] = x[1:10,]+0.4
#' res = uhclust(data=x)
#' res$groups
#'
#' @export
uhclust <- function(md = NULL, data = NULL, alpha = 0.05, rep = 15, plot = TRUE) {




    if (is.null(md)) {
        # Computing data matrix if one is not provided
        if (is.null(data)) {
            stop("No data provided")
        }
        md <- as.matrix(dist(data) ^ 2)
        if (is.null(rownames(data))) {
            rownames(md) <- 1:(dim(data)[1])
        }
    }

    if (class(md) != "matrix") {
        stop("md is not of class matrix")
    }

  if (is.null(rownames(md))) {
    rownames(md) <- 1:(dim(md)[1])
  }




    labels <- rownames(md)

    n <- dim(md)[1]
    grupo <- list()
    grupo[[1]] <- labels
    clusters <- list()
    pvalues <- vector()
    alphas <- vector()   # corrected alphas coresponding to each entry in pvalues
    taux <- vector() # Logical, whether uclust should should be applied
    taux[1] <- TRUE

    merge <- matrix(rep(0), ncol = 2, nrow = (n - 1))   # cladogram structure (see class uclust)
    height <- vector()                         # corresponding branch lengths (see class uclust)
    ind.grupo <- 1

    ########################################################################################
    ########### Test for the first division ################################################
    ########################################################################################

    ht <- suppressMessages(uclust(md, alpha = alpha))
    pvalues[n - 1] <- ht$ishomoResult$p.MaxTest
    alphas[n - 1] <- alpha
    height[n - 1] <- 1
    if (ht$p.value > alpha) {
        # if whole sample is homogeneous
        groups <- rep(1, n)
        ans <- list("homogeneous", ht$p.value, alpha, groups)
        names(ans) <- c("homogeneous", "pvalue", "alpha", "groups")

        message(
            paste(
                "\t Hierarachical significance clustering:\n\n",
                "Homogeneity p-value = ",
                round(ht$p.value, digits = 4),
                "\t alpha = ",
                alpha,
                "\n Accept H0: the sample is homogeneous."
            )
        )
    }
    else {
        #otherwise, do first division
        #organize first component
        ind.aux1 <- ind.grupo + 1

        ######## first branch of division #######
        if (length(ht$cluster1) == 1) {
            # if group has size one
            merge[n - 1, 1] <- -ht$cluster1
        }else {
            grupo[[ind.aux1]] <- labels[ht$cluster1]
            merge[n - 1, 1] <- n - ind.aux1                       # only makes sense if groups are larger than 3
            taux[ind.aux1] <- (length(grupo[[ind.aux1]]) > 3)  # logical, should this subgroup be tested (more than 3 elements)
            height[n - ind.aux1] <- 2
        }
        # repetes procedure for second component
        ind.aux2 <- ind.grupo + length(ht$cluster1)
        grupo[[ind.aux2]] <- labels[ht$cluster2]

        ######## second branch of division #######
        if (length(ht$cluster2) == 1) {
            merge[n - 1, 2] <- -ht$cluster2
        }
        else {
            merge[n - 1, 2] <- n - ind.aux2
            taux[ind.aux2] <- (length(grupo[[ind.aux2]]) > 3)
            height[n - ind.aux2] <- 2
        }

        ########################################################################################
        ##############  continues algorithm for the other divisions  ###########################
        ########################################################################################

        if (ht$p.value < alpha) {
            ind.grupo <- 2
        }
        while (ind.grupo < n) {
            ### If division should be performes (not homogeneous and n>3) ####
            if (taux[ind.grupo]) {
                indx.esses <- match(grupo[[ind.grupo]], labels) ## finds indices corresponding to group labels
                alpha_star <- alpha * (length(indx.esses) - 1) / (n - 1) #  Kimes alpha correction (2017)
                ht <- suppressMessages(uclust(md[indx.esses, indx.esses], alpha = alpha_star))

                pvalues[n - ind.grupo] <- ht$ishomoResult$p.MaxTest
                alphas[n - ind.grupo] <- alpha_star

                if (ht$p.value > ht$alpha_corrected) {
                    #if homogeneous
                    merge[n - ind.grupo, 1] <- -indx.esses[ht$cluster1[1]]
                    merge[n - ind.grupo, 2] <- n - (ind.grupo + 1)
                    grupo[[ind.grupo + 1]] <- labels[indx.esses[ht$cluster1[-1]]]
                    taux[ind.grupo + 1] <- FALSE
                    height[n - (ind.grupo)] <- 0.1   ####test
                    height[n - (ind.grupo + 1)] <- 0.1
                    ind.grupo <- ind.grupo + 1
                }
                else {
                    #if group not homogeneous
                    ######## first division branch #######
                    ind.aux1 <- (ind.grupo + 1) #organize group 1

                    if (length(ht$cluster1) == 1) {
                        #does group 1 have only one component?
                        merge[n - ind.grupo, 1] <- -indx.esses[ht$cluster1]
                    }
                    else {
                        #if group one has more than one component
                        merge[n - ind.grupo, 1] <- n - ind.aux1
                        height[n - ind.aux1] <- height[n - ind.grupo] + 1
                        grupo[[ind.aux1]] <- labels[indx.esses[ht$cluster1]]
                        taux[ind.aux1] <- (length(grupo[[ind.aux1]]) > 3)  #logical, should this subgroup be tested (more than 3 elements)
                    }
                    ######## second division branch #######
                    ind.aux2 <- ind.grupo + length(ht$cluster1) #organize group 2
                    if (length(ht$cluster2) == 1) {
                        #does group 2 have only one component?
                        merge[n - ind.grupo, 2] <- -indx.esses[ht$cluster2]
                    }
                    else {
                        merge[n - ind.grupo, 2] <- n - ind.aux2
                        height[n - ind.aux2] <- height[n - ind.grupo] + 1
                        grupo[[ind.aux2]] <- labels[indx.esses[ht$cluster2]]
                        taux[ind.aux2] <- (length(grupo[[ind.aux2]]) > 3)  ##logical, should this subgroup be tested (more than 3 elements)
                    }
                    ind.grupo <- ind.grupo + 1
                }
                ########   If the test should not be performed  #########

            }
            else {
                ng <- length(grupo[[ind.grupo]])
                while (ng > 1) {
                    indx.esses <- match(grupo[[ind.grupo]], labels)
                    if (ng == 2) {
                        height[n - (ind.grupo)] <- 0.1
                        merge[n - ind.grupo, 1] <- -match(grupo[[ind.grupo]][1], labels)
                        merge[n - ind.grupo, 2] <- -match(grupo[[ind.grupo]][2], labels)
                        taux[ind.grupo] <- FALSE
                        ind.grupo <- ind.grupo + 1
                        ng <- 1
                    }
                    else {
                        merge[n - ind.grupo, 1] <- -match(grupo[[ind.grupo]][1], labels)
                        merge[n - ind.grupo, 2] <- n - (ind.grupo + 1)
                        height[n - (ind.grupo)] <- 0.1
                        height[n - (ind.grupo + 1)] <- 0.1
                        grupo[[ind.grupo + 1]] <- grupo[[ind.grupo]][-1]
                        ind.grupo <- ind.grupo + 1
                        ng <- ng - 1
                        taux[ind.grupo] <- FALSE
                    }
                }
            }
        }

        ########################################################################################
        ##############  Organize answer and plot         #######################################
        ########################################################################################
        hh <- max(height)
        height <- (height != 0.1) * (1 + hh - height) + (height == 0.1) * 0.3
        order <- c(t(merge))
        order <- -order[order < 0]


        ans <- list(merge, height, labels, order, pvalues, alphas)
        names(ans) <- c("merge", "height", "labels", "order", "Pvalues", "alpha")
        class(ans) <- "hclust"
        hcd <- as.dendrogram(ans)
        groups <- dendextend::cutree_1h.dendrogram(hcd, h = 0.5)
        ans[[7]] <- groups
        names(ans) <- c("merge", "height", "labels", "order", "Pvalues", "alpha", "groups")

        if (plot)
            plot_uhclust(ans)

        ans <- as.hclust(ans)
    }

    invisible(ans)
}


###################################################################
# plot function
##################################################################
#' Plot function for the result of uhclust
#'
#' This function plots the p-value annotated dendrogram resulting from \code{uhclust}
#'
#' @param uhclust Result from \code{uhclust}
#' @param pvalues_cex Graphical parameter for p-value font size.
#' @param pvalues_dx Graphical parameter for p-value position shift on x axis.
#' @param pvalues_dy Graphical parameter for p-value position shift on y axis.
#' @param print_pvalues Logical. Should the p-values be printed?
#'
#'
#' @examples
#'
#' x = matrix(rnorm(100000),nrow=50)
#' x[1:35,] = x[1:35,]+0.7
#' x[1:15,] = x[1:15,]+0.4
#' res = uhclust(data=x, plot=FALSE)
#' plot_uhclust(res)
#'
#' @import dendextend
#' @export
plot_uhclust <- function(uhclust, pvalues_cex = 0.8, pvalues_dx = 2, pvalues_dy = 0.08, print_pvalues = TRUE) {
    dend <- as.dendrogram(uhclust, center = T)



    #Cluster colors
    c <- uhclust$groups
    cols <- sample(rainbow(max(c), s = 0.5, start = 0.3, end = 1))
    dendextend::labels_colors(dend) <- cols[c[labels(dend)]]

    #Plot P-values
    xy <- dend %>% get_nodes_xy()
    xy <-  get_nodes_xy(dend)
     is_internal_node <- is.na(dend %>% get_nodes_attr("leaf"))
    xy <- xy[is_internal_node, ]

    node_Lab <- (paste("p=", format(uhclust$Pvalues, digits = 2), " (", format(uhclust$alpha, digits = 2), ")", sep = ""))[length(uhclust$Pvalues):1]
    node_Lab[is.na(uhclust$Pvalues)[length(uhclust$Pvalues):1]] <- ""
    plot(dend, edge.root = T, edgePar = list(col = "grey", lwd = 2), axes = FALSE)
    if (print_pvalues)
        text(xy[, 1] + pvalues_dx, xy[, 2] + pvalues_dy, labels = node_Lab, col = "darkcyan", cex = pvalues_cex)
}
