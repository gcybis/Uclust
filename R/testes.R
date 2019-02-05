# x=matrix(rnorm(100000),nrow=100)  #creating homogeneous dataset
# is_homo(data=x)
#
#
#
#
# md=as.matrix(dist(x)^2)
# is_homo(md)
# ans=uclust(md)
#
#
# y=x
# x[1,]=x[1,]+0.3
# is_homo(data=y)
#
# x[1:30,]=x[1:30,]+0.3
# md=as.matrix(dist(x)^2)
# is_homo(md)
# is_homo(data=x)
#
#
#
#
#
# #### Para o help do is_homo
#
# x=matrix(rnorm(1000000),nrow=100)  #creating homogeneous Gaussian dataset
# res=is_homo(data=x)
#
# x[1:30,]=x[1:30,]+0.15   #Heterogeneous dataset (first 30 samples have different mean)
# res= is_homo(data=x)
#
# md=as.matrix(dist(x)^2)   #squared Euclidean distances for the same data
# res= is_homo(md)
#
# # Multidimensional sacling plot of distance matrix
# fit <- cmdscale(md, eig = TRUE, k = 2)
# x <- fit$points[, 1]
# y <- fit$points[, 2]
# plot(x,y, main=paste("Homogeneity test: p-value =",res$p.MaxTest))
#
#
#' #
#'
#' @examples
#' x = matrix(rnorm(100000),nrow=50)  #creating homogeneous Gaussian dataset
#' res = uhclust(data=x)
#'
#'
#' x[1:30,] = x[1:30,]+0.7   #Heterogeneous dataset
#' x[1:10,] = x[1:10,]+0.4
#' res = uhclust(data=x)
#' res$groups
#'



