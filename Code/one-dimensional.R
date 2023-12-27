#####################################################Jones-Henderson estimator#####################
##This Code File includes an implementation of the Jones-Henderson estimator, as well as the code
##for some visualizations of the Jones-Henderson estimator


##Function that calculates the kernel function (bivariate Gaussian copula density) over a 
##specific data-point (x,y) and with the bandwidth h.
##Input: x,y: numerics that the kernel is evaluated at
##        h: a numeric specifying the bandwidth
##Output: numeric result of the kernel function being calculated at a specific point and with a specific
##        bandwidth
gaussian_cop_kernel <- function(x, y, h) {
  b <- 1-h^2
  1/(sqrt(1-b^2))*exp((-(b^2*(qnorm(x))^2-2*b*qnorm(x)*qnorm(y)+b^2*(qnorm(y))^2)) / (2*(1-b^2)))
}

####One dimensional gaussian copula kernel estimation
##This function estimates a density with the Jones-Henderson estimator using the previously defined helper
##function. It calculates the estimate on the grid x, by using the data-points y and the bandwidth h.
##The function defaults an equally spaced grid of 101 values if no grid is specified and uses the rule of thumb
##bandwidth, presented by Jones and Henderson, if no bandwidth was provided.
##Inputs: x: numeric: grid on which the densities should be estimated.
##        y: numeric: data-points, that are used to estimate the density.
##        k: numeric: number of grid-points. If not specified, it will take the length of x.
##        h: numeric: bandwidth parameter. If not specified, Jones-Henderson's rule of thumb will be applied.
##Outputs: A list of results containing ind_est and avg_est:
##        ind_est: nxm Matrix with n the number of data-points and m the number of grid points.
##                 Kernel-values for each data values.
##        avg_est: nx2 matrix with the first column specifying the grid point and the second column
##                the estimated density value

gauss_cop_est_1d <- function(x = NULL, y, k = NULL, h = NULL){
  n <- length(y)
  #generate 101 equally spaced grid points
  if(is.null(x)) {
    x <- seq(from = 0, to = 1, by = 0.01)
  }
  #number of grid points
  if(is.null(k)) {
    k <- length(y)
  }
  #Bandwidth - calculates the optimal bandwidth according to the "rule of thumb" by Jones and Henderson
  #for the data-values y
   if(is.null(h)) {
     mu <- mean(qnorm(y))
     sigma <- sd(qnorm(y))
     h <- sigma*((2*mu^2*sigma^2+3*(1-sigma^2)^2)^(-1/5))*n^{-1/5}
    }
  #Matrix with rows for each data-value and columns for each grid-value
  Kgaussian <- matrix(0, nrow = n, ncol = k)
  for(i in 1:n) {
    for(j in 1:k) {
      #Kernel for each data-value evaluated at each grid-value
      Kgaussian[i, j] <-  gaussian_cop_kernel(x[j], y[i], h)
    }
  }
  #Mean over the data-values
  fhat <- colMeans(Kgaussian)
  #Store in list results$x: grid-values, results$y: estimated values,
  #results$individual: Kgaussian Matrix - Kernel for each data value
  results_m <- cbind(x, fhat)
  results <- list(ind_est = Kgaussian, avg_est = results_m)
  
  #class results
  class (results) <- c('list', 'GC')
  
  #print results
  results
}
 


 ## Function on some data-points:
 y <- c(0.1, 0.16, 0.25, 0.34, 0.42, 0.5, 0.58, 0.66, 0.75, 0.84, 0.9)
 k <- 101
 xx <- seq(from = 0, to = 1, by = 0.01) #grid
 #with specified bandwidth
 h <- sqrt(0.1)
 result2 <- gauss_cop_est_1d(x = xx, y = y, k = 101, h = h)
 #with "optimal" bandwidth
 result2 <- gauss_cop_est_1d(x = xx, y = y, k = 101)


# Create the plot from the thesis
plot(result2$avg_est[, "x"], result2$avg_est[, "fhat"], type = 'l', col = 'black', xlab = 'X', ylab = 'Kernel Density', ylim=c(0,5), xlim=c(0,1))
for (i in 1:length(y)) {
  lines(xx, result2$ind_est[i,], col = 'lightgrey', lty = 2)
}

