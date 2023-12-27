#####################################################Multivariate Gaussian copula kernel estimator#####################
##This Code File includes an implementation of the Multivariate Gaussian copula product kernel estimator.
##It is not used in the thesis itself but provided for the sake of completeness.

##Function that calculates the kernel function (bivariate Gaussian copula density) over a 
##specific data-point (x,y) and with the bandwidth h.
##Input: x,y: numerics, that the kernel is evaluated at
##       h a numeric specifying the bandwidth
##Output: numeric result of the kernel function being calculated at a specific point and with a specific
##        bandwidth
gaussian_cop_kernel <- function(x, y, h) {
  b <- 1-h^2
  1/(sqrt(1-b^2))*exp((-(b^2*(qnorm(x))^2-2*b*qnorm(x)*qnorm(y)+b^2*(qnorm(y))^2)) / (2*(1-b^2)))
}

####Multi-dimensional gaussian copula kernel estimation
##This function estimates a density with the Gaussian copula product estimator using the previously defined helper
##function. It calculates the estimate on the grid x, by using the data-points y and the bandwidth h.
##The function defaults an equally spaced grid of k values in each dimension, if no grid is specified. 
##It uses the same bandwidth for all dimensions.
##Inputs: x: kxd matrix with k the number of grid-points and d the dimension of the data-points: grid on 
##           which the densities should be estimated
##        y: nxd matrix: data-points, that are used to estimate the density
##        k: numeric: number of grid-points. If not specified, it will take the nuumber of rows of x.
##        h: numeric: bandwidth parameter. It will be used as the bandwidth in all dimensions.
##Outputs: A list of results containing ind_est and avg_est:
##        ind_est: nxk Matrix with n the number of data-points and k the number of grid points
##                Kernel-values for each data values.
##        avg_est: nx(d+1) matrix with the first d columns specifying each grid point and the last column
##                the estimated density value

gauss_cop_est_nd <- function(x=NULL, y, k=NULL, h){
  # dimension of y
  n <- nrow(y)
  d <- ncol(y)
  #Bandwidth
  bw_each <- rep(h, length.out = d)
  print(bw_each)
  
  #grid, in case x is specified as NULL
  if(is.null(x)) {
    grid_d <- vector(mode = "list", length = d)
  for(i in 1:d) {
    grid_d[[i]] <- seq(0.01, 0.99, length.out = k)
  }
  x <- do.call(expand.grid, grid_d)
  }
  
  
    #number of grid points
    k <- nrow(x)
  
  #Matrix with a row for each data-point and a column for each grid-point
  Kgaussian <- matrix(0, nrow = n, ncol = k)
  vals_d <- rep(0, length.out = d)
  
  for(i in 1:n) { 
    for(j in 1:k) {
      for(l in 1:d) {
      #for the index of each data-point and each grid-point, calculate the univariate Kernels for each dim;
      #then take the product of diff. dimensions; the result is again a density value for each data-point and grid-value
      vals_d[l] <- gaussian_cop_kernel(x[j,l], y[i,l], bw_each[l])
      Kgaussian[i, j] <- prod(vals_d)
      }
  }
  }
  #Mean over the data-points 
  fhat <- colMeans(Kgaussian)
  
  #Store in list results$x: grid-values, results$y: estimated values,
  #results$individual: Kgaussian Matrix - Kernel for each data value
  #avg_est is gives the grid-points by the first d dimensions and the estimated values in the last column
  results_m <- cbind(x, fhat)
  results <- list(ind_est = Kgaussian, avg_est = results_m)
  
  #class results
  class (results) <- c('list', 'GC')
  
  results
}


###Test on copula
# Test for a gaussian copula
set.seed(123456789)
library(copula)

# Create a Gaussian copula with rho = 0.3
copula <- normalCopula(param = 0.3, dim = 5)

# Generate random data from the copula
copula_data <- rCopula(1000, copula)
k <- 3
h <- (1-sqrt(0.3))
gauss_cop_est_nd(y = copula_data, h=h, k=k)

