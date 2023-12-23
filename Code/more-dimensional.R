##Function for the used kernel
gaussian_cop_kernel <- function(x, y, h) {
  b <- 1-h^2
  1/(sqrt(1-b^2))*exp((-(b^2*(qnorm(x))^2-2*b*qnorm(x)*qnorm(y)+b^2*(qnorm(y))^2)) / (2*(1-b^2)))
}

#function to estimate points with gaussian copula kernels
Gaussian_est_nd <- function(x = NULL, y, k = NULL, h = NULL, bw_each_l = NULL){
  ######Missing: Input checks
  # n is dimension of y
  n <- nrow(y)
  d <- ncol(y)
  # scheme to generate k grid points for d dimensions
  if(is.null(x)) {
    n_d <- rep(0, length.out = d)
    grid_d <- vector(mode = "list", length = d)
    for(i in 1:d) {
      n_d[i] <- k^(1/d)
      grid_d[[i]] <- seq(0.00001, 0.99999, length.out = n_d[i])
    }
    x <- do.call(expand.grid, grid_d)
  }
  
 # if(is.null(k)) {
    #number of grid points
    k <- nrow(x)
  #}
  
  #Bandwidth
  if(is.null(h)) {
    bw_each <- vector(mode = "numeric", length = d)
    #same bandwidth for each dimension
    if(bw_each_l == FALSE) {
      mu <- mean(qnorm(y[,1]))
      sigma <- sd(qnorm(y[,1]))
      h <- sigma*((2*mu^2*sigma^2+3*(1-sigma^2)^2)^(-1/5))*n^{-1/5}
      bw_each <- rep(h, length.out = d)
    }
    #different bandwidths for each dimension
    else {
      for(i in 1:d) {
        mu <- mean(qnorm(y[,i]))
        sigma <- sd(qnorm(y[,i]))
        h <- sigma*((2*mu^2*sigma^2+3*(1-sigma^2)^2)^(-1/5))*n^{-1/5}
        bw_each[i] <- h
      }
    }
  }
  else {
    bw_each <- rep(h, length.out = d)
  }
  
  #Matrix with a row for each data-point and a column for each grid-point
  Kgaussian <- matrix(0, nrow = n, ncol = k)
  vals_d <- rep(0, length.out = 3)
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
  results_m <- cbind(x, fhat)
  results <- list(ind_est = Kgaussian, avg_est = results_m)
  
  #class results?
  class (results) <- c('list', 'GC')
  
  results
}

y1 <-c(0.1, 0.16, 0.25, 0.34, 0.42, 0.5, 0.58, 0.66, 0.75, 0.84, 0.9)
y2 <- c(0.1, 0.16, 0.25, 0.34, 0.42, 0.5, 0.58, 0.66, 0.75, 0.84, 0.9)
y3 <-c(0.4, 0.26, 0.35, 0.84, 0.22, 0.45, 0.5, 0.36, 0.89, 0.24, 0.7)
y4 <-c(0.2, 0.26, 0.21, 0.14, 0.45, 0.9, 0.7, 0.12, 0.54, 0.95, 0.5)
y <- cbind(y1, y2, y3, y4)
h <- sqrt(0.1)
k <- 10000


#n_x1 <- 10  # Number of points along the x-axis
#n_x2 <- 10  # Number of points along the y-axis
#n_x3 <- 10
#n_x4 <- 10

# Create a 2D grid on the unit square
#x1_grid <- seq(0.0001, 0.9999, length.out = n_x1)
#x2_grid <- seq(0.0001, 0.9999, length.out = n_x2)
#x3_grid <- seq(0.0001, 0.9999, length.out = n_x2)
#x4_grid <- seq(0.0001, 0.9999, length.out = n_x2)
#x <- expand.grid(x1 = x1_grid, x2 = x2_grid, x3 = x3_grid, x4 = x4_grid)

Gaussian_est_nd(y = y, k=k, bw_each_l = FALSE)





###Test on copula
# Test for a gaussian copula
set.seed(123456789)
library(copula)

# Create a Gaussian copula with rho = 0.5
copula <- normalCopula(param = 0.3, dim = 5)

# Generate random data from the copula
copula_data <- rCopula(1000, copula)
k <- 1000
h <- (1-sqrt(0.3))
Gaussian_est_nd(y = copula_data, k=k, bw_each_l = FALSE)
