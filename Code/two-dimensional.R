###########################Gaussian copula kernel estimation in 2 dimensions

##Function for the used kernel
gaussian_cop_kernel <- function(x, y, h) {
  b <- 1-h^2
  1/(sqrt(1-b^2))*exp((-(b^2*(qnorm(x))^2-2*b*qnorm(x)*qnorm(y)+b^2*(qnorm(y))^2)) / (2*(1-b^2)))
}

library(tidyverse)
Gaussian_est_2d <- function(x = NULL, y, k = NULL, h = NULL, bw_each_l = NULL){
  ######Missing: Input checks
  # n is dimension of y
  n <- nrow(y)
  # scheme to generate 1000 grid points on the unit square
  if(is.null(x)) {
    # Create a 2D grid on the unit square
    x1_grid <- seq(0.01, 0.99, length.out = k)
    x2_grid <- seq(0.01, 0.99, length.out = k)
    x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
    }
  
  # number of grid points
  p <- nrow(x)
  
  # Bandwidth
  if(is.null(h)) {
    bw_each <- vector(mode = "numeric", length = 2)
    # same bandwidth for each dimension
    if(bw_each_l == FALSE) {
    mu <- mean(qnorm(y[,1]))
    sigma <- sd(qnorm(y[,1]))
    h <- sigma*((2*mu^2*sigma^2+3*(1-sigma^2)^2)^(-1/5))*n^{-1/5}
    bw_each <- c(h, h)
    }
    # different bandwidths for each dimension
    else {
      for(i in 1:2) {
        mu <- mean(qnorm(y[,i]))
        sigma <- sd(qnorm(y[,i]))
        h <- sigma*((2*mu^2*sigma^2+3*(1-sigma^2)^2)^(-1/5))*n^{-1/5}
        bw_each[i] <- h
      }
    }
  }
  else {
    bw_each <- c(h, h)
  }
  print(bw_each)
  
  # Matrix with a row for each data-point and a column for each grid-point
  Kgaussian <- matrix(0, nrow = n, ncol = p)
  for(i in 1:n) {
    for(j in 1:p) {
      # for the index of each data-point and each grid-point, calculate a product of the 
      # univariate Kernels. The result is a density value for each data-point and grid-value
      Kgaussian[i, j] <- prod(gaussian_cop_kernel(x[j,1], y[i,1], bw_each[1]), 
                              gaussian_cop_kernel(x[j,2], y[i,2], bw_each[2]))
    }
  }
  
  # Mean over the data-points
  fhat <- colMeans(Kgaussian)
  
  # Store in list results$x: grid-values, results$y: estimated values,
  # results$individual: Kgaussian Matrix - Kernel for each data value
  results_m <- cbind(x, fhat)
  # Make Matrix with each grid point as a cell
  vals_grid <- as.data.frame.matrix(results_m) %>%
    pivot_wider(names_from = x2, values_from = fhat) %>%
    column_to_rownames(var = "x1") %>%
    as.matrix()
  results <- list(ind_est = Kgaussian, avg_est = results_m, vals_grid = vals_grid)
  
  # class results?
  class (results) <- c('list', 'GC')
  
  results
}


##Data to test on
#y1 <- c(0.1, 0.16, 0.25, 0.34, 0.42, 0.5, 0.58, 0.66, 0.75, 0.84, 0.9)
#y2 <-c(0.1, 0.16, 0.25, 0.34, 0.42, 0.5, 0.58, 0.66, 0.75, 0.84, 0.9)
#y <- cbind(y1, y2)
#h <- sqrt(0.1)
#k <- 100
#n_x <- 10 
#n_y <- 10
# Create a 2D grid on the unit square
#x1_grid <- seq(0.0001, 0.9999, length.out = n_x)
#x2_grid <- seq(0.0001, 0.9999, length.out = n_y)
#x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
#result4 <- Gaussian_est_2d(x = x, y = y, k = k, h = h)


##Otther data:
#y1 <- runif(11)
#y2 <- runif(11)
#y <- cbind(y1, y2)
#h <- sqrt(0.1)
#k <- 100
#result <- Gaussian_est_2d(x = x, y = y, k = 100, h = h)



######################################################################
# Test for a gaussian copula
set.seed(123456789)
library(copula)

# Create a Gaussian copula
copula <- normalCopula(param = 0.9, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(1000, copula)

# Create a 2D grid on the unit square
#x1_grid <- seq(0.000000001, 0.999999999, length.out = 30)
#x2_grid <- seq(0.000000001, 0.999999999, length.out = 30)
x1_grid <- seq(0.01, 0.99, length.out = 30)
x2_grid <- seq(0.01, 0.99, length.out = 30)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
#h <- (1-sqrt(0.3))
#h <- 0.1

# Apply KDE to estimate the density with fixed bandwidth
kde_result <- Gaussian_est_2d(x = x, y = copula_data, h = 0.3)

#Apply KDE to estimate the density with optimal bandwidth for one dimension
kde_result <- Gaussian_est_2d(x = x, y = copula_data, bw_each_l = FALSE)

#Apply KDE to estimate the density with optimal bandwidth for each dimension
kde_result <- Gaussian_est_2d(x = x, y = copula_data, bw_each_l = TRUE)




########Plots of estimated density#############
#Surface Plot
persp(x = x1_grid, y = x2_grid, z = kde_result$vals_grid, col = "lightgrey", xlim = c(0,1), ticktype = "detailed", theta = -30, phi = 20, border = "blue")
#Contourplot
contour(x = x1_grid, y = x2_grid, z = kde_result$vals_grid)
#Normal contour plot?
#kde_transf_1 <- qnorm(x1_grid)
#kde_transf_2 <- qnorm(x2_grid)

########Plots of true density#############
#Scatterplot
plot(copula_data)
#Surface plot
persp(copula, dCopula, col = "purple", theta = -30, phi = 20)
#Contourplot
contour(copula, dCopula)
#Normal contour plot
#...........




######################################################################
# Test for a frank copula
set.seed(123456789)
library(copula)

# Create a frank copula
copula <- frankCopula(param = 0.3, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(1000, copula)

# Create a 2D grid on the unit square
#x1_grid <- seq(0.000000001, 0.999999999, length.out = 30)
#x2_grid <- seq(0.000000001, 0.999999999, length.out = 30)
x1_grid <- seq(0.01, 0.99, length.out = 30)
x2_grid <- seq(0.01, 0.99, length.out = 30)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
#h <- (1-sqrt(0.7))
#h <- 0.1

# Apply KDE to estimate the density with fixed bandwidth
kde_result <- Gaussian_est_2d(x = x, y = copula_data, h = h)

#Apply KDE to estimate the density with optimal bandwidth for one dimension
kde_result <- Gaussian_est_2d(x = x, y = copula_data, bw_each_l = FALSE)

#Apply KDE to estimate the density with optimal bandwidth for each dimension
kde_result <- Gaussian_est_2d(x = x, y = copula_data, bw_each_l = TRUE)




########Plots of estimated density#############
#Surface Plot
persp(x = x1_grid, y = x2_grid, z = kde_result$vals_grid, col = "lightgrey", xlim = c(0,1), ticktype = "detailed", theta = -30, phi = 20, border = "blue")
#Contourplot
contour(x = x1_grid, y = x2_grid, z = kde_result$vals_grid)
#Normal contour plot??
#kde_transf_1 <- qnorm(kde_result$avg_est[,1])
#kde_transf_2 <- qnorm(kde_result$avg_est[,2])


########Plots of true density#############
#Scatterplot
plot(copula_data)
#Surface plot
persp(copula, dCopula, col = "purple", theta = -30, phi = 20)
#Contourplot
contour(copula, dCopula)
#Normal contour plot
#...........



####################################Other packages for visualization##################################
##Plot surface
#library(plotly)
#plot_ly(x = seq(from = 0, to = 1, by = 0.00001), y = seq(from = 0, to = 1, by = 0.00001), z = as.matrix(kde_result$avg_est)) %>%
# add_surface() # %>%

#library(rgl)
#plot3d(x = kde_result$avg_est[,1], y = kde_result$avg_est[,2], z = as.matrix(kde_result$avg_est), xlim = c(0,1), ylim = c(0,1))

#library(plot3D)
#surf3D(x = as.matrix(kde_result$avg_est), y = as.matrix(kde_result$avg_est), z = as.matrix(kde_result$avg_est), colvar = as.matrix(kde_result$avg_est), xlim = c(0,1), ylim = c(0,1))
#persp( x = x1_grid, y = x2_grid, z = vals_grid, col = "purple")

