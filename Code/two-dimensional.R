#####################################################Bivariate Gaussian copula kernel estimator#####################
##This Code File includes an implementation of the Bivariate Gaussian copula kernel estimator, as well as the code
##for some visualizations of the estimator, that are used in the Chapter "Simulations and estimator performance"

#load libraries for this file
library(tidyverse)
library(copula)
library(VineCopula)

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


####Two dimensional gaussian copula kernel estimation
##This function estimates a density with the Gaussian copula product estimator using the previously defined helper
##function. It calculates the estimate on the grid x, by using the data-points y and the bandwidth h.
##The function defaults an equally spaced grid of kxk values if no grid is specified. It uses the same bandwidth for
##both dimensions.
##Inputs: x: kx2 matrix with k the number of grid-points: grid on which the densities should be estimated
##        y: nx2 matrix: data-points, that are used to estimate the density
##        k: numeric: number of grid-points. If not specified, it will take the number of rows of x.
##        h: numeric: bandwidth parameter. It will be used as the bandwidth in both dimensions.
##Outputs: A list of results containing ind_est, avg_est and vals_grid:
##        ind_est: nxk Matrix with n the number of data-points and k the number of grid points
##                Kernel-values for each data values.
##        avg_est: nx3 matrix with the first two columns specifying each grid point and the third column
##                the estimated density value
##        vals_grid: kxk matrix with the modified version of avg_est. Each cell represents a grid-point
##                    and contains the density estimate for this specific grid-point.

 gauss_cop_est_2d <- function(x = NULL, y, k = NULL, h = NULL){
  # n is dimension of y
  n <- nrow(y)
  # if not specified, generate grid points on the unit square
  if(is.null(x)) {
    x1_grid <- seq(0.01, 0.99, length.out = k)
    x2_grid <- seq(0.01, 0.99, length.out = k)
    x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
    }
  
  # number of grid points
  p <- nrow(x)
  
  # Bandwidth
    bw_each <- c(h, h)
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
  
  ##Results include individual estimates of the kernels (ind_est), averaged estimates (avg_est)
  ##with the first two columns specifying the grid-points and the last column the density-estimates.
  ##vals_grid specifies a kxk matrix with each cell representing the estimate of the designated grid-point
  results <- list(ind_est = Kgaussian, avg_est = results_m, vals_grid = vals_grid)
  
  # class results
  class (results) <- c('list', 'GC')
  
  results
}




######################################################################
###################################### Test for a gaussian copula
############Tau = 0.3
set.seed(123456789)
 
 
param <- BiCopTau2Par(family = 1, tau = 0.3)

# Create a Gaussian copula
copula <- normalCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
density_cop <- dCopula(as.matrix(x), copula)

# Get best case bandwidth
ph <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
mse <- vector(mode = "numeric", length = length(ph))
for(i in 1:length(ph)) {
  kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[i])
  mse[i] <- sum((density_cop - kde_result$avg_est[,3])^2)/1225
}
mse
ph[which.min(mse)] #-> 0.3

# Apply KDE to estimate the density with fixed bandwidth
kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[which.min(mse)])
#kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = 0.2)


#############Surface plot
z <- kde_result$vals_grid
kde_result$vals_grid<-ifelse(kde_result$vals_grid > 12, 12, kde_result$vals_grid)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (kde_result$vals_grid[-1, -1] + kde_result$vals_grid[-1, -ncol(z)] + kde_result$vals_grid[-nrow(z), -1] + kde_result$vals_grid[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = kde_result$vals_grid, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)

# ###############Contourplot
# contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
# title(xlab = "u1", line = 2.5, cex.lab= 1.5)
# title(ylab = "u2", line = 2.5, cex.lab= 1.5)
# 

###############Normal contour plot
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*kde_result$vals_grid, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)


############Tau = 0.7
set.seed(123456789)
 
 
param <- BiCopTau2Par(family = 1, tau = 0.7)

# Create a Gaussian copula
copula <- normalCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
density_cop <- dCopula(as.matrix(x), copula)

# Get best case bandwidth
ph <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
mse <- vector(mode = "numeric", length = length(ph))
for(i in 1:length(ph)) {
  kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[i])
  mse[i] <- sum((density_cop - kde_result$avg_est[,3])^2)/1225
}
mse
ph[which.min(mse)] #->0.15

# Apply KDE to estimate the density with fixed bandwidth
kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[which.min(mse)])


#################Scatterplot
#plot(copula_data, xlab = "", ylab = "", cex.axis = 1.2)
#title(xlab = "u1", line = 2.5, cex.lab= 1.5)
#title(ylab = "u2", line = 2.5, cex.lab= 1.5)


#############Surface plot
z <- kde_result$vals_grid
kde_result$vals_grid<-ifelse(kde_result$vals_grid > 12, 12, kde_result$vals_grid)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (kde_result$vals_grid[-1, -1] + kde_result$vals_grid[-1, -ncol(z)] + kde_result$vals_grid[-nrow(z), -1] + kde_result$vals_grid[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = kde_result$vals_grid, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)

# ###############Contourplot
# contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
# title(xlab = "u1", line = 2.5, cex.lab= 1.5)
# title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*kde_result$vals_grid, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)




######################################################################
###################################### Test for a Clayton copula
############Tau = 0.3
set.seed(123456789)
 
 
param <- BiCopTau2Par(family = 3, tau = 0.3)

# Create a Clayton copula
copula <- claytonCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
density_cop <- dCopula(as.matrix(x), copula)

# Get best case bandwidth
ph <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
mse <- vector(mode = "numeric", length = length(ph))
for(i in 1:length(ph)) {
  kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[i])
  mse[i] <- sum((density_cop - kde_result$avg_est[,3])^2)/1225
}
mse
ph[which.min(mse)] # -> 0.2

# Apply KDE to estimate the density with fixed bandwidth
kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[which.min(mse)])

#############Surface plot
z <- kde_result$vals_grid
kde_result$vals_grid<-ifelse(kde_result$vals_grid > 12, 12, kde_result$vals_grid)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (kde_result$vals_grid[-1, -1] + kde_result$vals_grid[-1, -ncol(z)] + kde_result$vals_grid[-nrow(z), -1] + kde_result$vals_grid[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = kde_result$vals_grid, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)

# ###############Contourplot
# contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
# title(xlab = "u1", line = 2.5, cex.lab= 1.5)
# title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*kde_result$vals_grid, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)


############Tau = 0.7
set.seed(123456789)
 
 
param <- BiCopTau2Par(family = 3, tau = 0.7)
# Create a Clayton copula
copula <- claytonCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
density_cop <- dCopula(as.matrix(x), copula)

# Get best case bandwidth
ph <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
mse <- vector(mode = "numeric", length = length(ph))
for(i in 1:length(ph)) {
  kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[i])
  mse[i] <- sum((density_cop - kde_result$avg_est[,3])^2)/1225
}
mse
ph[which.min(mse)] #-> 0.05

# Apply KDE to estimate the density with fixed bandwidth
kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[which.min(mse)])


#############Surface plot
z <- kde_result$vals_grid
kde_result$vals_grid<-ifelse(kde_result$vals_grid > 12, 12, kde_result$vals_grid)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (kde_result$vals_grid[-1, -1] + kde_result$vals_grid[-1, -ncol(z)] + kde_result$vals_grid[-nrow(z), -1] + kde_result$vals_grid[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = kde_result$vals_grid, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)

# ###############Contourplot
# contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
# title(xlab = "u1", line = 2.5, cex.lab= 1.5)
# title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
# Create a 2D grid on the unit square
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*kde_result$vals_grid, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)


##############################We also try clayton with tau = 0.7 for a different bandwidth:

kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = 0.2)


#############Surface plot
z <- kde_result$vals_grid
kde_result$vals_grid<-ifelse(kde_result$vals_grid > 12, 12, kde_result$vals_grid)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (kde_result$vals_grid[-1, -1] + kde_result$vals_grid[-1, -ncol(z)] + kde_result$vals_grid[-nrow(z), -1] + kde_result$vals_grid[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = kde_result$vals_grid, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)

# ###############Contourplot
# contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
# title(xlab = "u1", line = 2.5, cex.lab= 1.5)
# title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
# Create a 2D grid on the unit square
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*kde_result$vals_grid, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)












######################################################################
###################################### Test for a Frank copula
############Tau = 0.3
set.seed(123456789)
 
 
param <- BiCopTau2Par(family = 5, tau = 0.3)

# Create a Frank copula
copula <- frankCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
density_cop <- dCopula(as.matrix(x), copula)

# Get best case bandwidth
ph <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
mse <- vector(mode = "numeric", length = length(ph))
for(i in 1:length(ph)) {
  kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[i])
  mse[i] <- sum((density_cop - kde_result$avg_est[,3])^2)/1225
}
mse
ph[which.min(mse)] #-> 0.35

# Apply KDE to estimate the density with fixed bandwidth
kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[which.min(mse)])

#############Surface plot
z <- kde_result$vals_grid
kde_result$vals_grid<-ifelse(kde_result$vals_grid > 12, 12, kde_result$vals_grid)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (kde_result$vals_grid[-1, -1] + kde_result$vals_grid[-1, -ncol(z)] + kde_result$vals_grid[-nrow(z), -1] + kde_result$vals_grid[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = kde_result$vals_grid, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)

# ###############Contourplot
# contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
# title(xlab = "u1", line = 2.5, cex.lab= 1.5)
# title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*kde_result$vals_grid, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)


############Tau = 0.7
set.seed(123456789)
 
 
param <- BiCopTau2Par(family = 5, tau = 0.7)

# Create a Frank copula
copula <- frankCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
density_cop <- dCopula(as.matrix(x), copula)

# Get best case bandwidth
ph <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
mse <- vector(mode = "numeric", length = length(ph))
for(i in 1:length(ph)) {
  kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[i])
  mse[i] <- sum((density_cop - kde_result$avg_est[,3])^2)/1225
}
mse
ph[which.min(mse)] #-> 0.2

# Apply KDE to estimate the density with fixed bandwidth
kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[which.min(mse)])


#############Surface plot
z <- kde_result$vals_grid
kde_result$vals_grid<-ifelse(kde_result$vals_grid > 12, 12, kde_result$vals_grid)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (kde_result$vals_grid[-1, -1] + kde_result$vals_grid[-1, -ncol(z)] + kde_result$vals_grid[-nrow(z), -1] + kde_result$vals_grid[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = kde_result$vals_grid, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)

# ###############Contourplot
# contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
# title(xlab = "u1", line = 2.5, cex.lab= 1.5)
# title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*kde_result$vals_grid, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)











######################################################################
###################################### Test for a Gumbel copula
############Tau = 0.3
set.seed(123456789)
 
 
param <- BiCopTau2Par(family = 4, tau = 0.3)

# Create a Gumbel copula
copula <- gumbelCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
density_cop <- dCopula(as.matrix(x), copula)

# Get best case bandwidth
ph <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
mse <- vector(mode = "numeric", length = length(ph))
for(i in 1:length(ph)) {
  kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[i])
  mse[i] <- sum((density_cop - kde_result$avg_est[,3])^2)/1225
}
mse
ph[which.min(mse)] #-> 0.4

# Apply KDE to estimate the density with fixed bandwidth
kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[which.min(mse)])


#############Surface plot
z <- kde_result$vals_grid
kde_result$vals_grid<-ifelse(kde_result$vals_grid > 12, 12, kde_result$vals_grid)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (kde_result$vals_grid[-1, -1] + kde_result$vals_grid[-1, -ncol(z)] + kde_result$vals_grid[-nrow(z), -1] + kde_result$vals_grid[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = kde_result$vals_grid, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)

# ###############Contourplot
# contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
# title(xlab = "u1", line = 2.5, cex.lab= 1.5)
# title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*kde_result$vals_grid, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)


############Tau = 0.7
set.seed(123456789)
 
 
param <- BiCopTau2Par(family = 4, tau = 0.7)

# Create a Gumbel copula
copula <- gumbelCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
density_cop <- dCopula(as.matrix(x), copula)

# Get best case bandwidth
ph <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
mse <- vector(mode = "numeric", length = length(ph))
for(i in 1:length(ph)) {
  kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[i])
  mse[i] <- sum((density_cop - kde_result$avg_est[,3])^2)/1225
}
mse
ph[which.min(mse)] #-> 0.1

# Apply KDE to estimate the density with fixed bandwidth
kde_result <-  gauss_cop_est_2d(x = x, y = copula_data, h = ph[which.min(mse)])


#############Surface plot
z <- kde_result$vals_grid
kde_result$vals_grid<-ifelse(kde_result$vals_grid > 12, 12, kde_result$vals_grid)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (kde_result$vals_grid[-1, -1] + kde_result$vals_grid[-1, -ncol(z)] + kde_result$vals_grid[-nrow(z), -1] + kde_result$vals_grid[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = kde_result$vals_grid, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)

# ###############Contourplot
# contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
# title(xlab = "u1", line = 2.5, cex.lab= 1.5)
# title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*kde_result$vals_grid, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)

