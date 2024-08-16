####Checking performance of new estimator


### basic estimator

#load libraries for this file
library(tidyverse)
library(copula)
library(VineCopula)
library(rvinecopulib)

##Function for the used kernel
gaussian_cop_kernel_new <- function(x, y, h, rho = 0) {
  b <- 1-h^2
  k1 <- dbicop(cbind(x[, 1], y[, 1]), "gauss", 0, b)
  k2 <- dbicop(cbind(x[, 2], y[, 2]), "gauss", 0, b)
  
  # bisheriger kernel war k1 * k2, produkt von zwei bedingten dichten wie bei
  # unabhaengigkeit. jetzt multiplizieren wir noch einen term daran, der
  # abhaengigkeit im kernel erzeugt. die hbicop transformation stellt sicher
  # dass der kernel insgesamt eine dichte bleibt
  h1 <- hbicop(cbind(x[, 1], y[, 1]), cond_var = 2, "gauss", 0, b)
  h2 <- hbicop(cbind(x[, 2], y[, 2]), cond_var = 2, "gauss", 0, b)
  
  k12 <- dbicop(cbind(h1, h2), "gauss", 0, rho)
  k1 * k2 * k12
}

Gaussian_est_2d_new <- function(x = NULL, y, k = NULL, h = NULL, bw_each_l = NULL, indep = TRUE){
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
  if(indep == TRUE) {
    rho <- 0
  }
  else {
    rho <- cor(qnorm(y))[1, 2]
  }
  
  # for(i in 1:n) {
  for(j in 1:p) {
    # for the index of each data-point and each grid-point, calculate a product of the 
    # univariate Kernels. The result is a density value for each data-point and grid-value
    Kgaussian[, j] <- gaussian_cop_kernel_new(x[j, ], y, bw_each[1], rho)
  }
  # }
  
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
  
  # class results
  class (results) <- c('list', 'GC')
  
  results
}

##function surface plot
surface_plot <- function(copula, cutoff = 12, nbcol = 50, est, ph = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45),
                         x1_grid = seq(0.01, 0.99, length.out = 20), x2_grid = seq(0.01, 0.99, length.out = 20)) {
  x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
  density_cop <- dCopula(as.matrix(x), copula)
  copula_data <- rCopula(500, copula)
  
  
  if(est == "true") {
    vals_grid_t <- cbind(x, density_cop)
    vals_grid <- as.data.frame.matrix(vals_grid_t) %>%
      pivot_wider(names_from = x2, values_from = density_cop) %>%
      column_to_rownames(var = "x1") %>%
      as.matrix()
  }
  else if (est == "basic") {
    mse <- vector(mode = "numeric", length = length(ph))
    for(i in seq_along(ph)) {
      kde_result <-  Gaussian_est_2d_new(x = x, y = copula_data, h = ph[i], indep = TRUE)
      mse[i] <- sum((density_cop - kde_result$avg_est[,3])^2)/1225
    }
    print(paste0("optimal MSE:", ph[which.min(mse)]))
    vals_grid <-  Gaussian_est_2d_new(x = x, y = copula_data, h = ph[which.min(mse)], indep = TRUE)$vals_grid
  }
  else if (est == "new") {
    mse <- vector(mode = "numeric", length = length(ph))
    for(i in seq_along(ph)) {
      kde_result <-  Gaussian_est_2d_new(x = x, y = copula_data, h = ph[i], indep = FALSE)
      mse[i] <- sum((density_cop - kde_result$avg_est[,3])^2)/1225
    }
    print(paste0("optimal MSE:", ph[which.min(mse)]))
    vals_grid <-  Gaussian_est_2d_new(x = x, y = copula_data, h = ph[which.min(mse)], indep = FALSE)$vals_grid
  }
  
  
  z <- vals_grid
  #cut off values at 12
  vals_grid<-ifelse(vals_grid > cutoff, cutoff, vals_grid)
  jet.colors <- colorRampPalette( c("white", "blue"))
  nbcol <- 50
  color <- jet.colors(nbcol)
  ## Compute the z-value at the facet centres
  zfacet <-  (vals_grid[-1, -1] + vals_grid[-1, -ncol(z)] + vals_grid[-nrow(z), -1] + vals_grid[-nrow(z), -ncol(z)])/4
  ## Recode facet z-values into color indices
  facetcol <- cut(zfacet, 50)
  persp(x = x1_grid, y = x2_grid, z = vals_grid, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
        ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
  mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)
}




##function normal plot
normal_plot <- function(copula, est = NULL, nbcol = 50, x1_grid = seq(-3,3, length.out = 50), 
                        x2_grid = seq(-3,3, length.out = 50), ph = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45), levels) {
  x <- expand.grid(x1 = pnorm(x1_grid), x2 = pnorm(x2_grid))
  density_cop <- dCopula(as.matrix(x), copula)
  
  vals_grid_t <- cbind(x, density_cop)
  vals_grid <- as.data.frame.matrix(vals_grid_t) %>%
    pivot_wider(names_from = x2, values_from = density_cop) %>%
    column_to_rownames(var = "x1") %>%
    as.matrix()

  # Actual contour plot
  contour(x = seq(-3,3, length.out = 50), y = seq(-3,3, length.out = 50), 
          z = tcrossprod(dnorm(seq(-3,3, length.out = 50)), dnorm(seq(-3,3, length.out = 50)))*vals_grid, 
          xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8, levels = levels)
  title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
  title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)
  
 if(est == "basic") {
   # Get best case bandwidth
   mse <- vector(mode = "numeric", length = length(ph))
   for(i in 1:length(ph)) {
     kde_result <-  Gaussian_est_2d_new(x = x, y = copula_data, h = ph[i], indep = TRUE)
     mse[i] <- sum((density_cop - kde_result$avg_est[,3])^2)/1225
   }
   print(paste0("optimal MSE:", ph[which.min(mse)]))
   vals_grid <-  Gaussian_est_2d_new(x = x, y = copula_data, h = ph[which.min(mse)], indep = TRUE)$vals_grid
   
   contour(x = seq(-3,3, length.out = 50), y = seq(-3,3, length.out = 50), 
           z = tcrossprod(dnorm(seq(-3,3, length.out = 50)), dnorm(seq(-3,3, length.out = 50)))*vals_grid, 
           xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8, add = TRUE, col = "lightblue", levels = levels)
   title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
   title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)
 }
  
  else if(est == "new") {
    # Get best case bandwidth
    mse <- vector(mode = "numeric", length = length(ph))
    for(i in 1:length(ph)) {
      kde_result <-  Gaussian_est_2d_new(x = x, y = copula_data, h = ph[i], indep = FALSE)
      mse[i] <- sum((density_cop - kde_result$avg_est[,3])^2)/1225
    }
    print(paste0("optimal MSE:", ph[which.min(mse)]))
    vals_grid <-  Gaussian_est_2d_new(x = x, y = copula_data, h = ph[which.min(mse)], indep = FALSE)$vals_grid
    
    contour(x = seq(-3,3, length.out = 50), y = seq(-3,3, length.out = 50), 
            z = tcrossprod(dnorm(seq(-3,3, length.out = 50)), dnorm(seq(-3,3, length.out = 50)))*vals_grid, 
            xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8, add = TRUE, col = "lightblue", levels = levels)
    title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
    title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)
  }
}

levels <- seq(from = 0.05, to = 0.4, by = 0.05)

############################################# Frank 0.7##########################
################################################ True density
#Get parameter for copula from specific kendalls tau through https://search.r-project.org/CRAN/refmans/CopulaCenR/html/tau_copula.html
param <- BiCopTau2Par(family = 5, tau = 0.7)
# Create a frank copula
copula <- frankCopula(param = param, dim = 2)
copula_data <- rCopula(500, copula)

#plots
surface_plot(copula = copula, est = "true")
surface_plot(copula = copula,  est = "basic")
surface_plot(copula = copula, est = "new")


normal_plot(copula = copula, est = "basic", levels = levels)
normal_plot(copula = copula, est = "new", levels = levels)




############################################# Clayton 0.7##########################
################################################ True density
#Get parameter for copula from specific kendalls tau through https://search.r-project.org/CRAN/refmans/CopulaCenR/html/tau_copula.html
param <- BiCopTau2Par(family = 3, tau = 0.7)
# Create a clayton copula
copula <- claytonCopula(param = param, dim = 2)
copula_data <- rCopula(500, copula)

#plots
surface_plot(copula = copula, est = "true")
surface_plot(copula = copula,  est = "basic")
surface_plot(copula = copula, est = "new")


normal_plot(copula = copula, est = "basic", levels = levels)
normal_plot(copula = copula, est = "new", levels = levels)



############################################# Gumbel 0.7##########################
################################################ True density
set.seed(123456789)
#Get parameter for copula from specific kendalls tau through https://search.r-project.org/CRAN/refmans/VineCopula/html/BiCopTau2Par.html
param <- BiCopTau2Par(family = 4, tau = 0.7)
# Create a gumbel copula
copula <- gumbelCopula(param = param, dim = 2)
copula_data <- rCopula(500, copula)

#plots
surface_plot(copula = copula, est = "true")
surface_plot(copula = copula,  est = "basic")
surface_plot(copula = copula, est = "new")


normal_plot(copula = copula, est = "basic", levels = levels)
normal_plot(copula = copula, est = "new", levels = levels)


