##Function for the used kernel
gaussian_cop_kernel <- function(x, y, h) {
  b <- 1-h^2
  1/(sqrt(1-b^2))*exp((-(b^2*(qnorm(x))^2-2*b*qnorm(x)*qnorm(y)+b^2*(qnorm(y))^2)) / (2*(1-b^2)))
}

rot_bw <- function(y) {
  n <- length(y)
  mu <- mean(qnorm(y))
  sigma <- sd(qnorm(y))
  h <- sigma*(((2*mu^2*sigma^2)+(3*(1-sigma^2)^2))^(-1/5))*n^{-1/5}
  h
}

####one dimensional gaussian kernel estimation
Gaussian_est_1d <- function(x = NULL, y, k = NULL, h = NULL){
  n <- length(y)
  #generate 101 equally spaced grid points
  if(is.null(x)) {
    x <- seq(from = 0, to = 1, by = 0.01)
  }
  #number of grid points
  if(is.null(k)) {
    k <- length(y)
  }
  #Bandwidth - rule of thumb implementation
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
  
  #class results?
  class (results) <- c('list', 'GC')
  
  #print results
  results
}
 
 ##Test von Paper
 y <- c(0.1, 0.16, 0.25, 0.34, 0.42, 0.5, 0.58, 0.66, 0.75, 0.84, 0.9)
 h <- sqrt(0.1)
 k <- 101
 xx <- seq(from = 0, to = 1, by = 0.01)
 #mit gewählter bandwidth
 result2 <- Gaussian_est_1d(x = xx, y = y, k = 101, h = h)
 #mit "optimal" bandwidth
 result2 <- Gaussian_est_1d(x = xx, y = y, k = 101)


# Create the plot from the paper
plot(result2$avg_est[, "x"], result2$avg_est[, "fhat"], type = 'l', col = 'black', xlab = 'X', ylab = 'Kernel Density', ylim=c(0,5), xlim=c(0,1))
for (i in 1:length(y)) {
  lines(xx, result2$ind_est[i,], col = 'lightgrey', lty = 2)
}

