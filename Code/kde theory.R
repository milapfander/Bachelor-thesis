#####################################################KDE Theory Chapter#####################
##This Code File includes the code to all of the plots shown in the kde theory Chapter

##########################One-dimensional
set.seed(123)
data <- c(-0.2,  0.2,  0.4,  0.6,  1.0,  1.4,  1.6,  0.3,  2.0)

##Bandwidth at 0.2 -> undersmoothing
kde <- density(data, bw = 0.05, adjust = 1)
plot(kde, lwd = 2, xlim = c(-1, 3), ylim =c(0,10), main = "",  xlab = "", ylab = "", cex.axis = 1.2)
points(x = data, y = c(0, 0, 0, 0, 0, 0, 0, 0, 0))
for(i in 1:length(data)) {
  lines(density(data[i], bw = kde$bw), col = "darkgrey")
}
title(xlab = "x", line = 2.5, cex.lab= 1.5)
title(ylab = "Density", line = 2.5, cex.lab= 1.5)

##Bandwidth at 0.3 -> oversmoothing
kde <- density(data, bw = 0.3, adjust = 1)
plot(kde, lwd = 2,  ylim =c(0,2), main = "",  xlab = "", ylab = "", cex.axis = 1.2)
points(x = data, y = c(0, 0, 0, 0, 0, 0, 0, 0, 0))
for(i in 1:length(data)) {
  lines(density(data[i], bw = kde$bw), col = "darkgrey")
}
title(xlab = "x", line = 2.5, cex.lab= 1.5)
title(ylab = "Density", line = 2.5, cex.lab= 1.5)

##Overplotting
kde <- density(data, bw = 1.2, adjust = 1)
plot(kde, lwd = 2, xlim = c(-1,3), ylim = c(0,1), main = "",  xlab = "", ylab = "", cex.axis = 1.2)
points(x = data, y = c(0, 0, 0, 0, 0, 0, 0, 0, 0))
for(i in 1:length(data)) {
  lines(density(data[i], bw = kde$bw), col = "darkgrey")
}
title(xlab = "x", line = 2.5, cex.lab= 1.5)
title(ylab = "Density", line = 2.5, cex.lab= 1.5)


####################Two dimensional

##We use the package ks for the estimates. This package uses a standard multivariate normal as a Kernel. By setting
##the mixed terms in H to zero, the estimator is equivalent to the product of two standard normal distributions as Kernels

test <- as.matrix(cbind(c(-4, -3, -4, -3.5, 2, 3, 4, 3.5, 4), c(-4,  -3.4, -3, -4, 0.7, 1.5, 0, 2, 4)))
#Define grid
x1_grid <- seq(-8, 8, length.out = 30)
x2_grid <- seq(-8, 8, length.out = 30)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)

###good bandwidth
k <- ks::Hpi.diag(x = test)
fhat <- ks::kde(test, eval.points = x, H = matrix(c(0.7, 0, 0, 0.7), nrow = 2))$estimate

###Plotting the result with the same method as used in Copula theory chapter
z <- cbind(x, fhat)
# Make Matrix with each grid point as a cell
z <- as.data.frame.matrix(z) %>%
  pivot_wider(names_from = x2, values_from = fhat) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
res <- graphics::persp(x = x1_grid, y = x2_grid, z = z, col = color[facetcol], expand = 0.618, cex.axis = 0.7, cex.lab = 1,
                       theta = -30, phi = 20, border = "black", xlab = "\n \n x1", ylab = "\n \n x2", zlab = "", ticktype = "detailed")
mtext("                 Density", side = 2, line = 1, cex = 1)


#####Too large bandwidth
est <-  ks::kde(test,eval.points = x, H= matrix(c(15, 0, 0, 15), nrow = 2))
fhat <- ks::kde(test,eval.points = x, H= matrix(c(15, 0, 0, 15), nrow = 2))$estimate
###Plotting the result with the same method as used in Copula theory chapter
z <- cbind(x, fhat)
# Make Matrix with each grid point as a cell
z <- as.data.frame.matrix(z) %>%
  pivot_wider(names_from = x2, values_from = fhat) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
res <- graphics::persp(x = x1_grid, y = x2_grid, z = z, col = color[facetcol], expand = 0.618, cex.axis = 0.7, cex.lab = 1,
                       theta = -30, phi = 20, border = "black", xlab = "\n \n x1", ylab = "\n \n x2", zlab = "", ticktype = "detailed")
mtext("                 Density", side = 2, line = 1, cex = 1)



#####Too small bandwidth
fhat <- ks::kde(test,eval.points = x, H= matrix(c(0.15, 0, 0, 0.15), nrow = 2))$estimate
###Plotting the result with the same method as used in Copula theory chapter
z <- cbind(x, fhat)
# Make Matrix with each grid point as a cell
z <- as.data.frame.matrix(z) %>%
  pivot_wider(names_from = x2, values_from = fhat) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
res <- persp(x = x1_grid, y = x2_grid, z = z, col = color[facetcol], expand = 0.618, cex.axis = 0.7, cex.lab = 1,
              theta = -30, phi = 20, border = "black", xlab = "\n \n x1", ylab = "\n \n x2", zlab = "", ticktype = "detailed")
mtext("                 Density", side = 2, line = 1, cex = 1)


