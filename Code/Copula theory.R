#####################################################Copula Theory Chapter#####################
## This Code File includes the code for all of the plots shown in the Copula theory chapters, and
## certain calculations, such as of the maximum values of the copula densities.
###############################Visualization chapter
set.seed(123456789)
library(copula)
library(VineCopula)
library(tidyverse)



#Get parameter for copula from specific kendalls tau through https://search.r-project.org/CRAN/refmans/CopulaCenR/html/tau_copula.html
param <- BiCopTau2Par(family = 5, tau = 0.7)
# Create a frank copula
copula <- frankCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)

##True density at that points
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()

#Get maximum density
max(density_cop)

####################Plots of true density#############
#################Scatterplot
plot(copula_data, xlab = "", ylab = "", cex.axis = 1.2)
title(xlab = "u1", line = 2.5, cex.lab= 1.5)
title(ylab = "u2", line = 2.5, cex.lab= 1.5)


#############Surface plot
z <- vals_grid_f
vals_grid_f<-ifelse(vals_grid_f > 12, 12, vals_grid_f)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (vals_grid_f[-1, -1] + vals_grid_f[-1, -ncol(z)] + vals_grid_f[-nrow(z), -1] + vals_grid_f[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = vals_grid_f, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)


###############Contourplot
contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
title(xlab = "u1", line = 2.5, cex.lab= 1.5)
title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 50)
x2_grid <- seq(0.01, 0.99, length.out = 50)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)

##True density at that points
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
vals_grid_true <- cbind(x, density_cop)
vals_grid_true <- as.data.frame.matrix(vals_grid_true) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
# Actual contour plot
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*vals_grid_true, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)











###############################Some examples chapters

########################################Independence Copula
set.seed(123456789)

# Create an Independence copula
copula <- indepCopula(dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)

##True density at that points
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()

#Get maximum density
max(density_cop)


####################Plots of true density#############
#################Scatterplot
plot(copula_data, xlab = "", ylab = "", cex.axis = 1.2)
title(xlab = "u1", line = 2.5, cex.lab= 1.5)
title(ylab = "u2", line = 2.5, cex.lab= 1.5)


#############Surface plot
z <- vals_grid_f
vals_grid_f<-ifelse(vals_grid_f > 12, 12, vals_grid_f)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (vals_grid_f[-1, -1] + vals_grid_f[-1, -ncol(z)] + vals_grid_f[-nrow(z), -1] + vals_grid_f[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = vals_grid_f, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)


###############Normal contour plot
# Create a 2D grid on the unit square
x1_grid <- seq(0.0000001, 0.99999999, length.out = 50)
x2_grid <- seq(0.0000001, 0.99999999, length.out = 50)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)

##True density at that point
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
vals_grid_true <- cbind(x, density_cop)
vals_grid_true <- as.data.frame.matrix(vals_grid_true) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
# Actual plot
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*vals_grid_true, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8, xlim = c(-3,3), ylim = c(-3,3))
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)










########################################Gaussian Copula
###################Tau = 0.3
set.seed(123456789)


#Get parameter for copula from specific kendalls tau through https://search.r-project.org/CRAN/refmans/VineCopula/html/BiCopTau2Par.html

param <- BiCopTau2Par(family = 1, tau = 0.3)
# Create a gaussian copula
copula <- normalCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)

##True density at that points
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()

#Get maximum density
max(density_cop)


####################Plots of true density#############
#################Scatterplot
plot(copula_data, xlab = "", ylab = "", cex.axis = 1.2)
title(xlab = "u1", line = 2.5, cex.lab= 1.5)
title(ylab = "u2", line = 2.5, cex.lab= 1.5)


#############Surface plot
z <- vals_grid_f
vals_grid_f<-ifelse(vals_grid_f > 12, 12, vals_grid_f)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (vals_grid_f[-1, -1] + vals_grid_f[-1, -ncol(z)] + vals_grid_f[-nrow(z), -1] + vals_grid_f[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = vals_grid_f, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)


###############Normal contour plot
# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 50)
x2_grid <- seq(0.01, 0.99, length.out = 50)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that points
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
vals_grid_true <- cbind(x, density_cop)
vals_grid_true <- as.data.frame.matrix(vals_grid_true) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
# Actual plot
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*vals_grid_true, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)



###################Tau = 0.7
set.seed(123456789)


#Get parameter for copula from specific kendalls tau through https://search.r-project.org/CRAN/refmans/VineCopula/html/BiCopTau2Par.html

param <- BiCopTau2Par(family = 1, tau = 0.7)
# Create a Gaussian copula
copula <- normalCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that points
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()

#Get maximum density
max(density_cop)


####################Plots of true density#############
#################Scatterplot
plot(copula_data, xlab = "", ylab = "", cex.axis = 1.2)
title(xlab = "u1", line = 2.5, cex.lab= 1.5)
title(ylab = "u2", line = 2.5, cex.lab= 1.5)


#############Surface plot
z <- vals_grid_f
vals_grid_f<-ifelse(vals_grid_f > 12, 12, vals_grid_f)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (vals_grid_f[-1, -1] + vals_grid_f[-1, -ncol(z)] + vals_grid_f[-nrow(z), -1] + vals_grid_f[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = vals_grid_f, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)


###############Normal contour plot
# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 50)
x2_grid <- seq(0.01, 0.99, length.out = 50)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that point
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
vals_grid_true <- cbind(x, density_cop)
vals_grid_true <- as.data.frame.matrix(vals_grid_true) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)

#Actual plot
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*vals_grid_true, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)



















########################################Clayton Copula
###################Tau = 0.3
set.seed(123456789)


#Get parameter for copula from specific kendalls tau through https://search.r-project.org/CRAN/refmans/VineCopula/html/BiCopTau2Par.html

param <- BiCopTau2Par(family = 3, tau = 0.3)
# Create a clayton copula
copula <- claytonCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)

##True density at that points
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()

#Get maximum density
max(density_cop)


####################Plots of true density#############
#################Scatterplot
plot(copula_data, xlab = "", ylab = "", cex.axis = 1.2)
title(xlab = "u1", line = 2.5, cex.lab= 1.5)
title(ylab = "u2", line = 2.5, cex.lab= 1.5)


#############Surface plot
z <- vals_grid_f
vals_grid_f<-ifelse(vals_grid_f > 12, 12, vals_grid_f)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (vals_grid_f[-1, -1] + vals_grid_f[-1, -ncol(z)] + vals_grid_f[-nrow(z), -1] + vals_grid_f[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = vals_grid_f, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)



###############Normal contour plot
# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 50)
x2_grid <- seq(0.01, 0.99, length.out = 50)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that points
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
vals_grid_true <- cbind(x, density_cop)
vals_grid_true <- as.data.frame.matrix(vals_grid_true) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
# Actual plot
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*vals_grid_true, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)



###################Tau = 0.7
set.seed(123456789)


#Get parameter for copula from specific kendalls tau through https://search.r-project.org/CRAN/refmans/VineCopula/html/BiCopTau2Par.html

param <- BiCopTau2Par(family = 3, tau = 0.7)
# Create a Clayton copula
copula <- claytonCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that points
density_cop <- dCopula(as.matrix(x), copula)
range(density_cop)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()

#Get maximum density
max(density_cop)


####################Plots of true density#############
#################Scatterplot
plot(copula_data, xlab = "", ylab = "", cex.axis = 1.2)
title(xlab = "u1", line = 2.5, cex.lab= 1.5)
title(ylab = "u2", line = 2.5, cex.lab= 1.5)


#############Surface plot
z <- vals_grid_f
vals_grid_f<-ifelse(vals_grid_f > 12, 12, vals_grid_f)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (vals_grid_f[-1, -1] + vals_grid_f[-1, -ncol(z)] + vals_grid_f[-nrow(z), -1] + vals_grid_f[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = vals_grid_f, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)

###############Contourplot
#contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
#title(xlab = "u1", line = 2.5, cex.lab= 1.5)
#title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 50)
x2_grid <- seq(0.01, 0.99, length.out = 50)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that point
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
vals_grid_true <- cbind(x, density_cop)
vals_grid_true <- as.data.frame.matrix(vals_grid_true) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
# Actual plot
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*vals_grid_true, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)









########################################Frank Copula
###################Tau = 0.3
set.seed(123456789)


#Get parameter for copula from specific kendalls tau through https://search.r-project.org/CRAN/refmans/VineCopula/html/BiCopTau2Par.html

param <- BiCopTau2Par(family = 5, tau = 0.3)
# Create a frank copula
copula <- frankCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that point
density_cop <- dCopula(as.matrix(x), copula)
range(density_cop)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()

#Get maximum density
max(density_cop)


####################Plots of true density#############
#################Scatterplot
plot(copula_data, xlab = "", ylab = "", cex.axis = 1.2)
title(xlab = "u1", line = 2.5, cex.lab= 1.5)
title(ylab = "u2", line = 2.5, cex.lab= 1.5)


#############Surface plot
z <- vals_grid_f
vals_grid_f<-ifelse(vals_grid_f > 12, 12, vals_grid_f)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (vals_grid_f[-1, -1] + vals_grid_f[-1, -ncol(z)] + vals_grid_f[-nrow(z), -1] + vals_grid_f[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = vals_grid_f, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)
#persp(copula, dCopula, n.grid = 30, col = "lightblue", border = "black", theta = -30, phi = 20, xlab = "\n u1", ylab = "\n u2", zlab = "\n Density")

###############Contourplot
#contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
#title(xlab = "u1", line = 2.5, cex.lab= 1.5)
#title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
###############Normal contour plot
# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 50)
x2_grid <- seq(0.01, 0.99, length.out = 50)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that point
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
vals_grid_true <- cbind(x, density_cop)
vals_grid_true <- as.data.frame.matrix(vals_grid_true) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
# Actual plot
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*vals_grid_true, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)



###################Tau = 0.7
set.seed(123456789)


#Get parameter for copula from specific kendalls tau through https://search.r-project.org/CRAN/refmans/VineCopula/html/BiCopTau2Par.html

param <- BiCopTau2Par(family = 5, tau = 0.7)
# Create a frank copula
copula <- frankCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that point
density_cop <- dCopula(as.matrix(x), copula)
range(density_cop)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()

#Get maximum density
max(density_cop)


####################Plots of true density#############
#################Scatterplot
plot(copula_data, xlab = "", ylab = "", cex.axis = 1.2)
title(xlab = "u1", line = 2.5, cex.lab= 1.5)
title(ylab = "u2", line = 2.5, cex.lab= 1.5)


#############Surface plot
z <- vals_grid_f
vals_grid_f<-ifelse(vals_grid_f > 12, 12, vals_grid_f)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (vals_grid_f[-1, -1] + vals_grid_f[-1, -ncol(z)] + vals_grid_f[-nrow(z), -1] + vals_grid_f[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = vals_grid_f, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)

###############Contourplot
#contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
#title(xlab = "u1", line = 2.5, cex.lab= 1.5)
#title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 50)
x2_grid <- seq(0.01, 0.99, length.out = 50)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that point
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
vals_grid_true <- cbind(x, density_cop)
vals_grid_true <- as.data.frame.matrix(vals_grid_true) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)

# Actual plot
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*vals_grid_true, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)








########################################Gumbel Copula
###################Tau = 0.3
set.seed(123456789)


#Get parameter for copula from specific kendalls tau through https://search.r-project.org/CRAN/refmans/VineCopula/html/BiCopTau2Par.html

param <- BiCopTau2Par(family = 4, tau = 0.3)
# Create a gumbel copula
copula <- gumbelCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that point
density_cop <- dCopula(as.matrix(x), copula)
range(density_cop)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()

#Get maximum density
max(density_cop)


####################Plots of true density#############
#################Scatterplot
plot(copula_data, xlab = "", ylab = "", cex.axis = 1.2)
title(xlab = "u1", line = 2.5, cex.lab= 1.5)
title(ylab = "u2", line = 2.5, cex.lab= 1.5)


#############Surface plot
z <- vals_grid_f
vals_grid_f<-ifelse(vals_grid_f > 12, 12, vals_grid_f)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (vals_grid_f[-1, -1] + vals_grid_f[-1, -ncol(z)] + vals_grid_f[-nrow(z), -1] + vals_grid_f[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = vals_grid_f, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)

###############Contourplot
#contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
#title(xlab = "u1", line = 2.5, cex.lab= 1.5)
#title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 50)
x2_grid <- seq(0.01, 0.99, length.out = 50)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that point
density_cop <- dCopula(as.matrix(x), copula)
range(density_cop)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
vals_grid_true <- cbind(x, density_cop)
vals_grid_true <- as.data.frame.matrix(vals_grid_true) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)
# Actual plot
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*vals_grid_true, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)



###################Tau = 0.7
set.seed(123456789)


#Get parameter for copula from specific kendalls tau through https://search.r-project.org/CRAN/refmans/VineCopula/html/BiCopTau2Par.html

param <- BiCopTau2Par(family = 4, tau = 0.7)
# Create a gumbel copula
copula <- gumbelCopula(param = param, dim = 2)

# Generate random data from the copula
copula_data <- rCopula(500, copula)

# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 20)
x2_grid <- seq(0.01, 0.99, length.out = 20)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that point
density_cop <- dCopula(as.matrix(x), copula)
#Values on the grid for plotting
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()

#Get maximum density
max(density_cop)
x[which.max(density_cop),]
##-> maximum is in the top right tail and not on the bottom left as assumed from the plot -> perspective of the plot
density_cop[1]
##-> value at the other tail

####################Plots of true density#############
#################Scatterplot
plot(copula_data, xlab = "", ylab = "", cex.axis = 1.2)
title(xlab = "u1", line = 2.5, cex.lab= 1.5)
title(ylab = "u2", line = 2.5, cex.lab= 1.5)


#############Surface plot
z <- vals_grid_f
vals_grid_f<-ifelse(vals_grid_f > 12, 12, vals_grid_f)
jet.colors <- colorRampPalette( c("white", "blue"))
nbcol <- 50
color <- jet.colors(nbcol)
## Compute the z-value at the facet centres
zfacet <-  (vals_grid_f[-1, -1] + vals_grid_f[-1, -ncol(z)] + vals_grid_f[-nrow(z), -1] + vals_grid_f[-nrow(z), -ncol(z)])/4
## Recode facet z-values into color indices
facetcol <- cut(zfacet, 50)
persp(x = x1_grid, y = x2_grid, z = vals_grid_f, xlim = c(0,1), zlim = c(0,8), col = color[facetcol], expand = 0.618, cex.axis = 1.3, cex.lab = 1.5,
      ticktype = "detail", theta = -30, phi = 20, border = "black", xlab = "\n \n u1", ylab = "\n \n u2", zlab = "")
mtext("            c(u1,u2)", side = 2, line = 1, cex = 1.5)

###############Contourplot
#contour(copula, dCopula, xlab = "", ylab = "", cex.axis = 1.2)
#title(xlab = "u1", line = 2.5, cex.lab= 1.5)
#title(ylab = "u2", line = 2.5, cex.lab= 1.5)


###############Normal contour plot
# Create a 2D grid on the unit square
x1_grid <- seq(0.01, 0.99, length.out = 50)
x2_grid <- seq(0.01, 0.99, length.out = 50)
x <- expand.grid(x1 = x1_grid, x2 = x2_grid)
##True density at that point
density_cop <- dCopula(as.matrix(x), copula)
vals_grid_t <- cbind(x, density_cop)
vals_grid_f <- as.data.frame.matrix(vals_grid_t) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
vals_grid_true <- cbind(x, density_cop)
vals_grid_true <- as.data.frame.matrix(vals_grid_true) %>%
  pivot_wider(names_from = x2, values_from = density_cop) %>%
  column_to_rownames(var = "x1") %>%
  as.matrix()
kde_transf_1 <- qnorm(x1_grid)
kde_transf_2 <- qnorm(x2_grid)

#Actual plot
contour(x = kde_transf_1, y = kde_transf_2, z = tcrossprod(dnorm(kde_transf_1), dnorm(kde_transf_2))*vals_grid_true, 
        xlab = "", ylab = "", cex.axis = 1.2, vfont=c("sans serif", "bold italic"), labcex=0.8)
title(xlab = expression(Phi^-1 * (u1)), line = 2.7, cex.lab= 1.5)
title(ylab = expression(Phi^-1 * (u2)), line = 2.2, cex.lab= 1.5)






