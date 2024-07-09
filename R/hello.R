# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}
#
#
# s_i <- rnorm(10)  # Random vector for s_i
# y_i <- rnorm(10)  # Random vector for y_i
# mu_tte <- 0
# sd_tte <- 1
# mu_y <- rnorm(10)
# sigma_y_sq <- 2
# mu_r <- rnorm(4)
# sigma_r <- diag(rep(1, 4))
# sample_J <- 5
# x_min <- -1
# # MC_int_xy_cpp(s_i, y_i, mu_tte, sd_tte, mu_y, sigma_y_sq, mu_r, sigma_r, sample_J, x_min, obs = FALSE)
# # MC_int_xy_all_cpp(s_i, y_i, mu_tte, sd_tte, mu_y, sigma_y_sq, mu_r, sigma_r, sample_J, x_min, obs = FALSE)
#
# libs <- c("mvtnorm","RcppArmadillo","rbenchmark","bayesm","parallel")
# if (sum(!(libs %in% .packages(all.available = TRUE))) > 0) {
#   install.packages(libs[!(libs %in% .packages(all.available = TRUE))])
# }
#
# for (i in 1:length(libs)) {
#   library(libs[i],character.only = TRUE,quietly=TRUE)
# }
# set.seed(42)
# sigma <- rwishart(10,diag(4))$IW
# means <- rnorm(4)
# X     <- rmvnorm(500000, means, sigma)
# X     <- X[which(X[,1] > -1),]

# dmvnorm_arma(matrix(X[1,], ncol = 4),means,sigma,log=TRUE)
#
# library(MCEMCRJointCPP)
# cpp_res = rtruncnorm_cpp(10000,rep(0,10000),rep(1,10000),rep(0,10000),rep(0.01,10000))
# library(truncnorm)
# res = rtruncnorm(10000,0,1,0,0.01)
# mean(cpp_res)
# mean(res)


