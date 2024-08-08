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
# save(list = ls(), file = "/Users/yixiang/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Ibrahim/MCEM_CRJoint/tmp_res/tmp_res.RData")
# load("/Users/yixiang/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Ibrahim/MCEM_CRJoint/tmp_res/tmp_res.RData")
# Q_function_cpp(
#   Nobs, nobs, E_r_obs, E_r_rT_obs, Estep_2_obs, yobs, Xobs, Xtte_obs, visittime_obs, tobs, new_id_obs,
#   Ncen, ncen, E_r_cen, E_r_rT_cen, Estep_2_cen, ycen, Xcen, Xtte_cen, visittime_cen, E_g0_ti, E_ti, E_ti_sq, new_id_cen,
#   beta_tte, sigma_tte_sq, mu_r, Sigma_r, beta_y, sigma_y_sq,
#   E_Delta_cen
# )
# cpp_res = rtruncnorm_cpp(10000,rep(0,10000),rep(1,10000),rep(0,10000),rep(0.01,10000))
# cpp_res_2 = rtruncnorm_function(10000,rep(0,10000),rep(1,10000),rep(0,10000),rep(0.01,10000))
#
# library(microbenchmark)
#
# n <- 10000
# a <- rep(0, n)
# b <- rep(1, n)
# mu <- rep(0, n)
# sigma <- rep(0.01, n)
#
# # Benchmark the functions
# benchmark_results <- microbenchmark(
#   cpp_res = rtruncnorm_cpp(n, a, b, mu, sigma),
#   cpp_res_2 = rtruncnorm_function(n, a, b, mu, sigma),
#   times = 100  # Number of times to run each function for benchmarking
# )
#
# # Print the results
# print(benchmark_results)
