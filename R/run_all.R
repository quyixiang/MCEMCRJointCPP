# simulation_data <- CRsimulation(
#   fmla.tte = as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ 0 + Y0SCALE),
#   fmla.long = as.formula(PCHG ~ 0 + Y0SCALE),
#   beta.tte = c(0.2), normal.tte = TRUE, sd.tte = 0.2,
#   beta.y = c(-0.5), sd.y = 0.1, beta.y.cure = c(-0.2), sd.y.cured = 0.2,
#   cured.rate = 0.4, cured.mean = c(0, -0.2), cured.sd = c(0.2, 0.2),
#   cured.corr = matrix(
#     c(1, 0.5, 0.5, 1),
#     nrow = 2
#   ),
#   randeff.mean = c(0.5, 0, -0.5, 0.5), randeff.sd = rep(0.2, 4),
#   randeff.corr = matrix(
#     c(
#       1, -0.4, -0.2, -0.3,
#       -0.4, 1, 0.5, 0.20,
#       -0.2, 0.5, 1, 0.2,
#       -0.3, 0.20, 0.2, 1
#     ),
#     nrow = 4
#   ), n = 200, censor.parameter = 0.5, time.interval = 0.1, seed = 5
# )
#
# longdat <- simulation_data[["longdat"]]
# survdat <- simulation_data[["survdat"]]
# real_para <- simulation_data[["simulation.para"]]
#
# Sigma_r_init <- diag(c(0.2, 0.2, 0.2, 0.2)) %*% matrix(
#   c(
#     1, -0.4, -0.2, -0.3,
#     -0.4, 1, 0.5, 0.20,
#     -0.2, 0.5, 1, 0.2,
#     -0.3, 0.20, 0.2, 1
#   ),
#   nrow = 4
# ) %*% diag(c(0.2, 0.2, 0.2, 0.2))
#
# Sigma_r_cure_init <- diag(c(0.2, 0.2)) %*% matrix(
#   c(1,0.5,0.5,1), nrow = 2
# ) %*% diag(c(0.2, 0.2))
#
# init = list(
#   mu_r = c(0.5, 0, -0.5, 0.5),
#   Sigma_r = Sigma_r_init,
#   beta_y = c(-0.5),
#   sigma_y_sq = 0.01,
#   beta_tte = c(0.2),
#   sigma_tte_sq = 0.04,
#   pi_c = 0.4,
#   beta_cure = c(-0.2),
#   mu_r_cure = c(0, -0.2),
#   Sigma_r_cure = Sigma_r_cure_init,
#   sigma_y_cure_sq = 0.04
# )
# fmla.tte = as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ 0 + Y0SCALE)
# fmla.long = as.formula(PCHG ~ 0 + Y0SCALE)
# longdat.time <- "visittime"
# id.indicator <- "id"
# infer_CRJoint_MLE(survdat, longdat, init, fmla.tte, fmla.long, longdat.time, id.indicator)

infer_CRJoint_MLE <- function(survdat, longdat, init, fmla.tte, fmla.long, longdat.time, id.indicator, maxIter = 5000, tol = 8e-4, no_cure = FALSE) {
  data.list <- load_onearm_data(
    survdat = survdat,
    longdat = longdat,
    fmla.tte = fmla.tte,
    fmla.long = fmla.long,
    longdat.time = longdat.time,
    id.indicator = id.indicator
  )


  data.list[["tcen"]] <- log(data.list[["tcen"]])
  data.list[["tobs"]] <- log(data.list[["tobs"]])


  cure_res <- MCEM_cureJoint(data.list, initial = init, maxIter = maxIter, tol = tol, no_cure = no_cure)
}
