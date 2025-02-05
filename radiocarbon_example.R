# This script relies on the R package baydem available at:
#
# https://github.com/eehh-stanford/baydem
#
# The specific commit of baydem used for this script is:
#
# 8be9cc53bbedbbe0ba6f362e537ecdb152bc3c5f
#
# This version can be installed in R using the following command:
#
# devtools::install_github("eehh-stanford/baydem",ref="8be9cc53bbedbbe0ba6f362e537ecdb152bc3c5f")

# Clear the workspace
rm(list = ls())

library(baydem)
library(Bchron)
library(rcarbon)

# There are many excellent tools to can be used to calibrate a calendar date.
# Before plotting the illustration plots, compare baydem's calibration against
# rcarbon's to ensure they give the same results.

# The date to be calibrated, about 833 AD
t_m   <- 1950 - 830 # The uncalibrated radiocarbon date in years before present
sig_t_m <- 20       # The uncertainty of the measurement (20 "years")

# Do the calibration in rcarbon
rcarbonCalib <- rcarbon::calibrate(x=t_m,errors=sig_t_m,calCurves='intcal20',F14C=T,eps=1e-6)

# Now, do the calibration in baydem using the same minimum and maximum for the
# calendar dates (tau)
tau_min <- 1950 - max(rcarbonCalib$grids[[1]]$calBP)
tau_max <- 1950 - min(rcarbonCalib$grids[[1]]$calBP)
tau_mid <- (tau_min + tau_max)/2
dtau <- 1
tau <- seq(tau_min,tau_max,by=dtau)

# The calibration dataframe, which has three columns:
# yearBP	uncalYearBP	uncalYearBPError
calibDf = baydem::load_calib_curve("intcal20")

# The reference decay rate of carbon
kappa <- 8033

# The fraction modern and associated uncertainty
phi_m <- exp(-t_m/kappa)
sig_m <- sig_t_m * phi_m / kappa

# Calculate the likelihood using baydem
lik <- as.vector(calc_meas_matrix(tau,phi_m,sig_m,calibDf,T,F))

# Calculate the posterior density assuming a uniform prior
f_prior_unif <- rep(1,length(tau)) / (tau_max-tau_min)

f_posterior_unif <- f_prior_unif * lik
f_posterior_unif <- f_posterior_unif / sum(f_posterior_unif*dtau)

# Do the actual check to ensure rcarbon and baydem give the same results
# (within an appropriate tolerance)
# The rcarbon result is not necessarily evenly spaced. Where a calendar date is
# missing, set the density to zero
f_rcarbon <- rep(NA,length(tau))
for(k in 1:length(tau)) {
  ind_k <- which(rcarbonCalib$grids[[1]]$calBP == 1950 - tau[k])
  if(length(ind_k) == 0) {
    f_rcarbon[k] <- 0
  } else if(length(ind_k) == 1) {
    f_rcarbon[k] <- rcarbonCalib$grids[[1]]$PrDens[ind_k]
  } else {
    stop('There should be 0 or 1 matches to tau[k] in rcarbonCalib$grids[[1]]$calBP')
  }
}

testthat::expect_equal(
  f_posterior_unif,
  f_rcarbon,
  tol=1e-6
)

# Calculate the posterior density assuming a triangular prior

f_prior_tria <- (tau_mid - tau_min)-abs(tau-tau_mid)
f_prior_tria <- f_prior_tria / sum(dtau*f_prior_tria)

f_posterior_tria <- f_prior_tria * lik
f_posterior_tria <- f_posterior_tria / sum(f_posterior_tria*dtau)



fileName <- here::here("single_date_calibration.pdf")
### FIGURE 1: Calibration of a single date (using two priors)
pdf(fileName, width = 10, height = 6)
  plot(tau,f_prior_unif,col='black',lwd=3,type='l',xlab='Calendar Date [AD]',ylab='Density',lty=3,ylim=c(0,max(f_posterior_tria)))
  lines(tau,f_prior_tria,col='grey',lwd=3,lty=3)
  lines(tau,f_posterior_unif,col='black',lwd=3)
  lines(tau,f_posterior_tria,col='grey',lwd=3)
dev.off()
