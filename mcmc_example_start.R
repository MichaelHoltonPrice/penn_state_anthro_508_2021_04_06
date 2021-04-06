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
rm(list=ls(all=TRUE))


interpolate_cal_curve <- function(tau_AD,calib_df) {
  # tau_AD   The calendar data in years AD (more precisely, 1950 - years BP)
  # calib_df The data frame containing the calibration curve
  #  (see baydem::load_calib_curve)
  tau_BP <- 1950 - tau_AD
  tau_curve <- rev(calib_df$year_BP)

  # Calculate the fraction modern value
  phi_c_curve <- exp(-rev(calib_df$uncal_year_BP) / 8033)

  # Calculate the uncertainty for the fraction modern value
  sig_c_curve <- rev(calib_df$uncal_year_BP_error) * phi_c_curve / 8033

  # Interpolate the calibration curve and uncertainty at tau_BP to yield phi_c
  # and sig_c, the fraction modern and uncertainty of the calibration curve at
  # the input calendar date
  phi_c <- stats::approx(tau_curve, phi_c_curve, tau_BP)
  phi_c <- phi_c$y
  sig_c <- stats::approx(tau_curve, sig_c_curve, tau_BP)
  sig_c <- sig_c$y
  return(list(phi_c=phi_c,sig_c=sig_c))
}

calib_df <- load_calib_curve("intcal20")
baydem::vis_calib_curve(600,950,calib_df,xlab="Calendar Date [AD]",ylab="Fraction Modern")
phi_c_and_sig_c <- interpolate_cal_curve(800,calib_df)
points(800,phi_c_and_sig_c$phi_c,col="red",pch=19)
lines(c(800,800),phi_c_and_sig_c$phi_c+c(-1,1)*phi_c_and_sig_c$sig_c,col="red")

calc_likelihood <- function(tau,yrc_m,sig_yrc_m,calib_df) {
  # For the calendar date, tau, calculate the corresponding fraction modern
  # value and uncertainty using the calibration curve.
  phi_c_and_sig_c <- interpolate_cal_curve(tau,calib_df)

  phi_m <- exp(-yrc_m/8033)
  sig_m <- phi_m * sig_yrc_m / 8033

  # The overall uncertainty
  sig <- sqrt(sig_m^2 + phi_c_and_sig_c$sig_c^2)

  # The likelihood is the probability density of a normal distribution with
  # mean phi_c and standard deviation sig.
  return(dnorm(phi_m,phi_c_and_sig_c$phi_c,sig))
}

yrc_m     <- 1250
sig_yrc_m <-   20

tau_plot <- seq(600,950,by=1)
lik_plot <- rep(NA,length(tau_plot))
for(i in 1:length(tau_plot)) {
  lik_plot[i] <- calc_likelihood(tau_plot[i],yrc_m,sig_yrc_m,calib_df)
}
plot(tau_plot,lik_plot,xlab="Calendar Date [AD]",ylab="Likelihood",type="l")


# For a single date with no stratigraphic information, we can easily calculate
# posterior probability densities as a function of calendar date. That is not
# always possible if we wish to do more sophisticated things, such as use
# stratigraphic or other information to constrain a calibration. Let's build an
# MCMC sampler to sample from the posterior for the previous example (even
# though, in this case, we don't actually need to). We will use the "Swiss Army
# Knife" of MCMC algorithms, the Metropolis-Hastings algorithm. For a more
# sophisticated implementation of Metropolis-Hastings sampling (for which the
# "temperature" of sampling can be set) see enneal::do_mh_sampling_at_temp:
#
# https://github.com/MichaelHoltonPrice/enneal/blob/master/R/enneal.R

# (1) Set a starting parameter "vector" value
# Iterate over samples with a for loop, doing the following each time
# (2) Propose a new parameter vector (this requires defining a proposal
#     distribution)
# (3) Calculate the ratio of probabilities for the current and proposed
#     parameter vectors
# (4) Use an accept/reject step to choose which parameter vector to keep

#prop_sd <- 1
#tau_prop <- tau + rnorm(1,sd=prop_sd)
#alpha <- min(1,fprop/f)
