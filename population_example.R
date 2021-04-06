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

# Create the sampling grid
tau_min  <- 600  # AD
tau_max  <- 1200 # AD
dtau     <- 1
tau      <- seq(tau_min,tau_max,by=dtau)

# Create three notional Gaussian mixture models for the Bayesian inference
# illustrations.
#              pi1 pi2  mu1  mu2 s1 s2
th1 <- cbind(c(1.0,0.0, 800, 800,40,40))
th2 <- cbind(c(0.5,0.5, 900,1000,30,40))
th3 <- cbind(c(0.0,1.0,1000,1000,40,40))

# Calculate the probability densities
f1 <- baydem::calc_gauss_mix_pdf(th1,tau,tau_min=tau_min,tau_max=tau_max)
f2 <- baydem::calc_gauss_mix_pdf(th2,tau,tau_min=tau_min,tau_max=tau_max)
f3 <- baydem::calc_gauss_mix_pdf(th3,tau,tau_min=tau_min,tau_max=tau_max)

# Set the priors
prior1 = .4
prior2 = .2
prior3 = .4

# Set the likelihoods (in practice -- i.e., when not for an illustration -- 
# these likelihoods would be calculated from observed data)
lik1 <- .01
lik2 <- 1
lik3 <- .02

# Calculate the posteriors
post1 <- lik1 * prior1 / (lik1 * prior1 + lik2*prior2 + lik3*prior3)
post2 <- lik2 * prior2 / (lik1 * prior1 + lik2*prior2 + lik3*prior3)
post3 <- lik3 * prior3 / (lik1 * prior1 + lik2*prior2 + lik3*prior3)


# Set up variables for plotting
th <- 3 # Thickness of lines in plots
sz <- 14 # Character size in output PDF files
ppi <- 600 # points per inch for plots

segmentLength <- 4
secondaryColor <- "gray62"
plotPdf <- T
plotWidth <- 1600
plotHeight <- 1600
plotUnits <- 'px'
plotRes <- 400
plotMaxF <- .01

tau_range <- c(tau_min,tau_max)
# Create the actual plots
pdf(file="bayesian_update_illustration_th1.pdf",pointsize=sz)
par(mar=rep(0,4))
plot(tau,f1,col='black',lwd=th,lty=1,type='l',ylim=c(0,.01),xlim=tau_range,xlab='Calendar Date [AD]',ylab='Density',xaxt='n', ann=F,yaxt='n')
dev.off()

pdf(file="bayesian_update_illustration_th2.pdf",pointsize=sz)
par(mar=rep(0,4))
plot(tau,f2,col='black',lwd=th,lty=1,type='l',ylim=c(0,.01),xlim=tau_range,xlab='Calendar Date [AD]',ylab='Density',xaxt='n', ann=F,yaxt='n')
dev.off()

pdf(file="bayesian_update_illustration_th3.pdf",pointsize=sz)
par(mar=rep(0,4))
plot(tau,f3,col='black',lwd=th,lty=1,type='l',ylim=c(0,.01),xlim=tau_range,xlab='Calendar Date [AD]',ylab='Density',xaxt='n', ann=F,yaxt='n')
dev.off()

fprior <- prior1 * f1 + prior2 * f2 + prior3 * f3
pdf(file="bayesian_update_illustration_prior.pdf",pointsize=sz)
par(mar=rep(0,4))
plot(tau,fprior,col='black',lwd=th,lty=1,type='l',ylim=c(0,.01),xlim=tau_range,xlab='Calendar Date [AD]',ylab='Density',xaxt='n', ann=F,yaxt='n')
text(900,.008,'Mean across Priors',cex=3)
dev.off()

fpost <- post1 * f1 + post2 * f2 + post3 * f3
pdf(file="bayesian_update_illustration_posterior.pdf",pointsize=sz)
par(mar=rep(0,4))
plot(tau,fpost,col='black',lwd=th,lty=1,type='l',ylim=c(0,.01),xlim=tau_range,xlab='Calendar Date [AD]',ylab='Density',xaxt='n', ann=F,yaxt='n')
text(900,.008,'Mean across Posteriors',cex=3)
dev.off()

pdf(file="bayesian_update_illustration_blank.pdf",pointsize=sz)
par(mar=rep(0,4))
plot(NULL,NULL,col='red',lwd=th,lty=1,type='l',ylim=c(0,.01),xlim=tau_range,xlab='Calendar Date [AD]',ylab='Density',xaxt='n', ann=F,yaxt='n',frame.plot=F)
dev.off()
