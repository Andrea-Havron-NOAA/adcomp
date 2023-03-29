library(TMB)

compile("tmb_examples/tweedie_osa.cpp")
dyn.load(dynlib("tmb_examples/tweedie_osa"))

#Simulate data
set.seed(123)
n <- 100
data <- list(y = tweedie::rtweedie(n, mu = 2, 
                                   phi = 1.2, 
                                   power = 1.25))
parameters <- list(mu=1, ln_phi=0, logit_p = 0)


obj <- MakeADFun(data, parameters, DLL="tweedie_osa")
opt <- nlminb(obj$par, obj$fn, obj$gr)
## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
rep <- obj$report()
idx.0 <- data$y==0 
#evaluate 0 observations with a Poisson where lambda is defined based on the properties of the Compound Poisson-Gamma 
res.1 <- oneStepPredict(obj, observation.name = "y", 
                        data.term.indicator = "keep",
                        method = "cdf", seed = 123, 
                        discrete = TRUE)
res.2 <- oneStepPredict(obj, observation.name = "y", 
                        data.term.indicator = "keep",
                        method = "oneStepGeneric", seed = 123, 
                        discrete = FALSE)
#nlcdf.lower.raw <- ifelse(idx.0, res.1$nlcdf.lower, res.2$nlcdf.lower) 
nlcdf.lower <- ifelse(idx.0, res.1$nlcdf.lower, 
                      -log(
                           exp(-res.1$nlcdf.lower) +
                           exp(-res.2$nlcdf.lower)
                        )) 
nlcdf.upper <- ifelse(idx.0, res.1$nlcdf.upper, res.2$nlcdf.upper)

Fx.raw <- 1 / ( 1 + exp(nlcdf.lower - nlcdf.upper) )
Fx <- Fx.raw
set.seed(123)
Fx[idx.0] <- runif(sum(idx.0), 0, 
                          max(Fx[idx.0], na.rm = TRUE))
qqnorm(qnorm(Fx));abline(0,1)

#compare to analytical quantile residuals - not exact
cdf <- tweedie::ptweedie(data$y, mu = opt$par["mu"], 
                         phi = summary(sdr, "report")[1,1], 
                         power = summary(sdr, "report")[2,1])
r <- cdf
set.seed(123)
r[idx.0] <- runif(sum(idx.0), 0, max(cdf[idx.0]))

plot(r[idx.0], Fx[idx.0]);abline(0,1)
plot(r[!idx.0], Fx[!idx.0]);abline(0,1)
plot(qnorm(Fx), qnorm(r));abline(0,1)
