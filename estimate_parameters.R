
## Load libraries
library(truncnorm) # for truncated normal pdf
library(xtable) # for exporting parameter estimates

## Set seed
set.seed(20211117)

#### Import data
igr <- read.csv("data/incubator growth rate data - Sheet1.csv")
# clean up and get rid of individuals in dried leaves or killed by parasites
igr <- igr[-c(4, 46, 68, 80, 115, 138, 144, 145),]

# Temperature response function
getTPCval <- function(Topt, Gmax, a, b, Temp){  
	lnGR = Gmax * exp(-exp((b * (Temp - Topt)) - 6) - a * (Temp - Topt)^2)
	return(lnGR)
}

## Create functions that returns the negative log-likelihood
nll <- function(param, data, pdf, pos, logi, r.censor, big.T){
    pos.regex <- paste(pos, collapse = "|") # collapse "pos" to regex
    pos.tr <- grep(pos.regex, names(param)) # parameters to transform
    param.tr1 <- replace(param, pos.tr, exp(unlist(param[pos.tr]))) # transform (exp)
    logi.regex <- paste(logi, collapse = "|") # collapse "logi" to regex
    logi.tr <- grep(logi.regex, names(param)) # parameters to transform
    param.tr2 <- replace(param.tr1, logi.tr, exp(unlist(param[logi.tr])) / (1 + exp(unlist(param[logi.tr])))) # transform (logit)
    likestor <- dim(data)
    for(i in 1:length(data)){
        likestor[i] <-  log(do.call(pdf, append(list(x = data[i], obs.length = r.censor[i], Temp = big.T[i]),
                                                param.tr2)))
    }
    likestor[is.infinite(likestor)] <- -10000 # Punish 0 likelihoods
    return(-(sum(likestor)))
}

## Probability function
probdf <- function(x, obs.length, Temp, lambda1, lambda2, prob, Topt, Gmax, a.par, b.par, sd.par, maxTPC = 900){
    meanTPC <-  getTPCval(Topt, Gmax, a.par, b.par, Temp)
    if(meanTPC > maxTPC){
        return(0)
    } else{
        int0 <- integrate(f = function(H){
            exp(log(dtruncnorm(x, a = 0, b = maxTPC, mean = (obs.length - H) * meanTPC, sd = sd.par)) + dexp(H, rate = lambda1, log = TRUE))}, lower = 0, upper = obs.length)$value + ifelse(x > 0, 0, exp(-lambda1 * obs.length))
        int1 <- integrate(f = function(H){
            exp(log(dtruncnorm(x, a = 0, b = maxTPC, mean = H * meanTPC, sd = sd.par)) + dexp(H, rate = lambda2, log = TRUE))}, lower = 0, upper = obs.length)$value
        int1b <- dtruncnorm(x, a = 0, b = maxTPC, mean = obs.length * meanTPC, sd = sd.par) * (exp(-lambda2 * obs.length))
        marg <- int0 * prob + (int1 + int1b) * (1 - prob)
        return(marg)
    }
}

## Starting values for parameters
param <- list(lambda1 = log(.12),
              lambda2 = log(.04),
              prob = log(.27 / (1 - .27)),
              Topt = log(29.2),
              Gmax = log(60),
              a.par = log(0.005),
              b.par = log(.74),
              sd.par = log(40))
### Parameters that must be positive
positives <- c("lambda1", "lambda2", "Topt", "Gmax", "a.par", "b.par", "sd.par")
### Parameters on logit scale
logis <- c("prob")

## Fit model using optim
## Method = Broyden–Fletcher–Goldfarb–Shanno algorithm (unconstrained non-linear optimization; approximates Hessian)
fitdata.1 <- optim(param,
                   fn = nll,
                   data = igr[, "growth.mm2"] / igr[, "mass.start"],
                   r.censor = igr[, "duration.num"],
                   big.T = igr[, "Temperature.C"],
                   method = c("BFGS"),
                   hessian = TRUE,
                   pdf = probdf,
                   pos = positives,
                   logi = logis,
                   control = list(maxit = 10000, trace = TRUE))
### Summarize fit
fitdata.1

## Parameter estimates (untransformed)
estdata.1 <- c(exp(fitdata.1$par[1]),
               exp(fitdata.1$par[2]),
               exp(fitdata.1$par[3]) / (1 + exp(fitdata.1$par[3])),
               exp(fitdata.1$par[4]),
               exp(fitdata.1$par[5]),
               exp(fitdata.1$par[6]),
               exp(fitdata.1$par[7]),
               exp(fitdata.1$par[8]))
### Print
estdata.1

### Save results
write.table(t(estdata.1), file="output/estimates.csv", sep = ",", quote = FALSE, row.names = FALSE)

## ## OPTIONAL
## ## Estimate the standard error of the estimates by taking the inverse of the Hessian (producing covar matrix)
## fisher_info<-solve(fitdata.1$hessian)
## prop_sigma<-sqrt(diag(fisher_info))
## interval<-data.frame(value=fitdata.1$par, standarderror = prop_sigma)

## ## Format table of estimates for export to tex file
## names(estdata.1)[3] <- c("Pr. of initially molting")
## estdata.export <- as.table(estdata.1[c(3,1:2, 4:8)])
## ### Export
## print(xtable(estdata.export, digits = 4), type = "latex", file = "Estimate table.tex")
