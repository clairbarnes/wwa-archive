suppressMessages({
    library(extRemes)
    library(plyr) 
})

#=======================================================================================================================
## METHODS FOR THE ANALYSIS OF NONSTRATIONARY EXTREMES

# function to construct covariate matrix at required covariate value
event_qcov <- function(mdl, covariate) { make.qcov(mdl) + (1-make.qcov(mdl)) * covariate }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to get parameters of fitted model (with correct names)
model_pars <- function(mdl) {
    
    par_est = ci(mdl, 0.05, "parameter")[,2]
    
    if(!"sigma0" %in% names(par_est)) {
        
        if("phi0" %in% names(par_est)) { par_est["scale"] = exp(par_est["phi0"]) } 
        if("log.scale" %in% names(par_est))  par_est["sigma0"] = exp(par_est["log.scale"]) 
        if("scale" %in% names(par_est)) par_est["sigma0"] = par_est["scale"]
        
    }
    
    if(!"sigma1" %in% names(par_est)) {
        
        if("phi1" %in% names(par_est)) { par_est["sigma1"] = (par_est["phi1"]) } else { par_est["sigma1"] = NA }
        
    }
    
    return(par_est[c("mu0", "mu1", "sigma0", "sigma1", "shape")])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to get parameters of stationary GEV at given covariate level
sgev_pars <- function(mdl, covariate, burn.in = 499) {
    
    pars <- strip(mdl, burn.in = burn.in)
    cov <- event_qcov(mdl, covariate)[1:length(pars)]
    cpars <- pars * cov
    
    loc <- sum(cpars[grep("mu", names(cpars))])
    
    if(any(grepl("phi", names(cpars)))) {
        scale <- exp(sum(cpars[grep("phi", names(cpars))]))
    } else {
        if(any(grepl("sigma", names(cpars)))) {
            scale <- sum(cpars[grep("sigma", names(cpars))])
        } else {
            scale <- sum(cpars[grep("scale", names(cpars))])
        }
    }
    shape <- sum(cpars[grep("shape", names(cpars))])
    
    list("loc" = loc, "scale" = scale, "shape" = shape)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to extract return levels at given covariate (faster & more reliable than inbuilt methods)
return_level <- function(mdl, rp, covariate) {
    
    pars <- sgev_pars(mdl, covariate = covariate)
    rl <- rlevd(rp, loc = pars$loc, scale = pars$scale, shape = pars$shape)
    return(rl)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# bootstrapped confidence intervals
rl.boot.ci <- function(mdl, x, covariate_value, ci = 95, nsamp = 5000, seed = 1) {
    
    set.seed(seed)
    
    n <- length(mdl$x)
    
    # get quantiles
    if (ci > 1) { ci = ci/100 }
    qq = c((1-ci) / 2, 1 - (1-ci) / 2)
    
    # bootstrap the return level
    boot_sample <- sapply(1:nsamp, function(i) {
        
        # resample the source data, run the model again using the bootstrapped data
        boot_df <- mdl$cov.data[sample(1:n,n,replace = T),]
        boot_fit <- suppressWarnings(update(mdl, data = boot_df))
        
        # currently ignoring warnings, but should handle them properly at some point
        # provided methods to calculate return level can be slow, esp for Bayesian 
        return_level(boot_fit, x, covariate_value)
                        
        # possible alternative method - although this seems even slower
        # pars = lapply(findpars(boot_fit, qcov = event_qcov(boot_fit, covariate_value)), "sum")
        # rl = rlevd(x, loc = pars$location, scale = pars$scale, shape = pars$shape)
    })  
    # extract quantiles & return
    return(t(apply(boot_sample, 1, quantile, qq)))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# quickly compute change in intensity

delta_I <- function(mdl, rp, cov1, cov2) {
    unname(return_level(mdl, rp, cov1) - return_level(mdl, rp, cov2))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# quickly compute probability ratio

prob_ratio <- function(mdl, value, cov1, cov2) {
    ep1 = 1-pextRemes(mdl, q = value, qcov = event_qcov(mdl, cov1))
    ep2 = 1-pextRemes(mdl, q = value, qcov = event_qcov(mdl, cov2))
    return(unname(ep1/ep2))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~