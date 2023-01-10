suppressMessages({
    library(extRemes)
    library(plyr) 
    library(beepr)
})

#=======================================================================================================================
## METHODS FOR THE ANALYSIS OF NONSTATIONARY EXTREMES

# function to construct covariate matrix at required covariate value
event_qcov <- function(mdl, covariate) {
    
    if(mdl$type == 'normal_fixeddisp') {
        return(t(as.matrix(c(mu = 1, sigma = 1, alpha = covariate * unname(mdl$results$par["alpha"]), threshold = 0.925))))
    }
    
    if(mdl$type == 'GEV_fixeddisp') {
        return(t(as.matrix(c(mu = 1, sigma = 1, xi = 1, alpha = covariate * unname(mdl$results$par["alpha"]), threshold = 0.925))))
    } else {
        return(make.qcov(mdl) + (1-make.qcov(mdl)) * covariate )
    }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to get parameters of fitted model (with correct names)
model_pars <- function(mdl) {
    
    if((mdl$type %in% c('GEV_fixeddisp', 'normal_fixeddisp')) || (mdl$method == "MLE")) {
        
        return(mdl$results$par)
    } else {
        par_est = ci(mdl, 0.05, "parameter")[,2]
        
        if(!"sigma0" %in% names(par_est)) {
            
            if("phi0" %in% names(par_est)) { par_est["scale"] = exp(par_est["phi0"]) } 
            if("log.scale" %in% names(par_est))  par_est["sigma0"] = exp(par_est["log.scale"]) 
            if("scale" %in% names(par_est)) par_est["sigma0"] = par_est["scale"]
            
        }
        
        if(!"sigma1" %in% names(par_est)) {
            if("phi1" %in% names(par_est)) { par_est["sigma1"] = (par_est["phi1"]) } else { par_est["sigma1"] = NA }
        }
        
        if(!"mu0" %in% names(par_est)) {
            par_est["mu0"] <- unname(par_est["location"])
            par_est["mu1"] <- NA
        }
        return(par_est[c("mu0", "mu1", "sigma0", "sigma1", "shape")])
    } 
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to get parameters of stationary GEV at given covariate level
sgev_pars <- function(mdl, covariate, burn.in = 499) {
    
    if(mdl$type == 'normal_fixeddisp') {
        
        pars <- mdl$results$par
        sf <- exp(pars["alpha"] * covariate / pars["mu0"])   # scaling factor
        loc = pars["mu0"] * sf
        scale = pars["sigma0"] * sf
        shape = NA
        
    } else if(mdl$type == 'GEV_fixeddisp') {
        
        pars <- mdl$results$par
        sf <- exp(pars["alpha"] * covariate / pars["mu0"])   # scaling factor
        loc = pars["mu0"] * sf
        scale = pars["sigma0"] * sf
        shape = pars["shape"]
        
    } else {
        
        pars <- strip(mdl, burn.in = burn.in)
        cov <- event_qcov(mdl, covariate)[1:length(pars)]
        cpars <- pars * cov
        
        if("location" %in% names(cpars)) {
            loc <- cpars["location"]
        } else {
            loc <- sum(cpars[grep("mu", names(cpars))])
        }
        
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
    }
        
    return(list("loc" = unname(loc), "scale" = unname(scale), "shape" = unname(shape)))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to extract return periods at given covariate (allows compatibility with fixed-disp methods)

return_period <- function(mdl, value, covariate) {
    
    if(mdl$type == 'normal_fixeddisp') {
        
        pars <- sgev_pars(mdl, covariate)
        rp <- 1/(1-pnorm(value, mean = pars$loc, sd = pars$scale))
        
    } else if(mdl$type == 'GEV_fixeddisp') {
        
        pars <- sgev_pars(mdl, covariate)
        rp <- 1/(1-pevd(value, loc = pars$loc, scale = pars$scale, shape = pars$shape))

    } else {
        rp <- 1/(1-pextRemes(mdl, value))
    }
        
    return(rp)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to extract return levels at given covariate (faster & more reliable than inbuilt methods)
return_level <- function(mdl, rp, covariate) {
    
    if(mdl$type == 'normal_fixeddisp') {
        pars <- sgev_pars(mdl, covariate = covariate)
        rl <- qnorm(1-1/rp, mean = pars$loc, sd = pars$scale)
        return(rl)
    } else { 
        pars <- sgev_pars(mdl, covariate = covariate)
        rl <- rlevd(rp, loc = pars$loc, scale = pars$scale, shape = pars$shape)
        return(rl)
    }

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# bootstrapped confidence intervals
rl.boot.ci <- function(mdl, x, covariate_value, ci = 95, nsamp = 5000, seed = 1) {
    
    # needs updating to accommodate the fixed-dispersion model
    set.seed(seed)
    
    n <- length(mdl$x)
    
    # get quantiles
    if (ci > 1) { ci = ci/100 }
    qq = c((1-ci) / 2, 1 - (1-ci) / 2)
    
    # bootstrap the return level
    boot_sample <- sapply(1:nsamp, function(i) {
        
        # resample the source data, run the model again using the bootstrapped data
        boot_df <- mdl$cov.data[sample(1:n,n,replace = T),]
        
        if(mdl$type == "normal_fixeddisp") {
            boot_fit <- fnorm_fixeddisp(mdl$var.name, mdl$cov.name, data = boot_df)
        } else if(mdl$type == "GEV_fixeddisp") {
            boot_fit <- fevd_fixeddisp(mdl$var.name, mdl$cov.name, data = boot_df)
        } else {
            boot_fit <- suppressWarnings(update(mdl, data = boot_df))
        }
        
        # currently ignoring warnings, but should handle them properly at some point
        # provided methods to calculate return level can be slow, esp for Bayesian 
        tryCatch(return_level(boot_fit, x, event_gmst), return(NA))
          
        # possible alternative method - although this seems even slower
        # pars = lapply(findpars(boot_fit, qcov = event_qcov(boot_fit, covariate_value)), "sum")
        # rl = rlevd(x, loc = pars$location, scale = pars$scale, shape = pars$shape)
    })  
    # extract quantiles & return
    return(t(apply(boot_sample, 1, quantile, qq)))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# quickly compute change in intensity

delta_I <- function(mdl, rp, cov1, cov2, rel = F) {
    
    rl1 <- return_level(mdl, rp, cov1)
    rl2 <- return_level(mdl, rp, cov2)
    if(rel) {
        return(unname((rl1 - rl2) / rl2) * 100)
    } else {
        return(unname(rl1 - rl2))
    }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# quickly compute probability ratio

prob_ratio <- function(mdl, value, cov1, cov2) {
    
    if(mdl$type == 'normal_fixeddisp') {
        pars <- mdl$results$par
        ep1 = 1-pnorm(value, mean = pars["mu0"] * exp(pars["alpha"] * cov1 / pars["mu0"]), sd = pars["sigma0"] * exp(pars["alpha"] * cov1 / pars["mu0"]))
        ep2 = 1-pnorm(value, mean = pars["mu0"] * exp(pars["alpha"] * cov2 / pars["mu0"]), sd = pars["sigma0"] * exp(pars["alpha"] * cov2 / pars["mu0"]))
    } else if(mdl$type == 'GEV_fixeddisp') {
        pars <- mdl$results$par
        ep1 = 1-pevd(value, loc = pars["mu0"] * exp(pars["alpha"] * cov1 / pars["mu0"]), scale = pars["sigma0"] * exp(pars["alpha"] * cov1 / pars["mu0"]), shape = pars["xi"])
        ep2 = 1-pevd(value, loc = pars["mu0"] * exp(pars["alpha"] * cov2 / pars["mu0"]), scale = pars["sigma0"] * exp(pars["alpha"] * cov2 / pars["mu0"]), shape = pars["xi"])
    } else {
        ep1 = 1-pextRemes(mdl, q = value, qcov = event_qcov(mdl, cov1))
        ep2 = 1-pextRemes(mdl, q = value, qcov = event_qcov(mdl, cov2))
    }
    return(unname(ep1/ep2))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nonstationary GEV with fixed dispersion (as fitted by climate explorer)

gev_fixeddisp <- function(pars = c(mu0, sigma0, shape, alpha), covariate, x) {
    
    loc = pars["mu0"] * exp(pars["alpha"] * covariate / pars["mu0"])
    scale = pars["sigma0"] * exp(pars["alpha"] * covariate / pars["mu0"])
    shape = pars["shape"]
    
    # return negative log-likelihood to be minimised
    return(-sum(devd(x, loc = loc, scale = scale, shape = shape, log = T)))
}

fevd_fixeddisp <- function(x, covariate, data, method = "MLE", ...) {
    
    res <- list("results" = optim(par = c(mu0 = mean(data[,x]), sigma0 = sd(data[,x]), shape = 0, alpha = 0), gev_fixeddisp, 
                                  covariate = data[,covariate], x = data[,x]))
    res[["type"]] <- "GEV_fixeddisp"
    res[["x"]] <- data[,x]
    res[["cov.data"]] <- data
    res[["cov.name"]] <- covariate
    res[["var.name"]] <- x
    
    return(res)
}

# location, shape & scale parameters for fixed-dispersion model
fd_lss <- function(mdl, covariate = NA) {
    
    cov_list <- mdl$cov.data[,mdl$cov.name]
    pars <- mdl$results$par
    
    sf <- exp(pars["alpha"] * cov_list / pars["mu0"])   # scaling factor
    loc = pars["mu0"] * sf
    scale = pars["sigma0"] * sf
    
    if(mdl$type == 'normal_fixeddisp') { lss <- list(loc = loc, scale = scale) } 
    else if(mdl$type == 'GEV_fixeddisp') { lss <- list(loc = loc, scale = scale, shape = unname(pars["shape"])) }
        
    if(is.na(covariate)) {
        return(lss)
    } else {
        return(lapply(lss, "[", which.min(abs(cov_list - covariate))))
    }
}

trans_fd <- function(mdl) {
    
    if(mdl$type == 'normal_fixeddisp') {
        
    } else if(mdl$type == 'GEV_fixeddisp') {
        pars <- fd_lss(mdl)
    return(log(1 + (mdl$x - pars$loc) * pars$shape / pars$scale) / pars$shape)
    } else {
        return(trans(mdl))
    }
    
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nonstationary normal distribution with fixed dispersion (as fitted by climate explorer)

norm_fixeddisp <- function(pars = c(mu0, sigma0, alpha), covariate, x) {
    
    loc = pars["mu0"] * exp(pars["alpha"] * covariate / pars["mu0"])
    scale = pars["sigma0"] * exp(pars["alpha"] * covariate / pars["mu0"])
    
    # return negative log-likelihood to be minimised
    return(-sum(dnorm(x, mean = loc, sd = scale, log = T)))
}

fnorm_fixeddisp <- function(x, covariate, data, method = "MLE", ...) {
    
    res <- list("results" = optim(par = c(mu0 = mean(data[,x]), sigma0 = sd(data[,x]), alpha = 0), norm_fixeddisp, 
                                  covariate = data[,covariate], x = data[,x]))
    res[["type"]] <- "normal_fixeddisp"
    res[["x"]] <- data[,x]
    res[["cov.data"]] <- data
    res[["cov.name"]] <- covariate
    res[["var.name"]] <- x
    
    return(res)
}