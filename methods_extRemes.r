suppressMessages({
    library(extRemes)
    library(MASS)
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

# function to get parameters of nonstationary GEV at given covariate level
sgev_pars <- function(mdl, covariate, burn.in = 499) {
    
    if(grepl("fixeddisp", mdl$type)) {
        
        pars <- mdl$results$par
        sf <- exp(pars["alpha"] * covariate / pars["mu0"])   # scaling factor
        
        loc = pars["mu0"] * sf
        scale = pars["sigma0"] * sf
        if(any(grepl("shape", names(mdl$results$par)))) { shape = pars["shape"] } else { shape = NA }
        
    } else if(grepl("shiftscale", mdl$type)) {
        
        pars <- mdl$results$par
        
        loc = pars["mu0"] + pars["alpha"] * covariate
        scale = pars["sigma0"] + pars["beta"] * covariate
        if(any(grepl("shape", names(mdl$results$par)))) { shape = pars["shape"] } else { shape = NA }
        
    } else if(grepl("shift", mdl$type)) {
        
        pars <- mdl$results$par
        
        loc = pars["mu0"] + pars["alpha"] * covariate
        scale = rep(pars["sigma0"], length(covariate))
        if(any(grepl("shape", names(mdl$results$par)))) { shape = pars["shape"] } else { shape = NA }
        
    } else if(grepl("scale", mdl$type)) {
        
        pars <- mdl$results$par
        
        loc = rep(pars["mu0"], length(covariate))
        scale = pars["sigma0"] + pars["alpha"] * covariate
        if(any(grepl("shape", names(mdl$results$par)))) { shape = pars["shape"] } else { shape = NA }
        
    } else if(grepl("shape", mdl$type)) {
        
        pars <- mdl$results$par
        
        loc = rep(pars["mu0"], length(covariate))
        scale = pars["sigma0"] + pars["alpha"] * covariate
        shape = pars["shape"] + pars["alpha"] * covariate
        
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

return_period <- function(mdl, value, covariate, lower = F) {
    
    pars <- sgev_pars(mdl, covariate)
    
    # if(grepl("lnorm", mdl$type)) { value <- log(value) }
    
    if(grepl("gamma", mdl$type)) {
        p <- pgamma((value - pars$loc)/pars$scale, shape = pars$shape, lower.tail = lower) 
    } else if(grepl("lnorm", mdl$type)) {
        p <- plnorm(value, mean = pars$loc, sd = pars$scale, lower.tail = lower)
    } else if(grepl("norm", mdl$type)) {
        p <- pnorm(value, mean = pars$loc, sd = pars$scale, lower.tail = lower)
    } else if(mdl$type == 'GEV_fixeddisp') {
        p <- pevd(value, loc = pars$loc, scale = pars$scale, shape = pars$shape, lower.tail = lower)
    } else {
        p <- pextRemes(mdl, value, lower.tail = lower)
    }
    
    return(1/p)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to extract return levels at given covariate (faster & more reliable than inbuilt methods)
return_level <- function(mdl, rp, covariate, lower = F) {
    
    if(grepl("gamma", mdl$type)) {
        pars <- sgev_pars(mdl, covariate = covariate)
        rl <- qgamma(1/rp, shape = pars$shape, lower.tail = lower) * pars$scale + pars$loc
    } else if(grepl("norm", mdl$type)) {
        pars <- sgev_pars(mdl, covariate = covariate)
        rl <- qnorm(1/rp, mean = pars$loc, sd = pars$scale, lower.tail = lower)
    } else if(grepl("fixeddisp", mdl$type)) {
        rl <- sapply(covariate, function(st) {
            pars <- sgev_pars(mdl, covariate = st)
            unname(rlevd(rp, loc = pars$loc, scale = pars$scale, shape = pars$shape))
        })
    } else { 
        if(lower) { print("Check rlevd behaviour for lower tails") }
        pars <- sgev_pars(mdl, covariate = covariate)
        rl <- rlevd(rp, loc = pars$loc, scale = pars$scale, shape = pars$shape)
    }
    
    if(grepl("lnorm", mdl$type)) rl <- exp(rl)
    return(rl)
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

delta_I <- function(mdl, rp, cov1, cov2, rel = F, lower = F) {
    
    rl1 <- return_level(mdl, rp, cov1, lower = lower)
    rl2 <- return_level(mdl, rp, cov2, lower = lower)

    if(rel) {
        return(unname((rl1 - rl2) / rl2) * 100)
    } else {
        return(unname(rl1 - rl2))
    }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# quickly compute probability ratio

prob_ratio <- function(mdl, value, cov1, cov2, lower = F) {
    
    pars1 <- sgev_pars(mdl, cov1)
    pars2 <- sgev_pars(mdl, cov2)
    
    if(grepl("gamma", mdl$type)) {
        p1 = pgamma((value - pars1$loc)/pars1$scale, shape = pars1$shape, lower.tail = lower)
        p2 = pgamma((value - pars2$loc)/pars2$scale, shape = pars2$shape, lower.tail = lower)
    } else if(grepl("lnorm", mdl$type)) {
        p1 = plnorm(value, mean = pars1$loc, sd = pars1$scale, lower.tail = lower)
        p2 = plnorm(value, mean = pars2$loc, sd = pars2$scale, lower.tail = lower)
    } else if(grepl("norm", mdl$type)) {
        p1 = pnorm(value, mean = pars1$loc, sd = pars1$scale, lower.tail = lower)
        p2 = pnorm(value, mean = pars2$loc, sd = pars2$scale, lower.tail = lower)
    } else if(mdl$type == 'GEV_fixeddisp') {
        p1 = pevd(value, loc = pars1$loc, scale = pars1$scale, shape = pars1$shape, lower.tail = lower)
        p2 = pevd(value, loc = pars2$loc, scale = pars2$scale, shape = pars2$shape, lower.tail = lower)
    } else {
        p1 = pextRemes(mdl, q = value, qcov = event_qcov(mdl, cov1), lower.tail = lower)
        p2 = pextRemes(mdl, q = value, qcov = event_qcov(mdl, cov2), lower.tail = lower)
    }
    return(unname(p1/p2))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fit_results <- function(mdl, event_value, cov1, cov2, lower = F, dI_rel = F) {
        
    # method to compute useful results from fitted model and output as dataframe
    rp_event <- return_period(mdl, event_value, cov1, lower = lower)
    rp_alt = return_period(mdl, event_value, cov2, lower = lower)
    pr <- rp_alt / rp_event
    
    dI <- delta_I(mdl, rp_event, cov1, cov2, rel = dI_rel, lower = lower)
    
    res_list <- list(mdl = mdl$type, converged = mdl$results$convergence, event_value = event_value, 
                             alpha = unname(mdl$results$par["alpha"]), "loglik" = mdl$results$value,
                             rp_event = rp_event, rp_alt = rp_alt, pr = pr, delta_I = dI)
    
    if(grepl("shiftscale", mdl$type)) {res_list <- append(res_list, list(beta = unname(mdl$results$par["alpha"])), after = 3)}
    
    return(as.data.frame(res_list))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nonstationary GEV with fixed dispersion (as fitted by climate explorer)                                            ####
    
gev_fixeddisp <- function(pars = c(mu0, sigma0, shape, alpha), covariate, x) {
    
    loc = pars["mu0"] * exp(pars["alpha"] * covariate / pars["mu0"])
    scale = pars["sigma0"] * exp(pars["alpha"] * covariate / pars["mu0"])
    shape = pars["shape"]
    
    # return negative log-likelihood to be minimised
    return(-sum(devd(x, loc = loc, scale = scale, shape = shape, log = T)))
}

fevd_fixeddisp <- function(x, covariate, data, method = "MLE", solver = "L-BFGS-B", ...) {
    
    res <- list("results" = optim(par = c(mu0 = mean(data[,x]), sigma0 = sd(data[,x]), shape = 0, alpha = 0), gev_fixeddisp, 
                                  covariate = data[,covariate], x = data[,x], method = solver, ...))
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
    
    if(mdl$type == 'gamma_fixeddisp') { lss <- list(loc = loc, scale = scale, shape = unname(pars["shape"])) } 
    else if(mdl$type == 'normal_fixeddisp') { lss <- list(loc = loc, scale = scale) } 
    else if(mdl$type == 'GEV_fixeddisp') { lss <- list(loc = loc, scale = scale, shape = unname(pars["shape"])) }
        
    if(is.na(covariate)) {
        return(lss)
    } else {
        return(lapply(lss, "[", which.min(abs(cov_list - covariate))))
    }
}

stransf <- function(mdl, covariate = NA, lower = F) {
    
    # method to transform values to equivalent in stationary distribution
    
    # parameters of nonstationary distribution
    ns_pars <- sgev_pars(mdl, mdl$cov.data[,mdl$cov.name])
    
    # parameters of stationary distribution (set to standard form if not provided)
    if(is.na(covariate)) {
        s_pars <- list(loc = 0, scale = 1, shape = 1)
    } else {
        s_pars <- sgev_pars(mdl, covariate)
    }
    
    if(grepl("gamma", mdl$type)) {
        
        # use inverse probability transform
        pit <- sort(pgamma((mdl$x - ns_pars$loc) / ns_pars$scale, shape = ns_pars$shape, lower.tail = lower), decreasing = T)
        svalue <- qgamma(pit, s_pars$shape, lower.tail = lower) * s_pars$scale + s_pars$loc
        return(svalue)
        
    } else if(grepl("norm", mdl$type)) {
        
        # use inverse probability transform
        pit <- sort(pnorm(mdl$x, mean = ns_pars$loc, sd = ns_pars$scale, lower.tail = lower), decreasing = T)
        svalue <- qnorm(pit, mean = s_pars$loc, sd = s_pars$scale, lower.tail = lower)
        return(svalue)
        
    } else if(mdl$type == 'GEV_fixeddisp') {
        
        # use inverse probability transform
        pit <- sort(sapply(1:length(mdl$x), function(i) {
            pars <- sgev_pars(mdl, mdl$cov.data[i,mdl$cov.name])
            pevd(mdl$x[i], loc = pars$loc, scale = pars$scale, shape = pars$shape)
        }), decreasing = T)
        svalue <- qevd(pit, loc = s_pars$loc, scale = s_pars$scale, shape = s_pars$shape, lower.tail = lower)
        return(svalue)
        
    } else {
        return(trans(mdl))
    }
    
}

#=======================================================================================================================
# FITTING MODELS                                                                                                    ####

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nonstationary normal distribution with fixed dispersion (as fitted by climate explorer)

norm_fixeddisp <- function(pars = c(mu0, sigma0, alpha), covariate, x) {
    
    loc = pars["mu0"] * exp(pars["alpha"] * covariate / pars["mu0"])
    scale = pars["sigma0"] * exp(pars["alpha"] * covariate / pars["mu0"])
    
    # return negative log-likelihood to be minimised
    return(-sum(dnorm(x, mean = loc, sd = scale, log = T)))
}
      
norm_shift <- function(pars = c(mu0, sigma0, alpha), covariate, x) {
    
    loc = pars["mu0"] + pars["alpha"] * covariate
    scale = pars["sigma0"]
    
    # return negative log-likelihood to be minimised
    return(-sum(dnorm(x, mean = loc, sd = scale, log = T)))
}

norm_scale <- function(pars = c(mu0, sigma0, alpha), covariate, x) {
    
    loc = pars["mu0"]
    scale = pars["sigma0"] + pars["alpha"] * covariate
    
    # return negative log-likelihood to be minimised
    return(-sum(dnorm(x, mean = loc, sd = scale, log = T)))
}
    
    
norm_shiftscale <- function(pars = c(mu0, sigma0, alpha, beta), covariate, x) {
    
    loc = pars["mu0"] + pars["alpha"] * covariate
    scale = pars["sigma0"] + pars["beta"] * covariate
    
    # return negative log-likelihood to be minimised
    return(-sum(dnorm(x, mean = loc, sd = scale, log = T)))
}
    
    
fnorm <- function(x, covariate, data, type = "shift", method = "MLE", optim.method = "L-BFGS-B", init = NA, ...) {
    
    mtype <- paste0("norm_", type)
    fun <- get(mtype)
    
    if(is.na(init[1])) { init <- c(mu0 = mean(data[,x]), sigma0 = sd(data[,x]), alpha = 0) }
    
    if((type == "shiftscale") & !all(grepl("beta", init))) { init <- c(init, beta = 0) }
        
    # need to sort out a better way to estimate starting parameters
    res <- list("results" = optim(par = init, fun, covariate = data[,covariate], x = data[,x]), method = optim.method, lower = c(-Inf, 0, -Inf), ...)
    res[["type"]] <- mtype
    res[["dist"]] <- "norm"
    res[["x"]] <- data[,x]
    res[["cov.data"]] <- data
    res[["cov.name"]] <- covariate
    res[["var.name"]] <- x
    
    return(res)
}

    
       
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nonstationary gamma distribution (will always have fixed dispersion)

gamma_fixeddisp <- function(pars = c(rate, shape, alpha), covariate, x) {
    
    # gamma distribution has fixed coefficient of variation, so shouldn't be necessary to fit location + scale: instead, let rate vary with covariate
    ns_rate = pars["rate"] + pars["alpha"] * covariate
    shape = pars["shape"]
    
    # return negative log-likelihood to be minimised
    return(-sum(dgamma(x, shape = shape, rate = ns_rate, log = T)))
}

    
    
fgamma <- function(x, covariate, data, method = "MLE", type = "fixeddisp", optim.method = "Nelder-Mead", init = NA, ...) {
        
    mtype <- "gamma_fixeddisp"
    
    if(is.na(init[1])) init <- c(fitdistr(data[,x], "gamma")$estimate, "alpha" = 0)
            
    # need to sort out a better way to estimate starting parameters
    res <- list("results" = suppressWarnings(optim(par = init, gamma_fixeddisp, covariate = data[,covariate], x = data[,x], method = optim.method)),
                fittype = method,
                method = optim.method, 
                type = mtype,
                x = data[,x],
                cov.data = data,
                cov.name = covariate,
                var.name = x,
                ...)
    
    return(res)
}

    
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nonstationary normal distribution with fixed dispersion (as fitted by climate explorer)

lnorm_fixeddisp <- function(pars = c(mu0, sigma0, shape, alpha), covariate, x) {
    
    loc = pars["mu0"] * exp(pars["alpha"] * covariate / pars["mu0"])
    scale = pars["sigma0"] * exp(pars["alpha"] * covariate / pars["mu0"])
    
    # return negative log-likelihood to be minimised
    return(-sum(dlnorm(x, mean = loc, sd = scale, log = T)))
}
      
lnorm_shift <- function(pars = c(mu0, sigma0, shape, alpha), covariate, x) {
    
    loc = pars["mu0"] + pars["alpha"] * covariate
    scale = pars["sigma0"]
    
    # return negative log-likelihood to be minimised
    return(-sum(dlnorm(x, mean = loc, sd = scale, log = T)))
}

lnorm_scale <- function(pars = c(mu0, sigma0, shape, alpha), covariate, x) {
    
    loc = pars["mu0"]
    scale = pars["sigma0"] + pars["alpha"] * covariate
    
    # return negative log-likelihood to be minimised
    return(-sum(dlnorm(x, mean = loc, sd = scale, log = T)))
}
    
    
lnorm_shiftscale <- function(pars = c(mu0, sigma0, shape, alpha, beta), covariate, x) {
    
    loc = pars["mu0"] + pars["alpha"] * covariate
    scale = pars["sigma0"] + pars["beta"] * covariate
    
    # return negative log-likelihood to be minimised
    return(-sum(dlnorm(x, mean = loc, sd = scale, log = T)))
}
    
    
flnorm <- function(x, covariate, data, type = "shift", method = "MLE", optim.method = "Nelder-Mead", init = NA, ...) {
    
    mtype <- paste0("lnorm_", type)
    fun <- get(mtype)
    
    if(is.na(init)) { init <- c(mu0 = mean(log(data[,x])), sigma0 = sd(log(data[,x])), alpha = 0) }
    
    if((type == "shiftscale") & !all(grepl("beta", init))) { init <- c(init, beta = 0) }
    
    # need to sort out a better way to estimate starting parameters
    res <- list("results" = optim(par = init, fun, covariate = data[,covariate], x = data[,x]), method = optim.method, ...)
    res[["type"]] <- mtype
    res[["x"]] <- data[,x]
    res[["cov.data"]] <- data
    res[["cov.name"]] <- covariate
    res[["var.name"]] <- x
    
    return(res)
}
    
    
#=======================================================================================================================
# PLOTTING METHODS                                                                                                  ####
    
plot_gmsttrend <- function(mdl, event_year = "2022", lower = F, xlab = "GMST", ylab = "value", legend_pos = "topright", ...) {
    
    event_year = toString(event_year)
    
    covariate <- mdl$cov.data[,mdl$cov.name]
       
    plot(covariate, mdl$x, pch = 20, xlab = xlab, ylab = ylab, ...)
    points(df[event_year, c(mdl$cov.name, mdl$var.name)], col = "red", lwd = 3)
    lines(covariate, sgev_pars(mdl, covariate)$loc, lwd = 2)
    lines(covariate, return_level(mdl, 6, covariate, lower = lower), col = "blue", lwd = 2)
    lines(covariate, return_level(mdl, 40, covariate, lower = lower), col = "blue", lwd = 1) 
    
    legend(legend_pos, legend = c("location", "1-in-6-year event", "1-in-40-year event"), lty = 1, col = c("black", "blue", "blue"), lwd = c(2,2,1))
}
    

plot_returnperiods <- function(mdl, cov1, cov2, event_value, lower = F, ylim = NA, pch = 20, ylab = "value", legend_pos = "topright", main = "", ...) {
    
    rp_x <- unique(c(seq(1.1,2,0.1), seq(2,100,1), seq(100,1000,10), seq(100,1000,100), seq(1000,10000,1000))) # return periods at which to calculate values
    rp_th <- 1/seq(1,0,length.out = length(mdl$x)+2)[2:(length(mdl$x)+1)]                                      # theoretical return periods
    
    if(is.na(ylim[1])) { ylim <- range(pretty(mdl$x)) }
    if(min(ylim) <= 0) { ylim <- c(0.01, max(ylim)) }
    
    # prep axes
    plot(0,type = "n", xlim = c(1,10000), ylim = ylim, log = "x", xlab = "Return period (years)", ylab = ylab, main = main)
    
    # return period curves
    lines(rp_x, return_level(mdl, rp_x, cov1, lower = lower), lwd = 2, col = "firebrick")
    lines(rp_x, return_level(mdl, rp_x, cov2, lower = lower), lwd = 2, col = "blue")
    
    # expected return periods vs return levels transformed to stationarity at that covariate value
    points(rp_th, stransf(mdl, cov1, lower = lower), col = "firebrick", pch = pch)
    points(rp_th, stransf(mdl, cov2, lower = lower), col = "blue", pch = pch)
    
    # horizontal line showing observed event
    abline(h = event_value, col = "magenta")
    suppressWarnings(rug(return_period(mdl, event_value, cov1, lower = lower), lwd = 3, col = "firebrick"))
    suppressWarnings(rug(return_period(mdl, event_value, cov2, lower = lower), lwd = 3, col = "blue"))
    
    # abline(v = return_period(mdl, event_value, cov1, lower = lower))
    
    legend(legend_pos, legend = c("2022 GMST", "2022 GMST -1.2", "Observed event"), col = c("firebrick", "blue", "magenta"), lty = 1, pch = c(20,20,NA), bty = "n")
}