# METHODS USED IN WWA RAPID ATTRIBUTION STUDIES

suppressWarnings(suppressMessages({
    library(extRemes)
    library(abind)
    library(shape)
}))

###################################################################################################################
# Support functions

prep_window <- function(rc = c(1,1), w = 4, h = 4, ...) { options(repr.plot.width = rc[2]*w, repr.plot.height = rc[1]*h, repr.plot.res = 200); par(mfrow = rc, pch = 20, ...) }
load_ts <- function(fnm, col.names) { read.csv(fnm, comment.char = "#", sep = " ", header = F, col.names = col.names) }

###################################################################################################################

#' Compute nonstationary log-likelihood
#'
#' @param pars Vector of named parameters
#' @param cov1 Vector of covariates determining the nonstationary element
#' @param x Vector of values at which to compute the lof-likelihood
#' @param dist String determining parametric form of distribution (currently only 'norm' and 'gev' are accepted)
#' @param fittype String determining type of model to be fitted (currently only 'shift' and 'fixeddisp' are accepted): 'fixeddisp' is equivalent to the 'scale fit' implemented by the Climate Explorer
#' @return Scalar log-likelihood
#' @export
#'
ns_loglik <- function(pars, cov1, x, dist, fittype) {

    # compute nonstationary location & scale
    if(fittype == "fixeddisp") {
        const = exp((pars["alpha"] * cov1) / pars["mu0"])
        loc = pars["mu0"] * const
        scale = pars["sigma0"] * const
    } else if(fittype == "shift") {
        loc = pars["mu0"] + pars["alpha"] * cov1
        scale = pars["sigma0"]
    } else {
        print(paste(fittype, "not implemented"))
        return()
    }
    
    # constrain variance to be strictly positive
    if(any(scale <= 0)) return(NA)
        
    # return negative log-likelihood
    if(dist == "norm") {
        return(-sum(dnorm(x, mean = loc, sd = scale, log = T)))
    } else if(dist == "gev") {
        shape = pars["shape"]
        return(-sum(devd(x, loc = loc, scale = scale, shape = shape, log = T)))
    } else {
        print(paste(dist, "not implemented"))
        return()
    }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Estimate the parameters of a nonstationary distribution using max likelihood estimation
#'
#' @export
#'   
fit_ns <- function(dist, type = "fixeddisp", data, varnm, covnm_1, lower = F, mintemps = F, ev = NA, ...) {
    
    # currently only works for distributions fully specified by mean & sd: only tested for normal, lognormal
    if(! dist %in% c("norm", "gev")) {
        print("Not yet implemented: use norm or gev")
        return()
    }
    
    # if looking at minimum temperatures (or minima of negative values generally), need to flip data for stable model fitting
    x <- data[,varnm]; if(mintemps) x <- -x
    cov1 <- data[,covnm_1]
    
    # fit model with appropriate number of parameters, pad if necessary
    init <- c("mu0" = mean(x), "sigma0" = sd(x), "alpha" = 0)
    if(dist == "gev") init <- c(init, "shape" = 0)
    fitted <- suppressWarnings(optim(par = init, ns_loglik, cov1 = cov1, x = x, dist = dist, fittype = type))
    
    # if looking at minimum temperatures (or minima of negative values generally), so trend & location parameters have been flipped. This may cause some confusion so may have to modify later!
    if(mintemps) {
        fitted[["NOTE"]] <- "NB: model parameters are estimated for negative values"
        fitted$par["mu0"] <- -fitted$par["mu0"]
        fitted$par["alpha"] <- -fitted$par["alpha"]
        x <- -x
    }
            
    # attach assorted useful information
    fitted[["dist"]] <- dist
    fitted[["type"]] <- type
    fitted[["varnm"]] <- varnm
    fitted[["covnm_1"]] <- covnm_1
    fitted[["data"]] <- data
    fitted[["x"]] <- x
    fitted[["cov1"]] <- data[,covnm_1]
    
    fitted[["lower"]] <- lower               # saves having to specify every time later on
    fitted[["mintemps"]] <- mintemps         # look at maxima of 0-temps, rather than minima of observed temps
    
    if(is.na(ev)) { ev <- x[length(x)] } # event value: assume that event of interest is most recent, unless told otherwise (used in later plotting functions)
    fitted[["ev"]] <- ev

    return(fitted)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Re-estimate model parameters with new data
#'
#' @param mdl Fitted attribution model, as output by ns_fit()
#' @param new_df New data frame to be used in estimating model parameters. Must contain same column names as original data.
#' @export
#'
refit <- function(mdl, new_df) {
    fit_ns(dist = mdl$dist, type = mdl$type, data = new_df, varnm = mdl$varnm, covnm_1 = mdl$covnm_1, lower = mdl$lower, mintemps = mdl$mintemps, ev = mdl$ev)
}
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Get parameters of fitted distribution
#'
#' @export
#'   
ns_pars <- function(mdl, fixed_cov = NA) {
    
    # if no covariate value given, evaluate at all covariate values
    if(is.na(fixed_cov[1])) fixed_cov <- mdl$cov1
    
    # calculate the nonstationary parameter values
    pars <- mdl$par
    if(mdl$type == "fixeddisp") {
        
        ns_const = exp((pars["alpha"] * fixed_cov) / pars["mu0"])
        loc = pars["mu0"] * ns_const
        scale = pars["sigma0"] * ns_const
        
    } else if(mdl$type == "shift") {
        loc = pars["mu0"] + pars["alpha"] * fixed_cov
        scale = rep(pars["sigma0"], length(fixed_cov))

        
    } else {
        print(paste(mdl$type,"not implemented"))
        return()
    }
        
    # return the list of named parameters: location, scale, shape (if applicable)
    if("shape" %in% names(pars)) {
        return(lapply(list("loc" = loc, "scale" = scale, "shape" = rep(pars["shape"], length(fixed_cov))), unname))
    } else {
        return(lapply(list("loc" = loc, "scale" = scale), unname))
    }
    
}


    
#' Check sensitivity of parameter estimates to individual years by leave-one-out refitting
#' @param mdl Fitted attribution model, as output by ns_fit()
#' @rp Scalar return period to be used to compute relative change in intensity
#' @ev Scalar event value, to be used to compute return period & probability ratio
#' @cov Scalar value of covariate in factual world
#' @cov_cf Scalar value of covariate in counterfactual world
#'
#' @export
#'   
loo_pars <- function(mdl, check_year, rp = 100, ev, cov, cov_cf, plot = T, relative = T) {
    
    if(missing(ev)) { ev <- mdl$ev }
    if(missing(check_year)) { check_year <- length(mdl$x) }
    
    mdl_res <- c(mdl$par, 
                 "dI" = int_change(mdl, rp = rp, cov, cov_cf, relative = F),
                 "dI%" = int_change(mdl, rp = rp, cov, cov_cf, relative = T),
                 "rp" = return_period(mdl, ev, cov),
                 "pr" = prob_ratio(mdl, ev = ev, cov, cov_cf))
    
    mdl_df <- mdl$data
    
    loo_res <- sapply(1:nrow(mdl_df), function(i) {
        fit_i <- refit(mdl, mdl_df[-i,])
        
        c(fit_i$par, 
          "dI" = int_change(fit_i, rp = rp, cov, cov_cf, relative = F),
          "dI%" = int_change(fit_i, rp = rp, cov, cov_cf, relative = T),
          "rp" = return_period(fit_i, ev, cov),
          "pr" = prob_ratio(fit_i, ev = ev, cov, cov_cf))
    })
    loo_res <- cbind("all_years" = mdl_res, loo_res)
    if(!relative) loo_res <- loo_res[!grepl("%", rownames(loo_res)),]
    
    if(plot) {
            
        prep_window(c(1,nrow(loo_res)), w = 2)
        if(all(!is.finite(loo_res["pr",]))) { loo_res["pr",] <- 10e6 }
        
        invisible(sapply(row.names(loo_res), function(i) {
            boxplot(loo_res[i,-1], main = i, pch = 20)
            points(loo_res[i,1], pch = 21, bg = "darkgoldenrod2", cex = 1.4)
            points(loo_res[i,check_year+1], pch = 21, bg = "magenta", cex = 1.4)   # value if specified event omitted
            if(i == "pr") {
                abline(h = 1, lty = 2)
            } else {
                abline(h = 0, lty = 2)
            }
        }))
    } else {
        return(loo_res)
    }
}
    

###################################################################################################################

#' Map values from nonstationary distribution to stationary uniform distribution on (0,1)
#'
#' @export
#'   
map_to_u <- function(mdl, x, fixed_cov = NA) {
    
    pars <- ns_pars(mdl, fixed_cov = fixed_cov)
    if(missing(x)) x <- mdl$x
    
    # retrieve the actual fitted model parameters if they were flipped for fitting
    if(mdl$mintemps) {
        pars$loc <- -pars$loc
        x = -x
        mdl$lower <- !mdl$lower # also have to look at the opposite tail
    }
    
    # get exceedance probability
    if(mdl$dist == "norm") {
        pit <- pnorm(x, mean = pars$loc, sd = pars$scale, lower.tail = mdl$lower)
    } else if(mdl$dist == "gev") {
        pit <- sapply(1:length(pars$loc), function(i) pevd(x[i], loc = pars$loc[i], scale = pars$scale[i], shape = pars$shape[i], lower.tail = mdl$lower))
    } else {
        return(NULL)
    }
    return(pit)
}
                      
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Return period
#'
#' @export
#'   
return_period <- function(mdl, x, fixed_cov = NA) {
    
    1 / map_to_u(mdl, x, fixed_cov)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Probability ratio
#'
#' @export
#'   
prob_ratio <- function(mdl, ev, cov, cov_cf) {
    
    if(missing(ev)) ev <- mdl$ev
    
    ep_f <- map_to_u(mdl, ev, fixed_cov = cov)
    ep_cf <- map_to_u(mdl, ev, fixed_cov = cov_cf)
    
    ep_f / ep_cf
}


###################################################################################################################

#' Convert from exceedance probability to value under specific distribution
#'
#' @export
#'   
map_from_u <- function(u, mdl, fixed_cov = NA) {
    
    pars <- ns_pars(mdl, fixed_cov = fixed_cov)
    
    # retrieve the actual fitted model parameters if they were flipped for fitting
    if(mdl$mintemps) {
        pars$loc <- -pars$loc
        mdl$lower <- !mdl$lower       # also have to look at the opposite tail
    }
        
    # map quantile onto stationary distribution
    if(mdl$dist == "norm") {
        erl <- sapply(1:length(u), function(j) {
            qnorm(u[j], mean = pars$loc, sd = pars$scale, lower.tail = mdl$lower)
        })
    } else if(mdl$dist == "gev") {
        erl <- sapply(1:length(u), function(j) {
            sapply(1:length(pars$loc), function(i) {
                qevd(u[j], loc = pars$loc[i], scale = pars$scale[i], shape = pars$shape[i], lower.tail = mdl$lower)
            })
        })
    } else {
        return(NULL)
    }
                      
    # if parameters flipped for fitting, flip 'em back
    if(mdl$mintemps) erl <- -erl
                      
    return(erl)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Effective return level
#'
#' @export
#'   
eff_return_level <- function(rp, mdl, fixed_cov = NA) {
    
    map_from_u(1/rp, mdl, fixed_cov = fixed_cov)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Effective return level
#'
#' @export
#'   
int_change <- function(mdl, rp = NA, cov, cov_cf, relative = F) {
    
    if(is.na(rp)) {
        if(relative) {
            cat("Return period needed to calculate relative change")
            return(NA)
        } else {
            rp <- 10
        }
    }
    # if return period is less than 1, assume it's an exceedance probability and convert to a return period
    if(rp < 1) { rp <- 1/rp } else { rp <- rp }
    
    # get effective return levels
    rl <- eff_return_level(rp, mdl, fixed_cov = cov)
    rl_cf <- eff_return_level(rp, mdl, fixed_cov = cov_cf)
    
    # if variable is logged, convert to real values first
    if(substr(mdl$varnm, 1, 5) == "log10") {
        rl <- 10^rl
        rl_cf <- 10^rl_cf
    } else if (substr(mdl$varnm, 1, 3) == "log"){
        rl <- exp(rl)
        rl_cf <- exp(rl_cf)
    }
    
    if(relative) {
        (rl - rl_cf) / rl_cf * 100
    } else {
        rl - rl_cf
    }
}


###################################################################################################################

#' Return level plot
#'
#' @export
#'   
plot_returnlevels <- function(mdl, cov, cov_cf, ev, ylim = NA, pch = 20, ylab = NA, legend_pos = "topright", main = "", xlim = c(1,10000), legend_labels = c("Present climate", "Counterfactual climate"), seed = 42, nsamp = 500, ...) {
    
    x <- mdl$x
    if(missing(ev)) { ev <- mdl$ev }
    
    rp_x <- unique(c(seq(1.1,2,0.1), seq(2,100,1), seq(100,1000,10), seq(100,1000,100), seq(1000,10000,1000)))     # return periods at which to calculate values for curves
    rp_th <- 1/seq(1,0,length.out = length(x)+2)[2:(length(x)+1)]                                                  # quantiles to map against observations to check fit

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # calculate return periods & return levels
    
    rl_curve_pres <- map_from_u(1/rp_x, mdl, fixed_cov = cov)
    rl_curve_cf <- map_from_u(1/rp_x, mdl, fixed_cov = cov_cf)
    
    rl_obs_pres <- map_from_u(map_to_u(mdl), mdl, fixed_cov = cov)
    rl_obs_cf <- map_from_u(map_to_u(mdl), mdl, fixed_cov = cov_cf)
    
    rp_event_pres <- 1/map_to_u(mdl, ev, fixed_cov = cov)
    rp_event_cf <- 1/map_to_u(mdl, ev, fixed_cov = cov_cf)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # prep axes
    
    if(is.na(ylim[1])) { ylim <- range(pretty(c(x, rl_curve_pres, rl_curve_cf))) }
    if(is.na(ylab)) {ylab <- mdl$varnm}
    # if((substr(mdl$varnm,1,3) == "log") & (ylim[1] <= 0)) { ylim[1] <- 0.01 }
    
    # plot
    plot(0,type = "n", xlim = xlim, ylim = ylim, log = "x", xlab = "", ylab = "", main = main)
    mtext("Return period (years)", side = 1, line = 2.5, cex = par("cex"))
    mtext(ylab, side = 2, line = 2.5, cex = par("cex"))
    
    legend(legend_pos, legend = c(legend_labels, "Observed event"), col = c("firebrick", "blue", "magenta"), lty = 1, pch = c(pch,pch,NA), bty = "n")

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # return period curves
    lines(rp_x, rl_curve_pres, lwd = 2, col = "firebrick", lty = 1)       # present climate
    lines(rp_x, rl_curve_cf, lwd = 2, col = "blue", lty = 1)              # counterfactual
    
    # expected return periods vs return levels transformed to stationarity at that covariate value
    points(rp_th, sort(rl_obs_pres, decreasing = mdl$lower), col = "firebrick", pch = pch)      # present
    points(rp_th, sort(rl_obs_cf, decreasing = mdl$lower), col = "blue", pch = pch)             # counterfactual
    
    # horizontal line showing observed event, plus ticks showing return periods
    abline(h = ev, col = "magenta", lty = 2)
    suppressWarnings(rug(rp_event_pres, lwd = 3, col = "firebrick"))   # present
    suppressWarnings(rug(rp_event_cf, lwd = 3, col = "blue"))          # counterfactual
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Add confidence intervals to return periods
    
    if(!is.na(nsamp)) {
        x_ci <- c(5,10,20,50,100,200,500,1000,2000,5000,10000)
        set.seed(seed)
        
        mdl_df <- setNames(data.frame(mdl$x, mdl$cov1), c(mdl$varnm, mdl$covnm_1)) 
        boot_res <- sapply(1:nsamp, function(i) {
            boot_df <- mdl_df[sample(1:nrow(mdl_df), nrow(mdl_df), replace = T),]
            tryCatch({
                boot_mdl <- fit_ns(mdl$dist, mdl$type, boot_df, varnm = mdl$varnm, covnm_1 = mdl$covnm_1, lower = mdl$lower, mintemps = mdl$mintemps, ev = ev)
                c(map_from_u(1/x_ci, boot_mdl, fixed_cov = cov), map_from_u(1/x_ci, boot_mdl, fixed_cov = cov_cf))
            }, error = function(cond) {return(rep(NA, length(x_ci)*2))})
        })
        est_ci <- apply(boot_res, 1, quantile, c(0.025, 0.975), na.rm = T)
        
        # lines bounding confidence intervals
        # matplot(x_ci, t(est_ci[,1:length(x_ci)]), type = "l", lty = 1, lwd = 2, col = adjustcolor("firebrick", alpha = 0.3), add = T)
        # matplot(x_ci, t(est_ci[,-(1:length(x_ci))]), type = "l", lty = 1, lwd = 2, col = adjustcolor("blue", alpha = 0.3), add = T)
        
        # shaded region for confidence intervals
        polygon(x = c(x_ci, rev(x_ci)), y = c(est_ci[1,1:length(x_ci)], rev(est_ci[2,1:length(x_ci)])), density = NULL, border = NA, col = adjustcolor("firebrick", alpha = 0.1))
        polygon(x = c(x_ci, rev(x_ci)), y = c(est_ci[1,-(1:length(x_ci))], rev(est_ci[2,-(1:length(x_ci))])), density = NULL, border = NA, col = adjustcolor("blue", alpha = 0.1))
    }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' GMST trend plot
#'
#' @export
#'   
plot_gmsttrend <- function(mdl, cov, cov_cf, ev, ylim = NA, ylab = NA, legend_pos = "topleft", main = "", seed = 42, nsamp = 1000, xlab = "GMST anomaly") {

    if(is.na(ylab)) { ylab <- mdl$varnm}
    if(is.na(ylim[1])) { ylim <- range(pretty(mdl$x)) }
    if(missing(ev)) ev <- mdl$ev
    
    plot(mdl$cov1, mdl$x, pch = 20, main = main, xlab = "", ylab = "", ylim = ylim, xlim = range(c(mdl$cov1, cov, cov_cf)))
    mtext(xlab, side = 1, line = 2.5, cex = par("cex"))
    mtext(ylab, side = 2, line = 2.5, cex = par("cex"))
    
    points(cov, ev, col = "magenta", lwd = 2, pch = 0)
    
    # trend lines
    lines(mdl$cov1, ns_pars(mdl)$loc, lwd = 3, col = "black", lty = 1)
    lines(mdl$cov1, eff_return_level(6, mdl), col = "blue", lwd = 2, lty = 1)
    lines(mdl$cov1, eff_return_level(40, mdl), col = "blue", lwd = 1, lty = 1)
    
    # get confidence interval for mu'
    mdl_df <- mdl$data
    set.seed(seed)
    mu_ci <- apply(sapply(1:nsamp, function(i) {
        boot_df <- mdl_df[sample(1:nrow(mdl_df), nrow(mdl_df), replace = T),]
        tryCatch({
            boot_mdl <- fit_ns(mdl$dist, mdl$type, boot_df, varnm = mdl$varnm, covnm_1 = mdl$covnm_1, lower = mdl$lower, mintemps = mdl$mintemps, ev = ev)
            c("mu_ev" = ns_pars(boot_mdl, fixed_cov = cov)$loc,
              "mu_cf" = ns_pars(boot_mdl, fixed_cov = cov_cf)$loc)
        }, error = function(cond) {return(rep(NA, 2))})
    }), 1, quantile, c(0.025, 0.975), na.rm = T)
    
    # confidence interval & markers for mu' at factual & counterfactual covariates
    lines(rep(cov, 3), c(ns_pars(mdl, fixed_cov = cov)$loc, mu_ci[,"mu_ev"]), col = "black", lwd = 2, type = "o", pch = "_")
    lines(rep(cov_cf, 3), c(ns_pars(mdl, fixed_cov = cov_cf)$loc, mu_ci[,"mu_cf"]), col = "black", lwd = 2, type = "o", pch = "_")
    
    # add legend
    legend(legend_pos, legend = c("location", "1-in-6-year event", "1-in-40-year event"), lty = 1, col = c("black", "blue", "blue"), lwd = c(2,2,1))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Trend plot
#'
#' @export
#'   
plot_trend <- function(mdl, ev, ev_year, ylab = NA, legend_pos = "topleft", main = "", ylim = NA, ...) {
    
    if(is.na(ylab)) {ylab <- mdl$varnm}
    if(is.na(ylim[1])) { ylim <- range(pretty(mdl$x)) }
    if(missing(ev)) { ev <- mdl$ev }
    if(missing(ev_year)) { ev_year <- mdl$data$year[which.min(abs(mdl$x - ev))] }
    
    plot(mdl$data$year, mdl$x, type = "s", lwd = 2, col = adjustcolor("black", alpha = 0.5), xlab = "Year", ylab = ylab, main = main, ylim = ylim, ...)
    
    lines(mdl$data$year, ns_pars(mdl)$loc, col = "black", lwd = 2)
    matplot(mdl$data$year, eff_return_level(c(6,40), mdl), type = "l", lty = 1, add = T, col = "blue", lwd = c(2,1))
    
    points(ev_year, ev, col = "magenta", lwd = 2, pch = 0)
    
    # add legend
    legend(legend_pos, legend = c("location", "1-in-6-year event", "1-in-40-year event"), lty = 1, col = c("black", "blue", "blue"), lwd = c(2,2,1))
}


###################################################################################################################

#' Extract results from fitted model
#'
#' @export
#'   
mdl_ests <- function(mdl, cov, cov_cf, ev, rp = NA) {

    pars <- mdl$par
    disp <- unname(pars["sigma0"] / pars["mu0"])

    if(is.na(rp)) rp <- return_period(mdl, ev, fixed_cov = cov)

    if(is.finite(rp)) {
        c(pars,
          "disp" = disp, 
          "event_magnitude" = ev, 
          "return_period" = rp, 
          "PR" = prob_ratio(mdl, ev, cov, cov_cf),
          "dI_abs" = int_change(mdl, rp, cov, cov_cf, relative = F),
          "dI_rel" = int_change(mdl, rp, cov, cov_cf, relative = T))
    } else {
        return(rep(NA, 10))
    }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Confidence intervals for parameters from fitted model
#'
#' @export
#'   
# wrapper function to get bootstrapped confidence intervals for spreadsheet
boot_ci <- function(mdl, cov, cov_cf, ev = NA, rp = NA, seed = 42, nsamp = 500, dp = 5) {
    
    # get best estimate from the observed data
    if(is.na(ev)) ev  <- ev
    mdl_res <- mdl_ests(mdl, cov, cov_cf, ev, rp = rp)
    
    # get bootstrap sample
    set.seed(seed)    
    boot_res <- sapply(1:nsamp, function(i) {
        boot_df <- mdl$data[sample(1:nrow(mdl$data), replace = T),]
        tryCatch({
            boot_mdl <- fit_ns(mdl$dist, mdl$type, boot_df, varnm = mdl$varnm, covnm_1 = mdl$covnm_1, lower = mdl$lower, mintemps = mdl$mintemps, ev = ev)
            mdl_ests(boot_mdl, cov, cov_cf, ev, rp = rp)
        },
        error = function(cond) {return(rep(NA, 10))})
    })
    boot_qq <- t(rbind("bestimate" = mdl_res, apply(boot_res, 1, quantile, c(0.025, 0.975), na.rm = T)))
    if(!is.na(dp)) boot_qq <- round(boot_qq, dp)
    return(boot_qq)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Carry out attribution for climate model
#'
#' @export
#'   
cmodel_results <- function(mdl, rp = 10, cov_pres, cov_pi = NA, cov_fut = NA, nsamp = 5, seed = 42, di_relative = NA, y_start = 1979, y_now = 2023, y_fut = 2050) {
    
    set.seed(seed)
    
    # fill in missing parameters
    if(is.na(cov_pi)) cov_pi <- cov_pres - 1.2
    if(is.na(cov_fut)) cov_fut <- cov_pres + 0.8
    if(is.na(di_relative)) di_relative <- mdl$type == "fixeddisp"
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # FIT MODEL to three different subsets: for evaluation, attribution & projection
    
    df <- mdl$data
    
    # fit models
    mdl_eval <- fit_ns(mdl$dist, mdl$type, df[df$year >= y_start & df$year <= y_now,], varnm = mdl$varnm, covnm_1 = mdl$covnm_1, lower = mdl$lower, mintemps = mdl$mintemps, ev = mdl$ev)
    mdl_attr <- fit_ns(mdl$dist, mdl$type, df[df$year <= y_now,], varnm = mdl$varnm, covnm_1 = mdl$covnm_1, lower = mdl$lower, mintemps = mdl$mintemps, ev = mdl$ev)
    mdl_proj <- fit_ns(mdl$dist, mdl$type, df[df$year <= y_fut,], varnm = mdl$varnm, covnm_1 = mdl$covnm_1, lower = mdl$lower, mintemps = mdl$mintemps, ev = mdl$ev)
    
    # get return level to use for analysis
    event_rl <- eff_return_level(rp = rp, mdl_attr, fixed_cov = cov_pres)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Bootstrap each set of model results
    
    if(di_relative) { di_cnm <- "dI_rel" } else { di_cnm <- "dI_abs" }
    if(mdl$type == "fixeddisp") { key_par <- "disp" } else { key_par <- "sigma0" }
    if(mdl$dist == "gev") { key_par <- c(key_par, "shape") }
    
    ci_eval <- boot_ci(mdl_eval, cov = cov_pres, cov_cf = cov_pi, ev = event_rl, rp = rp, nsamp = nsamp)[c(key_par),,drop = F]
    ci_attr <- boot_ci(mdl_attr, cov = cov_pres, cov_cf = cov_pi, ev = event_rl, rp = rp, nsamp = nsamp)[c("PR", di_cnm),]
    ci_proj <- boot_ci(mdl_proj, cov = cov_pres, cov_cf = cov_fut, ev = event_rl, rp = rp, nsamp = nsamp)[c("PR", di_cnm),]
                   
    # invert future projections
    ci_proj["PR",] <- 1/ci_proj["PR",c(1,3,2)]
    ci_proj[di_cnm,] <- -ci_proj[di_cnm, c(1,3,2)]
    
    # flatten & rename
    ci_eval <- unlist(lapply(rownames(ci_eval), function(cnm) setNames(ci_eval[cnm,], paste("eval", gsub("_", "-", cnm), c("est", "lower", "upper"), sep = "_"))))
    ci_attr <- unlist(lapply(rownames(ci_attr), function(cnm) setNames(ci_attr[cnm,], paste("attr", gsub("_", "-", cnm), c("est", "lower", "upper"), sep = "_"))))
    ci_proj <- unlist(lapply(rownames(ci_proj), function(cnm) setNames(ci_proj[cnm,], paste("proj", gsub("_", "-", cnm), c("est", "lower", "upper"), sep = "_"))))
                         
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    res <- t(data.frame(c(ci_eval, "rp_value" = event_rl, ci_attr, ci_proj)))     
    rownames(res) <- paste0(mdl$varnm, " ~ ", mdl$covnm_1, " (rp", rp,")")
    return(res)
}