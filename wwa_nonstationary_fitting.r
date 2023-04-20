library("copula")
library("plyr")

###################################################################################################################
# Couple of quick support functions for repetetive tasks

prep_window <- function(rc = c(1,1), w = 4, h = 4) { options(repr.plot.width = rc[2]*w, repr.plot.height = rc[1]*h, repr.plot.res = 200); par(mfrow = rc, pch = 20) }
load_ts <- function(fnm, col.names) { read.csv(fnm, comment.char = "#", sep = " ", header = F, col.names = col.names) }

###################################################################################################################
# nonstationary log-likelihood function with two parameters

ns_mle_2cov <- function(pars = c(mu0, sigma0, alpha, beta), cov1, cov2, x, dist, fittype) {
    
    if(fittype == "fixeddisp") {
        
        const = exp((pars["alpha"] * cov1 + pars["beta"] * cov2) / pars["mu0"])
        loc = pars["mu0"] * const
        scale = pars["sigma0"] * const
        
    } else if(fittype == "shift") {
        loc = pars["mu0"] + pars["alpha"] * cov1 + pars["beta"] * cov2
        scale = pars["sigma0"]
        
    } else {
        print(paste(fittype,"not implemented"))
        return()
    }
    
    # return negative log-likelihood to be minimised
    ddist <- get(paste0("d", dist))
    return(-sum(ddist(x, mean = loc, sd = scale, log = T)))
}


# nonstationary log-likelihood function with one parameter
ns_mle <- function(pars = c(mu0, sigma0, alpha), cov1, x, dist, fittype) {

    # using separate function to avoid optimising a redundant parameter, which can lead to fitting issues
    
    if(fittype == "fixeddisp") {
        
        const = exp((pars["alpha"] * cov1) / pars["mu0"])
        loc = pars["mu0"] * const
        scale = pars["sigma0"] * const
        
    } else if(fittype == "shift") {
        loc = pars["mu0"] + pars["alpha"] * cov1
        scale = pars["sigma0"]
        
    } else {
        print(paste(fittype,"not implemented"))
        return()
    }
    
    # return negative log-likelihood to be minimised
    ddist <- get(paste0("d", dist))
    return(-sum(ddist(x, mean = loc, sd = scale, log = T)))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# wrapper to fit nonstationary fixed-dispersion model

fit_ns <- function(dist, type = "fixeddisp", data, varnm, covnm_1, covnm_2 = NA, lower = F, event_index = NA, ...) {
    
    # currently only works for distributions fully specified by mean & sd: only tested for normal, lognormal
    if(! dist %in% c("norm", "lnorm")) {
        print("Not yet implemented: use norm or lnorm")
        return()
    }
    
    x <- data[,varnm]
    if(is.na(event_index)) { event_index <- length(x) } # assume that year of interest is most recent, unless told otherwise
    
    # fit model with appropriate number of parameters, pad if necessary
    if(is.na(covnm_2)) {
        init <- c("mu0" = mean(x), "sigma0" = sd(x), "alpha" = 0)
        fitted <- suppressWarnings(optim(par = init, ns_mle, cov1 = data[,covnm_1], x = x, dist = dist, fittype = type, ...))
        fitted[["par"]]["beta"] <- 0
        mdl_call <- paste0(varnm," ~ ",covnm_1)
    } else {
        init <- c("mu0" = mean(x), "sigma0" = sd(x), "alpha" = 0, "beta" = 0)
        fitted <- suppressWarnings(optim(par = init, ns_mle_2cov, cov1 = data[,covnm_1], cov2 = data[,covnm_2], x = x, dist = dist, fittype = type, ...))
        mdl_call <- paste0(varnm," ~ ",covnm_1, " + ", covnm_2)
    }
        
    # estimate & return parameters & assorted useful information
    fitted[["dist"]] <- dist
    fitted[["type"]] <- type
    fitted[["lower"]] <- lower
    fitted[["call"]] <- mdl_call
    fitted[["varnm"]] <- varnm
    fitted[["covnm_1"]] <- covnm_1
    fitted[["covnm_2"]] <- covnm_2
    fitted[["ev_idx"]] <- event_index
    fitted[["x"]] <- x
    fitted[["cov1"]] <- data[,covnm_1]
    if(is.na(covnm_2)) {
        fitted[["cov2"]] <- 0
    } else {
        fitted[["cov2"]] <- data[,covnm_2]
    }

    return(fitted)
}

###################################################################################################################

# get nonstationary parameters from fitted model
ns_pars <- function(mdl, cov1 = NA, cov2 = 0) {
    
    pars <- mdl$par
    if(is.na(cov1[1])) cov1 <- mdl$cov1
    if(is.na(cov2[1])) cov2 <- mdl$cov2
    
    if(mdl$type == "fixeddisp") {
        
        ns_const = exp((pars["alpha"] * cov1 + pars["beta"] * cov2) / pars["mu0"])
        loc = pars["mu0"] * ns_const
        scale = pars["sigma0"] * ns_const
        
    } else if(mdl$type == "shift") {
        loc = pars["mu0"] + pars["alpha"] * cov1 + pars["beta"] * cov2
        scale = pars["sigma0"]
        
    } else {
        print(paste(mdl$type,"not implemented"))
        return()
    }
    return(list("loc" = loc, "scale" = scale))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# convert from nonstationary distribution to stationary uniform (PIT)
map_to_u <- function(mdl, x, cov1 = NA, cov2 = 0) {
    
    pars <- ns_pars(mdl, cov1 = cov1, cov2 = cov2)
    if(missing(x)) x <- mdl$x
    pit <- get(paste0("p", mdl$dist))(x, mean = pars$loc, sd = pars$scale, lower.tail = mdl$lower)
    
    return(pit)
}


# convert from exceedance probability to value conditioned on stationary distribution with specified covariates
map_from_u <- function(u, mdl, cov1, cov2 = 0) {
    
    pars <- ns_pars(mdl, cov1 = cov1, cov2 = cov2)
        
    u_transformed <- get(paste0("q", mdl$dist))(u, mean = pars$loc, sd = pars$scale, lower.tail = mdl$lower)
    return(u_transformed)
}


# wrapper for probability ratio
prob_ratio <- function(mdl, ev, cov1, cov1_cf, cov2 = 0, cov2_cf = 0) {
    
    ep_f <- map_to_u(mdl, ev, cov1 = cov1, cov2 = cov2)
    ep_cf <- map_to_u(mdl, ev, cov1 = cov1_cf, cov2 = cov2_cf)
    
    ep_f / ep_cf
}


# wrapper for change in intensity
Delta_I <- function(mdl, rp, cov1, cov1_cf, cov2 = 0, cov2_cf = 0, relative = T) {
    
    rl <- map_from_u(1/rp, mdl, cov1 = cov1, cov2 = cov2)
    rl_cf <- map_from_u(1/rp, mdl, cov1 = cov1_cf, cov2 = cov2_cf)
    
    if(substr(mdl$varnm, 1, 5) == "log10") {
        rl <- 10^rl
        rl_cf <- 10^rl_cf
    }
    
    if(relative) {
        (rl - rl_cf) / rl_cf * 100
    } else {
        rl - rl_cf
    }
}

###################################################################################################################

# compute joint exceedances over regular grid for easy plotting
copula_mesh <- function(mdl_x, mdl_y, copula, cov1, cov2 = 0, xrange, yrange, n = 32) {
    
    # compute joint exceedances over regular grid for easy plotting
    
    if(missing(xrange)) xrange <- range(pretty(mdl_x$x))
    if(missing(yrange)) yrange <- range(pretty(mdl_y$x))
    
    # define the regular mesh for plotting
    x_mesh <- seq(xrange[1], xrange[2],length.out = n)
    y_mesh <- seq(yrange[1], yrange[2], length.out = n)
    
    # convert the regular mesh to U space
    x_umesh <- map_to_u(mdl_x, x_mesh, cov1 = cov1, cov2 = cov2)
    y_umesh <- map_to_u(mdl_y, y_mesh, cov1 = cov1, cov2 = cov2)
    
    return(list("x" = x_mesh, "y" = y_mesh, "z" = sapply(y_umesh, function(y) sapply(x_umesh, function(x) pCopula(cbind(x,y), copula)))))
}
                                                                                     
###################################################################################################################
                                                                                     
# Get return periods, return levels etc from a single fitted model                                  
model_res <- function(mdl_x, x, cov1_hist, cov2_hist = 0, dI_rel = F) {
    
    u_x <- map_to_u(mdl_x)[mdl_x$x == x]
    uhist_x <- map_to_u(mdl_x, x, cov1 = cov1_hist, cov2 = cov2_hist)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # return levels
    
    rl_hist_x = map_from_u(u_x, mdl_x, cov1 = cov1_hist, cov2 = cov2_hist)
    
    if(dI_rel) { dI_x = (x - rl_hist_x) / rl_hist_x * 100 } else { dI_x = x - rl_hist_x }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # return periods 
    
    rp_2022_x = 1/u_x
    rp_hist_x = 1/uhist_x
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if(mdl_x$type == "fixeddisp") { vpar_x <- mdl_x$par["sigma0"] / mdl_x$par["mu0"] } else { vpar_x <- mdl_x$par["sigma0"] }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    res <- data.frame("dist" = c(mdl_x$dist),
                      "fit_type" = c(mdl_x$type),
                      # "converged" = c(mdl_x$convergence),
                      "mu0" = c(mdl_x$par["mu0"]),
                      "sigma0" = c(mdl_x$par["sigma0"]),
                      "alpha" = c(mdl_x$par["alpha"]),
                      "beta" = c(mdl_x$par["beta"]),
                      "var/disp" = c(vpar_x),
                      "rl_obs" = c(x),
                      "rl_cf" = c(rl_hist_x),
                      "delta_I" = c(dI_x),
                      "delta_I_rel" = c(dI_rel),
                      "rp_obs" = c(rp_2022_x),
                      "rp_cf" = c(rp_hist_x),
                      "prob_ratio" = c(rp_hist_x / rp_2022_x),
                      row.names = c(mdl_x$call))
    return(res)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                     
# Get return periods, return levels etc from joint model (currently limited to stationary copulas)
jmodel_res <- function(mdl_x, mdl_y, copula, x, y, cov1_hist, cov2_hist = 0, dI_x_rel = F, dI_y_rel = F) {
    
    res_x <- model_res(mdl_x, x, cov1_hist = cov1_hist, cov2_hist = cov2_hist, dI_rel = dI_x_rel)
    res_y <- model_res(mdl_y, y, cov1_hist = cov1_hist, cov2_hist = cov2_hist, dI_rel = dI_y_rel)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # return periods for joint model
    
    u_x <- map_to_u(mdl_x, x)[mdl_x$x == x]
    u_y <- map_to_u(mdl_y, y)[mdl_y$x == y]
        
    uhist_x <- map_to_u(mdl_x, x, cov1 = cov1_hist, cov2 = cov2_hist)
    uhist_y <- map_to_u(mdl_y, y, cov1 = cov1_hist, cov2 = cov2_hist)
    
    rp_2022_joint = 1/pCopula(cbind(u_x, u_y), copula)
    rp_hist_joint = 1/pCopula(cbind(uhist_x, uhist_y), copula)
    
    res_joint <- data.frame("dist" = c(class(cfit)),
                            "fit_type" = c("stationary"),
                            # "convergence" = c(NA),
                            "mu0" = c(NA),
                            "sigma0" = c(NA),
                            "alpha" = c(NA),
                            "beta" = c(NA),
                            "var/disp" = c(getSigma(copula)[1,2]),
                            "rl_obs" = c(NA),
                            "rl_cf" = c(NA),
                            "delta_I" = c(NA),
                            "delta_I_rel" = c(NA),
                            "rp_obs" = c(rp_2022_joint),
                            "rp_cf" = c(rp_hist_joint),
                            "prob_ratio" = c(rp_hist_joint / rp_2022_joint),
                           row.names = "joint")
    res <- rbind(res_x, res_y, res_joint)
    return(res)
}
                                                                                     
###################################################################################################################
                                                                                     
# Return period plots                                                                                  
plot_returnperiods <- function(mdl, cov1, cov1_cf, cov2 = 0, cov2_cf = 0, event_value = NA, ylim = NA, pch = 20, ylab = NA, legend_pos = "topright", main = "", xlim = c(1,10000),
                               add = F, col_f = "firebrick", col_cf = "blue", lty = 1, seed = 1, nsamp = 1000, ...) {
    
    x <- mdl$x
    rp_x <- unique(c(seq(1.1,2,0.1), seq(2,100,1), seq(100,1000,10), seq(100,1000,100), seq(1000,10000,1000))) # return periods at which to calculate values
    rp_th <- 1/seq(1,0,length.out = length(x)+2)[2:(length(x)+1)]                                      # theoretical return periods
    if(is.na(event_value)) { event_value <- x[mdl$ev_idx] }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    rl_curve_pres <- map_from_u(1/rp_x, mdl, cov1 = cov1, cov2 = cov2)
    rl_curve_cf <- map_from_u(1/rp_x, mdl, cov1 = cov1_cf, cov2 = cov2_cf)
    
    rl_obs_pres <- map_from_u(map_to_u(mdl), mdl, cov1 = cov1, cov2 = cov2)
    rl_obs_cf <- map_from_u(map_to_u(mdl), mdl, cov1 = cov1_cf, cov2 = cov2_cf)
    
    rp_event_pres <- 1/map_to_u(mdl, event_value, cov1 = cov1, cov2 = cov2)
    rp_event_cf <- 1/map_to_u(mdl, event_value, cov1 = cov1_cf, cov2 = cov2_cf)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # prep axes
    if(is.na(ylim[1])) { ylim <- range(pretty(x)) }
    if(is.na(ylab)) {ylab <- mdl$varnm}
    if((substr(mdl$varnm,1,3) == "log") & (ylim[1] <= 0)) { ylim[1] <- 0.01 }
    
    if(!add) {
        plot(0,type = "n", xlim = xlim, ylim = ylim, log = "x", xlab = "", ylab = "", main = main, ...)
    
        mtext("Return period (years)", side = 1, line = 2.5, cex = par("cex"))
        mtext(ylab, side = 2, line = 2.5, cex = par("cex"))
    }
    
    # return period curves
    lines(rp_x, rl_curve_pres, lwd = 2, col = col_f, lty = lty)          # present
    lines(rp_x, rl_curve_cf, lwd = 2, col = col_cf, lty = lty)         # counterfactual
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Add confidence intervals to return periods
    
    x_ci <- c(5,10,20,50,100,200,500,1000,2000,5000,10000)
    set.seed(seed)
    
    mdl_df <- setNames(data.frame(mdl$x, mdl$cov1, mdl$cov2), c(mdl$varnm, mdl$covnm_1, mdl$covnm_2)) 
    boot_res <- sapply(1:nsamp, function(i) {
        boot_df <- mdl_df[sample(1:nrow(mdl_df), nrow(mdl_df), replace = T),]
        boot_mdl <- fit_ns(mdl$dist, mdl$type, boot_df, varnm = mdl$varnm, covnm_1 = mdl$covnm_1, covnm_2 = mdl$covnm_2, lower = mdl$lower)
        c(map_from_u(1/x_ci, boot_mdl, cov1 = cov1, cov2 = cov2), map_from_u(1/x_ci, boot_mdl, cov1 = cov1_cf, cov2 = cov2_cf))
    })
    boot_ci <- apply(boot_res, 1, quantile, c(0.025, 0.975), na.rm = T)
    
    matplot(x_ci, t(boot_ci[,1:length(x_ci)]), type = "l", lty = 1, col = adjustcolor(col_f, alpha = 0.5), add = T)
    matplot(x_ci, t(boot_ci[,-(1:length(x_ci))]), type = "l", lty = 1, col = adjustcolor(col_cf, alpha = 0.5), add = T)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # expected return periods vs return levels transformed to stationarity at that covariate value
    points(rp_th, sort(rl_obs_pres, decreasing = mdl$lower), col = col_f, pch = pch)     # present
    points(rp_th, sort(rl_obs_cf, decreasing = mdl$lower), col = col_cf, pch = pch)          # counterfactual
    
    # horizontal line showing observed event
    abline(h = event_value, col = "magenta", lty = 2)
    suppressWarnings(rug(rp_event_pres, lwd = 3, col = col_f))   # present
    suppressWarnings(rug(rp_event_cf, lwd = 3, col = col_cf))          # counterfactual
            
    legend(legend_pos, legend = c("2022 GMST", "2022 GMST -1.2", "Observed event"), col = c(col_f, col_cf, "magenta"), lty = lty, pch = c(pch,pch,NA), bty = "n")
}                                                                              

         
# Plots of trend in GMST
plot_gmsttrend <- function(mdl, cov1, cov1_cf, cov2 = 0, ylim = NA, ylab = NA, legend_pos = "topleft", main = "", seed = 1, nsamp = 1000, add = F, col_main = "black", col_rl = "blue", lty = 1, ...) {
    
    if(is.na(ylab)) { ylab <- mdl$varnm}
    if(is.na(ylim[1])) { ylim <- range(pretty(mdl$x)) }
    
    # assumes that first covariate is GMST
    if(!add) {
        plot(mdl$cov1, mdl$x, pch = 20, main = main, xlab = "", ylab = "", ylim = ylim, xlim = range(c(mdl$cov1, cov1, cov1_cf)))
        mtext("GMST anomaly", side = 1, line = 2.5, cex = par("cex"))
        mtext(ylab, side = 2, line = 2.5, cex = par("cex"))
    }
    
    points(mdl$cov1[mdl$ev_idx], mdl$x[mdl$ev_idx], col = "magenta", lwd = 2, pch = 0)
    
    # trend lines
    lines(mdl$cov1, ns_pars(mdl, cov2 = cov2)$loc, lwd = 3, col = col_main, lty = lty)
    lines(mdl$cov1, map_from_u(1/6, mdl, cov1 = mdl$cov1, cov2 = cov2), col = col_rl, lwd = 2, lty = lty)
    lines(mdl$cov1, map_from_u(1/40, mdl, cov1 = mdl$cov1, cov2 = cov2), col = col_rl, lwd = 1, lty = lty)
    
    # get confidence interval for mu'
    mdl_df <- setNames(data.frame(mdl$x, mdl$cov1, mdl$cov2), c(mdl$varnm, mdl$covnm_1, mdl$covnm_2))
    set.seed(seed)
    mu_ci <- apply(sapply(1:nsamp, function(i) {
        boot_df <- mdl_df[sample(1:nrow(mdl_df), nrow(mdl_df), replace = T),]
        boot_mdl <- fit_ns(mdl$dist, mdl$type, boot_df, varnm = mdl$varnm, covnm_1 = mdl$covnm_1, covnm_2 = mdl$covnm_2, lower = mdl$lower)
        c("mu_ev" = ns_pars(boot_mdl, cov1 = cov1, cov2 = cov2)$loc,
          "mu_cf" = ns_pars(boot_mdl, cov1 = cov1_cf, cov2 = cov2)$loc)
    }), 1, quantile, c(0.025, 0.975))
    
    # confidence interval & markers for mu' at present & counterfactual covariates
    lines(rep(cov1, 2), mu_ci[,"mu_ev.mu0"], col = col_main, lwd = 2, pch = "_")
    lines(rep(cov1_cf, 2), mu_ci[,"mu_cf.mu0"], col = col_main, lwd = 2)
    
    points(rep(cov1, 3), c(mu_ci[,"mu_ev.mu0"], ns_pars(mdl, cov1 = cov1, cov2 = cov2)$loc), col = col_main, lwd = 3, pch = "_")
    points(rep(cov1_cf, 3), c(mu_ci[,"mu_cf.mu0"], ns_pars(mdl, cov1 = cov1_cf, cov2 = cov2)$loc), col = col_main, lwd = 3, pch = "_")
    
    # add legend
    legend(legend_pos, legend = c("location", "1-in-6-year event", "1-in-40-year event"), lty = 1, col = c(col_main, col_rl, col_rl), lwd = c(2,2,1))
}