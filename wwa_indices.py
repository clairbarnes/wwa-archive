
## MODULE TO COMPUTE MROE COMPLEX CLIMATE INDICES - TO BE USED ALONGSIDE WWA.PY

import sys; sys.path.append('/home/clair/wwa'); from wwa import *

from xclim.indices.stats import fit
from scipy.stats import norm, gamma

import lmoments3.distr                     # conda activate xclim; pip install git+https://github.com/OpenHydrology/lmoments3.git
glo = getattr(lmoments3.distr, "glo")      # generalised logistic distribution - needed to fit SPEI

###############################################################################################################
## COMPUTE INDICES

def get_spi(pr, months = range(1,13), calibration_period = slice("1980", "2010")):
    
    ## pr must be a DataArray with a time dimension labelled as a DateTime object

    ## NB CURRENTLY NOTHING IN PLACE TO HANDLE DRY PERIODS IN THE TIME SERIES - USE WITH CAUTION!
    
    # create empty list to hold results
    spi_monthly = []
    
    # fit each calendar month in turn
    for m in months:
        
        # select current calendar month
        pr_m = pr.sel(time = pr.time.dt.month == m)
        
        # estimate parameters over calibration period (PWM seems to give more robust fit, so using xclim fit method - also applies over whole spatial array)
        pr_cal = pr_m.sel(time = calibration_period).copy()
        pars = fit(pr_cal.dropna("time", "all"), dist = "gamma", method = "PWM")
        
        # normalise all values using fitted parameters (again, apply over whole spatial array)
        spi_m = xr.apply_ufunc(lambda pr, dparams : norm.ppf(gamma.cdf(pr, *dparams)), pr_m, pars, 
                               input_core_dims = [["time"],["dparams"]], output_core_dims = [["time"]], vectorize = True).assign_coords(time = pr_m.time)
        
        # replace +ve (-ve) infinite values with finite maximum (minimum) in each grid cell
        spi_finite = xr.concat([spi_m.where(~np.isinf(spi_m)), 
                                xr.ones_like(spi_m).where(spi_m == np.inf) * spi_m.where(spi_m < np.inf).max("time"),
                                xr.ones_like(spi_m).where(spi_m == -np.inf) * spi_m.where(spi_m > -np.inf).min("time")], "new").sum("new").copy()
        
        # append to list of monthly SPI for concatenation
        spi_monthly.append(spi_finite)
        
    # concatenate monthly fitted values, reorder & relabel
    spi = xr.concat(spi_monthly, "time").sortby("time").rename("spi")
    
    # clean out existing attributes and replace with new variable description
    for k in list(spi.attrs.keys()): del spi.attrs[k]
    spi = spi.assign_attrs(notes = "Calibrated against "+str(calibration_period.start)+"-"+str(calibration_period.stop))
    spi = spi.where(~np.isnan(pr)).dropna("time", "all")
    
    return spi



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# CDF of generalised logistic function used with PWM parametrisation (replicates that used in R)
def hosking_cdf(x, k, loc, scale):
    
    x_scaled = (x-loc) / scale
    
    if k == 0:
        y = x_scaled
    else:
        # if scaled value is greater than 1, replace with 0.99999 (otherwise log is NA):
        # this occurs when observed value falls above range of calibration distribution
        x_scaled = k*x_scaled
        x_scaled[x_scaled >= 1] = 0.99999
        y = -np.log(1 - x_scaled) / k
        
    Fx = 1/(1+np.exp(-y))
            
    return Fx


def get_spei(eff_pr, months = range(1,13), calibration_period = slice("1980", "2010")):
    
    ## eff_pr must be a DataArray with a time dimension labelled as a DateTime object

    ## NB CURRENTLY NOTHING IN PLACE TO HANDLE DRY PERIODS IN THE TIME SERIES - USE WITH CAUTION!
    
    # create empty list to hold results
    spei_monthly = []
    
    # fit each calendar month in turn
    for m in months:
        
        # select current calendar month
        epr_m = eff_pr.sel(time = eff_pr.time.dt.month == m)
        
        # estimate parameters over calibration period (PWM seems to give more robust fit, so using xclim fit method - also applies over whole spatial array)
        epr_cal = epr_m.sel(time = calibration_period).copy()
        
        # fit parameters for whole map (has to be done manually due to bug in xclim)
        pars = xr.apply_ufunc(lambda x : np.asarray(list(glo.lmom_fit(x.copy()).values())), epr_cal,
                              input_core_dims = [["time"]], output_core_dims = [["dparams"]], vectorize = True).assign_coords(dparams = ["k", "loc", "scale"])
        
        # running without dry-month normalisation for now - check if needed for this dataset
        spei_m = xr.apply_ufunc(lambda pr, dparams : norm.ppf(hosking_cdf(pr, *dparams)), epr_m, pars, 
                                input_core_dims=[["time"],["dparams"]], output_core_dims=[["time"]], vectorize = True).assign_coords(time = epr_m.time)
        
        # replace +ve (-ve) infinite values with finite maximum (minimum) in each grid cell
        spei_finite = xr.concat([spei_m.where(~np.isinf(spei_m)), 
                                 xr.ones_like(spei_m).where(spei_m == np.inf) * spei_m.where(spei_m < np.inf).max("time"),
                                 xr.ones_like(spei_m).where(spei_m == -np.inf) * spei_m.where(spei_m > -np.inf).min("time")], "new").sum("new").copy()
        
        spei_monthly.append(spei_finite)
        
    # concatenate monthly fitted values, reorder & relabel
    spei = xr.concat(spei_monthly, "time").sortby("time").rename("spei")
    
    # clean out existing attributes and replace with new variable description
    for k in list(spei.attrs.keys()): del spei.attrs[k]
    spei = spei.assign_attrs(notes = "Calibrated against "+str(calibration_period.start)+"-"+str(calibration_period.stop))
    spei = spei.where(~np.isnan(eff_pr)).dropna("time", "all")
    
    return spei