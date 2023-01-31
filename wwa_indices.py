
## MODULE TO COMPUTE MROE COMPLEX CLIMATE INDICES - TO BE USED ALONGSIDE WWA.PY

import sys; sys.path.append('/home/clair/wwa'); from wwa import *

from xclim.indices.stats import fit
from scipy.stats import norm, gamma

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