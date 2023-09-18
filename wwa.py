
## IMPORT NECESSARY MODULES

import xarray as xr; xr.set_options(keep_attrs = True)
import pandas as pd
import numpy as np

from scipy.stats import norm, gamma, lognorm, gaussian_kde as kde

from functools import reduce

from xclim.core.calendar import convert_calendar
from xclim.core.units import convert_units_to
from xclim.indices._conversion import potential_evapotranspiration

import os; os.environ['PROJ_LIB'] = '/home/clair/miniconda3/envs/wwa/share/proj'                         # fixes error message on import of cartopy etc

import cartopy
import geopandas as gpd
import regionmask
from geopy.geocoders import Nominatim

# needed when converting regionmask into polygon
from shapely.geometry import Polygon
# from xrspatial.experimental import polygonize

import re
import glob
from dateutil.relativedelta import relativedelta
from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Rectangle
matplotlib.rcParams['savefig.bbox'] = "tight"    # always save with tight bounding box
matplotlib.rcParams["savefig.facecolor"] = "w"   # always save with white (rather than transparent) background

import warnings
warnings.filterwarnings("ignore", message = "facecolor will have no effect.+")                           # warning about change to Cartopy plotting defaults
warnings.filterwarnings("ignore", message = "__len__ for multi-part geometries is deprecated.+")         # warning about change to Shapely defaults
warnings.filterwarnings("ignore", message = ".+Results from 'centroid' are likely incorrect.+")            # warning against using centroids without reprojecting

from  IPython.display import clear_output

####################################################################################################################

def load_ts(fnm, names, **kwargs): return(pd.read_csv(fnm, comment = "#", sep = " ", header = None, names = names, **kwargs))

####################################################################################################################
## GOODNESS OF FIT

def qqplot(ts, ax = None, dist = norm, marker = ".", ax_labels = True, **kwargs): 
    
    ts = ts[np.isfinite(ts)] 
    x = np.arange(0,1,1/(len(ts)+1))[1:]
    fitted = dist.fit(ts)
    
    if ax is None:
        fig, ax = plt.subplots(figsize = (5,5), dpi = 100, facecolor = "w")
        
    ax.scatter(dist.ppf(x, *fitted), sorted(ts), marker = marker, **kwargs)
    
    vmin = min([ts.min(), dist.ppf(x, *fitted).min()])
    vmax = max([ts.max(), dist.ppf(x, *fitted).max()])
    ax.plot((vmin, vmax), (vmin, vmax), ls = "--", color = "k")
    
    if ax_labels:
        ax.set_xlabel("Fitted"); ax.set_ylabel("Observed")
        

###############################################################################################################
## METHODS FOR EXTRACTING USEFUL INFORMATION FROM RESULTS FILES

def clean_line(l): return re.sub(" +", " ", re.sub("\\.\\.\\.", "", re.sub("<.+?>", " ", l)))

def read_results(fnm): 
    
    # initialise a couple of empty variables
    mu_prime = []; sigma_prime = []; disp = [None, None, None]; y_start = ""; y_end = ""; fitted_rp = ""
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # identify type of analysis carried out
    if "attr-now" in fnm:
        atype = "present"
    elif "attr-fut" in fnm:
        atype = "future"
    elif "val" in fnm:
        atype = "validation"
    else:
        atype = "?"
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # loop over all lines and extract sections of interest
    for line in open(fnm, "r").read().splitlines():
        
        # options used to run analysis (mainly for checking purposes)
        if "scripturl02" in line: opts = {s.split("=")[0].lower() : s.split("=")[1] for s in line.split("&")}
        if "covariate_description" in line: covariate_file = re.sub(".+:: ", "", line)
        if "covariate_file" in line: covariate_file = re.sub(".+:: ", "", line)
        if ">N:<" in line: N = re.sub("\D+", "", line)
        
        # parameter estimates & ranges
        if "mu;':" in line: mu_prime.append(clean_line(line).split(" ")[2:])
        if "sigma;':" in line: sigma_prime.append(clean_line(line).split(" ")[2:])
        if "sigma;/&mu;:" in line: disp = clean_line(line).split(" ")[2:]
        if "alpha;:" in line: alpha = clean_line(line).split(" ")[2:]
        
        if "return period event" in line: mag = clean_line(line).split(" ")[6]
        
        # estimated properties of distributions
        if "atr1" in line: rp = clean_line(line).split(" ")[-4:-1]
        if "atra" in line: pr = clean_line(line).split(" ")[-4:-1]
        if "change in intensity" in line: DeltaI = clean_line(line).split(" ")[-4:-1]
        
    # convert options into more useful information
    y = opts["year"]
    if "begin" in opts.keys(): y_start = opts["begin"]
    if "end" in opts.keys(): y_end = opts["end"]
    if "biasrt" in opts.keys(): fitted_rp = opts["biasrt"]
        
    covariate_matched = (covariate_file.split(" ")[2] in opts["station"]) and (covariate_file.split(" ")[4] in opts["station"])
    
    mu_prime = [m[1:] for m in mu_prime if m[0] == opts["year"]][0]
    sigma_prime = [s[1:] for s in sigma_prime if s[0] == opts["year"]][0]
    
    lower_tail = "changesign" in opts.keys()
    
    # will probably also need to confirm whether shape parameter was constrained
    # should also add whether data was transformed
            
    # return DataFrame of key information (some for confirmation, some to input into sheet)
    return pd.DataFrame({"dataset" : opts["station"],
                         "covariate" : covariate_file,
                         "covariate_matched" : covariate_matched,
                         "attr_type" : atype,
                         "distribution" : opts["fit"], 
                         "fit_type" : opts["assume"],
                         "lower_tail" : lower_tail,
                         "include_event" : opts["includelast"] == "on",
                         "fitted_years" : y_start+"-"+y_end,
                         "N" : N,
                         "event_year" : opts["year"],
                         "vs_gmst" : opts["cov1"],
                         "return_time" : fitted_rp,
                         "sigma_est" : sigma_prime[0],
                         "sigma_lower" : sigma_prime[1],
                         "sigma_upper" : sigma_prime[2],
                         "disp_est" : disp[0],
                         "disp_lower" : disp[1],
                         "disp_upper" : disp[2],
                         "alpha_est" : alpha[0],
                         "alpha_lower" : alpha[1],
                         "alpha_upper" : alpha[2],
                         "event_magnitude" : mag,
                         "rp_est" : rp[0],
                         "rp_lower" : rp[1],
                         "rp_upper" : rp[2],
                         "gmst_now" : str(-float(opts["cov1"])),
                         "pr_est" : pr[0],
                         "pr_lower" : pr[1],
                         "pr_upper" : pr[2],
                         "DI_est" : DeltaI[0],
                         "DI_lower" : DeltaI[1],
                         "DI_upper" : DeltaI[2],
                        },
                        index = [re.sub("\\..+", "", re.sub(".+/", "", fnm))])


###############################################################################################################
## PLOTTING

def sc_xlabels(dates, ax = None):
    
    # method to add labels to seasonal cycle plot
    
    if ax is None: ax = plt.gca()
    
    labelticks = [i for i in range(366) if dates.dt.day[i] == 1]
    labels = [dates[i].dt.strftime("%b").values for i in range(366) if dates.dt.day[i] == 1]

    ax.set_xticks(labelticks)
    ax.set_xticklabels(labels)

    
# method to get DOY offset for years starting other than in January
def y_offset(months): return datetime(2020,months,1).timetuple().tm_yday
    
    
# quickly plot a fitted polynomial
def xyline(x, y, ax = None, npoly = 1, **plot_kwargs):
    
    if not ax: fig, ax = plt.subplots(figsize = (5,3))
    ax.plot(np.sort(x), np.poly1d(np.polyfit(x, y, npoly))(np.sort(x)), **plot_kwargs)
    
###############################################################################################################
## DROUGHT PLOTTING FUNCTIONS

def drought_map(di, ax = None, **kwargs):
    
    dc = xr.apply_ufunc(np.digitize, di, kwargs={'bins': [-np.inf, -2, -1.55, -1.25, -.75, -.5]})
    
    if ax is None:
        fig, ax = plt.subplots(ncols = 1, dpi = 100)
    
    drought_cmap = matplotlib.colors.ListedColormap(['darkred', 'red', 'orange', 'gold','yellow']); drought_cmap.set_over('honeydew')

    cbar = dc.plot(ax = ax, cmap = drought_cmap, norm = matplotlib.colors.BoundaryNorm(np.arange(0.5,6.5,1), drought_cmap.N), add_colorbar = False, **kwargs)
    return cbar  


def drought_colorbar(cbar, ax, location = "bottom", label = "Drought classification", increasing = True, **cbar_kwargs): 
    
    cbar = plt.colorbar(cbar, ax = ax, location = location, ticks = list(range(1,6)), extend = "max", label = label, **cbar_kwargs)
    if location == "bottom":
        cbar.ax.set_xticklabels(["D" + str(x) for x in range(4,-1,-1)])
        if increasing: cbar.ax.invert_xaxis()
    else:
        cbar.ax.set_yticklabels(["D" + str(x) for x in range(4,-1,-1)])
        if increasing: cbar.ax.invert_yaxis()

###############################################################################################################
## MISC

def wrap_lon(ds):
    
    # method to wrap longitude from (0,360) to (-180,180)
    
    if "longitude" in ds.coords:
        lon = "longitude"
        lat = "latitude"
    elif "lon" in ds.coords:
        lon = "lon"
        lat = "lat"
    else: 
        # can only wrap longitude
        return ds
    
    if ds[lon].max() > 180:
        ds[lon] = (ds[lon].dims, (((ds[lon].values + 180) % 360) - 180), ds[lon].attrs)
        
    if lon in ds.dims:
        ds = ds.reindex({ lon : np.sort(ds[lon]) })
        ds = ds.reindex({ lat : np.sort(ds[lat]) })
    return ds



def decode_times(ts):
    
    # Method to manually decode times
    
    inc = re.sub(" .+", "", ts.time.units)
    startdate = pd.Timestamp(re.sub(".+since ", "", ts.time.units)+' 00:00:00.000000').to_pydatetime()
    
    if inc == "years":
        new_times = [np.datetime64(startdate + relativedelta(years = i)) for i in range(len(ts.time))]
    elif inc == "months":
        new_times = [np.datetime64(startdate + relativedelta(months = i)) for i in range(len(ts.time))]
    else:
        print("TBD: " +inc)
        return
        
    ts = ts.assign_coords(time = new_times)
    
    return ts



def get_latlon(city):
    
    # retrieve lat & lon for given location
    location = Nominatim(user_agent="GetLoc").geocode(city)
    if location is None:
        return {"lon" : None, "lat" : None}
    else:
        return {"lon" : location.longitude, "lat" : location.latitude}


def normalised_seasonal_cycle(ts):
    
    ts = convert_calendar(ts, "default", align_on = "date")
    sc = ts.groupby("time.dayofyear").mean()
    return sc / sc.mean()


def eval_df(ens, region = None):
    # create an empty DataFrame to store evaluation results
    if fnm is None: 
        fnm = ens+"_model-eval.txt"
    else:
        fnm = ens+"_"+region+"_model-eval.txt"
    pd.DataFrame({"seasonal_cycle" : "?", "spatial_pattern" : "?"}, index = [cordex_model(fnm) for fnm in glob.glob("cordex/pr-spatial_"+ens+"_*")]).to_csv(fnm)
    
    
    
def nearest_px(x,y,da, xcoord = "longitude", ycoord = "latitude", return_map = False):
    
#     if xcoord is None:
#         if "lon" in da.dims:
#             xcoord = "lon"
#         elif "longitude" in da.dims:
#             xcoord = "longitude"
#         else: 
#             print("No x-coords identified")
#             return None
        
#     if ycoord is None:
#         if "lat" in da.dims:
#             xcoord = "lat"
#         elif "latitude" in da.dims:
#             xcoord = "latitude"
#         else: 
#             print("No y-coords identified")
#             return None
   
    # get squared distance from (x,y) to each point
    dist2 = (da[ycoord] - y)**2 + (da[xcoord] - x)**2
   
    # exclude any cells where the gridded data is NaN
    dist2 = dist2.where(~np.isnan(da))
   
    # also limit distance to closest two squares (in case there really is no data nearby)
    dist2 = dist2.where(dist2 <= 5.76e8)
    
    if return_map:
        closest_px = xr.ones_like(da).where(dist2 == dist2.min())
        return closest_px
    else:
        # return time series
        # find value in cell containing minimum distance
        # if multiple equidistant cells, will average over them
        val = da.where(dist2 == dist2.min()).mean([xcoord, ycoord])
        return val



###############################################################################################################
# def cx_csv(da, fnm = None, dataset = None):
    
#     # write CSV for easy import into Climate Explorer
    
#     rnm = da.run.values[0]
#     da = da.squeeze(drop = True)
#     fnm_string = da.name+"_"+re.sub(" ", "_", rnm)
    
#     if dataset is not None:
#         fnm_string = dataset+"_"+fnm_string
        
#     if fnm is None:
#         fnm = "ts/"+fnm_string
    
#     if "time" in da.dims:
#         da = da.assign_coords(time = da.time.dt.year).rename(time = "#time")
#     elif "year" in da.dims:
#         da = da.rename(year = "#time")
#     else:
#         print(da.dims)
#         return

#     # write to csv
#     fnm = re.sub(".txt", "", fnm)+".txt"
#     da.to_dataframe().to_csv(fnm, sep = " ")
    
#     # add a text string specifying the units (don't think format is correct here)
#     if "units" in da.attrs:
#         unit_string = "# "+da.name+" ["+da.units+"]"
#         unit_string = "# variable ["+da.units+"]"
#         ! echo "$unit_string" >> $fnm
    
#     # add a line specifying the model & variable name, to be used as filename when uploading
#     fnm_string = "# "+fnm_string
#     ! echo "$fnm_string" >> $fnm

###############################################################################################################
## LISTS OF CORDEX MODELS

gcm = {'CCCma-CanESM2' : "CanESM2",
       'CNRM-CERFACS-CNRM-CM5' : "CNRM-CM5",
       'CSIRO-QCCCE-CSIRO-Mk3-6-0' : "CSIRO-Mk3-6-0",
       'CSIRO-BOM-ACCESS1-0' : "ACCESS1-0",
       'CSIRO-BOM-ACCESS1-3' : "ACCESS1-3",
       'ICHEC-EC-EARTH' : "EC-EARTH",
       'IPSL-IPSL-CM5A-LR' : "IPSL-CM5A-LR",
       'IPSL-IPSL-CM5A-MR' : "IPSL-CM5A-MR",
       'MIROC-MIROC5' : "MIROC5",
       'MOHC-HadGEM2-ES' : "HadGEM2-ES",
       'MPI-M-MPI-ESM-LR' : "MPI-ESM-LR",
       'MPI-M-MPI-ESM-MR' : "MPI-ESM-MR",
       'NCC-NorESM1-M' : "NorESM1-M",
       'NOAA-GFDL-GFDL-ESM2M' : "GFDL-ESM2M"}

rcm = {'CCCma-CanRCM4' : "CanRCM4",
       'CLMcom-CCLM4-8-17' : "CCLM4-8-17",
       'CLMcom-CCLM4-8-17-CLM3-5' : "CCLM4-8-17",
       'CLMcom-ETH-COSMO-crCLIM-v1-1' : 'crCLIM-v1-1',
       'CLMcom-HZG-CCLM5-0-15' : "CCLM5-0-15",
       'CLMcom-KIT-CCLM5-0-15' : "CCLM5-0-15",
       'CNRM-ALADIN53' : "ALADIN53",
       'CNRM-ALADIN63' : "ALADIN63",
       'DMI-HIRHAM5' : "HIRHAM5",
       'GERICS-REMO2009' : "REMO2009",
       'GERICS-REMO2015' : "REMO2015",
       'ICTP-RegCM4-3' : "RegCM4-3",
       'ICTP-RegCM4-6' : "RegCM4-6",
       'ICTP-RegCM4-7' : "RegCM4-7",
       'IITM-RegCM4-4' : "RegCM4-4",
       'IPSL-WRF381P' : "WRF381P",
       'KNMI-RACMO22E' : "RACMO22E",
       'KNMI-RACMO22T' : "RACMO22T",
       'MOHC-HadREM3-GA7-05' : "HadREM3-GA7-05",
       'MPI-CSC-REMO2009' : "REMO2009",
       'RMIB-UGent-ALARO-0' : 'ALARO-0',
       'SMHI-RCA4' : "RCA4",
       'UCAN-WRF341I' : "WRF341I",
       'UHOH-WRF361H' : "WRF361H",
       'UNSW-WRF360J' : "WRF360J",
       'UNSW-WRF360K' : "WRF360K"}


# method to decode filename into model name
def cordex_model(fnm): return gcm[fnm.split("_")[2]]+"_"+fnm.split("_")[4][:-4]+"_"+rcm[fnm.split("_")[5]]


###############################################################################################################

def reshape_df(fnm, da):
    
    # method to load dataframe of fitted results & reshape into DataArray for plotting
    df = pd.read_csv(fnm, index_col = 0)
    fitted = xr.Dataset(data_vars = {vnm : xr.DataArray(np.array(df.loc[vnm]).reshape(*da.shape), dims = ["lat", "lon"]) for vnm in df.index},
                        coords = {"lat" : da.lat, "lon" : da.lon})
    
    return fitted



def merge_byindex(df_list): 
    
    # merge a list of dataframes by matching indices
    return reduce(lambda left, right: pd.merge(left, right, left_index = True, right_index = True, how = 'outer'), df_list)