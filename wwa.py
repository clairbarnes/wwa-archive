import xarray as xr; xr.set_options(keep_attrs = True)
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib; matplotlib.rcParams['savefig.bbox'] = "tight"    # always save with tight bounding box

###############################################################################################################
## METHODS FOR EXTRACTING USEFUL INFORMATION FROM RESULTS FILES

def clean_line(l): return re.sub(" +", " ", re.sub("\\.\\.\\.", "", re.sub("<.+?>", " ", l)))

def read_results(fnm): 
    
    # initialise a couple of empty variables
    mu_prime = []; sigma_prime = []; disp = [None, None, None]; y_start = ""; y_end = ""
    
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
                         "return_time" : opts["biasrt"],
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
## MISC

def decode_times(ts):
    
    # Method to manually decode times
    
    inc = re.sub(" .+", "", ts.time.units)
    startdate = pd.Timestamp(re.sub(".+since ", "", ts.time.units)+' 00:00:00.000000').to_pydatetime()
    
    if inc == "years":
        new_times = [np.datetime64(startdate + relativedelta(years = i)) for i in range(len(ts.time))]
    else:
        print("TBD: " +inc)
        return
        
    ts = ts.assign_coords(time = new_times)
    
    return ts

###############################################################################################################