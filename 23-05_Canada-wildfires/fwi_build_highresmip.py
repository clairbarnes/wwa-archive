from fwi import *
import xarray as xr
import numpy as np
import time
import sys
import warnings
import glob
import os


#def arden_buck(daT):
#    daP = daT
#    daP[daT > 0] = 611.21 * np.exp((18.678 - daT[daT > 0]/234.5)*(daT[daT > 0]/(257.14 + daT[daT > 0])))
#    daP[daT < 0] = 611.15 * np.exp((23.036 - daT[daT < 0]/333.7)*(daT[daT < 0]/(279.82 + daT[daT < 0])))
#    return daP


#def relative_humidity(huss, ps, T):
#    w  = huss/(1-huss)
#    ws = 0.622 * arden_buck(T) / ps
#    return 100 * w / ws


def main(directory, lat_ind, member):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        da = xr.open_dataset(directory + f'hurs/hurs_3hr_{member}.nc').hurs.load()
        time_axis = da.time
        H = da[:,lat_ind,:].to_numpy()
        
        da = xr.open_dataset(directory + f'tas/tas_3hr_{member}.nc').hurs.load().interp(time = time_axis)
        T = da[:,lat_ind,:].to_numpy()

        da = xr.open_dataset(directory + f'snw/snw_day_{member}.nc').snw.load().interp(time = time_axis)
        S = da[:,lat_ind,:].to_numpy()

        da = xr.open_dataset(directory + f'pr_24h/pr_24hr_{member}.nc').pr.load().interp(time = time_axis)
        da = da.isel(time = da.time.dt.year < 2051)
        P = da[:,lat_ind,:].to_numpy()

        da = xr.open_dataset(directory + f'sfcWind/sfcWind_day_{member}.nc').sfcWind.load().interp(time = time_axis)
        W = da[:,lat_ind,:].to_numpy()

        months = da.time.dt.month.to_numpy()
        days = da.time.dt.day.to_numpy()

        FFMC = np.zeros_like(T)
        DMC = np.zeros_like(T)
        DC = np.zeros_like(T)
        ISI = np.zeros_like(T)
        BUI = np.zeros_like(T)
        FWI = np.zeros_like(T)

        for lon_ind in range(len(da.lon)):

            ffmc, dmc, dc, isi, bui, fwi = calculate_fwi(months, days,
                                                         T[:, lon_ind],
                                                         P[:, lon_ind],
                                                         W[:, lon_ind],
                                                         H[:, lon_ind],
                                                         S[:, lon_ind])
            FFMC[:, lon_ind] = ffmc
            DMC[:, lon_ind] = dmc
            DC[:, lon_ind] = dc
            ISI[:, lon_ind] = isi
            BUI[:, lon_ind] = bui
            FWI[:, lon_ind] = fwi

        output = xr.Dataset(data_vars = {'ffmc': (['time','lat','lon','member'], FFMC[:,np.newaxis,:,np.newaxis]),
                                         'dmc':  (['time','lat','lon','member'],  DMC[:,np.newaxis,:,np.newaxis]),
                                         'dc':   (['time','lat','lon','member'],   DC[:,np.newaxis,:,np.newaxis]),
                                         'isi':  (['time','lat','lon','member'],  ISI[:,np.newaxis,:,np.newaxis]),
                                         'bui':  (['time','lat','lon','member'],  BUI[:,np.newaxis,:,np.newaxis]),
                                         'fwi':  (['time','lat','lon','member'],  FWI[:,np.newaxis,:,np.newaxis])},
                            coords = {'time': da.time.to_numpy(),
                                      'member': np.atleast_1d(member),
                                      'lat': da.lat,
                                      'lon': da.lon})

        output.to_netcdf('/rds/general/user/cb2714/ephemeral/highresmip/fwi_slices/'+
                         f'fwi_output_{member}_lat_{lat_ind}.nc')
    return


if __name__ == '__main__':
    directory = '/rds/general/user/cb2714/projects/wwa/ephemeral/canada_fwi/highresmip/'
    tail = '.nc'

    paths = glob.glob(directory + 'tas/tas_*' + tail)

    lats = xr.open_dataset(paths[0]).lat.load()

    full_length = len(lats) * len(paths)

    J = int(os.getenv('PBS_ARRAY_INDEX'))

    i, j = int(np.floor(J/len(lats))), J%len(lats)                        

    path = paths[i]
    member = paths[i].split("tas_3hr_")[1].split('.n')[0]
    lat_ind = j

    main(directory, lat_ind, member)