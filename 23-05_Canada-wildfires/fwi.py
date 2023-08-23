import math
import numpy as np
from itertools import groupby


class FWICLASS:
    '''
    Developed from:
    
    "Updated source code for calculating fire danger indices
     in the canadian forest fire weather index system"
    Y. Wang, K.R. Anderson, and R.M. Suddaby
    Information Report NOR-X-424
    
    Code fixed to include limit on humidity at 100, and case for 
    saturation of the drought code to 0

    '''
    def __init__(self, temp, rhum, wind, prcp):
        self.h = rhum
        self.t = temp
        self.w = wind
        self.p = prcp
    
    def FFMCcalc(self, ffmc0):
        mo = (147.2 * (101.0 - ffmc0)) / (59.5 + ffmc0)
        h = self.h
        if h > 100:
            h = 100

        if self.p > 0.5:
            rf = self.p - 0.5
            
            if mo > 150.0:
                mo = (mo + 42.5 * rf * math.exp(-100.0 / (251.0 - mo)) * (1.0 - math.exp(-6.93 / rf))) + (0.0015 * (mo - 150.0) ** 2) * math.sqrt(rf)
            elif mo <= 150.0:
                mo = mo + 42.5 * rf * math.exp(-100.0 / (251.0 - mo)) * (1.0 - math.exp(-6.93 / rf))
            if mo > 250.0:
                mo = 250.0
        
        ed = 0.942 * (h ** 0.679) + (11.0 * math.exp((h - 100.0) / 10.0)) + 0.18 * (21.1 - self.t) * (1.0 - 1.0 / math.exp(0.1150 * h))
        
        if mo < ed:
            ew = 0.618 * (h ** 0.753) + (10.0 * math.exp((h - 100.0) / 10.0)) + 0.18 * (21.1 - self.t) * (1.0 - 1.0 / math.exp(0.115 * h))
            
            if mo <= ew:
                kl = 0.424 * (1.0 - ((100.0 - h) / 100.0) ** 1.7) + (0.0694 * math.sqrt(self.w)) * (1.0 - ((100.0 - h) / 100.0) ** 8)
                kw = kl * (0.581 * math.exp(0.0365 * self.t))
                m = ew - (ew - mo) / 10.0 ** kw
            elif mo > ew:
                m = mo
        elif mo == ed:
            m = mo
        elif mo > ed:
            kl = 0.424 * (1.0 - (h / 100.0) ** 1.7) + (0.0694 * math.sqrt(self.w)) * (1.0 - (h / 100.0) ** 8)
            kw = kl * (0.581 * math.exp(0.0365 * self.t))
            m = ed + (mo - ed) / 10.0 ** kw
        
        ffmc = (59.5 * (250.0 - m)) / (147.2 + m)
        
        if ffmc > 101.0:
            ffmc = 101.0
        elif ffmc <= 0.0:
            ffmc = 0.0
        
        return ffmc
    
    
    def DMCcalc(self, dmc0, mth):
        el = [6.5, 7.5, 9.0, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8.0, 7.0, 6.0]
        t = self.t
        h = self.h
        
        if h > 100:
            h = 100
        
        if t < -1.1:
            t = -1.1
        
        rk = 1.894 * (t + 1.1) * (100.0 - h) * (el[mth - 1] * 0.0001)
        
        if self.p > 1.5:
            ra = self.p
            rw = 0.92 * ra - 1.27
            wmi = 20.0 + 280.0 / math.exp(0.023 * dmc0)
            
            if dmc0 <= 33.0:
                b = 100.0 / (0.5 + 0.3 * dmc0)
            elif dmc0 > 33.0:
                if dmc0 <= 65.0:
                    b = 14.0 - 1.3 * math.log(dmc0)
                elif dmc0 > 65.0:
                    b = 6.2 * math.log(dmc0) - 17.2
            
            wmr = wmi + (1000 * rw) / (48.77 + b * rw)
            pr = 43.43 * (5.6348 - math.log(wmr - 20.0))
        elif self.p <= 1.5:
            pr = dmc0
        
        if pr < 0.0:
            pr = 0.0
        
        dmc = pr + rk
        
        if dmc <= 1.0:
            dmc = 1.0
        
        return dmc
    
    def DCcalc(self, dc0, mth):
        fl = [-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6]
        t = self.t
        
        if t < -2.8:
            t = -2.8
        
        pe = (0.36 * (t + 2.8) + fl[mth-1])/2
        
        if pe <= 0.0:
            pe = 0.0
        
        if self.p > 2.8:
            ra = self.p
            rw = 0.83 * ra - 1.27
            smi = 800.0 * math.exp(-dc0 / 400.0)
            dr = dc0 - 400.0 * math.log(1.0 + ((3.937 * rw) / smi))
            
            if dr > 0.0:
                dc = dr + pe
            else:
                dc = 0 
        
        elif self.p <= 2.8:
            dc = dc0 + pe
            
        if dc > 1000:
            dc = 1000
        
        return dc
    
    def ISIcalc(self, ffmc):
        mo = 147.2 * (101.0 - ffmc) / (59.5 + ffmc)
        ff = 19.115 * math.exp(mo * -0.1386) * (1.0 + (mo ** 5.31) / 49300000.0)
        isi = ff * math.exp(0.05039 * self.w)
        return isi
    
    def BUIcalc(self, dmc, dc):
        if dmc <= 0.4 * dc:
            bui = (0.8 * dc * dmc) / (dmc + 0.4 * dc)
        else:
            bui = dmc - (1.0 - 0.8 * dc / (dmc + 0.4 * dc)) * (0.92 + (0.0114 * dmc) ** 1.7)
        
        if bui < 0.0:
            bui = 0.0
        
        return bui
    
    def FWIcalc(self, isi, bui):
        if bui <= 80.0:
            bb = 0.1 * isi * (0.626 * bui ** 0.809 + 2.0)
        else:
            bb = 0.1 * isi * (1000.0 / (25. + 108.64 / math.exp(0.023 * bui)))
        
        if bb <= 1.0:
            fwi = bb
        else:
            fwi = math.exp(2.72 * (0.434 * math.log(bb)) ** 0.647)
        
        return fwi
    


def was_winter_snowy(i, s, months):
    # Checks if snow < 0.1mm in 75% of days
    if months[i] <= 2:
        test = s[i-90:i+365-90][np.logical_or(months[i-90:i+365-90] == 1,
                                              months[i-90:i+365-90] == 2)] > 0.1
    else:
        test = s[i-365:i][np.logical_or(months[i-365:i] == 1,
                                        months[i-365:i] == 2)] > 0.1
    
    test = np.sum(test) / len(test) > 0.75
    return test


#def will_fall_be_snowy(i, s, months):
#    # Testing for if there will be 7 consecutive snowy days in fall
#    if months[i] < 9:
#        test = s[i:i+365][np.logical_and(months[i:i+365] <= 11, months[i:i+365] >= 11)] > 0.1
#    else:
#        test = s[i-100:i+100][np.logical_and(months[i-100:i+100] <= 11, months[i-100:i+100] >= 11)] > 0.1
#    
#    test = np.max([len(list(g)) for k, g in groupby(test) if k == True]) >= 7
#    return test    


def does_season_start(i, s, t, snowy_winter_test):
    if snowy_winter_test:
        test = (np.mean(s[i-3:i] < 0.001) == 1)
    else:
        test = (np.mean(t[i-3:i] > 12) == 1)
    return test


def does_season_end(i, s, t, snowy_fall_test):
#    if snowy_fall_test:
#        test = (np.mean(s[i:i+7] > 0.0001) == 1)
#    else:
    test = (np.mean(t[i-3:i] < 5) == 1)
    
    return test
    
    
def calculate_fwi(months, days, t, p, w, h, s):
    
    ffmc = np.zeros_like(t) * np.nan
    dmc  = np.zeros_like(t) * np.nan
    dc   = np.zeros_like(t) * np.nan
    isi  = np.zeros_like(t) * np.nan
    bui  = np.zeros_like(t) * np.nan
    fwi  = np.zeros_like(t) * np.nan

    fire_season = False
    first_day = True
    p_ow = 0 # Precip over winter
    fall_dc = 15 # For first year as fall DC unknown

    for i in range(len(days)):    
        snowy_winter_test = was_winter_snowy(i, s, months)
        
        if fire_season == False:
            if first_day == True:
                first_day = False

            p_ow += p[i]

            start_test = does_season_start(i, s, t, snowy_winter_test)

            if start_test == True:
                fire_season = True
                first_day = True

        if fire_season == True:
            
            if first_day == True:
                snowy_fall_test = True#will_fall_be_snowy(i, s, months)
                first_day = False
                
                ffmc0 = 85.0
                dmc0 = 6.0
                # Overwintering of drought code:
                if p_ow <= 200:
                    # Overwintering formula from Hanes et al (2020) - coefs from Yan Boulanger
                    moisture_equivalent = 1 * (800.0 * math.exp(-fall_dc / 400.0)) + 0.75 * (3.937 * p_ow)
                    dc0 = 400 * math.log(800 / moisture_equivalent)
                    if dc0 < 0:
                        dc0 = 0
                    if dc0 > 1000:
                        dc0 = 1000
                else:
                    dc0 = 15

            # Calculating FWI values:
            fwisystem = FWICLASS(t[i], h[i], w[i], p[i])
            ffmc[i]   = fwisystem.FFMCcalc(ffmc0)
            dmc[i]    = fwisystem.DMCcalc(dmc0, months[i])
            dc[i]     = fwisystem.DCcalc(dc0, months[i])
            isi[i]    = fwisystem.ISIcalc(ffmc[i])
            bui[i]    = fwisystem.BUIcalc(dmc[i], dc[i])
            fwi[i]    = fwisystem.FWIcalc(isi[i], bui[i])

            ffmc0 = ffmc[i]
            dmc0 = dmc[i]
            dc0 = dc[i]

            end_test = does_season_end(i, s, t, snowy_fall_test)
            
            # Adding in a hard limit on the fire season end (beginning of winter)
            #if months[i] == 12:
            #    if days[i] == 1:
            #        end_test = True
                    
            if end_test == True:
                fall_dc = dc0
                p_ow = 0
                first_day = True
                fire_season = False

    return ffmc, dmc, dc, isi, bui, fwi
    
    
def calculate_fwi_no_overwintering(months, days, T, P, W, H):
    
    if (((True in np.isnan(T)) or 
         (True in np.isnan(P))) or 
        ((True in np.isnan(W)) or 
         (True in np.isnan(H)))):
        ffmc = np.nan * np.zeros_like(T)
        dmc = np.nan * np.zeros_like(T)
        dc = np.nan * np.zeros_like(T)
        isi = np.nan * np.zeros_like(T)
        bui = np.nan * np.zeros_like(T)
        fwi = np.nan * np.zeros_like(T)
    else:
    
        ffmc0 = 85.0
        dmc0 = 6.0
        dc0 = 15.0

        ffmc = np.zeros_like(T)
        dmc = np.zeros_like(T)
        dc = np.zeros_like(T)
        isi = np.zeros_like(T)
        bui = np.zeros_like(T)
        fwi = np.zeros_like(T)

        for i in range(len(days)):
            fwisystem = FWICLASS(T[i], H[i], W[i], P[i])
            ffmc[i]   = fwisystem.FFMCcalc(ffmc0)
            dmc[i]    = fwisystem.DMCcalc(dmc0, months[i])
            dc[i]     = fwisystem.DCcalc(dc0, months[i])
            isi[i]    = fwisystem.ISIcalc(ffmc[i])
            bui[i]    = fwisystem.BUIcalc(dmc[i], dc[i])
            fwi[i]    = fwisystem.FWIcalc(isi[i], bui[i])

            ffmc0 = ffmc[i]
            dmc0 = dmc[i]
            dc0 = dc[i]
        
    return ffmc, dmc, dc, isi, bui, fwi