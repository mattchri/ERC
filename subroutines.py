#-----------------------------------------------------------------
#Module Containing all Subroutines
#required to run the SMIP code
#The first group of subroutines deal with extracting
#times and converting them to julian day or 
#collocating files together.
#The second group deals with collocation.
#
# 02/11/17, MC: upload initial version of the code to the repo
# 15/01/18, MC: updated code to accept Accumulation and Instantaneous types
#-----------------------------------------------------------------
from jdcal import *
import sys
import numpy as np
from netCDF4 import Dataset
import netCDF4
from datetime import datetime
import pdb
import numpy.ma as ma
#import iris
#import datetime
#iris.FUTURE.cell_datetime_objects = True

def doy_to_month(doy,year):
    #fetch julian day
    jday0 = gcal2jd(year,1,1)
    jday1 = (jday0[0],jday0[1]+doy-1)
    caldat = jd2gcal(*jday1)
    month = caldat[1]
    return month

def convert_time(year,month,dom):
    #fetch julian day
    jday = gcal2jd(year,month,dom)
    return jday

def jday_to_doy(jday):
    caldat = jd2gcal(*jday)
    year = caldat[0]
    jday0 = gcal2jd(year,1,1)
    doy = jday[1]-jday0[1]+1
    return doy

def extract_time_seviri(fname):
    year  = float(fname[38:42])
    month = float(fname[42:44])
    dom   = float(fname[44:46])
    hour  = float(fname[46:48])
    minute= float(fname[48:50])
    fday = (hour+minute/60.)/24.
    jday = gcal2jd(year,month,dom)
    jdayOUT = (jday[0],jday[1] + fday)     
    OUT = { "filename" : fname, "MONTH" : month, "YEAR" : year, "DOM" : dom, "HOUR" : hour, "MINUTE" : minute, "JDAY" : jdayOUT}
    return OUT

def extract_time_imerg(fname):
    year = float(fname[23:27])
    month = float(fname[27:29])
    dom   = float(fname[29:31])
    hour  = float(fname[33:35])
    minute = float(fname[35:37]) + 15.
    fday = (hour+minute/60.)/24.
    jday = gcal2jd(year,month,dom)
    jdayOUT = (jday[0],jday[1] + fday)
    OUT = { "filename" : fname, "MONTH" : month, "YEAR" : year, "DOM" : dom, "HOUR" : hour, "MINUTE" : minute, "JDAY" : jdayOUT}
    return OUT

def extract_time_ecmwf(fname):
    year = float(fname[4:8])
    month = float(fname[8:10])
    dom   = float(fname[10:12])
    step  = float(fname[12:14])
    hour  = float(fname[14:16])+step
    minute = step
    fday = (hour+minute/60.)/24.
    jday = gcal2jd(year,month,dom)
    jdayOUT = (jday[0],jday[1] + fday)
    OUT = { "filename" : fname, "MONTH" : month, "YEAR" : year, "DOM" : dom, "HOUR" : hour, "MINUTE" : minute, "JDAY" : jdayOUT}
    return OUT


#-----------------------------------------------------------------
#Read Accumulated (gafs) ECWMF variable from /badc/ on JASMIN
#Inputs
#variable name (e.g. 'TP' total precipitation)
#year [int]
#month [int]
#dom day of the month [int]
#Output
#ARR [accumulated field]
#
#Two types of files
#Accumulated: indxType=0, Instantaneous: indxType=1
#-----------------------------------------------------------------
def read_ecmwf_variable(vname,year,month,dom,ecmwffile):

    #print(vname)

    #Fetch ECWMF file based on variable name
    varNames = np.array(['TTR','STRD','TISR','TSR','SSRD','SSR','STR','TP','SSHF','p82.162','p83.162','p85.162'])
    varTypeIndx  = np.array([  0  ,  0   ,  0   ,  0  ,  0   ,  0   ,  0  , 0  ,  0   ,   1     ,   1     ,    1])
    indxType = (varTypeIndx[np.where( vname == varNames )])[0]
    inputfile = ecmwffile[indxType]

    #print(ecmwffile)
    #print('index = ',indxType)
    #print('inputfile = ',inputfile)

    #Read Downloaded ECWMF Data
    if(len(inputfile) > 1):
        ncfile = Dataset( inputfile, mode='r')
        time = ncfile.variables['time']
        dates = netCDF4.num2date(time[:], time.units, time.calendar)
        tmpTime = time[:]
        t_months = np.array([x.month for x in dates])
        t_days = np.array([x.day for x in dates])
        t_hours = np.array([x.hour for x in dates])
        t_index = (np.where( t_months == month))[0]
        
        #Now sort the index
        tmpsortID  = np.asarray(np.argsort(time[t_index]))
        sort_index = t_index[tmpsortID[:]]

        scaleF  = ncfile.variables[vname.lower()].scale_factor
        add_off = ncfile.variables[vname.lower()].add_offset
        fillVal = ncfile.variables[vname.lower()]._FillValue
        unitStr = ncfile.variables[vname.lower()].units
        ARR = (ncfile.variables[vname.lower()])[sort_index,:,:]


    # Read ECWMF Data from Jasmin
    if(len(inputfile) == 0):
        path_ecmwf = '/badc/ecmwf-era-interim/data/ga/fs/'+str(year).zfill(4)+'/'+str(month).zfill(2)+'/'+str(dom).zfill(2)+'/' 
        prefix = str(year).zfill(4)+str(month).zfill(2)+str(dom).zfill(2)
        ftimes = ['0003','0006','0009','0012','1203','1206','1209','1212']
        ARR = np.empty( [8,256,512] )
        for i in range(len(ftimes)):
            ecmwf_file = path_ecmwf+'gafs'+prefix+ftimes[i]+'.nc'
            ncfile = Dataset( ecmwf_file, mode='r')
            unitStr = ncfile.variables[vname].units
            fillVal=-999.
            ARR[i,:,:] = (ncfile.variables[vname][:])[0,0,:,:]


    #if(ma.is_masked(((ARR[:,0,0]))[0]) == False):
    if(indxType == 1):
        #print('Instantaneous: '+vname)
        #pdb.set_trace()
        tmpARR = ARR[:,:,:]
        aARR = instantaneous_ecmwf_3hr(tmpARR)

    #if(ma.is_masked(((ARR[:,0,0]))[0]) == True):
    if(indxType == 0):
        #print('Accumulated: '+vname)
        tmpARR = ARR[:,:,:]
        aARR = accumulate_ecmwf_3hr(tmpARR)


    #Daily means
    bARR = np.mean(aARR, axis=0)
    

    lat = ncfile.variables['latitude'][:]
    lon = ncfile.variables['longitude'][:]
    OUT = {'lat':lat, 'lon':lon, 'data':aARR, 'unit':unitStr+'/s', 'daily':bARR}

    return OUT

def accumulate_ecmwf_3hr(ARR):
    #Express quantity as instantaneous
    sz=ARR.shape
    tstep    = 3. * 3600. #3 hours converted to seconds
    aARR     = np.empty( [9,sz[1],sz[2]] )
    aARR[0,:,:] = ( ARR[0,:,:]-0.        ) / tstep
    aARR[1,:,:] = ( ARR[1,:,:]-ARR[0,:,:]) / tstep
    aARR[2,:,:] = ( ARR[2,:,:]-ARR[1,:,:]) / tstep
    aARR[3,:,:] = ( ARR[3,:,:]-ARR[2,:,:]) / tstep
    aARR[4,:,:] = ( ARR[4,:,:]-0.        ) / tstep
    aARR[5,:,:] = ( ARR[5,:,:]-ARR[4,:,:]) / tstep
    aARR[6,:,:] = ( ARR[6,:,:]-ARR[5,:,:]) / tstep
    aARR[7,:,:] = ( ARR[7,:,:]-ARR[6,:,:]) / tstep
    aARR[8,:,:] = (ARR[3,:,:]+ARR[7,:,:]) / (tstep*8.)
    return aARR    

def instantaneous_ecmwf_3hr(ARR):
    #Express quantity as instantaneous
    sz=ARR.shape
    INTERP_times = 1.5+np.array([0,3,6,9,12,15,18,21])
    ARR_times = np.array([3,9,15,21])
    aARR     = np.empty( [9,sz[1],sz[2]] )
    aARR[0,:,:] = ARR[0,:,:] #np.interp(INTERP_times[0],ARR_times,ARR)
    aARR[1,:,:] = ARR[0,:,:]
    aARR[2,:,:] = ARR[1,:,:]
    aARR[3,:,:] = ARR[1,:,:]
    aARR[4,:,:] = ARR[2,:,:] 
    aARR[5,:,:] = ARR[2,:,:] 
    aARR[6,:,:] = ARR[3,:,:] 
    aARR[7,:,:] = ARR[3,:,:] 
    aARR[8,:,:] = np.mean(ARR,axis=0)
    return aARR    


def process_ecmwf_variables(year,month,dom,inputfile):
    #Read Data
    STRD = read_ecmwf_variable('STRD',year,month,dom,inputfile)
    TISR = read_ecmwf_variable('TISR',year,month,dom,inputfile)
    TSR  = read_ecmwf_variable('TSR',year,month,dom ,inputfile)
    TTR  = read_ecmwf_variable('TTR',year,month,dom ,inputfile)
    SSRD = read_ecmwf_variable('SSRD',year,month,dom,inputfile)
    SSR  = read_ecmwf_variable('SSR',year,month,dom ,inputfile)
    STRD = read_ecmwf_variable('STRD',year,month,dom,inputfile)
    STR  = read_ecmwf_variable('STR',year,month,dom ,inputfile)
    TP   = read_ecmwf_variable('TP',year,month,dom  ,inputfile)
    SSHF = read_ecmwf_variable('SSHF',year,month,dom,inputfile)
    sz=SSHF['data'].shape
    DIV_DRY_STATIC_ENERGY_data = np.zeros( [sz[0],sz[1],sz[2]] )
    DIV_DRY_STATIC_ENERGY_daily = np.zeros( [sz[1],sz[2]] )

    if(len(inputfile) > 1):
        #NOTE THESE QUANTITIES ARE NOT ACCUMULATED!!!
        p82 = read_ecmwf_variable('p82.162',year,month,dom,inputfile)
        p83 = read_ecmwf_variable('p83.162',year,month,dom,inputfile)
        p85 = read_ecmwf_variable('p85.162',year,month,dom,inputfile)
        #pdb.set_trace()
        DIV_DRY_STATIC_ENERGY_data  = (p82['data']+p83['data']+p85['data']   )
        DIV_DRY_STATIC_ENERGY_daily = (p82['daily']+p83['daily']+p85['daily'])
        #pdb.set_trace()

    # Calculate Fluxes
    L = 2.5e6 #latent heat of vaporization at 0C (J/kg)
    boa_lwdn = STRD['data']
    boa_lwup = STRD['data'] - STR['data']
    boa_swdn = SSRD['data']
    boa_swup = SSRD['data'] - SSR['data']
    toa_lwup = 0. - TTR['data']
    toa_swdn = TISR['data']
    toa_swup = TISR['data'] - TSR['data']
    prate    = TP['data'] * 1000. * 3600. #(m/s * 3600. * 1000mm/m --> mm/hr)
    SH     = 0. - SSHF['data']  #not sure if this is right!!!!!
    LP       = L * prate *(1/3600.) #(J/kg * mm/hr * 1kg/m2 * 1hr/3600s --> W/m2)
    Rtoa = toa_swdn - toa_swup + 0.       - toa_lwup
    Rs   = boa_swdn - boa_swup + boa_lwdn - boa_lwup
    Ra   = Rtoa - Rs
    dFa  = Ra + LP + SH
    H    = DIV_DRY_STATIC_ENERGY_data
    Ea   = Ra + LP + SH - H
    ecmwf_flux_hourly = [ {'data':STRD['lat'], 'name':'lat',  'long':'latitude'},
                          {'data':STRD['lon'], 'name':'lon',  'long':'longitude'},
                          {'data':boa_lwdn, 'name':'boa_lwdn','long':'bottom of atmosphere downwelling longwave radiative flux'},
                          {'data':boa_lwup, 'name':'boa_lwup','long':'bottom of atmosphere upwelling longwave radiative flux'},
                          {'data':boa_swdn, 'name':'boa_swdn','long':'bottom of atmosphere downwelling shortwave radiative flux'},
                          {'data':boa_swup, 'name':'boa_swup','long':'bottom of atmosphere upwelling shortwave radiative flux'},
                          {'data':toa_lwup, 'name':'toa_lwup','long':'top of atmosphere upwelling longwave radiative flux'},
                          {'data':toa_swdn, 'name':'toa_swdn','long':'top of atmosphere downwelling shortwave radiative flux'},
                          {'data':toa_swup, 'name':'toa_swup','long':'top of atmosphere upwelling shortwave radiative flux'},
                          {'data':SH, 'name':'SH','long':'accumulated surface sensible heat flux'},
                          {'data':prate, 'name':'prate','long':'accumulated total precipitation rate'},
                          {'data':LP, 'name':'LP','long':'diabatic heating by precipitation'},
                          {'data':Rtoa, 'name':'Rtoa','long':'net top of atmopshere radiation flux'},
                          {'data':Rs, 'name':'Rs','long':'net surface radiation flux'},
                          {'data':Ra, 'name':'Ra','long':'net atmospheric column radiation flux'},
                          {'data':dFa, 'name':'dFa','long':'residual constrained dry static energy flux'},
                          {'data':H  , 'name':'H','long':'dry static energy computed from vertical integral of the column dry static energy'},
                          {'data':Ea , 'name':'Ea','long':'time rate of change of the energy content of an atmospheric column of unit horizontal area extending from the surface to the top of the atmosphere'}]



    boa_lwdn = STRD['daily']
    boa_lwup = STRD['daily'] - STR['daily']
    boa_swdn = SSRD['daily']
    boa_swup = SSRD['daily'] - SSR['daily']
    toa_lwup = 0. - TTR['daily']
    toa_swdn = TISR['daily']
    toa_swup = TISR['daily'] - TSR['daily']
    prate    = TP['daily'] * 1000. * 3600. #(m/s * 3600. * 1000mm/m --> mm/hr)
    SH     = 0. - SSHF['daily']
    LP       = L * prate *(1/3600.) #(J/kg * mm/hr * 1kg/m2 * 1hr/3600s --> W/m2)
    Rtoa = toa_swdn - toa_swup + 0.       - toa_lwup
    Rs   = boa_swdn - boa_swup + boa_lwdn - boa_lwup
    Ra   = Rtoa - Rs
    dFa  = Ra + LP + SH
    H = DIV_DRY_STATIC_ENERGY_daily
    Ea   = Ra + LP + SH - H
    ecmwf_flux_daily = [ {'data':STRD['lat'], 'name':'lat',  'long':'latitude'},
                          {'data':STRD['lon'], 'name':'lon',  'long':'longitude'},
                          {'data':boa_lwdn, 'name':'boa_lwdn','long':'bottom of atmosphere downwelling longwave radiative flux'},
                          {'data':boa_lwup, 'name':'boa_lwup','long':'bottom of atmosphere upwelling longwave radiative flux'},
                          {'data':boa_swdn, 'name':'boa_swdn','long':'bottom of atmosphere downwelling shortwave radiative flux'},
                          {'data':boa_swup, 'name':'boa_swup','long':'bottom of atmosphere upwelling shortwave radiative flux'},
                          {'data':toa_lwup, 'name':'toa_lwup','long':'top of atmosphere upwelling longwave radiative flux'},
                          {'data':toa_swdn, 'name':'toa_swdn','long':'top of atmosphere downwelling shortwave radiative flux'},
                          {'data':toa_swup, 'name':'toa_swup','long':'top of atmosphere upwelling shortwave radiative flux'},
                          {'data':SH, 'name':'SH','long':'accumulated surface sensible heat flux'},
                          {'data':prate, 'name':'prate','long':'accumulated total precipitation rate'},
                          {'data':LP, 'name':'LP','long':'diabatic heating by precipitation'},
                          {'data':Rtoa, 'name':'Rtoa','long':'net top of atmopshere radiation flux'},
                          {'data':Rs, 'name':'Rs','long':'net surface radiation flux'},
                          {'data':Ra, 'name':'Ra','long':'net atmospheric column radiation flux'},
                          {'data':dFa, 'name':'dFa','long':'residual constrained dry static energy flux'},
                          {'data':H  , 'name':'H','long':'dry static energy computed from vertical integral of the column dry static energy'},
                          {'data':Ea , 'name':'Ea','long':'time rate of change of the energy content of an atmospheric column of unit horizontal area extending from the surface to the top of the atmosphere'}]

    ecmwf_flux = {'hourly':ecmwf_flux_hourly,'daily':ecmwf_flux_daily}
    return ecmwf_flux
