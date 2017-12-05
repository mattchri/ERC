#-----------------------------------------------------------------
#Module Containing all Subroutines
#required to run the SMIP code
#The first group of subroutines deal with extracting
#times and converting them to julian day or 
#collocating files together.
#The second group deals with collocation.
#
# 02/11/17, MC: upload initial version of the code to the repo
#-----------------------------------------------------------------
from jdcal import *
import sys
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import pdb

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
#COLLOCATE_FILES
#Input the base SEVIRI files (merged output in this case) and
#all other file types (e.g. imerg, ECMWF, ect) to collocate to SEVIRI
#
#Input list of files for each data product
#Output list of collocated files from each product
#-----------------------------------------------------------------
def collocate_files(infiles,dt):
    #SEVIRI Base files Information
    BaseDict     = infiles[0] #by default zeroth element
    BaseFileInfo = BaseDict['files']
    BaseFiles = [item['filename'] for item in BaseFileInfo]
    BaseJDAY = [item['JDAY'] for item in BaseFileInfo]
    BasePath = BaseDict['path']

    #Define Matched Array
    nFileTypes = len(infiles)
    nFiles = len(BaseFiles)
    print '# Data Products: ',nFileTypes
    print '# Base (collocation) Files: ',nFiles
    matchedFiles = np.empty([nFileTypes,nFiles], dtype="S512")
    for iT in range(nFiles):
        matchedFiles[0,iT] = BasePath+BaseFiles[iT]

    i=1 #initialised value required 1st file type in for loop
    #Loop over each data product
    for iT in range(nFileTypes-1):
        TMPDict     = infiles[i] #by default zeroth element
        TMPFileInfo = TMPDict['files']
        TMPFiles    = [item['filename'] for item in TMPFileInfo]
        TMPJDAY     = [item['JDAY'] for item in TMPFileInfo]
        TMPJDAY1 = [iii[1] for iii in TMPJDAY]
        TMPPath  = TMPDict['path']
        y = np.asarray(TMPJDAY1)

        #Loop over each file
        for j in range(nFiles):
            x = float((BaseJDAY[j])[1]) #julian day temp var
            fid = np.where(abs(y-x) == min(abs(y-x)))
            z = int(fid[0])
            #print i,j,BaseFiles[j],'  ',TMPFiles[z]
            matchedFiles[ i, j ] = TMPPath + TMPFiles[z]

        i=i+1 #start with first data product (+1)
    return matchedFiles

#-----------------------------------------------------------------
#COLLOCATE_IMERG_SEVIRI
#Routine collocates the IMERG lat/lon to the SEVIRI lat/lon grid
#The top-level code produces a save file so this only should be run
#one time since the lat/lon fields are invariant for both sensors
#Input
# lonIM, latIM, lonSV, latSV
#Output
#index for longitude and latitude arrays
#-----------------------------------------------------------------
def collocate_imerg_seviri(lonIM,latIM,lonSV,latSV):
    print('Collocating IMERG to SEVIRI (~1hr processing time)')
    xDim = (lonSV.shape)[0]
    yDim = (latSV.shape)[1]
    LATidIM = np.zeros( (xDim, yDim) )
    LONidIM = np.zeros( (xDim, yDim) )
    for iX in range(xDim):
        print(iX)
        for iY in range(yDim):
            tmpLat = latSV[iX,iY]
            tmpLon = lonSV[iX,iY]
            if tmpLon > -999. and tmpLat > -999.:
                LONidIM[iX,iY] = int(((np.where(abs(lonIM-tmpLon) == min(abs(lonIM-tmpLon))) )[0])[0])
                LATidIM[iX,iY] = int(((np.where(abs(latIM-tmpLat) == min(abs(latIM-tmpLat))) )[0])[0])
            else:
                LONidIM[iX,iY] = -999.
                LATidIM[iX,iY] = -999.                   
    IMERG_ID = {'latid':LATidIM, 'lonid':LONidIM, 'xdim':xDim, 'ydim':yDim}
    return IMERG_ID


#Create NetCDF file
def make_netcdf_file(ncname,vName,vData,vLong,vStan,vUnit,vFill):
    print('Creating: '+ncname)
    sz = vData.shape
    nVars= sz[0]
    xDim = sz[1]
    yDim = sz[2]

    f = Dataset(ncname,'w', format='NETCDF4_CLASSIC')
    f.createDimension('xdim', xDim)
    f.createDimension('ydim', yDim)
    
    for iV in range(len(vName)):
        tmpstr = vName[iV]
        temp = f.createVariable(tmpstr, 'f4', ('xdim','ydim'), zlib=True)
        temp[:,:] = vData[iV,:,:]
        temp.standard_name  = vStan[iV]
        temp.long_name  = vLong[iV]
        temp.units = vUnit[iV]
        temp.FillValue = vFill[iV]

    #global attributes
    today = datetime.today()
    f.description = "Satellite Model Integrated Product"
    f.history = "Created " + today.strftime("%d/%m/%y")
    f.close()


#-----------------------------------------------------------------
#Read Accumulated (gafs) ECWMF variable from /badc/ on JASMIN
#Inputs
#variable name (e.g. 'TP' total precipitation)
#year [int]
#month [int]
#dom day of the month [int]
#Output
#ARR [accumulated field]
#-----------------------------------------------------------------
def read_ecmwf_variable(vname,year,month,dom,testfile):

    #Read Downloaded ECWMF Data
    if(len(testfile) > 1):
        ncfile = Dataset( testfile, mode='r')
        time = ncfile.variables['time'][:]
        #units: hours since 1900-01-01 00:00:0.0
        jdayST = gcal2jd( 1900, 1, 1 )
        jday   = (jdayST[0],jdayST[1]+time/24.)
        sort_index = np.argsort(jday[1])
        jday_sorted= (jdayST[0], (jday[1])[sort_index] )
        alljday = jday_sorted[1]

        #Find indices for corresponding input day
        jdayInput = gcal2jd(year,month,dom)
        injday = jdayInput[1]

        ecmwf_tID = np.where( alljday.astype(int) == int(injday))

        id_sorted = sort_index[ecmwf_tID]

        ncfile = Dataset( testfile, mode='r')
        scaleF  = ncfile.variables[vname.lower()].scale_factor
        add_off = ncfile.variables[vname.lower()].add_offset
        fillVal = ncfile.variables[vname.lower()]._FillValue
        unitStr = ncfile.variables[vname.lower()].units

        ARR = (ncfile.variables[vname.lower()])[id_sorted[:],:,:]

    # Read ECWMF Data from Jasmin
    if(len(testfile) == 0):
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


    #if(vname == 'ishf'):
    #    pdb.set_trace()

    if((np.asarray(ARR[:,0,0]))[1] == fillVal):
        print('Instantaneous: '+vname)
        aARR = instantaneous_ecmwf_3hr(ARR)
    if((np.asarray(ARR[:,0,0]))[1] != fillVal):
        print('Accumulated: '+vname)
        aARR = accumulate_ecmwf_3hr(ARR)


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
    aARR     = np.empty( [sz[0],sz[1],sz[2]] )
    aARR[0,:,:] = ( ARR[0,:,:]-0.        ) / tstep
    aARR[1,:,:] = ( ARR[1,:,:]-ARR[2,:,:]) / tstep
    aARR[2,:,:] = ( ARR[2,:,:]-ARR[3,:,:]) / tstep
    aARR[3,:,:] = ( ARR[3,:,:]-ARR[4,:,:]) / tstep
    aARR[4,:,:] = ( ARR[4,:,:]-0.        ) / tstep
    aARR[5,:,:] = ( ARR[5,:,:]-ARR[5,:,:]) / tstep
    aARR[6,:,:] = ( ARR[6,:,:]-ARR[6,:,:]) / tstep
    aARR[7,:,:] = ( ARR[7,:,:]-ARR[7,:,:]) / tstep
    return aARR    

def instantaneous_ecmwf_3hr(ARR):
    #Express quantity as instantaneous
    sz=ARR.shape
    aARR     = np.empty( [sz[0],sz[1],sz[2]] )
    aARR[0,:,:] = ARR[0,:,:]
    aARR[1,:,:] = ARR[0,:,:]
    aARR[2,:,:] = ARR[2,:,:]
    aARR[3,:,:] = ARR[2,:,:]
    aARR[4,:,:] = ARR[4,:,:]
    aARR[5,:,:] = ARR[4,:,:]
    aARR[6,:,:] = ARR[6,:,:]
    aARR[7,:,:] = ARR[6,:,:]
    return aARR    


def process_ecmwf_variables(year,month,dom,testfile):
    #Read Data
    STRD = read_ecmwf_variable('STRD',year,month,dom,testfile)
    TISR = read_ecmwf_variable('TISR',year,month,dom,testfile)
    TSR  = read_ecmwf_variable('TSR',year,month,dom ,testfile)
    TTR  = read_ecmwf_variable('TTR',year,month,dom ,testfile)
    SSRD = read_ecmwf_variable('SSRD',year,month,dom,testfile)
    SSR  = read_ecmwf_variable('SSR',year,month,dom ,testfile)
    STRD = read_ecmwf_variable('STRD',year,month,dom,testfile)
    STR  = read_ecmwf_variable('STR',year,month,dom ,testfile)
    TP   = read_ecmwf_variable('TP',year,month,dom  ,testfile)
    SSHF = read_ecmwf_variable('SSHF',year,month,dom,testfile)
    sz=SSHF['data'].shape
    DIV_DRY_STATIC_ENERGY_data = np.zeros( [sz[0],sz[1],sz[2]] )
    DIV_DRY_STATIC_ENERGY_daily = np.zeros( [sz[1],sz[2]] )

    if(len(testfile) > 1):
        #NOTE THESE QUANTITIES ARE NOT ACCUMULATED!!!
        p82 = read_ecmwf_variable('p82.162',year,month,dom,testfile)
        p83 = read_ecmwf_variable('p83.162',year,month,dom,testfile)
        p85 = read_ecmwf_variable('p85.162',year,month,dom,testfile)
        #pdb.set_trace()
        DIV_DRY_STATIC_ENERGY_data  = (p82['data']+p83['data']+p85['data']   )
        DIV_DRY_STATIC_ENERGY_daily = (p82['daily']+p83['daily']+p85['daily'])
        #pdb.set_trace()
        ISHF = read_ecmwf_variable('ishf',year,month,dom,testfile)
        ISHF_data = ISHF['data'] * 3. * 3600.
        ISHF_daily = ISHF['daily'] * 3. * 3600.

    # Calculate Fluxes
    L = 2.5e6 #latent heat of vaporization at 0C (J/kg)
    boa_lwdn = STRD['data']
    boa_lwup = STR['data'] - STRD['data']
    boa_swdn = SSRD['data']
    boa_swup = SSR['data'] - SSRD['data']
    toa_lwup = TTR['data']
    toa_swdn = TISR['data']
    toa_swup = TSR['data'] - TISR['data']
    prate    = TP['data'] * 1000. * 3600. #(m/s * 3600. * 1000mm/m --> mm/hr)
    #SH     = 0. - SSHF['data']  #not sure if this is right!!!!!
    SH     = 0. - ISHF_data
    LP       = L * prate *(1/3600.) #(J/kg * mm/hr * 1kg/m2 * 1hr/3600s --> W/m2)
    Rtoa = toa_swdn - toa_swup + 0.       - toa_lwup
    Rs   = boa_swdn - boa_swup + boa_lwdn - boa_lwup
    Ra   = Rtoa - Rs
    dFa  = Ra + LP + SH
    H    = DIV_DRY_STATIC_ENERGY_data
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
                          {'data':LP, 'name':'LP','long':'bottom of atmosphere downwelling longwave radiative flux'},
                          {'data':Rtoa, 'name':'Rtoa','long':'bottom of atmosphere downwelling longwave radiative flux'},
                          {'data':Rs, 'name':'Rs','long':'bottom of atmosphere downwelling longwave radiative flux'},
                          {'data':Ra, 'name':'Ra','long':'bottom of atmosphere downwelling longwave radiative flux'},
                          {'data':dFa, 'name':'dFa','long':'bottom of atmosphere downwelling longwave radiative flux'},
                          {'data':H  , 'name':'H','long':'bottom of atmosphere downwelling longwave radiative flux'}]

    boa_lwdn = STRD['daily']
    boa_lwup = STR['daily'] - STRD['daily']
    boa_swdn = SSRD['daily']
    boa_swup = SSR['daily'] - SSRD['daily']
    toa_lwup = TTR['daily']
    toa_swdn = TISR['daily']
    toa_swup = TSR['daily'] - TISR['daily']
    prate    = TP['daily'] * 1000. * 3600. #(m/s * 3600. * 1000mm/m --> mm/hr)
    #SH     = 0. - SSHF['daily']
    SH     = 0. - ISHF_daily
    LP       = L * prate *(1/3600.) #(J/kg * mm/hr * 1kg/m2 * 1hr/3600s --> W/m2)
    Rtoa = toa_swdn - toa_swup + 0.       - toa_lwup
    Rs   = boa_swdn - boa_swup + boa_lwdn - boa_lwup
    Ra   = Rtoa - Rs
    dFa  = Ra + LP + SH
    H = DIV_DRY_STATIC_ENERGY_daily
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
                          {'data':LP, 'name':'LP','long':'bottom of atmosphere downwelling longwave radiative flux'},
                          {'data':Rtoa, 'name':'Rtoa','long':'bottom of atmosphere downwelling longwave radiative flux'},
                          {'data':Rs, 'name':'Rs','long':'bottom of atmosphere downwelling longwave radiative flux'},
                          {'data':Ra, 'name':'Ra','long':'bottom of atmosphere downwelling longwave radiative flux'},
                          {'data':dFa, 'name':'dFa','long':'bottom of atmosphere downwelling longwave radiative flux'},
                          {'data':H  , 'name':'H','long':'bottom of atmosphere downwelling longwave radiative flux'}]

    ecmwf_flux = {'hourly':ecmwf_flux_hourly,'daily':ecmwf_flux_daily}
    return ecmwf_flux
