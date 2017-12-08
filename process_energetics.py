#CURRENTLY, THE INSTANTANEOUS FIELDS ARE NOT CORRECT BECAUSE
#ONLY 00:00 AND 12:00 ARE BEING USED TO CONSTRUCT THE MEAN
#CAUSING LARGE VARIATION ACROSS LONGITUDE.

#---------------------------------------------------------------
#TOP LEVEL CODE
#---------------------------------------------------------------
#PROCESS_ENERGETICS
#
#Calculate the atmospheric energy budget
#Ra + LP + SH = dFa
#
#Each term
#Ra: net radiative heating of the atmospheric column
# --> Ra = Rnet(TOA) - Rnet(BOA)
#       -->Rnet = SWDN - SWUP + LWDN - LWUP
#LP: heating of the atmospheric column by latent heat relase during precipitation
#SH: sensible heat transfer from the surface to the atmosphere
#dFa: horizontal divergence of energy out of the column by transpoft in the atmosphere
#
#Example
#python2.7 -i process_energetics.py
#
# 02/11/17, MC: upload initial version of the code to the repo
# 08/11/17, MC: included plotting code and updated averaging method
#               identified bug in instantaneous code but have not
#               corrected it.
#---------------------------------------------------------------
import sys
import os
from subroutines import *
from plotting import *

figpath = '/home/users/mchristensen/Desktop/ERC/'
figpath = '/group_workspaces/cems/cloud_ecv/public/temp_transfer/erc/'

#Compute Fluxes based on monthly, seasonal, and annual

years = [2015]
months = [1,2,3,4,5,6,7,8,9,10,11,12]
for tY in range(len(years)):
    year=years[tY]
    #Fetch File
    ecmwfFile ='/group_workspaces/cems2/nceo_generic/model_data/ERA_INTERIM/ERA_INTERIM_ENERGETICS_monthly_'+str(year).zfill(4)+'.nc'
    data_month = []
    for tM in range(len(months)):
        month=months[tM]
        dom=15
        print('Processing: ',str(months[tM]).zfill(2))
        flux_month = process_ecmwf_variables(year,month,dom,ecmwfFile)['daily']
        data_month.append(flux_month)
        
varNames = np.array([x['name'] for x in data_month[0]])
VARS = ['Ra','LP','SH','dFa','H','Ea','toa_swdn','toa_swup','toa_lwup','boa_lwup','boa_lwdn','boa_swdn','boa_swup']
VARS_PLOT = [r'$R_a$','LP','SH',r'$\Delta$'+r'$F_a$',r'$H$',r'$E_a$','toa_swdn','toa_swup','toa_lwup','boa_lwup','boa_lwdn','boa_swdn','boa_swup']
VARS_UNIT = r'W/m$^{2}$'
VARS_titl = [VARS_PLOT[0]+' ('+VARS_UNIT+')', VARS_PLOT[1]+' ('+VARS_UNIT+')',
             VARS_PLOT[2]+' ('+VARS_UNIT+')', VARS_PLOT[3]+' ('+VARS_UNIT+')',
             VARS_PLOT[4]+' ('+VARS_UNIT+')', VARS_PLOT[5]+' ('+VARS_UNIT+')',
             VARS_PLOT[6]+' ('+VARS_UNIT+')', VARS_PLOT[7]+' ('+VARS_UNIT+')',
             VARS_PLOT[8]+' ('+VARS_UNIT+')', VARS_PLOT[9]+' ('+VARS_UNIT+')',
             VARS_PLOT[10]+' ('+VARS_UNIT+')', VARS_PLOT[11]+' ('+VARS_UNIT+')',
             VARS_PLOT[12]+' ('+VARS_UNIT+')']
VARS_mins = [-400, 0   , -60, -200, -400, -100, 0  ,0  ,100, 100, 100, 0  ,0]
VARS_maxs = [-100, 300 ,  60,  200,  400,  100, 420,300,300, 500, 400, 300,200]
nvars = len(VARS)
lats = ((flux_month[0])['data'])[:]
lons = ((flux_month[1])['data'])[:]


#Fetch dimensions
ecmwf_flux = data_month[0]
index = ((np.where(varNames == 'Ra'))[0])[0]
tmpdata = (ecmwf_flux[index])['data']
xdim = (tmpdata.shape)[0]
ydim = (tmpdata.shape)[1]

#Monthly Mean
varData_month = np.zeros( [xdim, ydim, 12, nvars] )
for tM in range(len(months)):
    ecmwf_flux = data_month[tM]
    for iV in range(nvars):
        index = ((np.where(varNames == VARS[iV]))[0])[0]
        varData_month[:,:,tM,iV]=(ecmwf_flux[index])['data']

#Annual Mean
print('Averaging over Year')
varData_annual = np.zeros( [xdim, ydim, nvars] )
for iV in range(nvars):
    varData_annual[:,:,iV] = np.mean(varData_month[:,:,:,iV],axis=2)

#Seasonal Mean
print('Averaging over Season')
sid = np.array( [ [12,1,2], [3,4,5 ] , [6,7,8], [9,10,11] ] )
varData_season = np.zeros( [xdim, ydim, 4, nvars] )
for iS in range(4):
    for iV in range(nvars):
        varData_season[:,:,iS,iV] = np.mean(varData_month[:,:,sid[iS,:]-1,iV], axis=2)

ARR = np.mean( varData_annual, axis=1 )
#Plot Annual Meridional Distribution
plt.plot(lats, ARR[:,0], 'r--', 
         lats, ARR[:,1], 'b--', 
         lats, ARR[:,2], 'g--', 
         lats, ARR[:,3], 'y--')
plt.ylabel('Watts/Square Meter')
plt.xlabel('Latitude')
plt.legend((VARS_PLOT[0:3]),loc='upper left', shadow=True)
plt.title('Atmospheric Energy Budget')
#plt.show()
plt.savefig(figpath+'merid_constrained_annual.png')
plt.close()

#Plot Annual Meridional Distribution
plt.plot(lats, ARR[:,0], 'r--', 
         lats, ARR[:,1], 'b--', 
         lats, ARR[:,2], 'g--', 
         lats, ARR[:,3], 'y--',
         lats, ARR[:,4], 'm--',
         lats, ARR[:,5], 'b-')
plt.ylabel('Watts/Square Meter')
plt.xlabel('Latitude')
plt.legend((VARS_PLOT[0:6]),loc='upper left', shadow=True)
plt.title('Atmospheric Energy Budget')
#plt.show()
plt.savefig(figpath+'merid_annual.png')
plt.close()

#Geographic distribution
for iV in range(nvars):
    data = varData_annual[:,:,iV]
    minv = VARS_mins[iV]
    maxv = VARS_maxs[iV]
    titl = VARS_titl[iV]
    figtitl = figpath+'global_annual_'+VARS[iV]+'.png'
    plot_map(lons,lats,data,minv,maxv,titl,figtitl)
      
