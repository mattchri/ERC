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
#python2.7 -i process_energetics.py 2017 8 14 fv1.0
#
# 02/11/17, MC: upload initial version of the code to the repo
#---------------------------------------------------------------
import sys
import os
from netCDF4 import Dataset
from subroutines import *
import matplotlib.pyplot as plt
import matplotlib.text as text

#---------------------------------------------------------------
# input variables read from command line
#---------------------------------------------------------------
year  = int(sys.argv[1])
month = int(sys.argv[2])
dom   = int(sys.argv[3])
version   = sys.argv[4]

#---------------------------------------------------------------
#Paths
#---------------------------------------------------------------
outpath = '/group_workspaces/jasmin2/aopp/mchristensen/erc/energetics/'
testFile='/group_workspaces/cems2/nceo_generic/model_data/ERA_INTERIM/ERA_INTERIM_ENERGETICS_'+str(year).zfill(4)+str(month).zfill(2)+'.nc'
#---------------------------------------------------------------

#TP = read_ecmwf_variable('TP',year,month,dom,testFile)
#TP = read_ecmwf_variable('TP',year,month,dom,'')

#ecmwf_flux = ( process_ecmwf_variables(year,month,dom,'') )['daily']
ecmwf_flux = ( process_ecmwf_variables(year,month,dom,testFile) )['daily']
varNames = np.array([x['name'] for x in ecmwf_flux])

#Fetch array size
sz=(ecmwf_flux[10])['data'].shape
xDim=sz[0]
yDim=sz[1]
print(xDim,yDim)

VARS = ['Ra','LP','SH','dFa','H']
VARS_PLOT = [r'$R_a$','LP','SH',r'$\Delta$'+r'$F_a$',r'$H$']
ARR  = np.zeros( [xDim,len(VARS)] )
for i in range(len(VARS)):
    index = ((np.where(varNames == VARS[i]))[0])[0]
    tmp=(ecmwf_flux[index])['data']
    ARR[:,i] = np.mean(tmp, axis=1)

lat = ((ecmwf_flux[0])['data'])[:]
plt.plot(lat, ARR[:,0], 'r--', 
         lat, ARR[:,1], 'b--', 
         lat, ARR[:,2], 'g--', 
         lat, ARR[:,3], 'y--',
         lat, ARR[:,4], 'm--')
plt.ylabel('Watts/Square Meter')
plt.xlabel('Latitude')
plt.legend((VARS_PLOT),loc='upper left', shadow=True)
plt.title('Atmospheric Energy Budget')
plt.show()
#plt.savefig('/home/users/mchristensen/Desktop/test.png')
