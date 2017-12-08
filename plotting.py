import sys
import os
from netCDF4 import Dataset
from subroutines import *
import matplotlib.pyplot as plt
import matplotlib.text as text
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import (LATITUDE_FORMATTER,
                                   LONGITUDE_FORMATTER)
import numpy.ma as ma

def plot_map(lons,lats,data,minv,maxv,titl,figtitl):
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True)
    gl.xlocator = mticker.FixedLocator([-180, -60, 60, 180])
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    p = ax.pcolormesh(lons, lats, data, vmin=minv, vmax=maxv)
    clb=plt.colorbar(p, orientation="horizontal", pad=0.2)
    clb.ax.set_title(titl)
    #plt.show()
    plt.savefig(figtitl)
    plt.close()
