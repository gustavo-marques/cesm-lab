#!/usr/bin/env /home/gmarques/miniconda2/bin/python

# generates MOM6 diagnostics

import argparse
import xarray as xr
from netCDF4 import MFDataset, Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import m6plot
import warnings
import os

class MyError(Exception):
  """
  Class for error handling
  """
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Computes MOM6 diagnostics.
      ''',
  epilog='Written by Gustavo Marques.')

  parser.add_argument('-case_name', type=str, default='test',
      help='''Case name. Default is test.''')

  parser.add_argument('-geometry', type=str, default='ocean_geometry.nc',
      help='''The name of the ocean geometry file. Default is ocean_geometry.nc''')

  parser.add_argument('-surface_file', type=str, default='sfc_daily__*',
      help='''The name of the daily surface file. Default is sfc_daily__*.nc.''')

  parser.add_argument('-forcing_file', type=str, default='forcing_*',
      help='''The name of the forcing file with surface fluxes. Default is forcing_*.nc.''')

  parser.add_argument('-month_file', type=str, default='ocean_month__*',
      help='''The name of monthly mean file. Default is ocean_month__*''')

  parser.add_argument('-month_z_file', type=str, default='ocean_month_z__*',
      help='''The name of monthly mean file remapped to a z coordinate. Default is ocean_month_z__*''')

  parser.add_argument('-prog_file', type=str, default='prog__*',
      help='''The name of prognostic ocean file. Default is prog__*''')

  parser.add_argument('-ndays', type=int, default=2,
      help='''Number of days to skip when computing time series. Default is 2.''')

  parser.add_argument('-year_start', type=int, default=80,
      help='''Start year to compute averages. Default is 80.''')

  parser.add_argument('-year_end', type=int, default=100,
      help='''End year to compute averages. Default is 100.''')

  parser.add_argument('-to_netcdf', help='''Save data into a netCDF file.''', 
      action="store_true")
  
  parser.add_argument('-savefigs', help='''Save figures in a PNG format.''',
      action="store_true")

  optCmdLineArgs = parser.parse_args()
  global case_name
  case_name = optCmdLineArgs.case_name
  driver(optCmdLineArgs)

#-- This is where all the action happends, i.e., functions for each diagnostic are called.

def driver(args):
  os.system('mkdir PNG')
  os.system('mkdir ncfiles')
  # mom6 grid
  grd = MOM6grid(args.geometry)

  # extract mean surface latlon time series from forcing and surface files
  #mean_latlon_time_series(args, grd, ['SSS','SST','MLD_003','SSH','hfds','PRCmE','taux','tauy'])

  # FIXME: SSU and SSV need to be plotted on u and v points, instead of tracer points
  mean_latlon_plot(args,grd,['SSH','SSS','SST','MLD_003','SSU','SSV','hfds','PRCmE','taux','tauy'])

  return

# -- mean surface latlon time series from forcing and surface files
def mean_latlon_time_series(args, grd, variables):
  # check if stats have been done
  filename = str('ncfiles/%s_stats.nc' % (args.case_name))
  if os.path.isfile(filename):
    print("File {} already exists! Moving to the next one...\n".format(filename))
  else:
    nc1 = xr.open_mfdataset(args.surface_file)
    nc2 = xr.open_mfdataset(args.forcing_file)
    tm = len(nc1.time)

    #TODO check if time in nc1 and nc2 are the same, otherwise give error

    print("### About to extract spatial mean time series for following variables {} ### \n".format(variables))

    # construct a list with the units for each variable
    units = []
    for var in variables:
      if var in nc1.variables:
        units.append(nc1.variables[var].attrs['units'])
      elif var in nc2.variables:
        units.append(nc2.variables[var].attrs['units'])
      else:
        print 'Error!'

    # create Dataset
    dtime = nc1.time[0:tm:args.ndays].values
    ds = create_xarray_dataset(variables,units,dtime)
    # loop in time
    n = 0
    for t in range(0,tm,args.ndays):
      print ("==> ' + 'step # {} out of {}  ...\n".format(t+1,tm))
      for var in variables:
        if var in nc1.variables:
          data = nc1.variables[var][t,:].values
        else: 
          data = nc2.variables[var][t,:].values

        data = np.ma.masked_invalid(data)
        # get stats
        sMin, sMax, mean, std, rms = m6plot.myStats(data, grd.Ah)

        # update Dataset
        ds[var][0,n] = sMin; ds[var][1,n] = sMax; ds[var][2,n] = mean
        ds[var][3,n] = std; ds[var][4,n] = rms

      # increment n
      n += 1

    nc1.close(); nc2.close()

    if args.to_netcdf:
      # save in a netcdf file
      ds.to_netcdf('ncfiles/'+args.case_name+'_stats.nc')

    if args.savefigs:
      # save figures
      n = 0
      for var in variables:
        # TODO: this can be done in a more elegant way...
        plt.figure()
        f, ax = plt.subplots(5, sharex=True)
        ax[0].plot(ds['time'], ds[var][2,:])
        ax[0].set_title(r'%s, [%s]'%(var,units[n]))
        ax[0].set_ylabel('Mean')
        ax[1].plot(ds['time'], ds[var][0,:])
        ax[1].set_ylabel('Min.')
        ax[2].plot(ds['time'], ds[var][1,:])
        ax[2].set_ylabel('Max.')
        ax[3].plot(ds['time'], ds[var][3,:])
        ax[3].set_ylabel('Std')
        ax[4].plot(ds['time'], ds[var][4,:])
        ax[4].set_ylabel('Rms')
        ax[4].set_xlabel('Year')
        plt.savefig('PNG/%s_stats.png'%(var))
        plt.close()
        n += 1


  return

# -- time-mean latlon plot
def mean_latlon_plot(args, grd, variables):
  xr_mfopen_args = {'decode_times':False,
                  'decode_coords':False,
                  'data_vars':'minimal'}

  nc1 = xr.open_mfdataset(args.surface_file, **xr_mfopen_args)
  nc2 = xr.open_mfdataset(args.forcing_file, **xr_mfopen_args)

  if not nc1.time.attrs['calendar'] == 'NOLEAP':
    raise NameError('Only noleap calendars are supported at this moment!')

  # TODO: assign a new variable called time_years
  # convert time in years
  nc1['time'] = nc1.time/365.
  nc2['time'] = nc2.time/365.

  ti = args.year_start
  tf = args.year_end
 
  # TODO: check if data includes years between ti and tf 
 
  for var in variables:
    filename = str('PNG/%s.png' % (var))
    if os.path.isfile(filename):
      print (' \n' + '==> ' + '{} has been saved, moving to the next one ...\n' + ''.format(var))
    else:
      print("About to plot surface time-average for variable {}... \n".format(var))
      if var in nc1.variables:
        data = np.ma.masked_invalid(nc1[var].sel(time=slice(ti,tf)).mean('time').values)
        units = nc1[var].attrs['units']
      elif var in nc2.variables:
        data = np.ma.masked_invalid(nc2[var].sel(time=slice(ti,tf)).mean('time').values)
        units = nc2[var].attrs['units']
      else:
        raise NameError('Variable {} does not exist in {} or {}!'.format(var,args.surface_file,args.forcing_file))

      if args.savefigs:    
        #long_name = nc.variables[var].long_name
        m6plot.xyplot( data , grd.geolon, grd.geolat, area=grd.Ah,
          suptitle=args.case_name,
          title=r'%s, [%s] averaged over years %i-%i'%(var,units,args.year_start,args.year_end),
          extend='both',
          save=filename)
      else:
        m6plot.xyplot( data , grd.geolon, grd.geolat, area=grd.Ah,
          suptitle=args.case_name,
          title=r'%s, [%s] averaged over years %i-%i'%(var,units,args.year_start,args.year_end),
          extend='both',
          show=True)
     

  nc1.close(); nc2.close()
  return

# -- create a xarray Dataset given variables, units and time
def create_xarray_dataset(variables,units,time):
   ds = xr.Dataset() 
   # TODO: fix the time using pandas, # of years are limited, see link below
   # http://pandas-docs.github.io/pandas-docs-travis/timeseries.html#timestamp-limitations
   # It should be days since 0001-01-01 00:00:00
   dt=(time[1]-time[0]).days
   f=str(int(dt))+'D'
   ds.coords['time'] = pd.date_range(str(time[0]),freq=f, periods=len(time))
   print "ds.coords['time']",ds.coords['time']
   #ds.coords['time'].attrs = {'units': 'days since 1900-01-01 00:00:00'}
   stats = ['min','max','mean','std','rms']
   ds.coords['stats'] = stats
   tm = len(time)
   for i in range(len(variables)):
     ds[variables[i]] = (('stats', 'time'), np.zeros((len(stats),tm)))
     ds[variables[i]].attrs['units'] = units[i]

   return ds

# -- process log files and plot the data
def get_log_data(args):


  return

# -- return MOM6 grid object
def MOM6grid(grd_file):
    """
    grd = MOM6grid()
    """

    # create an empty class object
    class MOM6_grd:
        pass

    # open grid file
    #nc = Dataset(grd_file)
    nc = xr.open_dataset(grd_file)
    # fill grid object
    for var in nc.variables:
       #dummy = str("MOM6_grd.%s = nc.variables[var][:]"% (var))
       dummy = str("MOM6_grd.%s = nc.%s[:].values"% (var,var))
       exec(dummy)

    # close netcdf file
    nc.close()
    print('MOM6 grid successfully loaded... \n')
    return MOM6_grd

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()



