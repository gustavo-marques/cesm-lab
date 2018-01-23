#!/usr/bin/env python

# generates MOM6 diagnostics

import argparse
from netCDF4 import MFDataset, Dataset
import numpy as np
import matplotlib.pyplot as plt
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
  epilog='Written by Gustavo Marques, Dec. 2017.')

  parser.add_argument('-case_name', type=str, default='test',
      help='''Case name. Default is test.''')

  parser.add_argument('-geometry', type=str, default='ocean_geometry.nc',
      help='''The name of the ocean geometry file. Default is ocean_geometry.nc''')

  parser.add_argument('-surface_file', type=str, default='sfc_daily*',
      help='''The name of the daily surface file. Default is sfc_daily*.nc.''')

  parser.add_argument('-forcing_file', type=str, default='forcing*',
      help='''The name of the forcing file with surface fluxes. Default is forcing*.nc.''')

  parser.add_argument('-month_file', type=str, default='ocean_month__*',
      help='''The name of monthly mean file. Default is ocean_month__*''')

  parser.add_argument('-month_z_file', type=str, default='ocean_month_z__*',
      help='''The name of monthly mean file remapped to a z coordinate. Default is ocean_month_z__*''')

  parser.add_argument('-prog_file', type=str, default='prog__*',
      help='''The name of prognostic ocean file. Default is prog__*''')

  parser.add_argument('-year_start', type=int, default=80,
      help='''Start year to compute averages. Default is 80.''')

  parser.add_argument('-year_end', type=int, default=100,
      help='''End year to compute averages. Default is 100.''')

  optCmdLineArgs = parser.parse_args()
  global case_name
  case_name = optCmdLineArgs.case_name
  driver(optCmdLineArgs)

#-- This is where all the action happends, i.e., functions for each diagnostic are called.

def driver(args):
  os.system('mkdir PNG')
  # mom6 grid
  grd = MOM6grid()

  # daily surface data
  # snapshots and domain-averaged plots
  latlon_plot(args,args.surface_file,grd,['SSS','SST','MLD_003'])

  # FIXME: SSU and SSV need to be plotted on u and v points, instead of tracer points
  mean_latlon_plot(args,args.surface_file,grd,['SSH','SSS','SST','MLD_003','SSU','SSV'])

  # surface fluxes
  mean_latlon_plot(args,args.forcing_file,grd,['hfds','PRCmE','taux','tauy'])

  return

# -- time-mean latlon plot
def latlon_plot(args, ncfile, grd, variables):
  nc = MFDataset(ncfile)
  time = nc.variables['time'][:]
  tm = len(time)
  for var in nc.variables:
    if var in variables:
      print("### About to plot variable {} ### \n".format(var))
      # TODO, input clim via yaml,json
      if var == 'SSH':
        clim=[-2,2]
      elif var == 'SSS':
        clim=[30.75,38.0]
      elif var == 'SST':
        clim = [-1,30]
      elif var == 'MLD_003':
        clim = [0,2000]

      dtime=[]; dmin=[]; dmax=[]; dmean=[]; dstd=[]; drms=[]
      plot_stats=True
      for t in range(tm):
        filename = str('PNG/%s_%05d.png' % (var,t))
        if os.path.isfile(filename):
          print("File {} already exists! Moving to the next one...\n".format(filename))
          plot_stats=False
        else:
          print ("time index {} of {}".format(t, tm))
          data = nc.variables[var][t,:]
          units = nc.variables[var].units
          #TODO: convert days to date
          m6plot.xyplot( data , grd.geolon, grd.geolat, area=grd.Ah,
            suptitle=case_name,
            title=r'%s, [%s] - Day: %5.1f'%(var,units,time[t]),
            extend='both',
            clim=clim,
            save=filename)

          # get stats
          sMin, sMax, mean, std, rms = m6plot.myStats(data, grd.Ah)
          dmin.append(sMin); dmax.append(sMax); dmean.append(mean)
          dstd.append(std); drms.append(rms); dtime.append(time[t])

      if plot_stats:
        plt.figure()
        f, ax = plt.subplots(5, sharex=True)
        ax[0].plot(dtime, dmean)
        ax[0].set_title(r'%s, [%s]'%(var,units))
        ax[0].set_ylabel('Mean')
        ax[1].plot(dtime, dmin)
        ax[1].set_ylabel('Min.')
        ax[2].plot(dtime, dmax)
        ax[2].set_ylabel('Max.')
        ax[3].plot(dtime, dstd)
        ax[3].set_ylabel('Std')
        ax[4].plot(dtime, drms)
        ax[4].set_ylabel('Rms')
        ax[4].set_xlabel('Time [days]')
        plt.savefig('PNG/%s_stats.png'%(var))
        plt.close()

  nc.close()

  return

# -- time-mean latlon plot
def mean_latlon_plot(args, ncfile, grd, variables):
  nc = MFDataset(ncfile)
  time = nc.variables['time'][:]
  ti = np.nonzero(time<=args.year_start*365)[0][-1]
  tf = np.nonzero(time>=args.year_end*365)[0][0]
  print 'ti, tf', ti, tf
  for var in nc.variables:
    if var in variables:
      filename = str('PNG/%s.png' % (var))
      if os.path.isfile(filename):
        print (' \n' + '==> ' + 'File has been saved, moving to the next one ...\n' + '')
      else:
        print("About to plot surface time-average for variable {}... \n".format(var))
        data = nc.variables[var][ti:tf,:].mean(axis=0)
        units = nc.variables[var].units
        #long_name = nc.variables[var].long_name
        m6plot.xyplot( data , grd.geolon, grd.geolat, area=grd.Ah,
          suptitle=case_name,
          title=r'%s, [%s] averaged over years %i-%i'%(var,units,args.year_start,args.year_end),
          extend='both',
          save=filename)

  nc.close()
  return


# -- process log files and plot the data
def get_log_data(args):


  return

# -- return MOM6 grid object
def MOM6grid(grd_file='/glade/scratch/gmarques/g.c2b6.GNYF.T62_t061.control.001/run/ocean_geometry.nc'):
    """
    grd = MOM6grid()
    """

    # create an empty class object
    class MOM6_grd:
        pass

    # open grid file
    nc = Dataset(grd_file)

    # fill grid object
    for var in nc.variables:
       dummy = str("MOM6_grd.%s = nc.variables[var][:]"% (var))
       exec(dummy)

    # close netcdf file
    nc.close()
    print('MOM6 grid successfully loaded... \n')
    return MOM6_grd

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()



