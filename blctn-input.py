#!/usr/bin/env python

'''
description:    create sea surface temperatures and sea ice concentrations
                input files for Blue Action experiments 1 and 2
license:        APACHE 2.0
author:         Ronald van Haren, NLeSC (r.vanharen@esciencecenter.nl)
'''

from netCDF4 import Dataset
from netCDF4 import date2num as nc_date2num
from netCDF4 import num2date as nc_num2date
import numpy as np
import os
import time
import sys
import subprocess
import argparse


def main(startyear, endyear, overwrite):
    '''
    Main function, call experiment 1 and experiment 2 functions

    Args:
        startyear (int):        First year to use in calculations
        endyear (int):          Stop calculations in the beginning of this year
        overwrite (bool):       Force overwriting exsisting files if True
    '''
    exp1(startyear, endyear, overwrite)
    exp2(startyear, endyear, overwrite)
    exp3_4(startyear, endyear, 3, overwrite)
    exp3_4(startyear, endyear, 4, overwrite)


def create_dimensions_netcdf(ncfile, dtobj, lat, lon):
    '''
    Define dimensions in netCDF file

    Args:
        ncfile:     netCDF file handle
        dtobj:      python datetime object
        lat:        latitude [degrees]
        lon:        longitude [degrees]

    Returns:
        ncfile:     netCDF file handle
    '''
    # description of the file
    ncfile.description = 'Blue Action'
    ncfile.history = 'Created ' + time.ctime(time.time())
    # create time dimension
    timevar = ncfile.createDimension('time', None)
    # netcdf time variable UTC
    timevar = ncfile.createVariable('time', 'float32', ('time',),
                                    zlib=True)
    timevar.units = 'days since 1850-01-01 00:00:00'
    timevar.calendar = 'gregorian'
    timevar.standard_name = 'time'
    timevar.long_name = 'Time'
    timevar.axis = 'T'
    # convert dtobj to num
    timevar[:] = nc_date2num(dtobj,
                             units=ncfile['time'].units,
                             calendar=ncfile['time'].calendar)
    # write lon/lat variables
    lonvar = ncfile.createDimension('longitude', len(lon))
    lonvar = ncfile.createVariable('longitude', 'float32', ('longitude',))
    lonvar.units = 'degrees_east'
    lonvar.axis = 'X'
    lonvar.standard_name = 'longitude'
    lonvar.long_name = 'Longitude'
    lonvar[:] = lon
    latvar = ncfile.createDimension('latitude', len(lat))
    latvar = ncfile.createVariable('latitude', 'float32', ('latitude',))
    latvar.units = 'degrees_north'
    latvar.axis = 'Y'
    latvar.standard_name = 'latitude'
    latvar.long_name = 'Latitude'
    latvar[:] = lat
    return ncfile


def create_variable_sic(ncfile, sic):
    '''
    Add sea ice concentraction[time, lat, lon] to netCDF file

    Args:
        ncfile:     netCDF file handle
        sic:        sea ice concentration in [frac]

    Returns:
        ncfile:     netCDF file handle
    '''
    sicvar = ncfile.createVariable('sic', 'float32',
                                   ('time', 'latitude', 'longitude',),
                                   zlib=True, fill_value=-1e30)
    sicvar.standar_name = 'sea_ice_area_fration'
    sicvar.long_name = 'SIC'
    sicvar.units = '1'
    sicvar.cell_methods = 'time: lat: lon: mean'
    sicvar[:] = sic
    return ncfile


def create_variable_sst(ncfile, sst):
    '''
    Add sea surface temperatures[time, lat, lon] to netCDF file

    Args:
        ncfile:     netCDF file handle
        sst:        sea surface temperatures in [K]

    Returns:
        ncfile:     netCDF file handle
    '''
    sstvar = ncfile.createVariable('sst', 'float32',
                                   ('time', 'latitude', 'longitude',),
                                   fill_value=-1e30, zlib=True)
    sstvar.standar_name = 'sea_surface_temperature'
    sstvar.long_name = 'SST'
    sstvar.units = 'K'
    sstvar.cell_methods = 'time: lat: lon: mean'
    sstvar[:] = sst
    return ncfile


def sst_sic_adjustment(sst, siconc, units_sst='K', units_sic='frac'):
    '''
    Perform SST and SIC adjustments as per instructions:
        1. Set minimum SST to -1.8 degC
        2. Set SST to -1.8 degC if SIC>0.9
        3. If SST>5 degC, set SIC to 0
        4. If SIC<0.9, we calculate SSTmax, where
           SSTmax=9.328(*0.729-SIC^3)-1.8. If SST>SSTmax,
           reduce SIC so that SST=SSTmax.

    Args:
        sst:        sea surface temperatures in [degC] of [K]
        siconc:     sea ice concentration in [frac] of [perc]
        units_sst:  units of sst input, [degC] or [K]
        units_sic:  units of siconc input [frac] of [perc]

    Returns:
        sst:        adjusted sea surface temperatures in [K]
        siconc:     adjusted sea ice concentration in [frac]
    '''
    # get sst mask
    sst_mask = np.ma.getmask(sst)
    # convert units for calculation
    if (units_sst == 'K'):
        # convert to degC for calculations
        sst = sst - 273.15
    elif (units_sst == 'C'):
        pass
    else:
        print('Unknown units for SST: ' + str(units_sst))
        sys.exit()
    if (units_sic == 'perc'):
        # convert to fraction
        siconc = siconc/100.
    elif (units_sic == 'frac'):
        pass
    else:
        print('Unknown units for SIC: ' + str(units_sic))
        sys.exit()
    # set minimum SST to -1.8 degC
    np.ma.MaskedArray.clip(sst, -1.8, None, out=sst)
    # set maximum SIC to 0.999478 (IFS segfaults when SIC=1)
    np.ma.MaskedArray.clip(siconc, None, 0.999478, out=siconc)
    # unshare mask
    sst.unshare_mask()
    siconc.unshare_mask()
    # set SST to -1.8degC if sic>0.9
    sst[siconc > 0.9] = -1.8
    # If SST>5 degC, set SIC to 0
    siconc[sst > 5] = 0
    # If SIC<90%, we calculate SSTmax (SSTmax=9.328*(0.729-SIC^3)-1.8)
    sstmax = np.zeros(np.shape(siconc))
    idx_sic = (siconc < 0.9) & (siconc > 0)
    sstmax[idx_sic] = (9.328 * (0.729-(siconc[idx_sic])**3)-1.8)
    # If SST > SSTmax,reduce the SIC, so that SST=SSTmax
    idx = (sst > sstmax) & (siconc < 0.9) & (siconc > 0)
    siconc[idx] = (0.729 - ((sst[idx] + 1.8)/9.328))**(1/3)
    # set maximum SIC to 0.999478 (IFS segfaults when SIC=1)
    np.ma.MaskedArray.clip(siconc, None, 0.999478, out=siconc)
    # convert SST to K
    sst = sst + 273.15
    # reapply original sst mask
    sst = np.ma.masked_where(sst_mask, sst)
    return sst, siconc


def exp1(startyear, endyear, overwrite=False):
    '''
    create input files for Blue Action experiment 1

    Args:
        startyear (int):      First year to use in calculations
        endyear (int):        Stop calculations in the beginning of this year
        overwrite (bool):     Force overwriting existing files
    '''
    # basestring input files
    siconc_bs = ("siconc_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-" +
                 "HadISST-2-2-0-0-0_gn_")
    tos_bs = ("tos_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-" +
              "HadISST-2-2-0-0-0_gn_")
    # basestring output files
    tos_bsout = 'HadISST2_prelim_0to360_alldays_sst_'
    siconc_bsout = 'HadISST2_prelim_0to360_alldays_sic_'
    # create output directory if needed
    if not os.path.exists('exp1'):
        os.makedirs('exp1')
    # loop over all years
    for yr in range(int(startyear), int(endyear)):
        # input filenames
        timestr = str(yr) + '0101-' + str(yr) + '1231'
        filename_sic_in = os.path.join('siconc', siconc_bs + timestr + '.nc')
        filename_tos_in = os.path.join('tos', tos_bs + timestr + '.nc')
        # output filenames
        filename_sic_out = os.path.join('exp1', siconc_bsout + str(yr) + '.nc')
        filename_tos_out = os.path.join('exp1', tos_bsout + str(yr) + '.nc')
        # check if existing file can be used
        if (all(os.path.isfile(fl) for fl in
                [filename_sic_out, filename_tos_out]) and not overwrite):
            print('Keeping existing files: ' + filename_sic_out + ' and ' +
                  filename_tos_out)
            continue
        # open output netCDF files
        ncfile_tos_in = Dataset(filename_tos_in, 'r')
        ncfile_sic_in = Dataset(filename_sic_in, 'r')
        # open output netCDF files
        ncfile_tos = Dataset(filename_tos_out, 'w')
        ncfile_sic = Dataset(filename_sic_out, 'w')
        # convert time to datetimeobject
        dtobj = nc_num2date(ncfile_tos_in['time'][:],
                            units=ncfile_tos_in['time'].units,
                            calendar=ncfile_tos_in['time'].calendar)
        sst = ncfile_tos_in.variables['tos'][:]
        siconc = ncfile_sic_in.variables['siconc'][:]
        lat = ncfile_tos_in.variables['latitude'][:]
        lon = ncfile_tos_in.variables['longitude'][:]
        # adjust sst and sic per instructions
        sst, siconc = sst_sic_adjustment(sst, siconc,
                                         units_sst='C', units_sic='perc')
        # write sst
        ncfile_tos = create_dimensions_netcdf(ncfile_tos, dtobj, lat, lon)
        ncfile_tos = create_variable_sst(ncfile_tos, sst)
        # write sic
        ncfile_sic = create_dimensions_netcdf(ncfile_sic, dtobj, lat, lon)
        ncfile_sic = create_variable_sic(ncfile_sic, siconc)
        # close netCDF files
        ncfile_tos.close()
        ncfile_sic.close()
        ncfile_tos_in.close()
        ncfile_sic_in.close()

def exp3_4(startyear, endyear, expno, overwrite=False):
    '''
    create input files for Blue Action experiment 3 and/or 4

    Args:
        startyear (int):      First year to use in calculations
        endyear (int):        Stop calculations in the beginning of this year
        expno (int):          Experiment number (3 or 4)
        overwrite (bool):     Force overwriting existing files
    '''
    if expno not in [3, 4]:
        print('experiment number should be 3 or 4, returning...')
        return
    # basestring input files
    siconc_bs = ("siconc_input4MIPs_SSTsAndSeaIce_HighResMIP_MOHC-" +
                 "HadISST-2-2-0-0-0_gn_")
    tos_bs = ("tos_EXP" + str(expno) + "-HadISST-2-2-0-0-0_gn_")
    # basestring output files
    tos_bsout = 'HadISST2_prelim_0to360_alldays_sst_'
    siconc_bsout = 'HadISST2_prelim_0to360_alldays_sic_'
    # create output directory if needed
    if not os.path.exists('exp' + str(expno)):
        os.makedirs('exp' + str(expno))
    # loop over all years
    for yr in range(int(startyear), int(endyear)):
        # input filenames
        timestr = str(yr) + '0101-' + str(yr) + '1231'
        filename_sic_in = os.path.join('siconc', siconc_bs + timestr + '.nc')
        filename_tos_in = os.path.join('input', 'EXP' + str(expno), tos_bs + timestr + '.nc')
        # output filenames
        filename_sic_out = os.path.join('exp' + str(expno), siconc_bsout + str(yr) + '.nc')
        filename_tos_out = os.path.join('exp' + str(expno), tos_bsout + str(yr) + '.nc')
        # check if existing file can be used
        if (all(os.path.isfile(fl) for fl in
                [filename_sic_out, filename_tos_out]) and not overwrite):
            print('Keeping existing files: ' + filename_sic_out + ' and ' +
                  filename_tos_out)
            continue
        # open output netCDF files
        ncfile_tos_in = Dataset(filename_tos_in, 'r')
        ncfile_sic_in = Dataset(filename_sic_in, 'r')
        # open output netCDF files
        ncfile_tos = Dataset(filename_tos_out, 'w')
        ncfile_sic = Dataset(filename_sic_out, 'w')
        # convert time to datetimeobject
        dtobj = nc_num2date(ncfile_tos_in['time'][:],
                            units=ncfile_tos_in['time'].units,
                            calendar=ncfile_tos_in['time'].calendar)
        sst = ncfile_tos_in.variables['tos'][:]
        siconc = ncfile_sic_in.variables['siconc'][:]
        lat = ncfile_tos_in.variables['latitude'][:]
        lon = ncfile_tos_in.variables['longitude'][:]
        # adjust sst and sic per instructions
        sst, siconc = sst_sic_adjustment(sst, siconc,
                                         units_sst='C', units_sic='perc')
        # write sst
        ncfile_tos = create_dimensions_netcdf(ncfile_tos, dtobj, lat, lon)
        ncfile_tos = create_variable_sst(ncfile_tos, sst)
        # write sic
        ncfile_sic = create_dimensions_netcdf(ncfile_sic, dtobj, lat, lon)
        ncfile_sic = create_variable_sic(ncfile_sic, siconc)
        # close netCDF files
        ncfile_tos.close()
        ncfile_sic.close()
        ncfile_tos_in.close()
        ncfile_sic_in.close()


def exp2_climate(startyear, endyear, filename_sic_climate,
                 filename_sst_climate):
    '''
    Calculate daily climate for SST/SIC for Blue Action experiment 2

    Args:
        startyear (int):        First year to use in calculations
        endyear (int):          Stop calculations in the beginning of this year
        filename_sic_climate:   Path sea ice concentration climate netCDF file
        filename_sst_climate:   Path sea surface temperatures climate netCDF
                                file
    '''
    # basestring files
    sst_bs = 'HadISST2_prelim_0to360_alldays_sst_'
    siconc_bs = 'HadISST2_prelim_0to360_alldays_sic_'
    # calculate sea ice daily climate
    fn_sic = [os.path.join('exp1', siconc_bs + str(yr) + '.nc') for yr in
              range(int(startyear), int(endyear))]
    fn_sic_str = ' '.join(map(str, fn_sic))
    tmpfile = os.path.join('exp2', 'tmpfile.nc')
    command = 'cdo cat ' + fn_sic_str + ' ' + tmpfile
    subprocess.check_call(command, shell=True)
    command = 'cdo ydaymean ' + tmpfile + ' ' + filename_sic_climate
    subprocess.check_call(command, shell=True)
    os.remove(tmpfile)
    # calculate sst daily climate
    fn_sst = [os.path.join('exp1', sst_bs + str(yr) + '.nc') for yr in
              range(int(startyear), int(endyear))]
    fn_sst_str = ' '.join(map(str, fn_sst))
    command = 'cdo cat ' + fn_sst_str + ' ' + tmpfile
    subprocess.check_call(command, shell=True)
    command = 'cdo ydaymean ' + tmpfile + ' ' + filename_sst_climate
    subprocess.check_call(command, shell=True)
    os.remove(tmpfile)


def exp2(startyear, endyear, overwrite=False):
    '''
    Create input files for Blue Action experiment 2

    Args:
        startyear (int):      First year to use in calculations
        endyear (int):        Stop calculations in the beginning of this year
        overwrite (bool):     Force overwriting existing files
    '''
    # basestring files
    sst_bs = 'HadISST2_prelim_0to360_alldays_sst_'
    siconc_bs = 'HadISST2_prelim_0to360_alldays_sic_'
    # input filenames climate
    filename_sic_climate = os.path.join('exp2',
                                        'sic_climate_' + str(startyear) +
                                        '_' + str(endyear) + '.nc')
    filename_sst_climate = os.path.join('exp2',
                                        'sst_climate_' + str(startyear) +
                                        '_' + str(endyear) + '.nc')
    # create output directory if needed
    if not os.path.exists('exp2'):
        os.makedirs('exp2')
    # check if existing file can be used
    if (all(os.path.isfile(fl) for fl in
            [filename_sic_climate, filename_sst_climate]) and not overwrite):
        print('Keeping existing files: ' + filename_sic_climate +
              ' and ' + filename_sst_climate)
    else:
        exp2_climate(startyear, endyear,
                     filename_sic_climate, filename_sst_climate)
    # open sst and sic daily climate netCDF files (r)
    ncfile_sic_climate = Dataset(filename_sic_climate, 'r')
    ncfile_sst_climate = Dataset(filename_sst_climate, 'r')
    sst_clim = ncfile_sst_climate.variables['sst'][:]
    sic_clim = ncfile_sic_climate.variables['sic'][:]
    # loop over all years of exp 1 files
    for idx, yr in enumerate(range(int(startyear), int(endyear))):
        # input filenames
        filename_sic_in = os.path.join('exp1', siconc_bs + str(yr) + '.nc')
        filename_sst_in = os.path.join('exp1', sst_bs + str(yr) + '.nc')
        # output filenames
        filename_sic_out = os.path.join('exp2', siconc_bs + str(yr) + '.nc')
        filename_sst_out = os.path.join('exp2', sst_bs + str(yr) + '.nc')
        # check if existing file can be used
        if ((os.path.isfile(filename_sic_out) and
             (os.path.isfile(filename_sst_out)) and
             not overwrite)):
            print('Keeping existing files: ' + filename_sic_out +
                  ' and ' + filename_sst_out)
            continue
        # open input files (r)
        ncfile_sic_in = Dataset(filename_sic_in, 'r')
        ncfile_sst_in = Dataset(filename_sst_in, 'r')
        # open output files (rw)
        ncfile_sic_out = Dataset(filename_sic_out, 'w')
        ncfile_sst_out = Dataset(filename_sst_out, 'w')
        # convert time to datetimeobject
        dtobj = nc_num2date(ncfile_sst_in['time'][:],
                            units=ncfile_sst_in['time'].units,
                            calendar=ncfile_sst_in['time'].calendar)
        sst = ncfile_sst_in.variables['sst'][:]
        lat = ncfile_sst_in.variables['latitude'][:]
        lon = ncfile_sst_in.variables['longitude'][:]
        # get sst mask
        sst_mask = np.ma.getmask(sst)
        # use sic_clim for sic and sst_clim where sic>0
        sst = np.where(sic_clim[0:np.shape(sst)[0], :] > 0,
                       sst_clim[0:np.shape(sst)[0], :], sst)
        # reapply original sst mask
        sst = np.ma.masked_where(sst_mask, sst)
        sic = sic_clim[0:np.shape(sst)[0], :]  # make sure ndays are equal
        # redo the sst and sic adjustments as per instruction
        sst, sic = sst_sic_adjustment(sst, sic,
                                      units_sst='K', units_sic='frac')
        # write sst
        ncfile_sst_out = create_dimensions_netcdf(ncfile_sst_out, dtobj,
                                                  lat, lon)
        ncfile_sst_out = create_variable_sst(ncfile_sst_out, sst)
        # write sic
        ncfile_sic_out = create_dimensions_netcdf(ncfile_sic_out, dtobj,
                                                  lat, lon)
        ncfile_sic_out = create_variable_sic(ncfile_sic_out, sic)
        # close netCDF files
        ncfile_sic_out.close()
        ncfile_sst_out.close()
        ncfile_sic_in.close()
        ncfile_sst_in.close()
    # close daily climate netCDF files
    ncfile_sic_climate.close()
    ncfile_sst_climate.close()


if __name__ == "__main__":
    # define command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--startyear', type=int, default=1979,
                        help="startyear to use for processing files [default 1979]")
    parser.add_argument('-e', '--endyear', type=int, default=2016,
                        help="stop processing in the beginning of endyear [default 2016]")
    parser.add_argument('-f', '--force', action='store_true',
                        help="force overwriting existing files")
    # get arguments
    args = parser.parse_args()
    # call main()
    main(args.startyear, args.endyear, args.force)
