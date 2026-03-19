#!/usr/bin/env python3
#
# Original code by Gijsbert Tilstra (KNMI) and Maarten Sneep. 
# Adapted by MdG for
#___  ____/_  __ \    ___    |__  __/_________(_)___________ _
#__  __/  _  / / /    __  /| |_  /_ __  ___/_  /_  ___/  __ `/
#_  /___  / /_/ /     _  ___ |  __/ _  /   _  / / /__ / /_/ / 
#/_____/  \____/      /_/  |_/_/    /_/    /_/  \___/ \__,_/  

#
# @author Martin de Graaf
# @version 0.1
# @since 2022-11-11
# @copyright 2022 KNMI.
#   This software is the proprietary information of KNMI.
#   All rights Reserved
# TROPOMI NL L2 Data Processors
#
# @
"""
Created on Mar 10, 2026

@author: M. de Graaf
"""
import os, sys
import datetime
import numpy as np
import netCDF4
import scipy.interpolate
import sys

class s5p_aai_class:

    def find_var_recursively(self, grp, varname):

        """Find a variable in a netCDF file recursively in all groups.
        returns None if the variable isn't found."""

        if varname in list(grp.variables.keys()):
            return grp.variables[varname]

        grp_list = list(grp.groups.keys())
        for grp_name in grp_list:
            rval = self.find_var_recursively(grp.groups[grp_name], varname)
            if rval is not None:
                return rval

        return None

    def __init__(self,aaifile,index=None):
        
        if index is None:
            index = [-1]

        print("{0:%Y-%m-%d %H:%M:%S} Opening {1}".format(datetime.datetime.now(),os.path.basename(aaifile)))
        ref = netCDF4.Dataset(aaifile, mode = 'r')

        self.orbit_number = ref.orbit
        self.geojson = ref['METADATA/EOP_METADATA/om:featureOfInterest/eop:multiExtentOf/gml:surfaceMembers/gml:exterior/'].getncattr('gml:posList')

        #Extract the post processing offset to the AAI from the METADATA. This is the number that HAS BEEN SUBTRACTED from 
        #the AAI. I.e. the negative correction is the correction that was subtracted, raising the AAI as a result. 
        self.aai_add_offset_340380 = float(ref['METADATA/ALGORITHM_SETTINGS/'].getncattr('algo.pair_1.aai_add_offset'))
        self.aai_add_offset_354388 = float(ref['METADATA/ALGORITHM_SETTINGS/'].getncattr('algo.pair_2.aai_add_offset'))

        field_name = "PRODUCT"

        scanline=self.find_var_recursively(ref[field_name], "scanline")[:]

        if index is None:
            index = scanline
        if index[0] == -1:
            index = scanline
        if index[0] < -1: 
            print("s5p_aai_class: Indeces have to be positive, returning")
            return
        if np.max(index) > np.max(scanline):
            print("s5p_aai_class: Index maximum must be smaller than {0}, returning".format(np.max(scanline)))
            return
        
        self.vza = self.find_var_recursively(ref[field_name], "viewing_zenith_angle")[0,index]
        self.sza = self.find_var_recursively(ref[field_name], "solar_zenith_angle")[0,index]
        self.vaa = self.find_var_recursively(ref[field_name], "viewing_azimuth_angle")[0,index]
        self.saa = self.find_var_recursively(ref[field_name], "solar_azimuth_angle")[0,index]


        # compute derived geometry params
        self.mu0 = np.cos( np.deg2rad(self.sza))
        self.mu = np.cos( np.deg2rad(self.vza))
        phi = np.abs(self.saa - self.vaa)
        self.phi =np.where(phi>180., 360.-phi, phi)
        self.phi_rt = 180. -  self.phi             # DAK compatible"

        self.aai340380   =self.find_var_recursively(ref[field_name], "aerosol_index_340_380")[0,index]
        self.aai354388   =self.find_var_recursively(ref[field_name], "aerosol_index_354_388")[0,index]
        
        self.latitude    =self.find_var_recursively(ref[field_name], "latitude")[0,index] 
        self.longitude   =self.find_var_recursively(ref[field_name], "longitude")[0,index] 
        
        # compute the epoch for converting the tai93 time
        reference_time = self.find_var_recursively(ref[field_name], "time")[0]
        self.abs_time = datetime.datetime(2010,1,1,0,0,0) + datetime.timedelta(seconds=float(reference_time))
        #
        self.delta_time = self.find_var_recursively(ref[field_name], "delta_time")[0]
        time = []
        for i in range(self.delta_time.shape[0]):
            time.append(self.abs_time + datetime.timedelta(seconds=float(self.delta_time[i])/1000.) )
        time = np.array(time)
        self.time = time[index]
        #

        field_name = "PRODUCT/SUPPORT_DATA/GEOLOCATIONS"
        self.corner_latitude  = self.find_var_recursively(ref[field_name], "latitude_bounds" )[0,index]
        self.corner_longitude = self.find_var_recursively(ref[field_name], "longitude_bounds")[0,index]

        self.vza = self.find_var_recursively(ref[field_name], "viewing_zenith_angle")[0,index]
        self.sza = self.find_var_recursively(ref[field_name], "solar_zenith_angle")[0,index]
        self.vaa = self.find_var_recursively(ref[field_name], "viewing_azimuth_angle")[0,index]
        self.saa = self.find_var_recursively(ref[field_name], "solar_azimuth_angle")[0,index]

        
        return 

if __name__ == '__main__':
    
    main()
