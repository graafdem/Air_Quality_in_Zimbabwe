#!/usr/bin/env python3
#
# Original code by Gijsbert Tilstra (KNMI)
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

import os,sys
import datetime
import numpy as np
import netCDF4
"""
s5p_l1b_class

perband

Read L1b files to return the UV reflectance. Input are the s5p solar irradiance
l1b file for a day, and two radiance files for that day, containing the radiance
for band 3 and band 4 respectively. 
E.g. the call 

refl = s5p_l1b_class(irr_filename,rad1_filename,rad2_filename,index)

will return the reflectance and relevant UV L1B input in refl. A non-exhaustive list is

refl.reflectance
relf.irradiance
refl.radiance
refl.wavelength
refl.mu0
refl.sza
refl.longitude
refl.latitude

A wavelength array is returned for each footprint, which is the solar irradiance
nominal wavelength grid, on which the radiance is interpolated, the reflectance 
is calculated, and copied for each scanline. 

"""

class s5p_l1b_class:

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

    def read_rad(self,radfile,index,read=None):

        """
        read_tropomi_radiance_fast.py

        Opens and reads the requested TROPOMI Earth radiance NetCDF file 
        (radfile) and returns the contents. Index can be used to filter the
        number of pixels (currently scanlines) that are being read. (should
        be changed to scanlines x ground_pixels)

        Longitudes and latitudes are ordered in the following way:

           4----3
           |    |   (near the equator, for the sunlit orbit part)
           1----2

        The Earth radiance is only equipped with "nominal wavelength". The solar 
        irradiance is *also* equipped with a field "calibrated wavelength".

        This will not work because the "nominal wavelength" and "calibrated 
        wavelength" grids are not very compatible. It is in fact better to use 
        "nominal wavelength" for both. The "nominal wavelength" of the solar 
        irradiance is more or less the same as that of the Earth radiance 
        but a correction for Doppler shift has been applied here.

        A correction for differences in the Earth-Sun distance (esd) of the 
        radiance and solar irradiance measurements is *NOT* performed: These 
        were already brought back to 1.0 a.u.

        """

#         print(radfile)
        print("{0:%Y-%m-%d %H:%M:%S} Opening {1}".format(datetime.datetime.now(),os.path.basename(radfile)))
        # Read the Earth radiance and other variables from the radiance file.

        ref = netCDF4.Dataset(radfile, mode = 'r')
        #
        #RETRIEVE the band number from the "summary" topgroup (/) metadata field
        field_name = 'BAND{}_RADIANCE'.format(ref.summary[22])
        band=ref.summary[22]

        #READ the dimensions
        ground_pixel=self.find_var_recursively(ref[field_name], "ground_pixel")[:]
        scanline=self.find_var_recursively(ref[field_name], "scanline")[:]
        spectral_channel=self.find_var_recursively(ref[field_name], "spectral_channel")[:]
        
        if index is None:
            index = scanline
        if index[0] == -1:
            index = scanline
        if index[0] < -1: 
            print("Indeces have to be positive, returning")
            return
        if np.max(index) > np.max(scanline):
            print("Index maximum must be smaller than {0}".format(np.max(scanline)))
            return

        # compute the epoch for converting the tai93 time
        reference_time = self.find_var_recursively(ref[field_name], "time")[0]
        abs_time = datetime.datetime(2010,1,1,0,0,0) + datetime.timedelta(seconds=float(reference_time))
        #
        delta_time = self.find_var_recursively(ref[field_name], "delta_time")[0,index]
        time = []
        for i in range(delta_time.shape[0]):
            time.append(abs_time + datetime.timedelta(seconds=float(delta_time[i])/1000.) )
        time = np.array(time)
        self.reference_time = reference_time
        self.delta_time = delta_time
        self.abs_time = abs_time
        self.time=time

        #READ the radiance (the entire orbit if index was not given, or for the 
        #indices in index only 
        if read is not None:
            rad=self.find_var_recursively(ref[field_name], "radiance")[0,index]
#            print(np.shape(rad))

            #READ the wavelength
            wvl=self.find_var_recursively(ref[field_name], "nominal_wavelength")[0]
        else: 
            rad = 0
            wvl = 0    

        lat = self.find_var_recursively(ref[field_name], "latitude" )[0,index]
        lon = self.find_var_recursively(ref[field_name], "longitude")[0,index]

        clat = self.find_var_recursively(ref[field_name], "latitude_bounds" )[0,index]
        clon = self.find_var_recursively(ref[field_name], "longitude_bounds")[0,index]

        vza = self.find_var_recursively(ref[field_name], "viewing_zenith_angle")[0,index]
        sza = self.find_var_recursively(ref[field_name], "solar_zenith_angle")[0,index]
        vaa = self.find_var_recursively(ref[field_name], "viewing_azimuth_angle")[0,index]
        saa = self.find_var_recursively(ref[field_name], "solar_azimuth_angle")[0,index]

        raa = vaa - saa             # This definition is "180 - definition DAK"

        # compute derived geometry params
        mu0 = np.cos( np.deg2rad(sza))
        mu = np.cos( np.deg2rad(vza))
        phi = np.abs(saa - vaa)
#         print( np.shape(delta_time))
        phi =np.where(phi>180., 360.-phi, phi)
        phi_rt = 180. -  phi             # DAK compatible"

        ref.close()

        return time,ground_pixel,scanline,spectral_channel,lat,lon,clat,clon,\
               vza,sza,vaa,saa,raa,mu,mu0,phi_rt,wvl,rad


    def read_uv_irradiance(self,irrfile,band=None):

        # READ the irradiance in band 3

        # Open reference to netCDF file
        ref = netCDF4.Dataset(irrfile, "r")
        field_name = 'BAND{}_IRRADIANCE'.format(band)

        with ref:
            #                 print (field_name)
            irradiance = self.find_var_recursively(ref[field_name], "irradiance")[0,0]
#            print(np.shape(irradiance))
            solar_wavelength = self.find_var_recursively(ref[field_name], "nominal_wavelength")[0]

        return irradiance, solar_wavelength

    def __init__(self,irrfile,rad1file,read=1,index=None, 
                      wavelengths=None, wavelength_bandwidth=None,band=None):
        
        if wavelength_bandwidth is None:
            wavelength_bandwidth=1	

        if band is None:
            band=3
            
        print("{0:%Y-%m-%d %H:%M:%S} Opening {1}".format(datetime.datetime.now(),os.path.basename(irrfile)))

        if index is None:
            index = [-1]
        if index[0] < -1: 
            print("Indeces have to be positive, returning")
            return
            
#        print("{0:%Y-%m-%d %H:%M:%S} In {1}".format(datetime.datetime.now(),os.path.dirname(irrfile)))
#        print(os.path.basename(rad2file))

        if read is not None:

            #READ the UV irradiance and the radiance in band 3
            irradiance0, solar_wavelength0=self.read_uv_irradiance(irrfile,band=band)
#            print(np.shape(irradiance0))
            lambda_min = np.nanmin(solar_wavelength0)
            lambda_max = np.nanmax(solar_wavelength0)

            wav_index = np.where((wavelengths > lambda_min) & (wavelengths < lambda_max))[0]
#            print(lambda_min, lambda_max, wav_index)
            
            #
            time,ground_pixel,scanline,spectral_channel0,lat,lon,clat,clon,\
                   vza,sza,vaa,saa,raa,mu,mu0,phi_rt,wvl0,rad0 = self.read_rad(rad1file,index,read=1)

            if index[0] == -1 :
                index = scanline

            # Interpolate the radiance to the solar irradiance and calculate the Earth reflectance.
            reflect0 = np.nan	*rad0

            reflectance = np.zeros( (len(index),len(ground_pixel),len(wavelengths)) )  * np.nan
            counter = len(index)*len(ground_pixel)

            for iobs in range(counter):
                i = np.floor(iobs / len(ground_pixel))  #GET the scanline index
                i_idx = i.astype(int)
                j_gp = iobs % len(ground_pixel)         #GET the ground_pixel index
    #            print(iobs, i_idx, j_gp)

                #INTERPOLATE radiance to solar irr wavelength grid             
                radiance_new = np.interp(solar_wavelength0[j_gp,:], wvl0[j_gp,:], rad0[i_idx,j_gp,:])
                #
                #KEEP the interpolated reflectance only
                r = np.pi * radiance_new / irradiance0[j_gp,:] / mu0[i_idx,j_gp]
                #
    #
    #CREATE monochromatic reflectances bands by extracting triangular weighted reflectances 
    #around the selected wavelengths. Tent-shaped:
    #        
    #        ^
    #       / \
    #      /   \
    #_____/__|__\___
    # -dx/2  x  +dx/2  
    # 
                #GET the reflectance on the selected wavelength bands only
                for iwav,sel_wav in enumerate(wavelengths[wav_index]):
                    k_wav=wav_index[iwav]
    #                print("iwav,sel_wav,wav_index[iwav]",iwav,sel_wav,wav_index[iwav])
                    wav_range = np.where(np.abs(solar_wavelength0[j_gp,:] - sel_wav) < 0.5 * wavelength_bandwidth) 
                    dx = solar_wavelength0[j_gp,wav_range] - sel_wav        
                    
    #       IF dx LT 0 THEN wf = dx + 0.5D * wavelength_bandwidth ELSE $
    #           wf = 0.5D * wavelength_bandwidth - dx
    #       wf >= 0 
                    wf = -np.abs(dx) + 0.5*wavelength_bandwidth
                    reflectance[i_idx,j_gp,k_wav] = np.sum(wf * r[wav_range]) / np.sum(wf)
    #

        #
        #GET the dimensions of the final arrays. 
        #num_wavelengths = len(self.spectral_channel)
        num_footprints=len(index)*len(ground_pixel)

        self.latitude=lat
        self.longitude=lon

        self.corner_latitude =clat#.reshape(num_footprints,4)
        self.corner_longitude=clon#.reshape(num_footprints,4)
        
        self.sza=sza
        self.vza=vza
        self.saa=saa
        self.vaa=vaa
        self.mu0=mu0
        self.mu =mu
        self.phi=phi_rt

        self.num_scanlines=len(index)
        self.num_groundpixels=len(ground_pixel)
        #self.num_footprints=num_footprints
        #self.num_wavelengths=num_wavelengths
        self.reflectance = reflectance
       
        return

if __name__ == '__main__':

    main()
