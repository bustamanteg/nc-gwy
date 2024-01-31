#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 14:22:35 2022

@author: lfairgrievepark12

This function is used to convert .nc file types from GXSM into gwyddion readable files
with the useful information included in the Metadata.

This is set up to provide Topography, Freq shift = ADC1, Dissipation = ADC2,
cantilever amplitude = ADC4, and PID Phase = ADC6. It should be configureable though
if you're files contain different data

### NOTE ### The y coordinates gwyddion vs GXSM read files will be flipped along the
center line /xaxis. This is because each system uses a different coordinate system. Be
mindful - if you move to the gwyddion coordinates on GXSM you will go to the wrong spot
"""

from netCDF4 import Dataset
from gwyfile.objects import GwyContainer, GwyDataField
import numpy as np
import fire



def run(direct,prefix, cantsens, signalgain):
    '''
    direct = File directory
    prefix = File prefix (everything that comes before -M or -xp etc)
    cant sens = cantilever sensitivity in V/nm
    signal gain = signal amp used on preamp
    '''
    #direct = '/Users/lfairgrievepark12/20220327/'
    #prefix = 'FerroceneA_LH_6_scan041'
    # Make sure sensitvity is correct
    #cantsens = 0.120*5
    fileprefix = direct+prefix
    
    def datapull(filename, scalefactor):
        # This function pulls the xyz data from the nc file
        # Floatfield is the actual x,y,z data, weird indexingneeded to make it an nxn matrix
        # Multiplied by dz to scale correctly
        # data is in float32 as a masked array, need to turn into float64 and regular np array to read it
        # scale factor depends on data type
        file  = Dataset(filename, mode='r')
        dat = np.float64(np.array(file.variables['FloatField'][::10,::10][0][0]*file.variables['dz'][:]))*scalefactor
        return dat
    
    def rangepull(filename):
        # This function pulls the correct range and offset for the image
        # Should only be neccessary to do on one file in a set because all should have the same range
        file  = Dataset(filename, mode='r')
        xrange =  file.variables['rangex'][:]*1e-10
        yrange =  file.variables['rangey'][:]*1e-10
        xoff = file.variables['offsetx'][:]*1e-10-xrange/2
        yoff = file.variables['offsety'][:]*1e-10-xrange/2
        return (xrange, yrange, xoff, yoff)
    
    def metapull(filename):
        # This pulls the relevant metadata from the file
        # Again, should only be needed once because all files in set have same metadata
        file  = Dataset(filename, mode='r')
        
        # Date of scan and comment or byte encoded and in an array so this is necessary
        dateofscan=''
        for i in file.variables['dateofscan'][:-1]:
            dateofscan+=i.decode('ascii')
        
        comment=''
        for i in file.variables['comment'][:-1]:
            comment+=i.decode('ascii')
        
        meta ={
            'Date of Scan':dateofscan,
            'Comment':comment,
            'Bias':str(file.variables['sranger_mk2_hwi_bias'][:])+'V',    
            'CP':str(file.variables['sranger_mk2_hwi_z_servo_CP'][:]),
            'CI':str(file.variables['sranger_mk2_hwi_z_servo_CI'][:]),
            'Scan Speed':str(file.variables['sranger_mk2_hwi_scan_speed_x'][:]*.1)
           }
        for key, value in meta.items():
            print(key, ' : ', value)
        
        #metadata needs to be in the form of a GwyContainer
        meta = GwyContainer(meta)
        return meta
        
    
    Topo_r = datapull(fileprefix+'-M-Xp-Topo.nc',1e-10) #ang to m scale factor
    Topo_l = datapull(fileprefix+'-Xm-Topo.nc',1e-10)
    
    Freq_r = datapull(fileprefix+'-Xp-ADC1.nc',-1/0.02) # Labone gain scaling
    Freq_l = datapull(fileprefix+'-Xm-ADC1.nc',-1/0.02)
    
    Diss_r = datapull(fileprefix+'-Xp-ADC2.nc',1) #frac difference scaling
    Diss_r = (Diss_r-Diss_r.min())/Diss_r.min()
    Diss_l = datapull(fileprefix+'-Xm-ADC2.nc',1)
    Diss_r = (Diss_l-Diss_l.min())/Diss_l.min()
    
    Amp_r = datapull(fileprefix+'-Xp-ADC4.nc',1e-10/cantsens/signalgain) #V to nm to m, potential amplifier, labone gain of 110
    Amp_l = datapull(fileprefix+'-Xm-ADC4.nc',1e-10/cantsens/signalgain)
    
    Phase_r = datapull(fileprefix+'-Xp-ADC6.nc',1/0.01) # Labone scalign
    Phase_l = datapull(fileprefix+'-Xm-ADC6.nc',1/0.01)
    
    (xrange, yrange, xoff, yoff) = rangepull(fileprefix+'-M-Xp-Topo.nc')
    
    meta=  metapull(fileprefix+'-M-Xp-Topo.nc')
    
    # This is how the gwyddion file is actually built out
    obj = GwyContainer()
    obj['/0/data/title'] = 'Topography ->'
    obj['/0/data'] = GwyDataField(Topo_r, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='m')
    obj['/1/data/title'] = 'Topography <-'
    obj['/1/data'] = GwyDataField(Topo_l, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='m')
    
    obj['/2/data/title'] = 'Frequency Shift ->'
    obj['/2/data'] = GwyDataField(Freq_r, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='Hz')
    obj['/3/data/title'] = 'Frequency Shift <-'
    obj['/3/data'] = GwyDataField(Freq_l, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='Hz')
    
    obj['/4/data/title'] = 'Dissipation ->'
    obj['/4/data'] = GwyDataField(Diss_r, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='Hz')
    obj['/5/data/title'] = 'Dissipation <-'
    obj['/5/data'] = GwyDataField(Diss_l, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='Hz')
    
    obj['/6/data/title'] = 'Amplitude ->'
    obj['/6/data'] = GwyDataField(Amp_r, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='m')
    obj['/7/data/title'] = 'Amplitude <-'
    obj['/7/data'] = GwyDataField(Amp_l, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='m')
    
    obj['/8/data/title'] = 'PID Phase ->'
    obj['/8/data'] = GwyDataField(Phase_r, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='deg')
    obj['/9/data/title'] = 'PID Phase <-'
    obj['/9/data'] = GwyDataField(Phase_l, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='deg')
    
    # Populating metadat field
    for i in range(0,10):
        obj['/'+str(i)+'/meta'] = meta
        
    # Making trace for each channel visible on file load
    for i in [0,2,4,6,8]:
        obj['/'+str(i)+'/data/visible'] = True
    
    
    # These are things that can't be included in regular metadata (maybe if they are externally supplied to GXSM
    # they could be) These factors should be included as a "Comment" in the metadata
    '''
    # Not included
    res freq
    q factor
    Sens
    Lab one gains
    Freq shift
    '''
    
    obj.tofile(fileprefix+'.gwy')

if __name__ == '__main__':
    fire.Fire(run)