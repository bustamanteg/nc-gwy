#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Code that converts scans in the nc format, produced by gxsm to gwy format. Each channel is saved in the same file. 
'''
This code creates the following channels: 




0 Topography (right)
1 Topography (left)
2. Frequency shift (right)
3. Frequency shift (left)
4. Excitation amplitude (right) [V]            Voltage Amplitude aplied to the excitation system to mantain a constant oscillation of the cantilever. 
5. Excitation amplitude (left)  [V]             ''
6. Amplitude (right)                            Oscillation amplitude nm  
7. Amplitude (left )                            Oscillation amplitude nm
8. 1w amplitude  (right)                  If the voltage applied is modulated with an AC voltage V(t)=Aac*sin(wt) , an amplitude A_1w appears at f0+w/2pi
9. 1w amplitude (left)                   ''
10. 2w  amplitude
11. 2w  amplitude
12. PID Phase right 
13. PID Phase (left)
14. Excess excitation (right)  [Hz] or [eV]    Either by calculating the mean excitation amplitude of the scan without outliers A0_avg=1, or manually entered in the parameter Diss0 (A0_avg=0) 
15. Excess excitation (left)    [Hz] of [eV]    ''


E0_J=np.pi*A**2*k/Q
E0_eV=E0_J/1.602e-19

'''

# The dissipation channel can be saved in units of eV/cycle or Hz. 
# Change the user parameters to suit your needs. 

#To do: 

#Make an if to only save Diss_r (the excess dissipation with respect to the same image) or Diss_r0 (the excess dissipation with respect to a fixed value)
#Output a text file that sumarises the analysis. 

#--------------------------
#USER  PARAMETERS
#edit the following for your file

direct='/data/20231218_run24/'  #Directory where the data is stored
prefix='r23_brockley_lini026'        #File prefix
cantsens=0.03                  # V/nm Conversion factor of the deflection signal
signalgain=1                    # Preamplifier signal gain.
Q=29329                         # Cantilever Q factor
f0=155022                       # Hz Resonance frequency of the cantilever
s_VHz=0.02# V/Hz                # V/Hz Conversion factor from V to Hz in the frequency shift
Diss0=3.354                     # V Excitation amplitude when the sample forces are not present
A0_avg=0                        # Calculate the excitation A_0 from the same scan? 0= no, 1=yes.
ex_units='eV'
#--------------------------
#

from netCDF4 import Dataset
from gwyfile.objects import GwyContainer, GwyDataField
import numpy as np
import fire
import os



def run(direct,prefix, cantsens, signalgain, Q, f0, ex_units="eV"):
    '''
    direct = File directory
    prefix = File prefix (everything that comes before -M or -xp etc)
    cant sens = cantilever sensitivity in V/nm
    signal gain = signal amp used on preamp
    Q = Q factr
    f0 = res freq
    exc_units=output units for the excitation (dissipation) channel. 
    '''
 
    fileprefix = direct+prefix
    w0 = f0*2*np.pi
    
    def botqrtavg(nums): 
        # This function takes the mean of the middle range of the data given
        out = np.sort(nums.flatten())[int(0.02*len(nums)*len(nums[0])):int(0.05*len(nums)*len(nums[0]))].mean()
        return out
    
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
        print(file.variables['opt_zpiezo_av'][:])
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
    
    Freq_r = datapull(fileprefix+'-Xp-ADC1.nc',-1/s_VHz) # Frequency shift divided by the conversion factor Df (V)/s_V/Hz=Df(Hz)Labone gain scaling
    Freq_l = datapull(fileprefix+'-Xm-ADC1.nc',-1/s_VHz)
    
    Amp_r = datapull(fileprefix+'-Xp-ADC4.nc',1e-10/cantsens/signalgain) #V to nm to m, potential amplifier, labone gain of 10
    Amp_l = datapull(fileprefix+'-Xm-ADC4.nc',1e-10/cantsens/signalgain)
    
    #printing some characteristics of the amplitude
    print('the amplitude type is',type(Amp_l))
    print('the amplitude size is',np.size(Amp_l))
    print('botqrtavg(Amp_l)',botqrtavg(Amp_l))

    #Dissipation ------------------------------------------------------------------------------------------------------------------------
    #dissipation before converting
    Exc_rV = datapull(fileprefix+'-Xp-ADC2.nc',1)                   #Excitation amplitude in V, trace
    print('average excitation amplitude right',botqrtavg(Exc_rV), "[V]")   #Excitation amplitude in V, retrace
    Exc_lV = datapull(fileprefix+'-Xm-ADC2.nc',1)
    print('average excitation amplitude left', botqrtavg(Exc_lV), "[V]")
    
    #Creating a copy of the excitaion amplitude in V save excitation in V     
    Diss_r=Exc_rV
    Diss_l=Exc_lV
    print(' ')
    print('Average value of dissipation in V',botqrtavg((Exc_rV+Exc_lV)/2) )
    #print('Average value of dissipation in V',botqrtavg((Exc_rV+Exc_lV))/2 )
    


    if ex_units=="Hz":
         #frac difference scaling
        
        Diss_r = (Diss_r/botqrtavg(Diss_r)-1)*w0/Q  # only valid when most of the image has a baseline dissipation, eg. rings on nanoparticles

        Diss_l = (Diss_l/botqrtavg(Diss_l)-1)*w0/Q  # only valid when most of the image has a baseline dissipation, eg. rings on nanoparticles

        Diss_r0 = (Exc_rV/Diss0-1)*w0/Q             # Dissipation in Hz using Diss0, which is the dissipation when the sample is not present. Valid if Diss0 is correct. 

        Diss_l0 = (Exc_lV/Diss0-1)*w0/Q             # # Dissipation in Hz using Diss0, which is the dissipation when the sample is not present. Valid if Diss0 is correct. 
    
    if ex_units=="eV": 
        
        '''
        # Dissipation units in eV 
        Ecuation 2.27, pg 26 in 
        S. Morita, R. Wiesendanger, and E. Meyer, eds., Noncontact Atomic Force Microscopy.
        Springer, hardcover ed., 7 2002.

        L. Fairgireve-Park thesis eq. 2.7, 2.8
        '''
        k=40 # [N/m]
        Amanual=4.5e-9 #m
        print('Amplitude from comments',Amanual)
        A_from_scan=1/2*(botqrtavg(Amp_r)+botqrtavg(Amp_l)) #nm
        print('A_from_scan',A_from_scan) 


        A=A_from_scan
        E0_J=np.pi*A**2*k/Q
        E0_eV=E0_J/1.602e-19

        print('E0=',E0_J,'[J]')
        print('E0=',E0_eV,'[eV]')
        #Diss0=0.772 # V, from the comments of the scan it  is the dissipation at 0V, far away. 


        #calculating the dissipation as the average of the center of the values of the image
        Diss_r = (Diss_r/botqrtavg(Diss_r)-1)*E0_eV
        Diss_l = (Diss_l/botqrtavg(Diss_l)-1)*E0_eV

        # Calculate the dissipation with respect to the dissipation far away

        # This is the dissipation as defined in the ncafm book
        Diss_r0 = (Exc_rV/Diss0-1)*E0_eV

        Diss_l0 = (Exc_lV/Diss0-1)*E0_eV


        #using the dissipation measured far away from the sample with vias=0V
        
         
    print(" ")
    print("average dissipation",botqrtavg(Diss_l+Diss_r) )
    
    A1w_r=datapull(fileprefix+'-Xp-ADC5.nc',1/1000) # one omega amplitude
    A1w_l=datapull(fileprefix+'-Xm-ADC5.nc',1/1000)

    A2w_r=datapull(fileprefix+'-Xp-ADC6.nc',1/10000) # two omega amplitude
    A2w_l=datapull(fileprefix+'-Xm-ADC6.nc',1/10000)

    Phase_r = datapull(fileprefix+'-Xp-ADC0mITunnel.nc',1/0.01) # Labone scalign
    Phase_l = datapull(fileprefix+'-Xm-ADC0mITunnel.nc',1/0.01)
  
#    r23_brockley_lini013-Xm-ADC0mITunnel.nc
    Vcpd_r = datapull(fileprefix+'-Xp-ADC7.nc',1) # Labone scalign
    Vcpd_l = datapull(fileprefix+'-Xm-ADC7.nc',1)
  

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
    
#dissipation with respect the the average of the image (excluding outliers, rings)
    obj['/4/data/title'] = 'Excitation amplitude ->'
    obj['/4/data'] = GwyDataField(Exc_rV, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='V')
    obj['/5/data/title'] = 'Excitation amplitude <-'
    obj['/5/data'] = GwyDataField(Exc_rV, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='V')

    obj['/6/data/title'] = 'Amplitude ->'
    obj['/6/data'] = GwyDataField(Amp_r, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='m')
    
    obj['/7/data/title'] = 'Amplitude <-'
    obj['/7/data'] = GwyDataField(Amp_l, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='m')
    
    obj['/8/data/title']=' Amplitude 1w <-'
    obj['/8/data']=GwyDataField(A1w_r, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='V')

    obj['/9/data/title']=' Amplitude 1w ->'
    obj['/9/data']=GwyDataField(A1w_l, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='V')

    obj['/10/data/title']=' Amplitude 2w <-'
    obj['/10/data']=GwyDataField(A2w_r, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='V')

    obj['/11/data/title']=' Amplitude 2w ->'
    obj['/11/data']=GwyDataField(A2w_l, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='V')


    obj['/12/data/title'] = 'PID Phase ->'
    obj['/12/data'] = GwyDataField(Phase_r, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='deg')
 
    obj['/13/data/title'] = 'PID Phase <-'
    obj['/13/data'] = GwyDataField(Phase_l, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='deg')
    

    #saving the dissipation Diss_i0 with i=r,l to a channel
    
    if A0_avg==1:

        #---------------------------------------
    # Ifs for dissipation units
        if ex_units=='Hz':  
            #dissipation with respect the the average of the image (excluding outliers, rings)
            obj['/14/data/title'] = 'Excess excitation ->'
            obj['/14/data'] = GwyDataField(Diss_r, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='Hz')
            obj['/15/data/title'] = 'Excess excitation  <-'
            obj['/15/data'] = GwyDataField(Diss_l, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='Hz')

        if ex_units=='eV':
            #dissipation with respect the the average of the image (excluding outliers, rings)
            obj['/14/data/title'] = 'Excess excitation ->'
            obj['/14/data'] = GwyDataField(Diss_r, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='eV/cycle')
            obj['/15/data/title'] = 'Excess excitation  <-'
            obj['/15/data'] = GwyDataField(Diss_l, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='eV/cycle')
        # ----------------------------------------
    elif A0_avg==0:
        
        if ex_units=='Hz':
            #Excitacion with respect to an excitation far away from the sample Diss0, manually entered
            obj['/14/data/title'] = 'Excess excitation ->'
            obj['/14/data'] = GwyDataField(Diss_r0, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='Hz')
            obj['/15/data/title'] = 'Excess excitation <-'
            obj['/15/data'] = GwyDataField(Diss_l0, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='Hz')

        if ex_units=='eV':

            obj['/14/data/title'] = 'Excess excitation ->'
            obj['/14/data'] = GwyDataField(Diss_r0, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='eV/cycle')
            obj['/15/data/title'] = 'Excess excitation <-'
            obj['/15/data'] = GwyDataField(Diss_l0, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='eV/cycle')



    obj['/16/data/title'] = 'Vcpd ->'
    obj['/16/data'] = GwyDataField(Vcpd_r, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='V')
 
    obj['/17/data/title'] = 'Vcpd <-'
    obj['/17/data'] = GwyDataField(Vcpd_l, xreal=xrange, yreal=yrange, xoff=xoff, yoff=yoff, si_unit_xy='m',si_unit_z='V')
    
    
    # Populating metadat field
    for i in range(0,10):
        obj['/'+str(i)+'/meta'] = meta
        
    # Making trace for each channel visible on file load
    for i in [0,2,4,6,8,10]:
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
    #Check if the directory to save the new file exists and if not create one
    dir_to_save = direct+"converted_to_gwy"

    if os.path.isdir(dir_to_save):
        print(f"{dir_to_save} exists.")
    else:
        print(f"{dir_to_save} does not exist. Creating it. ")
        os.mkdir(dir_to_save)



    obj.tofile(dir_to_save+'/'+prefix+ex_units+'A0_'+str(Diss0)+'.gwy')



run(direct, prefix, cantsens, signalgain, Q, f0,"eV")
