# nc-gwy
 Code that converts scans in the nc format, produced by gxsm to gwy format. Each channel is saved in the same file. 
'''
Code originally written by Logan Fairgrieve-Park

Description of the scripts: 


# nc_to_gwy_kpfm.py

Script created for kpfm scans. 
The script takes all the files with the "prefix" and converts it into a single .gwy file. 
It also, produces a .txt report, and copies itself into the data output directory for future reference. 


This script reads several nc files and  creates the following channels: 



0. Topography (right)
1. Topography (left)
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

# nc_to_gwy_stable.py

This script is for regular topography, dissipation, frequency shift scanning. 

 Code that converts scans in the nc format, produced by gxsm to gwy format. Each channel is saved in the same file. 
'''
This code creates the following channels: 


1 Topography (right)
2 Topography (left)
3. Frequency shift (right)
4. Frequency shift (left)
5. Dissipation (right) This dissipation is calculated assuming the mean values of the scan represent E_0=\pi k A^2 /Q
6. Dissipation (left)  This dissipation is calculated assuming the mean values of the scan represent E_0=\pi k A^2 /Q
7. Amplitude (right)  Oscillation amplitude nm
8. Amplitude (left )  Oscillation amplitude nm
9. 1 \omega amplitude  (right) If the voltage applied is modulated, this is the amplitude of the first sideband
10. 1 \omega amplitude (left)
11. PID Phase right 
12. PID Phase (left)
13. Excitation amplitude (right) [V]            Voltage Amplitude aplied to the excitation system to mantain a constant oscillation of the cantilever. 
14. Excitation amplitude (left)  [V]
15. Excess excitation (right)  [Hz] or [eV]     With respect to a parameter Diss0 which is the dissipation of the cantilever only, (far away from the sample)
16. Excess excitation (left)    [Hz] of [eV]    With respect to a parameter Diss0 which is the dissipation of the cantilever only, (far away from the sample)

Channel 15 and 16 are calculated with 

Diss_r0 = (Diss_r/Diss0-1)*E0_eV
Diss_l0 = (Diss_l/Diss0-1)*E0_eV

E0_J=np.pi*A**2*k/Q
E0_eV=E0_J/1.602e-19

'''