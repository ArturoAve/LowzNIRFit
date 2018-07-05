# Fit with Snoopy the light-curve data for a bunch of SNe
# It works for low-z and RAISINs.
#
#     Use:
#
# i) Open this file in a Text Editor. Select the settings for the fitting.
# ii) Open a new iPython Notebook and type in a cell:
#
#     %run "/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/Codes/\
#Canopy/3_Canopy_Fitting_v1_10.py"   "/Users/arturo/Dropbox/Research/\
#Articulos/12_RAISINs/Data/RAISIN_1/Data/PS1/2018_01_30/DavidJones/\
#snoopy_fit_v1/data"
#
# where the last path is the location of the folder containing the snoopy files.
# NOTE: Do not put a slash ("/") at the end of the path.
#
# The plots and ".snpy" files of the fitted data will be saved in folder
# automatically created called like ~/OpticalNIR/Fit/.
# It will be created also a log file listing all the SNe with errors during
# the fitting process.

#############################################################################

#   Version 1.2
# May 29, 2017
# - K-correction computation by Monte-Carlo sampling.

#############################################################################

#       USER

import sys # To read arguments in command line

# old Sample = '%s'%sys.argv[1]
DirSnoopyFiles = '%s'%sys.argv[1]
# print DirSnoopyFiles

# What bands to fit, optical or optical+nir:
BandType = 'Optical'   # Options: (OpticalNIR, Optical)

# Fit specific observer-frame bands? For instance, fit BVR bands only?
FitSpecificBands = True    # True, False

if FitSpecificBands == True:
    # For the the low-z SNe Ia, it is highly recommended  to include observer
    # frame B,V bands always: my initial fit to determine T_Bmax comes from
    # fitting BV only.
    # NIR LOW-Z. NOTE: It will used a dictionary containing the values of
    # T_Bmax derived from the optical LCs.

	#      OPTICAL LOW-Z
    # SpecificBands = ['BANDI', 'Bs', 'B', 'VANDI', 'Vs', 'V', 'V0', 'V1',
    #                 'RANDI','Rs', 'r', 'r_s']   # BVR bands

	#      NIR LOW-Z
	# NOTE: It will used a dictionary containing the values of T_Bmax derived
	# from the optical LCs.
    # SpecificBands = ['Y', 'Ydw']

    #      RAISIN 1
    # SpecificBands = ['r_ps1', 'i_ps1', 'z_ps1']
    # SpecificBands = ['r_ps1', 'i_ps1', 'z_ps1', 'f125w', 'f160w']

    #      RAISIN 2
    SpecificBands = ['r_des', 'i_des', 'z_des']
    # SpecificBands = ['r_des', 'i_des', 'z_des', 'f125w', 'f160w']


# All the posible names for a given band:
# B bands: ['BANDI', 'Bs', 'B']
# V bands: ['VANDI', 'Vs', 'V', 'V0']
# R bands: ['RANDI', 'Rs', 'r', 'r_s']

# Restframe band sets. Choose the restframe filter sets to use (the options are
# either CfA (NIR 2MASS) or the CSP natural filters):
# WATCH OUT!!: The only valid option for now is "CSP" becase Snoopy only accept
# the snoopy natural filters as the options for NIR restframes.
RestframeFilterSet =  'CSP'   # (CfA, CSP), but "CSP" is the only valid option for now.
# The only filters that Snoopy accepts for restframes are:
# EBV_model:
# ['u','B','V','g','r','i','Y','J','H','K','Bs','Vs','Rs','Is','J_K','H_K']

# Compute and use k-corrections during the fitting to correct the data?
# The option "True" is the standard and most commmon way I will use. The special
# case of  "kcorrData = False" is when somehow I already have the photometry
# k-corrected and I just want to fit the data (then I don't need to
# kcorrect the data again).
kcorrData = True

# Use the stretch (or dm15) to stretch the SED template in time before
# computing k-corrections?
# For the low-z paper use 'Use_k_stretch = False', to avoid to apply this
# additional correction during the fitting and computation of the k-corrections.
# For anything else it is better to use 'Use_k_stretch = True'.
Use_k_stretch = False # (True, False)

# In k-corr UNCERTAINTIES computation, for dates with no data.
# Use interpolation? "False" is preferred option.
interpol_kcorrError = False # (True, False).
LCModel_kcorrError = False # (True, False). Use model? "False" is preferred option.

# Number of simulation of the photometry to compute k-corr uncertainties
# When NumSim = 100, it takes ~4 min, 30 sec to run for one SN.
# "NumSim" can be any positive integer number. If NumSim = 0 then it is not
# tried to compute the uncertainties in kcorrs at all.
NumSim = 0          # Suggested numbers: 0, 50, 60, 100


#############################################################################

#       AUTOMATIC

from snpy import *
import numpy as np
import glob # To read the files in my directory
import os # To use command line like instructions
from matplotlib import pyplot as plt
import random # To compute k-corr uncertainties
import json # To save the simulated mag and k-corr uncertainties.

cc = 299792.458  # Speed of light (km/s)

# Change working directory to the folder where the snoopy files are.
# cwd = os.getcwd()
# print cwd
os.chdir(DirSnoopyFiles)
# os.getcwd()

#- Reading the LC data file names with the snoopy format.
the_list = glob.glob('*snoopy.dat')

print "# %s SNe in the initial list to be fitted."%len(the_list)

DirSaveOutput = DirSnoopyFiles+'/'+BandType+'/Fit/'
#- "If the subdirectory does not exist then create it"
if not os.path.exists(DirSaveOutput): os.makedirs(DirSaveOutput)

#-----------------------------------------------
# Get the name of this PYTHON (this doesn't work on iPython notebooks).
# To print it in the output text files as reference
NotebookName = os.path.basename(__file__)
print '#', (NotebookName)

#-----------------------------------------------
# Get the current date and time
import datetime
# Read the time and date now
now = datetime.datetime.now()

#---------------------------------------------------------------------------
#      (OPTIONAL) Low-z SNe Metadata

# I use the metadata information just to WRITE IN THE TEXT FILE comparison
# between the value of zcmb in the metadata file vs the one retrieved from NED
# by Snoopy. Both values and its difference is written in the
# "Fitting_Settings_Verbose.txt" file.

# DirMetadata = '/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/\
# AndyLCComp_2018_02/'
#
# # Reading the metadata file
# infoSNe_data = np.genfromtxt(DirMetadata+'carrick_Flow_corrections_snnames_v1.txt',
#                             dtype=['S17', float,float, 'S40',float,float,
#                                    float,float,float,float,'S16',float ])
#
# # Create a dictionary: {snname: ra, dec, zhelio}
# infoSNe_dict = {infoSNe_data['f0'][i]: np.array( [ infoSNe_data['f5'][i]/cc ] )
#                 for i in range(len(infoSNe_data)) }

#-----------------------------------------------
#---- Creation of some text files -----
#- Open a text file to write the failures on it.

textfile_1 = open(DirSaveOutput + 'Fitting_Settings_Verbose.txt', 'w')
fail_list = open(DirSaveOutput+'Fitting_Failure.log', 'w')

fail_list.write('Failed SNe to be fitted: \n')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_Author = '# Data table created by: Arturo Avelino \n'
text_Date   = '# On date: %s \n'%text_timenow
text_script = '# Script used: %s \n'%NotebookName
text_line = '#'+'-'*79 + '\n'

textfile_1.write("#         Fitting settings and results \n")
# old. textfile_1.write("# Sample: %s \n"%Sample)
textfile_1.write("# BandType: %s \n"%BandType)
textfile_1.write("# Fit specific bands?: %s \n"%FitSpecificBands)
if FitSpecificBands:
    textfile_1.write("# What are the specific bands to fit?: %s \n"%SpecificBands)
textfile_1.write("# Choose the restframe filter sets to use (the options \
are either CfA (NIR 2MASS) or the CSP natural filters): %s \n"%RestframeFilterSet)
textfile_1.write("# Compute and use k-corrections during the fitting to \
correct the data? %s \n"%kcorrData)
textfile_1.write("# Use_k_stretch = %s \n"%Use_k_stretch)
textfile_1.write("# Use GLOEs interpolation for dates with no data to determine \
colors during the k-corr uncertainty computation? = %s, \n"%interpol_kcorrError)
textfile_1.write("# or use a light-curve model to interpolate the dates with no \
data to determine colors during the k-corr uncertainty computation? = \
\%s \n"%LCModel_kcorrError)
textfile_1.write("# Number of Monte-Carlo simulations to determine the k-corr \
errors: %r \n"%NumSim)
textfile_1.write("# Directory of the snoopy-format data to be fitted: \n")
textfile_1.write('# %s \n'%DirSnoopyFiles)

textfile_1.write(text_line)
textfile_1.write(text_Author); textfile_1.write(text_Date); textfile_1.write(text_script);
textfile_1.write(text_line)
#
#------- Loop over all data -------

countSN = 0 # Counter number of SNe fitted correctly
countSNFail = 0 # # Counter number of SNe failed during the fitting.

#------------------------------

for file in the_list:
    try:
        print " "
        print "\n==================== %s ===================\n"%file[0:14]
        print "%s"%file

        # cwd = os.getcwd()
        # print cwd

        s = get_sn(file)
        # s = get_sn("%s"%file)
        s.summary()
        # print s.summary()

        #- Creation of an array with the name of the filters of this SN.
        FilterNames_array = []
        for band in s.restbands:
            FilterNames_array += [band]

        #--- Setting the NIR natural Snoopy bands as the NIR restframes bands

        if 'J2m' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['J2m'] = 'J'
            elif RestframeFilterSet == 'CfA': s.restbands['J2m'] = 'J2m'
        if 'JANDI' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['JANDI'] = 'J'
            elif RestframeFilterSet == 'CfA': s.restbands['JANDI'] = 'J2m'
        if 'J' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['J'] = 'J'
            elif RestframeFilterSet == 'CfA': s.restbands['J'] = 'J2m'
        if 'Jrc1' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['Jrc1'] = 'J'
            elif RestframeFilterSet == 'CfA': s.restbands['Jrc1'] = 'J2m'
        if 'Jrc2' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['Jrc2'] = 'J'
            elif RestframeFilterSet == 'CfA': s.restbands['Jrc2'] = 'J2m'
        if 'Jdw' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['Jdw'] = 'J'
            elif RestframeFilterSet == 'CfA': s.restbands['Jdw'] = 'J2m'
        if 'J_K' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['J_K'] = 'J'
            elif RestframeFilterSet == 'CfA': s.restbands['J_K'] = 'J2m'


        if 'H2m' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['H2m'] = 'H'
            elif RestframeFilterSet == 'CfA': s.restbands['H2m'] = 'H2m'
        if 'HANDI' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['HANDI'] = 'H'
            elif RestframeFilterSet == 'CfA': s.restbands['HANDI'] = 'H2m'
        if 'H' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['H'] = 'H'
            elif RestframeFilterSet == 'CfA': s.restbands['H'] = 'H2m'
        if 'Hdw' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['Hdw'] = 'H'
            elif RestframeFilterSet == 'CfA': s.restbands['Hdw'] = 'H2m'
        if 'H_K' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['H_K'] = 'H'
            elif RestframeFilterSet == 'CfA': s.restbands['H_K'] = 'H2m'


        if 'Ks2m' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['Ks2m'] = 'K'
            elif RestframeFilterSet == 'CfA': s.restbands['Ks2m'] = 'Ks2m'
        if 'KANDI' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['KANDI'] = 'K'
            elif RestframeFilterSet == 'CfA': s.restbands['KANDI'] = 'Ks2m'
        if 'Kd' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['Kd'] = 'K'
            elif RestframeFilterSet == 'CfA': s.restbands['Kd'] = 'Ks2m'
        if 'K' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['K'] = 'K'
            elif RestframeFilterSet == 'CfA': s.restbands['K'] = 'Ks2m'


        if 'Y' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['Y'] = 'Y'
            elif RestframeFilterSet == 'CfA': s.restbands['Y'] = 'Y'
        if 'Ydw' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['Ydw'] = 'Y'
            elif RestframeFilterSet == 'CfA': s.restbands['Ydw'] = 'Y'

        #--- Setting all the optical bands to the same restframe bands --

        if 'Us' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['Us'] = 'u'
            elif RestframeFilterSet == 'CfA': s.restbands['Us'] = 'Us'
        if 'u_s' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['u_s'] = 'u'
            elif RestframeFilterSet == 'CfA': s.restbands['u_s'] = 'Us'
        if 'u' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['u'] = 'u'
            elif RestframeFilterSet == 'CfA': s.restbands['u'] = 'Us'


        if 'Bs' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['Bs'] = 'B'
            elif RestframeFilterSet == 'CfA': s.restbands['Bs'] = 'Bs'
        if 'BANDI' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['BANDI'] = 'B'
            elif RestframeFilterSet == 'CfA': s.restbands['BANDI'] = 'Bs'
        if 'B' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['B'] = 'B'
            elif RestframeFilterSet == 'CfA': s.restbands['B'] = 'Bs'


        if 'g' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['g'] = 'g'
            elif RestframeFilterSet == 'CfA': s.restbands['g'] = 'g'


        if 'Vs' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['Vs'] = 'V'
            elif RestframeFilterSet == 'CfA': s.restbands['Vs'] = 'Vs'
        if 'VANDI' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['VANDI'] = 'V'
            elif RestframeFilterSet == 'CfA': s.restbands['VANDI'] = 'Vs'
        if 'V' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['V'] = 'V'
            elif RestframeFilterSet == 'CfA': s.restbands['V'] = 'Vs'
        if 'V0' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['V0'] = 'V'
            elif RestframeFilterSet == 'CfA': s.restbands['V0'] = 'Vs'
        if 'V1' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['V1'] = 'V'
            elif RestframeFilterSet == 'CfA': s.restbands['V1'] = 'Vs'


        if 'Rs' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['Rs'] = 'r'
            elif RestframeFilterSet == 'CfA': s.restbands['Rs'] = 'Rs'
        if 'RANDI' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['RANDI'] = 'r'
            elif RestframeFilterSet == 'CfA': s.restbands['RANDI'] = 'Rs'
        if 'r' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['r'] = 'r'
            elif RestframeFilterSet == 'CfA': s.restbands['r'] = 'Rs'
        if 'r_s' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['r_s'] = 'r'
            elif RestframeFilterSet == 'CfA': s.restbands['r_s'] = 'Rs'


        if 'Is' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['Is'] = 'i'
            elif RestframeFilterSet == 'CfA': s.restbands['Is'] = 'Is'
        if 'IANDI' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['IANDI'] = 'i'
            elif RestframeFilterSet == 'CfA': s.restbands['IANDI'] = 'Is'
        if 'i' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['i'] = 'i'
            elif RestframeFilterSet == 'CfA': s.restbands['i'] = 'Is'
        if 'i_s' in FilterNames_array:
            if   RestframeFilterSet == 'CSP': s.restbands['i_s'] = 'i'
            elif RestframeFilterSet == 'CfA': s.restbands['i_s'] = 'Is'

        #--------------------------------------

        #- Creation of an array with the specific band names to fit for this SN:
        if FitSpecificBands == True:
       	    BandsToFit = []; BandsExcludedOfFit = []
       	    for band in s.data.keys():
          		if band in SpecificBands: BandsToFit += [band]
          		else: BandsExcludedOfFit += [band]

        elif FitSpecificBands == False:
            #- Creation of an array with the name the OPTICAL only and
            # NIR only filters.
            # Observer-frame NIR bands in Andy's compilation or RAISINs
            # (but using Snoopy names).
            All_NIR_bands = ['Y','Ydw','JANDI','J2m', 'J', 'Jrc1', 'Jrc2','Jdw',
                            'HANDI', 'H2m', 'H', 'Hdw', 'KANDI', 'Ks2m', 'K',
                            'f125w', 'f160w']
            OpticalBands = [] # List to put the optical-only bands
            NIRbands = [] # List to put the NIR-only bands
            for band in s.data.keys():
                if band not in All_NIR_bands: OpticalBands += [band]
                else: NIRbands += [band]

            # print "Optical bands:", OpticalBands

        #--- Find out if filters (Bs, Vs) or (B,V) are present in the photometry.
        # If so, do a quick fit (without k-corr) to find T_Bmax,
        # if (Bs, Vs) or (B,V) are -not- present, then fit all the LC data with "s.fit()" and
        # write the name of the SN file in the failure text file.
        #-------
        if ('Bs' and 'Vs') in FilterNames_array:
            # print ("Bands ('Bs','Vs') in: %s \n" % (s.name ))
            # To quickly find the TB_max. It is needed 2 bands necessarily
            s.fit(['Bs','Vs'], kcorr=0)
        #-------
        elif ('B' and 'V')  in FilterNames_array:
            # print ("Bands ('B','V') in: %s \n" % (s.name ))
            # To quickly find the TB_max. It is needed 2 bands necessarily
            s.fit(['B','V'], kcorr=0)
	    #-------
        elif ('B' and 'V0')  in FilterNames_array:
            # print ("Bands ('B','V') in: %s \n" % (s.name ))
            # To quickly find the TB_max. It is needed 2 bands necessarily
            s.fit(['B','V0'], kcorr=0)
	    #-------
        elif ('B' and 'V1')  in FilterNames_array:
            # print ("Bands ('B','V') in: %s \n" % (s.name ))
            # To quickly find the TB_max. It is needed 2 bands necessarily
            s.fit(['B','V1'], kcorr=0)
        #-------
        elif ('BANDI' and 'VANDI')  in FilterNames_array:
            # print ("Bands ('BANDI','VANDI') in: %s \n" % (s.name ))
            # To quickly find the TB_max. It is needed 2 bands necessarily
            s.fit(['BANDI','VANDI'], kcorr=0)
        #-------
        else: # When there is NOT the bands (B, V) nor (Bs, Vs) in the LC data
            # print ("No bands in %s \n" % (s.name ))
            print "No B,V observed-frame bands found, but don't worry."
            fail_list.write('%s: No (B, V) or (Bs, Vs) bands. \n' % file)
            # s.fit(OpticalBands) #  Fit all the data

        print "%s. Prefitting (B,V) with no kcorrections: done."%s.name
        #-----------------------

        # if s.dm15 > 1.7:
            # using the 91bg SED template
        #     s.k_version = '91bg'

        #-------------------------------------------------------------------
        # MAIN FITTING, either, Optical only, Optica+NIR, or specific band data:

        if FitSpecificBands == True: # Final fit all the data
            print "Bands to plot:", BandsToFit
            s.fit(BandsToFit, kcorr=kcorrData, k_stretch=Use_k_stretch,
                reset_kcorrs=True)
    	    #old.  print s.summary()
    	    print 'Fitted bands:', BandsToFit

        elif FitSpecificBands == False:
            if BandType == 'Optical': # Final fit all the data
                print "Bands to plot:", OpticalBands
                s.fit(OpticalBands, kcorr=kcorrData, k_stretch=Use_k_stretch,
                    reset_kcorrs=True)
                #old. print s.summary()
                print "Fitting optical bands only: done"
            elif BandType == 'OpticalNIR': # Final fit all the data
                print "Bands to plot:", s.data.keys()
                s.fit(kcorr=kcorrData, k_stretch=Use_k_stretch, reset_kcorrs=True)
                #old. print s.summary()
                print "Fitting optical+nir bands: done"

        print "%s. Fitting all the bands: done with no issues."%s.name
        #-------------------------------------------------------------------
	    #  Saving the snpy data

	    # Removing the words "_snoopy.dat" at the end of the name for each SN.
        RemoveExtensionName = len(file)-11
        NameDataFileToSave = file[0:RemoveExtensionName]

        s.save(DirSaveOutput+NameDataFileToSave+'_1stFit.snpy')
        print "%s. The '_1stFit.snpy' file created and saved."%s.name

        s.summary(textFileOut=DirSaveOutput+NameDataFileToSave+'_SummaryFit_.txt')
        print "%s. The '_SummaryFit_.txt' file created and saved."%s.name

        # TO CHECK:
        # WHY I OBTAIN AN ERROR IN THE FOLLOWING INSTRUCTIONS WHEN I USE THE
        # COMBINED OPTION: BandType = 'Optical' ,  FitSpecificBands = True
        """
        #--- Saving the kcorrections ---
        kcorr_original_dict = {} # To store the original kcorr values.
        for band in s.data.keys():
            kcorr_original_dict[band] = list(s.ks[band])

        #--- Saving the dictionaries created above above using JSON format ---
        # Dictionary of the original kcorr values (i.e., no simulated)
        with open(DirSaveOutput+NameDataFileToSave+'_kcorr_CSP.json', 'w') as outfile:
            json.dump(kcorr_original_dict, outfile, sort_keys=True, indent=4)
        print "%s. Kcorrs to CSP system saved to JSON file."%s.name
        """

        #-------------------------------------------------------------------

        #    Write the summary to the global text file
        #- Open a text file to write the failures on it.
        textfile_1.write('------------------------------------------------------------ \n')
        textfile_1.write('%s \n'%NameDataFileToSave)
        textfile_1.write('SN %s \n'%s.name)
        textfile_1.write('z_hel = %r    ra=%r    dec=%r \n'%(s.z, s.ra, s.decl))
        textfile_1.write('Data in the following bands: %s \n'%s.data.keys())
        textfile_1.write('Fit results: \n')

        if BandType == 'OpticalNIR':
            for band in s.data.keys():
                textfile_1.write('    Observed %s fit to restband %s \n'%(
                band, s.restbands[band]))
        elif BandType == 'Optical':
            if FitSpecificBands == True:
                for band in BandsToFit:
                    textfile_1.write('    Observed %s fit to restband %s \n'%(
                    band, s.restbands[band]))
            elif FitSpecificBands == False:
                for band in OpticalBands:
                    textfile_1.write('    Observed %s fit to restband %s \n'%(
                    band, s.restbands[band]))

        textfile_1.write(' \n')
        textfile_1.write('    EBVhost = %.5f +/- %.5f \n'%(s.EBVhost, s.e_EBVhost))
        textfile_1.write('    EBVgal = %.5f +/- %.7f \n'%(s.EBVgal, s.e_EBVgal))
        textfile_1.write('    Tmax = %.4f +/- %.4f \n'%(s.Tmax, s.e_Tmax))
        textfile_1.write('    Distance mu = %.4f +/- %.4f \n'%(s.DM, s.e_DM))
        textfile_1.write('    dm15 = %.4f +/- %.4f \n'%(s.dm15, s.e_dm15))
        textfile_1.write('    z_CMB (NED) = %.5f  \n'%s.zcmb)
        # textfile_1.write('    z_CMB(Andy) = %.5f  \n'%infoSNe_dict['%s'%(s.name)][0])
        # textfile_1.write('    z_CMB(Andy) - zCMB(NED)  = %.5f  \n'%(
        #                         infoSNe_dict['%s'%(s.name)][0] - s.zcmb))
        textfile_1.write(' \n')

        #-------------------------------------------------------------------

        #       PLOTTING

        print "%s. Preparing to plot the fit."%s.name
        plt.close() # Close any possible plot unfinished/leftover.
        # Plot fit in phase
        s.plot(epoch=True, outfile=DirSaveOutput+NameDataFileToSave+'_PlotFit.png')
        plt.close()
        print "%s. Plot fit: done."%s.name

        # Plot over
        s.plot(epoch=True, single=True, offset=1,
               outfile=DirSaveOutput+NameDataFileToSave+'_PlotOver.png')
        plt.close()
        print "%s. Plot over fit: done."%s.name

        # Plot filters
        s.plot_filters(fill=True, outfile=DirSaveOutput+NameDataFileToSave+'_Filters.png')
        plt.close()
        print "%s. Plot filters: done."%s.name

        # Plot kcorrs
        s.plot_kcorrs(outfile=DirSaveOutput+NameDataFileToSave+'_PlotKcorrs.png')
        plt.close()
        print "%s. Plot kcorrs: done."%s.name

        countSN = countSN + 1

        # --- Plotting 1 --->

        s.plot(epoch=True) # Epoch for time axis

        #- Determine the y location for the text info. It is going to be below
        # (1.5 mag) from the maximum of the last filter to be plotted
        if FitSpecificBands == True:
            yLoc = s.get_max(bands=BandsToFit[0])[1]+0.8
        elif FitSpecificBands == False:
            if BandType == 'OpticalNIR':
                yLoc = s.get_max(bands=s.filter_order[-1])[1]+0.8
            elif BandType == 'Optical':
                yLoc = s.get_max(bands=s.filter_order[-(len(NIRbands)+1)])[1]+0.8

        plt.text(0, yLoc, r"$\Delta$m15 = %.2f $\pm$ %.2f"%(s.dm15,s.e_dm15),
                     color='black', fontsize=10)
        plt.text(0, yLoc+0.5, r"$z_{\rm hel}$ = %.3f" %(s.z),
                    color='black', fontsize=10)
        plt.text(0, yLoc+1.0, r"$\mu$ = %.3f $\pm$ %.3f"%(s.DM,s.e_DM),
                    color='black', fontsize=10)
        plt.text(0, yLoc+1.5, r"$T_{\rm Bmax}$ = %.2f $\pm$ %.3f"%(s.Tmax,s.e_Tmax),
                    color='black', fontsize=10)
        plt.text(0, yLoc+2.0, r"E(B-V)$_{\rm host}$ = %.3f $\pm$ %.3f"%(
                    s.EBVhost,s.e_EBVhost),color='black', fontsize=10)
        plt.text(0, yLoc+2.5, r"E(B-V)$_{\rm MW}$ = %.3f $\pm$ %.4f "%(
                    s.EBVgal, s.e_EBVgal),color='black', fontsize=10)
        plt.savefig(DirSaveOutput+NameDataFileToSave+"_PlotFitText.png",
                    format='png', dpi=90)
        plt.close()

        print "%s. Plot fit with text: done."%s.name
        print '%s. All Plotting done.'%s.name

        # <--- Plotting 1 ---
        # print 'End plot 1'

        """
        # --- Plotting 2 --->
        s.plot(epoch=True, single=1) # Epoch for time axis
        #- Determine the y location for the text info. It is going to be below (0.6 mag) from
        #- the minimum of the first filter to be plotted
        yLoc = s.get_max(bands=s.filter_order[0])[1]+2.5

        plt.text(-7, yLoc, r"$\Delta$m15 = %.2f $\pm$ %.2f"%(s.dm15,s.e_dm15),
                     color='black', fontsize=10)
        plt.text(-7, yLoc+0.2, r"$z_{\rm hel}$ = %.3f" %(s.z),
                    color='black', fontsize=10)
        plt.text(-7, yLoc+0.6, r"$\mu$ = %.3f $\pm$ %.3f"%(s.DM,s.e_DM),
                    color='black', fontsize=10)
        plt.text(-7, yLoc+0.4, r"E(B-V)$_{\rm host}$ = %.3f $\pm$ %.3f"%(s.EBVhost,s.e_EBVhost),
                    color='black', fontsize=10)

        plt.savefig(DirSaveOutput+NameDataFileToSave+"_PlotFitOverlay.png", format='png')
        plt.close()
        # <--- Plotting 2 ---
        """

        #-----------------------------------------------------------------------
        #    k-corr from CSP to Cousin & 2MASS passband systems

        """
        # Compute the kcorrs using my function
        # s.kcorrArt(interp=0, use_model=1, use_stretch=Use_k_stretch)
        s.kcorrArt(interp=1, use_model=0, use_stretch=Use_k_stretch) # Preferred
        print "%s. Kcorr to Cousin & 2MASS system: done."%s.name

        # Saving the fitting with the kcorrs for 2MASS and Johnson/Cousing
        s.save(DirSaveOutput+NameDataFileToSave+'_StdFilt.snpy')
        print "%s. The '_StdFilt.snpy' file created and saved."%s.name

        #--- Save the new k-corrections to a JSON file

        # kcorr values dictionary to store all the data for all bands and all random draws.
        # Saving the k-corrections.
        kcorr_StdFilt_dict = {}
        for band in s.data.keys():
            kcorr_StdFilt_dict[band] = list(s.ks[band])

        # Saving the dictionaries created above above using JSON format
        # Dictionary of the original kcorr values (i.e., no simulated)

        with open(DirSaveOutput+NameDataFileToSave+'_kcorr_Cousin2MASS.json', 'w') as outfile:
            json.dump(kcorr_StdFilt_dict, outfile, sort_keys=True, indent=4)
        print "%s. Kcorrs to Cousin & 2MASS system saved to JSON file."%s.name

        print "%s. The restframe filters now and the best fitted parameters are:"%s.name
        s.summary()

        #-----------------------------------------------------------------------

        #      Compute k-correction uncertainties

        if NumSim == 0:
            print "%s. No kcorr uncertainties will be computed."%s.name
        else:
            print "%s. Start Monte-Carlo simulations to determine k-corr uncertainties."%s.name

            s = get_sn(DirSaveOutput+NameDataFileToSave+'_1stFit.snpy')
            print "%s. Reading the '_1stFit.snpy' file: OK"%s.name

            # Creating a dictionary with all the original magnitude data.
            errmag_fix_dict = {} # error_mag dictionary for all bands
            mag_dict = {} # mag dict to store all the simulated photometry
            for band in s.data.keys():
                errmag_fix_dict[band] = s.data[band].e_mag
                mag_dict[band+'_0'] = list(s.data[band].mag) # initializing this dict
            print "%s. Create initial mag and err_mag dictionaries: OK"%s.name

            # kcorr values dictionary to store all the data for all bands and all random draws.
            # Copying the original kcorr dictionary to initialize the dictionary
            kcorr_dict = {} # To add the simulated kcorr values.
            for band in s.data.keys():
                kcorr_dict[band+'_0'] = list(s.ks[band])
            print "%s. Create initial kcorr dictionary: OK"%s.name

            #----- MAIN LOOP -------
            # This cell may take few minutes to run for each SN

            print "%s. Simulating %s times the photometry and computing their \
k-corrections."%(s.name, NumSim)

            for j in range(NumSim): # Loop over simulations
                for band in s.data.keys(): # Loop over bands.
                    # print band

                    # Loop over photometry in a given band.
                    for i in range(len(s.data[band].mag)):
                        muInt = 0; sigmaInt=0; # initialize these values.
                        muInt = mag_dict[band+'_0'][i]
                        sigmaInt = errmag_fix_dict[band][i]

                        # Generate Gaussian random photometry and redefine it as the "actual".
                        mag_sim_int = random.gauss(muInt, sigmaInt)
                        s.data[band].mag[i] = mag_sim_int

                        # Checking that everything is ok: Print on the screen for a
                        # given band and phase
                        # if band == 'Bs' and i==0: # for low-z
                        # if band == 'r_ps1' and i==0: # for low-z
                        #     print muInt, '|' , sigmaInt , '|' , mag_sim_int , '|',
                        #     s.data[band].mag[i]

                # k-correct the new simulated photometry:
                s.kcorr(interp=interpol_kcorrError, use_model=LCModel_kcorrError,
                    use_stretch=Use_k_stretch)

                # Save the new kcorr values and simulated magnitudes to the
                # "kcorr_dict" and "mag_dict" dicts.
                for band2 in s.data.keys():
                    kcorr_dict[band2+'_%s'%(j+1)] = list(s.ks[band2])
                    mag_dict[band2+'_%s'%(j+1)] = list(s.data[band2].mag)
            # <<--- end main loop for k-corr uncertainties

            print '%s. k-corr uncertainties computed. Now save the simulation.'%s.name

            # Create a dictionary with the mean and standard deviation of the kcorr values
            # at a given MJD for a given band.
            meanStd_kcorr_dict = {}
            for band in s.data.keys(): # Loop over bands.
                meanStd_kcorr_dict[band] = {}
                for i in range(len(s.data[band].MJD)): # Loop over MJD for a given band.
                    MJD_int = s.data[band].MJD[i] # Define the MJD.
                    # Loop over simulations for a given MJD and band.
                    meanStd_kcorr_dict[band][str(MJD_int)] = [
                    mean(np.array([kcorr_dict[band+'_%s'%j][i] for j in range(NumSim+1) ])),
                    std(np.array([kcorr_dict[band+'_%s'%j][i] for j in range(NumSim+1) ])) ]


            #--- Saving the dictionaries created above above using JSON format ---

            # Dictionary of Monte-Carlo simulated k-corrs
            with open(DirSaveOutput+NameDataFileToSave+'_kcorr_sim.json', 'w') as outfile:
                json.dump(kcorr_dict, outfile, sort_keys=True, indent=4)

            # Dictionary of mean and standard deviation.
            with open(DirSaveOutput+NameDataFileToSave+'_kcorr_Mean_Std.json', 'w') as outfile:
                json.dump(meanStd_kcorr_dict, outfile, sort_keys=True, indent=4)

            # Dictionary of simulated photometry.
            with open(DirSaveOutput+NameDataFileToSave+'_mag_sim.json', 'w') as outfile:
                json.dump(mag_dict, outfile, sort_keys=True, indent=4)

            print "%s. Simulation to determine k-corr uncertainties saved to JSON files."%s.name
            """


        #-----------------------------------------------------------------------

        print '%s: all done with no issues.'%s.name

    except:
        countSNFail = countSNFail + 1
        # Removing the words "snoopy.dat" at the end of the name for each SN:
        RemoveExtensionName = len(file)-11
        NameDataFileToSave = file[0:RemoveExtensionName]
        # fail_list.write('Failed fitting on object: %s \n' % file)
        fail_list.write('%s: Unable to fit. \n'%NameDataFileToSave)
        print "%s. Failed in some part during running this code."%s.name

fail_list.close()
textfile_1.close()

print "   "
print "-- Sweet: %r SNe fitted and %r failed -- "%(countSN, countSNFail)
