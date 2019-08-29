#!/usr/bin/env python
# coding: utf-8

# # Distance modulus computation using the NIR template
# 
# #### I have to run this notebook 2 times: 
# 1) The first time to determine:
# - the apparent magnitudes (and their uncertainties, sigma_m) at t_Bmax infered from the template fit.
# - uncertainty in the photometric distance moduli of the SNe sample defined as: 
#     error_mu_photometric^2 = sigma_m^2
# where sigma_m is the uncertainty in the apparent magnitude at t_Bmax computed in this first run.
# - the absolute magnitude for each SN from AbsMag = appMag - mu_LCDM(z).
# - the average absolute magnitude at t_Bmax, mean_AbsMag, with its standard deviation of the SNe sample. Write down these values in (AverageAbsMag_atMax, err_AverageAbsMag_atMax).
# 
# Run this notebook until the end of section "Determining average Absolute magnitude of the sample".
# Set: AbsMagFromHisto = False
# 
# 2) The second time to determine: 
# - the photometric distance moduli of the SNe sample defined as mu_photometric = appMag - mean_AbsMag.
# - the intrinsic dispersion of the Hubble-diagram residual, sigma_intrinsic. 
# 
# Run the entire notebook
# Set: AbsMagFromHisto = True
#     
# #### Diverse
# 
# Relationship between index (a.k.a, row or line) of the "Template_phase_mu_tau_FromR.dat" and the day (= phase):
# day = (index - 71)/2

# In[1]:


import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import quad as intquad
from scipy.optimize import fmin as simplex

import os
import glob # To read the name of the files in a given directory

# To read arguments in command line
# Used in the ".py" version of this notebook.
import sys 

#--------------------------------------------------------60
code_created_by = 'Arturo_Avelino'
# On date: 2017.01.10 (yyyy.mm.dd)
code_name = '11_DistanceMu_HubbleDiagram.ipynb'
version_code = '0.3.11'
last_update = '2019.08.06'
#--------------------------------------------------------60


# In[2]:


##############################################################################80


# ---

# # User interface

# In[ ]:


# There are 27 arguments to set in the 'terminal' version. They are:

#  1.- BandName = sys.argv[1]. Options: (Y, J, H ,K)
#  2.- vpecFix = int(sys.argv[2]). Peculiar velocity (km/s). Options: (150, 250)
#  3.- AbsMagFromHisto = sys.argv[3] == 'True'. Mean Absolute magnitude 
#         determined from histogram of 'appMagTmax_s - mu_s'?:
#  4.- NotebookToPlotOnly = sys.argv[4] == 'True'.
#  5.- DirSaveOutput = sys.argv[5].
#  6.- DirAppMag = sys.argv[6].
#  7.- DirTemplate = sys.argv[7].
#  8.- HoFix = float(sys.argv[8]).
#  9.- zcmbUpperLim = float(sys.argv[9]). Redshift cutoff. 
#         In the low-z paper I use zcmbUpperLim = 0.04.
# 10.- Average_NIRAbsMag_TBmax = float(sys.argv[10]).
# 11.- error_Average_NIRAbsMag_TBmax = float(sys.argv[11]).
# 12.- l_kern = float(sys.argv[12]);.
# 13.- PhaseMinTemp = int(sys.argv[13]) # -8 days.
# 14.- PhaseMaxTemp = int(sys.argv[14]) # 30, 41, days.
# 15.- EBVhostMin = float(sys.argv[15]) # -0.4 # host galaxy.
# 16.- EBVhostMax = float(sys.argv[16]) # 0.4 # host galaxy..
# 17.- EBVMWLim = float(sys.argv[17]) # Milky-Way galaxy.
# 18.- dm15LowerLim = float(sys.argv[18]) # I assume 0.8.
# 19.- dm15UpperLim = float(sys.argv[19]). # 1.6
# 20.- Chi2dofPrint = sys.argv[20] == 'True'.
# 21.- deltamu_print = sys.argv[21] == 'True'.
# 22.- DirSNeWithCepheid = sys.argv[22].
# 23.- BandMax = sys.argv[23]. Options: ( Bmax , NIRmax , Bmax_GP , Snoopy , SALT2 ).
# 24.- PlotTotalMu = sys.argv[24] == 'True'. Plot the "total" distance modulus 
#         derived from the three distance modulus computed from each band?
# 25.- BandsCombination = sys.argv[25]. Options: ( AllBands , JH , YJH , JHK ,  YJHK )
# 26.- plot_raisins =  sys.argv[26] == 'True'.
# 27.- minimize_residuals = sys.argv[27] == 'True'.


# In[3]:



## Terminal or notebook version of this script?
ScriptVersion = 'terminal' # ( terminal , notebook ) 

#-----------------------------------------
#    Command line version

if ScriptVersion == 'terminal':
    
    # What band to fit:(Y, J, H ,K)
    BandName = sys.argv[1]
    
    # Peculiar velocity (km/s). Options: (150, 250)
    # This number must be an integer: it is used to name the output folder.
    vpecFix = int(sys.argv[2])

    # Mean Absolute magnitude determined from histogram of 'appMagTmax_s - mu_s'?:
    # First run the notebook with this option with setting the value to 
    # "False" in order to determine 
    # the mean abs mag, then once I get that number, run a second time 
    # the notebook with this option as "True"
    # to compute the final tables and results.
    # If "True" it is generated the final 'DistanceMu_Good_AfterCutoffs_Main_.txt' 
    # text file, otherwise if "False"
    # generate temporal text files.

    AbsMagFromHisto = sys.argv[3] == 'True'

    # Use all this notebook just to plot the Hubble diagram?
    # Useful to create the HD for other cases  different to the template method.
    # If 'False' then compute the distance modulus using the template method.
    NotebookToPlotOnly = sys.argv[4] == 'True'
    
    # Dir where the "DistanceMu_Good_AfterCutoffs_Main_.txt" files is located.
    # Here will be saved the output too.
    DirSaveOutput = sys.argv[5]
    
    ## Location of the apparent magnitude light curves:
    DirAppMag = sys.argv[6]

    ## Location of the NIR template to fit:
    DirTemplate = sys.argv[7]    

    ## Hubble constant:
    HoFix = float(sys.argv[8])
    
    # Redshift cutoff.
    # In the low-z paper we use zcmbUpperLim = 0.04.
    zcmbUpperLim = float(sys.argv[9])      
    
#--------------------------------------------------------60
#   Notebook version

elif ScriptVersion == 'notebook':
    
    # What band to fit:(Y, J, H ,K)
    BandName = 'J' 
    
    # Peculiar velocity (km/s). Options: (150, 250)
    # This number must be an integer: it is used to name the output folder.
    vpecFix = 150  

    # Mean Absolute magnitude determined from histogram of 'appMagTmax_s - mu_s'?:
    # First run the notebook with this option with setting the value to 
    # "False" in order to determine 
    # the mean abs mag, then once I get that number, run a second time 
    # the notebook with this option as "True"
    # to compute the final tables and results.
    # If "True" it is generated the final 'DistanceMu_Good_AfterCutoffs_Main_.txt' 
    # text file, otherwise if "False"
    # generate temporal text files.

    AbsMagFromHisto = True  

    # Use all this notebook just to plot the Hubble diagram?
    # Useful to create the HD for other cases  different to the template method.
    # If 'False' then compute the distance modulus using the template method.
    NotebookToPlotOnly = False # Notebook version.  
    
    ## Dir where I'll save the output, or  dir where the 
    ## "DistanceMu_Good_AfterCutoffs_Main_.txt" files is located.
    # DirSaveOutput = '/Users/arturo/Dropbox/Research/Articulos/14_RAISINs/Data/\
# raisin12/hubblediagram/SALT2_optical/v2/plots_HD/'
    
    DirSaveOutput = '/Users/arturo/Downloads/tmp/SNIRfit/'
    
    ##----automatic when ScriptVersion=='notebook' ----
    
    ## Location of the apparent magnitude light curves:
    # DirAppMag = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/\
# 10Compute/TheTemplates/%s_band/Std_filters/1_AllData_InitialFit/\
# AbsMag/AllSamples/'%BandName
    
    DirAppMag = '/Users/arturo/Downloads/tmp/SNIRfit/SNLCs/'

    ## Location of the NIR template to fit:
    # DirTemplate = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/\
# 10Compute/TheTemplates/%s_band/Std_filters/3_Template_FlatPrior/\
 #AllSamples_vpec_%s/z_gr_0/'%(BandName,vpecFix)
    
    DirTemplate = '/Users/arturo/Downloads/tmp/SNIRfit/'

    HoFix = 73.24 # Valued used by default in the paper
    # HoFix = 72.0 # 72 # 73.24  # Hubble constant (km/(s Mpc))
    # HoFix = 72.78 # value reported by Dhawan et al 2017.
    
    # Redshift cutoff.
    # In the low-z paper we use zcmbUpperLim = 0.04.
    zcmbUpperLim = 0.04
    
##############################################################################80

#--- Fixed values ---
OmMFix = 0.28 # Omega_Matter
OmLFix = 0.72 # Omega_Lambda
wFix = -1.0 # Dark energy EoS
c = 299792.458  # Speed of light (km/s)
# In the process to convert all the 'c' to 'cc'
cc = 299792.458  # Speed of light (km/s) 

#--- Uncertainty in z_CMB:---

# Used to plot the -theoretical- peculiar velocity uncertanty in the
# Hubble residual plot.

# Dan Scolnic gave me the value of err_cz_CMB = 150 km/s about the 
# collection of z_CMB values he provided me.
# err_zCMB_fix =  0.0005003461427972281 when err_cz_CMB = 150 km/s
err_zCMB_fix = 150.0/cc  # 150 has units of km/s

# Just to plot the peculiar velocity uncertainty curve in the 
# Hubble-diagram residual, 
# I was using,  err_zCMB_fix = 0.001, that is the average 
# error_zcmb in Andy's compilation.

#--------------------------------------------

#   SELECT THE APPARENT MAG DATA TO USE TO CONSTRUCT THE HUBBLE DIAGRAM
# The plots with cuts z>0 and z>0.01 are done automatically.

# KindOfData4HD = 'CfA'
# KindOfData4HD = 'CSP'
# KindOfData4HD = 'Others'
KindOfData4HD = 'AllSamples'

#--------------------------------------------

#   SELECT THE TEMPLATE TO USE TO COMPUTE THE DISTANCE MODULUS

# This will be also the data used to construct the Hubble diagram
# Selecting from '3_Template' folder

#- Use the template constructed from the subsample?:
# KindOfTemp = 'CfA'
# KindOfTemp = 'CSP'
# KindOfTemp = 'Others'
KindOfTemp = 'AllSamples'

KindOfTempSubgroup = 'z_gr_0'
# KindOfTempSubgroup = 'z_gr_001'
# KindOfTempSubgroup = 'AllGoodData'

# Indicate the technique used to compute the GP fitting in 
# the '1_AllData_InitialFit' step:
# TempPrior_Hyper = Using a template prior (computed from 
# the Moving Windows Average template) for the Gaussian Process 
# fitting, then determining the hyperparameters using all 
# the LCs simultaneously.
# FlatPrior= Assuming a flat prior at ~ -17 Abs mag, then 
# computing the hyperparameters for each LC independently.
# MWA = Moving windows average

# TempType = 'TempPrior_Hyper' 
TempType = 'FlatPrior' # I use this option for the paper
# TempType = 'MWA' 

# - Smoothed Moving windows average to construct the template
# TempType = 'MWA' # Moving windows average
# If TempType = 'MWA' then specify the template file to use
MWATempTypeFile = 'TempWeightedSmooth_Box7_Step05_Window21_Poly3.dat'

# Use a build-in normalized template?:
# Normalized = True
# It is, use a template that was constructed assuming 
# vPec=0 km/s during the GP fitting of the individual LCs, 
# and then with the option 'NormalizedTemp' in the 
# hierarchical Bayesian code.


# In[4]:


#-----------------------------------------------------------------------------80


# In[5]:


# Reset the values:
Average_NIRAbsMag_TBmax = None;
error_Average_NIRAbsMag_TBmax = None;

#--------------------------------------

if BandName == 'J':

    # These values have to be computed the first time I run this notebook.
    # NOTE: The values of "Average_NIRAbsMag_TBmax" and
    # "error_Average_NIRAbsMag_TBmax" do NOT depend on the value of "vpecFix".

    if ScriptVersion == 'notebook':

        # Values used in the low-z paper with the template normalized
        # at phase=0 days.
        #   J band | vpecFix = 150 km/s | 0 < z < 0.04
        Average_NIRAbsMag_TBmax = -18.339435  # 2018-09-28; 12:50 hrs.
        error_Average_NIRAbsMag_TBmax = 0.174834;

        # Fitting with a template normalized at phase=+15 days.
        #   J band | vpecFix = 150 km/s | 0 < z < 0.04
        # Average_NIRAbsMag_TBmax = -16.832741  # 2019-02-22; 10:48 hrs.
        # error_Average_NIRAbsMag_TBmax = 0.187574;

        # Length hyperparameter from the GP fit.
        l_kern = 7.0199;

        # PHASE RANGE IN DAYS OF TEMPLATE TO USE TO COMPUTE
        # THE DISTANCE MODULUS
        PhaseMinTemp = -8 # -8 days
        PhaseMaxTemp = 30 # 30, 41, days

        #-- EBV cutoff
        EBVhostMin = -0.4 # -0.4 # host galaxy
        EBVhostMax = 0.4 # 0.4 # host galaxy.
        EBVMWLim = 1.0 # Milky-Way galaxy

        #-- dm15 cutoff
        dm15LowerLim = 0.8 # I assume 0.79.
        dm15UpperLim = 1.6

    #------------------------------

    elif ScriptVersion == 'terminal':
        Average_NIRAbsMag_TBmax = float(sys.argv[10])
        error_Average_NIRAbsMag_TBmax = float(sys.argv[11])

        # Length hyperparameter from the GP fit.
        l_kern = float(sys.argv[12]);

        # PHASE RANGE IN DAYS OF TEMPLATE TO USE TO COMPUTE
        # THE DISTANCE MODULUS
        PhaseMinTemp = int(sys.argv[13]) # -8 days
        PhaseMaxTemp = int(sys.argv[14]) # 30, 41, days

        #-- EBV cutoff
        EBVhostMin = float(sys.argv[15]) # -0.4 # host galaxy
        EBVhostMax = float(sys.argv[16]) # 0.4 # host galaxy.
        EBVMWLim = float(sys.argv[17]) # Milky-Way galaxy

        #-- dm15 cutoff
        dm15LowerLim = float(sys.argv[18]) # I assume 0.79.
        dm15UpperLim = float(sys.argv[19])


    # Ignore data with a chi^2_dof larger than a given threshold:
    # Use the values ( chi2_dof_Max = 1e6; chi2_dof_Max_Label = 'chi_1e6')
    # the first time I run the notebook.
    chi2_dof_Max = 1e6; chi2_dof_Max_Label = 'chi_1e6'
    # chi2_dof_Max = 3.5; chi2_dof_Max_Label = 'chi3_5'

    # - Smoothed Moving windows average to construct the template
    # If TempType = 'MWA' then specify the template file to use
    MWATempTypeFile = 'TempWeightedSmooth_Box7_Step05_Window21_Poly3.dat'

###########################################################

if BandName == 'Y':

    # These values have to be computed the first time I run this notebook.
    # NOTE: The values of "Average_NIRAbsMag_TBmax" and
    # "error_Average_NIRAbsMag_TBmax" do NOT
    # depend on the value of "vpecFix".

    if ScriptVersion == 'notebook':
        #   Y band | vpecFix = 150 km/s | 0 < z < 0.04
        Average_NIRAbsMag_TBmax = -18.124377  # 2018-09-29; 12:59 hrs.
        error_Average_NIRAbsMag_TBmax = 0.152012;

        # Length hyperparameter from the GP fit.
        l_kern = 7.9002;

        # PHASE RANGE IN DAYS OF TEMPLATE TO USE TO COMPUTE THE DISTANCE MODULUS
        PhaseMinTemp = -8 # -8 days
        PhaseMaxTemp = 30 # 30, 41, days

        #-- EBV cutoff
        EBVhostMin = -0.4 # -0.4 # host galaxy
        EBVhostMax = 0.4 # 0.4 # host galaxy.
        EBVMWLim = 1.0 # Milky-Way galaxy

        #-- dm15 cutoff
        dm15LowerLim = 0.8 # I assume 0.8.
        dm15UpperLim = 1.6

    #------------------------------

    elif ScriptVersion == 'terminal':
        Average_NIRAbsMag_TBmax = float(sys.argv[10])
        error_Average_NIRAbsMag_TBmax = float(sys.argv[11])

        # Length hyperparameter from the GP fit.
        l_kern = float(sys.argv[12]);

        # PHASE RANGE IN DAYS OF TEMPLATE TO USE TO COMPUTE
        # THE DISTANCE MODULUS
        PhaseMinTemp = int(sys.argv[13]) # -8 days
        PhaseMaxTemp = int(sys.argv[14]) # 30, 41, days

        #-- EBV cutoff
        EBVhostMin = float(sys.argv[15]) # -0.4 # host galaxy
        EBVhostMax = float(sys.argv[16]) # 0.4 # host galaxy.
        EBVMWLim = float(sys.argv[17]) # Milky-Way galaxy

        #-- dm15 cutoff
        dm15LowerLim = float(sys.argv[18]) # I assume 0.79.
        dm15UpperLim = float(sys.argv[19])


    # Ignore data with a chi^2_dof larger than a given threshold:
    # Use the values ( chi2_dof_Max = 1e6; chi2_dof_Max_Label = 'chi_1e6')
    # the first time I run the notebook.
    # chi2_dof_Max = 3.5; chi2_dof_Max_Label = 'chi3_5'
    # chi2_dof_Max = 3; chi2_dof_Max_Label = 'chi3' # Low-z
    chi2_dof_Max = 10; chi2_dof_Max_Label = 'chi10' # RAISIN

    # - Smoothed Moving windows average to construct the template
    # If TempType = 'MWA' then specify the template file to use
    MWATempTypeFile = 'TempWeightedSmooth_Box7_Step05_Window21_Poly3.dat'

###########################################################

if BandName == 'H':

    # These values have to be computed the first time I run this notebook.
    # NOTE: The values of "Average_NIRAbsMag_TBmax"
    # and "error_Average_NIRAbsMag_TBmax" do NOT
    # depend on the value of "vpecFix".

    if ScriptVersion == 'notebook':
        #   H band | vpecFix = 150 km/s | 0 < z < 0.04
        Average_NIRAbsMag_TBmax = -18.181440  # 2018-09-29; 14:39 hrs.
        error_Average_NIRAbsMag_TBmax = 0.165906;

        # Length hyperparameter from the GP fit.
        l_kern = 9.8131;

        # PHASE RANGE IN DAYS OF TEMPLATE TO USE TO COMPUTE THE DISTANCE MODULUS
        PhaseMinTemp = -8 # -8 days
        PhaseMaxTemp = 30 # 30, 41, days

        #-- EBV cutoff
        EBVhostMin = -0.4 # -0.4 # host galaxy
        EBVhostMax = 0.4 # 0.4 # host galaxy.
        EBVMWLim = 1.0 # Milky-Way galaxy

        #-- dm15 cutoff
        dm15LowerLim = 0.8 # I assume 0.78.
        dm15UpperLim = 1.6

    #------------------------------

    elif ScriptVersion == 'terminal':
        Average_NIRAbsMag_TBmax = float(sys.argv[10])
        error_Average_NIRAbsMag_TBmax = float(sys.argv[11])

        # Length hyperparameter from the GP fit.
        l_kern = float(sys.argv[12]);

        # PHASE RANGE IN DAYS OF TEMPLATE TO USE TO COMPUTE
        # THE DISTANCE MODULUS
        PhaseMinTemp = int(sys.argv[13]) # -8 days
        PhaseMaxTemp = int(sys.argv[14]) # 30, 41, days

        #-- EBV cutoff
        EBVhostMin = float(sys.argv[15]) # -0.4 # host galaxy
        EBVhostMax = float(sys.argv[16]) # 0.4 # host galaxy.
        EBVMWLim = float(sys.argv[17]) # Milky-Way galaxy

        #-- dm15 cutoff
        dm15LowerLim = float(sys.argv[18]) # I assume 0.79.
        dm15UpperLim = float(sys.argv[19])


    # Ignore data with a chi^2_dof larger than a given threshold:
    # Use the values ( chi2_dof_Max = 1e6; chi2_dof_Max_Label = 'chi_1e6')
    # the first time I run the notebook.
    chi2_dof_Max = 1e6; chi2_dof_Max_Label = 'chi_1e6'
    # chi2_dof_Max = 3.5; chi2_dof_Max_Label = 'chi3_5'

    # - Smoothed Moving windows average to construct the template
    # If TempType = 'MWA' then specify the template file to use
    MWATempTypeFile = 'TempWeightedSmooth_Box7_Step05_Window21_Poly3.dat'

###########################################################

if BandName == 'K':

    # These values have to be computed the first time I run this notebook.
    # NOTE: The values of "Average_NIRAbsMag_TBmax" and
    # "error_Average_NIRAbsMag_TBmax" do NOT
    # depend on the value of "vpecFix".

    if ScriptVersion == 'notebook':
        #   K band | vpecFix = 150 km/s | 0 < z < 0.04
        Average_NIRAbsMag_TBmax = -18.349367  # 2018-09-29; 14:56 hrs.
        error_Average_NIRAbsMag_TBmax = 0.206866;

        # Length hyperparameter from the GP fit.
        l_kern = 8.1853;

        # PHASE RANGE IN DAYS OF TEMPLATE TO USE TO COMPUTE THE DISTANCE MODULUS
        PhaseMinTemp = -8 # -8 days
        PhaseMaxTemp = 30 # 30, 41, days

        #-- EBV cutoff
        EBVhostMin = -0.4 # -0.4 # host galaxy
        EBVhostMax = 0.4 # 0.4 # host galaxy.
        EBVMWLim = 1.0 # Milky-Way galaxy

        #-- dm15 cutoff
        dm15LowerLim = 0.8 # I assume 0.78.
        dm15UpperLim = 1.6

    #------------------------------

    elif ScriptVersion == 'terminal':
        Average_NIRAbsMag_TBmax = float(sys.argv[10])
        error_Average_NIRAbsMag_TBmax = float(sys.argv[11])

        # Length hyperparameter from the GP fit.
        l_kern = float(sys.argv[12]);

        # PHASE RANGE IN DAYS OF TEMPLATE TO USE TO COMPUTE
        # THE DISTANCE MODULUS
        PhaseMinTemp = int(sys.argv[13]) # -8 days
        PhaseMaxTemp = int(sys.argv[14]) # 30, 41, days

        #-- EBV cutoff
        EBVhostMin = float(sys.argv[15]) # -0.4 # host galaxy
        EBVhostMax = float(sys.argv[16]) # 0.4 # host galaxy.
        EBVMWLim = float(sys.argv[17]) # Milky-Way galaxy

        #-- dm15 cutoff
        dm15LowerLim = float(sys.argv[18]) # I assume 0.79.
        dm15UpperLim = float(sys.argv[19])


    # Ignore data with a chi^2_dof larger than a given threshold:
    # Use the values ( chi2_dof_Max = 1e6; chi2_dof_Max_Label = 'chi_1e6')
    # the first time I run the notebook.
    # chi2_dof_Max = 1e6; chi2_dof_Max_Label = 'chi_1e6'
    chi2_dof_Max = 4; chi2_dof_Max_Label = 'chi4'

    # - Smoothed Moving windows average to construct the template
    # If TempType = 'MWA' then specify the template file to use
    MWATempTypeFile = 'TempWeightedSmooth_Box7_Step05_Window21_Poly3.dat'

#########################################################################

# Warning in case that the mean abs mag was not defined:
if Average_NIRAbsMag_TBmax == None and  error_Average_NIRAbsMag_TBmax == None:
    print "# ERROR: The mean abs magnitude and its std dev have to be defined."

#########################################################################

# Text in the notes about the origin of (Average_NIRAbsMag_TBmax, error_Average_NIRAbsMag_TBmax)
if AbsMagFromHisto == True:
    Notes1 = 'Determined from the histogram of {appMagTmax_s - mu_LCDM_s}'
elif AbsMagFromHisto == False:
    Notes1 = 'Just a random guess about these values. It doesnt matter the value of this guess because this numbers are to determine the distance modulus but in this first step I m computing the appMag_TBmax only.'

#------------------

# About 'l_kern': Lenght scale 'l' of the kernel in the covariance matrix for the residual.
# this value is ignored in case of 'CovMat_MeanMu = False'
# The values I obtained when using the global marginal likehood function to compute
# the hyperparameter 'l' using all the LCs together.

#==============================================================

#             Method to determine the distance modulus
# (see description below)

Method = 7 # Options  = (1,2,3,4,5,6,7). In the paper I use method 7, before I was using 5.
# The best methods are 5 (unnormalized template) and 7 (normalized template).

# COVARIANCE MATRICES
# It is not necessary to modify the following 2 lines unless some very particular computation is needed.
# Using the covariance matrix of the photometry to determine the -mean- and -uncertainty-
# of the distance modulus?:
if Method==1 or Method==5 or Method==6 or Method==7: CovMat_MeanMu = False; CovMat_ErrorMu = True
elif Method==2 or Method==3 or Method==4: CovMat_MeanMu = True; CovMat_ErrorMu = True

# CovMat_MeanMu = False # (True, False). # In the paper I use 'False'
# CovMat_ErrorMu = False # (True, False). # In the paper I use 'True'

# Use the peculiar velocity covariance square matrix to determine the distance modulus?:
#OLD.  Use_CovMat_PecVel = True # (True, False).
# NOTE: When "Use_CovMat_PecVel == True" then in the total covariance matrix it is used the -uncertainty-
# in the mean template instead of the -population standard deviation- of the template.

#-------------

#-- Minimal number of data in LC:
MinNumOfDataInLC = 1

# Print the (RMS, WRMS, sigma_int) text in the residual plot?:
WRMS_label = True # Print the WRMS of the total sample (or total+subsamples)?
RMS_simple = True # or print instead the simple RMS only of total and subsamples?
WRMS_subsamples = True # Print the RMS, WRMS of each subsample?

#=================================================================
#( Fixed values usually)

# --- (FIX) Residual cutoff ---
# residualMax = 0.7; residualMax_Label = 'resid07'
# residualMax = 0.8; residualMax_Label = 'resid08'
# residualMax = 0.9; residualMax_Label = 'resid09'
# residualMax = 1; residualMax_Label = 'resid1'
# residualMax = 3; residualMax_Label = 'resid3'
# residualMax = 5; residualMax_Label = 'resid5'
residualMax = 20; residualMax_Label = 'resid20'

#---- (FIX) Limits in the plots for the Hubble diagram and residual  ----

# For the Hubble diagram:
xlimPlots = 0.0015, zcmbUpperLim+0.006 # For low-z only
ylimPlots = 29.8, 38 # For low-z only

# For the Hubble residual
ylimPlots_residual = -1.1, 1.1

#----------------

# Plotting range of the fitted LC figures.
x_RangePlots = -10, 60;

if ScriptVersion == 'notebook':
    # In the plots of the individual LC, print the chi2_dof value?:
    Chi2dofPrint = True

    # Print residual value:
    deltamu_print = False

elif ScriptVersion == 'terminal':
    # In the plots of the individual LC, print the chi2_dof value?:
    Chi2dofPrint = sys.argv[20] == 'True'

    # Print residual value:
    deltamu_print = sys.argv[21] == 'True'

#--------------------------------------------

#---- (FIX) Filter system ----
FilterSyst = 'Std_filters/'
# FilterSyst = 'CSP_filters/'


# In[6]:


#-----------------------------------------------------------------------------80


# In[7]:


# Description of the methods

""" 
- Method 1: 
	- No use of CovMatrix to determine the best estimate of distance mu (but still optional using 
    "CovMat_use == True" [default: "CovMat_use == False" ]),  BUT using it (i.e., "CovMat_use == True") 
    to determine its uncertainty [default: "CovMat_use == True" ].
	- using the uncertainty in the mean template instead of the population standard deviation
	- adding the uncertainty in the distance modulus due to the peculiar velocity uncertainty in the total 
    covariance matrix -before- to compute the distance-modulus uncertainty.
	- CovMat_MeanMu == False
	- CovMat_ErrorMu == True
	- Template = np.genfromtxt('Template_phase_mu_stdError_FromR.dat')
	- CovMatrix_Mu = CovMatResidualMu + Cov_appmag # Total covariance matrix -before- computing
	- CovMatrix_ErrorMu = CovMatCorrData_4ErrorMu + Cov_appmag + Cov_pecvel # Total covariance matrix -before- 
    computing
	- sigma2_mu = (1/Denominator_ErrorMu) # Variance mu
	

- Method 2: 
	- using the uncertainty in the mean template instead of the population standard deviation
	- mean mu: using the CovMat_MeanMu == True and ignoring Cov_pecvel to compute it.
	- error mu: adding in quadratures the uncertainty in the distance modulus due to the peculiar velocity 
    uncertainty with the uncertainty in the distance modulus from fitting the template -after- computing 
    the distance modulus
    - Apparently, this method gives exactly the same results than Method 3
    - CovMat_MeanMu == True
	- CovMat_ErrorMu == True
	- Template = np.genfromtxt('Template_phase_mu_stdError_FromR.dat')
	- CovMatrix_Mu = CovMatResidualMu + Cov_appmag # Total covariance matrix -before- computing
	- CovMatrix_ErrorMu = CovMatCorrData_4ErrorMu + Cov_appmag # Total covariance matrix -before- computing
	- sigma2_mu = (1/Denominator_ErrorMu) +  sigma2_pec(zcmbInt, err_zcmb, vpecFix) # Variance mu


- Method 3: 
	- using the uncertainty in the mean template instead of the population standard deviation
	- adding the uncertainty in the distance modulus due to the peculiar velocity uncertainty in the total 
    covariance matrix -before- to compute the mean distance modulus and its uncertainty.
    - Apparently, this method gives exactly the same results than Method 2
	- CovMat_MeanMu == True
	- CovMat_ErrorMu == True	
	- Template = np.genfromtxt('Template_phase_mu_stdError_FromR.dat')
	- CovMatrix_Mu = CovMatResidualMu + Cov_appmag + Cov_pecvel # Total covariance matrix -before- computing
	- CovMatrix_ErrorMu = CovMatCorrData_4ErrorMu + Cov_appmag + Cov_pecvel # Total covariance matrix -before- 
    computing
	- sigma2_mu = (1/Denominator_ErrorMu) # Variance mu
	
	
- Method 4: 
	- using the population standard deviation of the template
	- NOT including the uncertainty in the distance modulus due to the peculiar velocity during the fitting.
	= It produces very small error bars for the distance-modulus uncertainty
	- CovMat_MeanMu == True
	- CovMat_ErrorMu == True	
	- Template = np.genfromtxt('Template_phase_mu_tau_FromR.dat')
	- CovMatrix_Mu = CovMatResidualMu + Cov_appmag # Total covariance matrix -before- computing
	- CovMatrix_ErrorMu = CovMatCorrData_4ErrorMu + Cov_appmag # Total covariance matrix -before- computing
	- sigma2_mu = (1/Denominator_ErrorMu) # Variance mu
    
- Method 5: 
    - Equal to the method 1 with the difference of using the -population standard deviation- of the template 
    instead of its uncertainty. The idea of this is that during the fitting of the template to the app mag data 
    is that it gives more weight to the data around the first NIR peak than the data in the 2nd peak because the 
    population standard deviation of the template is larger around the 2nd peak, so that my fitting is not 
    screwed up by the fact that the 2nd peak seems to have an independent variability between the SNe, then 
    biasing my distance mu estimations.
	- CovMat_MeanMu == False
	- CovMat_ErrorMu == True
	- Template = np.genfromtxt('Template_phase_mu_tau_FromR.dat')
	- CovMatrix_Mu = CovMatResidualMu + Cov_appmag # Total covariance matrix -before- computing
	- CovMatrix_ErrorMu = CovMatCorrData_4ErrorMu + Cov_appmag + Cov_pecvel # Total covariance matrix -before- computing
	- sigma2_mu = (1/Denominator_ErrorMu) # Variance mu
	
    
- Method 6:
    - So far, equal to method 5 but now normalizing the template, then fit it to determine the apparent 
    magnitude at T_Bmax, then subtract the absolute magnitude at T_Bmax to determine the photometric
    distance modulus.
    - I take the regular template constructed in the absolute magnitude frame then I normalize it here in this 
    notebook, so that 
    at T_Bmax I set the magnitude to be equal zero.
    
- Method 7:
    - I used a normalized template that was normalized by construction in the hierarchical-Bayesian-R script.
    - The rest, same than Method 6.
    - The notebook has to be run twice. The first run is to determine the NIR apparent magnitude at t_Bmax only
      by fitting the template to the apparent-magnitude light curves, then using the relation, 
      AbsMag_s = appMag_TBmax - mu_LCDM, it is determined the NIR absolute magnitude at t_Bmax, AbsMag_s, for
      each SN. In this run it is created a value of the distance modulus but it is not relevant/reliable. 
      After the main loop, there is a section to determine <Average_NIRAbsMag_TBmax> from the collection of 
      {Average_NIRAbsMag_TBmax}
      
      In the second run is computed final distance modulus from: 
              mu_photo_Analytic = appMag_TBmax - <Average_NIRAbsMag_TBmax>
      In the second run the following quantities are again re-computed like in the first run: 
      (appMag_TBmax, error_appMagTBmax, mu_LCDM, sigma_muLCDM_vPec, AbsMagTBmax, error_AbsMagTBmax, chi2_dof),
      so, their values are -exactly- the same as those obtained during the first run. In the second run, the only 
      new quantities compared with run 1 are: (mu_photo_Analytic, stdDev_mu, mu_resid)

"""
0


# In[8]:


##############################################################################80


# ---

# # Automatic

# #### Get the name of this ipython notebook
# To print it in the output text files as reference.

# In[9]:


# %%javascript
# var kernel = IPython.notebook.kernel;
# var thename = window.document.getElementById("notebook_name").innerHTML;
# var command = "NotebookName = " + "'"+thename+".ipynb"+"'";
# kernel.execute(command);


# In[10]:


# if ScriptVersion == 'notebook': print(NotebookName)


# In[11]:


# Get the current date and time
import datetime 

# Read the time and date now
now = datetime.datetime.now()


# In[12]:


#- Function to convert from index to days (phase)

if TempType=='MWA': shiftNum = 24
else: shiftNum = 71
    
#- Function to convert from index to days (phase).
def index2day(index):
    day = (index-shiftNum)/2.
    return day

#- Function to convert from days (phase) to index.
def day2index(day):
    index = 2*day + shiftNum
    return index

if ScriptVersion == 'notebook':
    print 'Testing the functions (%s):'%TempType, index2day(105), ',', day2index(-11)

#--------------------------------
# Converting phases to indices

index_LowerPhase = day2index(PhaseMinTemp)
index_UpperPhase = day2index(PhaseMaxTemp)

LowerLabel = str(PhaseMinTemp)
UpperLabel = str(PhaseMaxTemp)

#--------------------------------
# OLD. Testing the functions (GP): 42 , 155
# Testing the functions (MWA): 40.5 , 2
# Testing the functions (FlatPrior): 17.0 , 49


# In[13]:


# Plotting Settings:

if BandName == 'J':
    # For overplotting the template to the fitted LC.
    # MinMagTempPlot = minimum value in the template plot
    if Method==7: MinMagTempPlot = 3.2
    else: MinMagTempPlot = -15
    
elif BandName == 'Y':
    # For overplotting the template to the fitted LC.
    # MinMagTempPlot = minimum value in the template plot
    if Method==7: MinMagTempPlot = 2
    else: MinMagTempPlot = -15.5
    
elif BandName == 'H':
    # For overplotting the template to the fitted LC.
    # MinMagTempPlot = minimum value in the template plot
    if Method==7: MinMagTempPlot = 2
    else: MinMagTempPlot = -15.5
    
elif BandName == 'K':
    # For overplotting the template to the fitted LC.
    # MinMagTempPlot = minimum value in the template plot
    if Method==7: MinMagTempPlot = 2
    else: MinMagTempPlot = -15.5


# In[14]:


#-----------------------------------------------------------------------------80


# #### Defining and creating directories

# In[15]:



# Defining the directories

DirMain_1 = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/TheTemplates/' 

DirMain = (DirMain_1+BandName+'_band/'+FilterSyst)  

if CovMat_MeanMu == True:
    CovMatLabel_MeanMu = 'CovMatMu'
elif CovMat_MeanMu == False:
    CovMatLabel_MeanMu = 'NoCovMatMu'
    
if CovMat_ErrorMu == True:
    CovMatLabel_ErrorMu = 'CovMatErrorMu'
elif CovMat_ErrorMu == False:
    CovMatLabel_ErrorMu = 'NoCovMatErrorMu'

#--------------------------
#   Dir Save Output
    
# OLD, OK:
# DirHubbleDiag = DirMain+'4_HubbleDiagram_%s/'%(TempType)
# DirSaveOutput = (DirHubbleDiag + KindOfData4HD +'/Templ_'+KindOfTemp+'_'+
#                 KindOfTempSubgroup+'/Phase'+LowerLabel+'_'+UpperLabel+'_'+
#                  residualMax_Label+'_'+
#                  chi2_dof_Max_Label+'_EBVh%r'%EBVhostMax+'_Method%r'%Method+
#                '_MinData%r'%MinNumOfDataInLC+'_vpec%r'%vpecFix+'/plots_HD/')
    
#- Force the creation of the directory to save the outputs.
#- "If the subdirectory does not exist then create it"
import os # To use command line like instructions
if not os.path.exists(DirSaveOutput): os.makedirs(DirSaveOutput)

#--------------------------

# Dir to save the plots of each fitted SNe
if NotebookToPlotOnly == False:
    if not os.path.exists(DirSaveOutput+'Plot_Fits/'): os.makedirs(DirSaveOutput+'Plot_Fits/')
    if not os.path.exists(DirSaveOutput+'Plot_Fits/AfterCutoffs_z001/'): os.makedirs(DirSaveOutput+'Plot_Fits/AfterCutoffs_z001/')
    if not os.path.exists(DirSaveOutput+'Plot_Fits/AfterCutoffs_z0/'): os.makedirs(DirSaveOutput+'Plot_Fits/AfterCutoffs_z0/')
    if not os.path.exists(DirSaveOutput+'Plot_Fits/OutCutoffs/'): os.makedirs(DirSaveOutput+'Plot_Fits/OutCutoffs/')

# Visualizing the directories
if ScriptVersion == 'notebook':
    print 'Dir to save the outputs:'
    print DirSaveOutput
    print ' '

    print 'Template to be used to fit the LC data:'
    print DirTemplate
    print ' '

    print 'Dir where the LCs in APPARENT mag are located:'
    print DirAppMag
    print ' '


# In[16]:


if NotebookToPlotOnly == False:
    
    import os # To use command line like instructions
    import glob # To read the files in my directory

    # Change the working directory where the data files are located
    os.chdir(DirAppMag)

    #--- Reading the data files in the app mag folder 
    if KindOfData4HD == 'CfA':
        list_SNe2 = glob.glob('*'+'CfA_'+BandName+'.txt')
    elif KindOfData4HD == 'CSP':
        list_SNe2 = glob.glob('*'+'CSP_'+BandName+'.txt')
    elif KindOfData4HD == 'Others':
        list_SNe2 = glob.glob('*'+'Others_'+BandName+'.txt')
    elif KindOfData4HD == 'AllSamples':
        list_SNe2 = glob.glob('*'+BandName+'.txt')

    print '# Number of -%s- SNe in the APPARENT mag list (%s):'%(
        KindOfData4HD,TempType), len(list_SNe2)

    # Number of SNe in APPARENT mag list: 132
    # Number of -AllSamples- SNe in the APPARENT mag list (MWA): 142
    # Number of -AllSamples- SNe in the APPARENT mag list (FlatPrior): 174


# In[ ]:





# In[17]:


if NotebookToPlotOnly == False:

    #     CREATE &/or READ THE NAME OF THE SNe FROM THE TEXT FILE

    # Change the working directory where the data files are located
    os.chdir(DirSaveOutput)
    # Reading the data files in that folder 
    Textfiles = glob.glob('ListSNe_AppMag*.txt')

    # Reset the variable
    list_SNe = np.array([0])

    if 'ListSNe_AppMag_Notes_.txt' in Textfiles:
        list_SNe = np.genfromtxt(DirSaveOutput+'ListSNe_AppMag_Notes_.txt', dtype=['S33', int])
        print '# Reading the existing < ListSNe_AppMag_Notes_.txt > file.'
        print "# %s SNe file names in this list and uploaded."%len(list_SNe)

    elif 'ListSNe_AppMag_.txt' in Textfiles:
        list_SNe = np.genfromtxt(DirSaveOutput+'ListSNe_AppMag_.txt', dtype=['S33', int])
        print '# Reading the existing < ListSNe_AppMag_.txt > file.'
        print "# %s SNe file names in this list and uploaded."%len(list_SNe)

    #--------------------------------------------------

    # if ListSNe_AppMag_.txt doesn't exist, then create it and read it.

    else: 
        # Create a list text file with the name of the SNe of the apparent mag LCs.
        list_SNe_file = open(DirSaveOutput+'ListSNe_AppMag_.txt', 'w')

        list_SNe_file.write('# SN list of LCs to be used to \
construct the Hubble diagram \n')
        list_SNe_file.write('# I m NOT applying any quality cutoff yet: These \
SNe are ALL those located in: \n')
        list_SNe_file.write('# %s \n'%DirAppMag)

        #------
        now = datetime.datetime.now() # Read the time and date right now
        text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd)")
        text_Date   = '# On date: %s \n'%text_timenow
        text_Author = '# Data table created by: Arturo Avelino \n'
        text_script = '# Script used: %s (version %s | last update: %s)\n'%(
                code_name, version_code, last_update)
        text_line = '#'+'-'*50 + '\n'

        list_SNe_file.write(text_line); 
        list_SNe_file.write(text_Author); list_SNe_file.write(text_Date);
        list_SNe_file.write(text_script);
        list_SNe_file.write(text_line);
        #------

        # Upload the list of repeated SNe between CfA and CSP
        DirRepeatedSNe = '/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/AndyLCComp_2018_02/MyNotes/'
        RepeatedSNeList = np.genfromtxt(DirRepeatedSNe+'RepeatedSNe_between_CfA_CSP.txt',
                                       dtype=['S10', int])

        # Upload the list of SNe to discard from the Hubble diagrams
        DirSNeWithIssues = '/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/AndyLCComp_2018_02/MyNotes/'
        SNeWithIssuesList = np.genfromtxt(DirRepeatedSNe+'SNeWithIssues.txt',
                                       dtype=['S30', 'S150'])

        for name in list_SNe2:
            snname_int =    '%s'%name[0:8]
            subsample_int = '%s'%name[-9:-6]
            # print snname_int, subsample_int

            if snname_int not in SNeWithIssuesList['f0']:
                if snname_int not in RepeatedSNeList['f0']:
                    # SNe to be NO commented:
                    list_SNe_file.write('%-32s  0 \n'%name)

                # Comment automatically the repeated SNe.
                elif snname_int in RepeatedSNeList['f0']:
                    for i1 in range(len(RepeatedSNeList['f0'])):
                        if snname_int == RepeatedSNeList['f0'][i1]:
                            if RepeatedSNeList['f1'][i1] == 1 and subsample_int =='CSP':
                                list_SNe_file.write('%-32s    1 ## Repeated. CfA selected. \n'%name)
                            elif RepeatedSNeList['f1'][i1] == 2 and subsample_int =='CfA':
                                list_SNe_file.write('%-32s    1 ## Repeated. CSP selected. \n'%name)
                            else:
                                list_SNe_file.write('%-32s  0 \n'%name)

            elif snname_int in SNeWithIssuesList['f0']:
                for i2 in range(len(SNeWithIssuesList['f0'])):
                    if snname_int == SNeWithIssuesList['f0'][i2]:
                        list_SNe_file.write('%-32s    1 ## %s \n'%(
                            name, SNeWithIssuesList['f1'][i2]))

        list_SNe_file.write(text_line)
        list_SNe_file.write('# %s SNe in total in this list. \n'%(len(list_SNe2)))
        list_SNe_file.close()

        #--------------------------------------

        # Read the list I've just created:

        list_SNe = np.genfromtxt(DirSaveOutput+'ListSNe_AppMag_.txt', dtype=['S33', int])
        print '# Reading: < ListSNe_AppMag_.txt >'

    print '# %s SNe in total in this list.'%len(list_SNe)
    print '# %s SNe masked in the list.'%sum(list_SNe['f1'])


# In[18]:


# Reading: <ListSNe_AppMag_.txt>
# Number of SNe in the list = 163
# Number of masked SNe in the list:  0


# ### Interpolating the template

# In[19]:


if NotebookToPlotOnly == False:

    from scipy.interpolate import interp1d # To interpolate

    if TempType == 'MWA':
        Template = np.genfromtxt(DirTemplate+MWATempTypeFile) 
        # Index of max abs mag in template. It doesn't need to be very precise; 
        # it is just to define the plot limit in y-axis
        indexMaxTempl = 17 # --> phase = -3.5
        indexLowerPhase = index_LowerPhase-1
        indexUpperPhase = index_UpperPhase
        #-- Interpolating the template
        MTemplInt       = interp1d(Template[indexLowerPhase:indexUpperPhase,0] , 
                                   Template[indexLowerPhase:indexUpperPhase,1] , 
                                   kind='linear')
        error_MTemplInt = interp1d(Template[indexLowerPhase:indexUpperPhase,0] , 
                                   np.sqrt(Template[indexLowerPhase:indexUpperPhase,3]) , 
                                   kind='linear')

    else:
        if Method==1 or Method==2 or Method==3: TemplateFile = 'Template_phase_mu_stdError_FromR.dat'
        elif Method==4 or Method==5 or Method==6: TemplateFile = 'Template_phase_mu_tau_FromR.dat'
        elif Method==7: TemplateFile = 'Template_phase_mu_tau_FromR_Norma.dat'

        Template = np.genfromtxt(DirTemplate+TemplateFile)

        # Index of max abs mag in template. It doesn't need to be very precise, 
        # it is just to define the plot limit in y-axis.
        indexMaxTempl = 64 # --> phase = -3.5

        #--- Interpolating the template ---

        indexLowerPhase = index_LowerPhase-1
        indexUpperPhase = index_UpperPhase


        if Method==6: # Normalizing the template
            Average_NIRAbsMag_TBmax = Template[(day2index(0)-1),1] # Abs mag at T_Bmax in the template.

            MTemplInt       = interp1d(Template[indexLowerPhase:indexUpperPhase,0] , 
                                       Template[indexLowerPhase:indexUpperPhase,1]-Average_NIRAbsMag_TBmax , 
                                       kind='linear') 

            TemplateError = np.genfromtxt(DirTemplate+'Template_phase_mu_stdError_FromR.dat')
            # Uncertainty int the abs mag at T_Bmax in the template:
            error_Average_NIRAbsMag_TBmax = TemplateError[(day2index(0)-1),2] 

        else:
            MTemplInt       = interp1d(Template[indexLowerPhase:indexUpperPhase,0] , 
                                       Template[indexLowerPhase:indexUpperPhase,1] , kind='linear')

        error_MTemplInt = interp1d(Template[indexLowerPhase:indexUpperPhase,0] , 
                                   Template[indexLowerPhase:indexUpperPhase,2] , kind='linear')
    
    if ScriptVersion == 'notebook':
        print '# Test (%s). Value of the template at phase = -1.2: (%.6f, %.6f). %s band'%(
            TempType, MTemplInt(-1.2), error_MTemplInt(-1.2), BandName)

    # Test (FlatPrior). Value of the template at phase = -1.2: (-0.066778, 0.016821). J band


# In[91]:


if NotebookToPlotOnly == False and ScriptVersion == 'notebook':
    
    # Checking the range of interpolation
    print Template[indexLowerPhase:indexUpperPhase,0]


# In[21]:


if NotebookToPlotOnly == False and ScriptVersion == 'notebook':
    # Checking the range of interpolation
    print Template[indexLowerPhase][0], Template[indexUpperPhase-1][0]
    # -8.0 31.0


# In[22]:


if NotebookToPlotOnly == False and ScriptVersion == 'notebook':
    
    # Checking the interpolation at T_Bmax
    print MTemplInt(0)
    # 0.0
    # -0.000032300235868e-05
    # -8.17052766133e-06


# In[23]:


if NotebookToPlotOnly == False and ScriptVersion == 'notebook':

    if Method==6 or Method==7:
        print 'Average_NIRAbsMag_TBmax = ', Average_NIRAbsMag_TBmax
        print 'error_Average_NIRAbsMag_TBmax', error_Average_NIRAbsMag_TBmax

# J band:
# Average_NIRAbsMag_TBmax =  -18.3805198211
# error_Average_NIRAbsMag_TBmax 0.0247540912503


# In[ ]:





# ## $\mu_{\rm \Lambda CDM}$

# In[24]:


# Inverse of the dimensionless Hubble parameter
def InvEHubblePar(z, OmM, wde):
    "Dimensionless Hubble parameter"
    InvEHubbleParInt = 1.0/(np.sqrt(OmM*(1.0+z)**3.0 + (1.0-OmM)*(1.+z)**(3.*(1.+wde))))
    return InvEHubbleParInt

# ---- The luminosity distance ----
def LumDistanceVec(z, OmM, wde, Ho):
    "Luminosity distance"
    LumDistanceVecInt = 0.
    LumDistanceVecInt = c*(1.+z)*intquad(InvEHubblePar, 0., z, args=(OmM, wde))[0]/Ho 
    return LumDistanceVecInt

# ---- Distance modulus scalar ----
def DistanceMu(z, OmM, wde, Ho):
    "Distance modulus"     
    DistanceMuInt = 5.0*np.log10(LumDistanceVec(z, OmM, wde, Ho)) + 25.0
    return DistanceMuInt

# ---- Distance modulus Vector ----
def DistanceMuVector(z, OmM, wde, Ho):
    "Distance modulus"     
    DistanceMuInt= []
    for i in range(len(z)):
        DistanceMuInt += [5.0*np.log10(LumDistanceVec(z[i], OmM, wde, Ho)) + 25.0] 
    return DistanceMuInt

#--------------------------------------------------

if ScriptVersion == 'notebook':
    ztest1 = 0.01
    print 'Checking that the functions work well:', DistanceMu(ztest1, OmMFix, wFix, HoFix)
    # Checking that the functions work well: 33.1141460988 # Ho=72
    # Checking that the functions work well: 33.0773926577 # Ho=73.24


# In[25]:


# Creation of an array of redshift vs theoretical distance modulus (for Pete Challis)
""" 
z1 = np.linspace(0.01, 0.7, 501)
DistMu_array = DistanceMuVector(z1, OmMFix, wFix, HoFix)

DirSaveOutput = '/Users/arturo/Dropbox/Research/Articulos/0_LCDM/'
textfile = open(DirSaveOutput+'DistanceModulus.txt', 'w')
textfile.write('# Cosmology: flat Universe, Om_matter = %r, Om_Lambda = %r, Ho = %r km/s/Mpc \n'%(OmMFix, 1-OmMFix, HoFix))
textfile.write('# z        Distance modulus  \n')

for i in range(len(z1)):
    textfile.write('%.6f   %.6f \n'%(z1[i], DistMu_array[i]))
    
textfile.close()
"""
0


# ### Uncertainty in distance modulus due to the peculiar-velocity uncertainty

# In[26]:


# sigma^2_mu from the peculiar velocity uncertainty
# This function is used to determine in the sections "Intrinsic dispersion" and "Optical RMS", to
# determine the intrinsic dispersion.

def sigma2_pec(zcmb, err_zcmb, vpec):
    sigma2_pecInt = ((5/(zcmb*np.log(10)))*np.sqrt((vpec/c)**2 + err_zcmb**2))**2
    return sigma2_pecInt

if ScriptVersion == 'notebook':
    # Test
    print sigma2_pec(0.0109942726, 0.0010420420, 150)
    # 0.0521250373543


# In[27]:


# Test: Uncertainty on distance mu due to uncertainty in pec velocity.
if ScriptVersion == 'notebook':
    print np.sqrt(sigma2_pec(0.005073733274739549,  0.00020317696892145927 , 0))


# #### Function to identify string or number

# In[28]:


# Function to identify if a string is an integer number or a letter.
# This will be used in the dictionary construction to properly read some SN names.

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

if ScriptVersion == 'notebook':
    # Tests
    print is_number('5'), is_number('e')
    # True False


# ### Cepheid distances

# In[29]:


if ScriptVersion == 'notebook':
    DirSNeWithCepheid = '/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/AndyLCComp_2018_02/MyNotes/'
else: DirSNeWithCepheid = sys.argv[22]

# From the Cepheid SNe list the only part that I use in the entire
# notebook is the first column, i.e., the SN Name column.
ListSNeCepheid = np.genfromtxt(DirSNeWithCepheid+'SNeWithCepheidDistances.txt', dtype=['S10', 
                                                float,float,float,float,float,float])
if ScriptVersion == 'notebook':
    print "# %s SNe with Cepheid distances in Andy's compilation. "%len(ListSNeCepheid['f0'])
    ListSNeCepheid['f0']


# In[30]:


"""
# 13 SNe with Cepheid distances in Andy's compilation. 
Out[23]:
array(['sn1981B', 'sn1998aq', 'sn2001el', 'sn2002fk', 'sn2003du',
       'sn2005cf', 'sn2007af', 'sn2007sr', 'sn2009ig', 'sn2011by',
       'sn2011fe', 'sn2012cg', 'sn2012fr'],
      dtype='|S10')
"""
0


# ### Function to determine $c z_{\rm cmb}$ from Cepheid $\mu$
# 
# These functions are NOT used during any part of the computations in the notebook; I have written here just as part of my library of functions.

# In[31]:


def cz_CMB_Cepheid(mu, err_mu, Ho, err_Ho):
    
    # Determine cz_cmb, z_cmb:
    cz_int =  Ho * 10**((mu-25)/5)  # cz_cmb
    zcmb_int = cz_int/c  # z_cmb
    
    # Determine error_cz_cmb, error_z_cmb:
    error_cz_int = cz_int * np.sqrt((err_Ho/Ho)**2 + (np.log(10)*err_mu/5)**2) # error_cz_cmb
    error_zcmb_int = error_cz_int/c  # error_z_cmb
    
    # (z_cmb, error_z_cmb, cz_cmb, error_cz_cmb)
    return zcmb_int, error_zcmb_int, cz_int, error_cz_int

    
if ScriptVersion == 'notebook':
    # TEST:
    print '# sn2003du (data from Riess+16):'
    print '#',cz_CMB_Cepheid(32.919, 0.063, 73.24, 1.74)

    print ''
    print '# sn2005cf (data from Riess+16):'
    print '#',cz_CMB_Cepheid(32.263, 0.102, 73.24, 1.74)
    
    print ''
    print '# sn2003hv (data from Riess+16):'
    print '#', cz_CMB_Cepheid(31.514,  0.150, 73.24, 0)

    # sn2003du (data from Riess+16):
    # (0.009369742003029419, 0.00035135265754374572, 2808.977985914033, 105.33287682987176)

    # sn2005cf (data from Riess+16):
    # (0.006926720004322207, 0.00036461514433908711, 2076.578415973525, 109.3088703454397)


# In[32]:


if ScriptVersion == 'notebook': print cz_CMB_Cepheid(31.587, 0.070, 73.24, 1.74)


# In[ ]:





# In[33]:


if NotebookToPlotOnly == False and ScriptVersion == 'notebook':
    
    HoNew = 72.0; err_HoNew = 1.74

    # err_HoNew = 1.74 comes from assuming the same uncertainty than the 
    # case of  Ho=73.24 +/- 1.74 km/s/Mpc  (Riess+2016).
    """ 
    print '# sn2003du (data from Riess+16):'
    print '#',cz_CMB_Cepheid(32.919, 0.063, HoNew, err_HoNew)
    """   
0


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[34]:


##############################################################################80


# ### Preeliminary writting of  text files

# In[35]:


if NotebookToPlotOnly == False:

    # Creation of the datatable text file of (z_CMB, dist_mu, dist_mu_error, delta_mu) 

    InfoSN_NumLinesSkip = 8 # Number info lines in LC files

    # Mean Absolute magnitude determined from histogram of 'appMagTmax_s - mu_s'?:
    # If so it is generated the final 'DistanceMu_Good_AfterCutoffs_Main_.txt' text file, otherwise
    # generate temporal text files.
    if  AbsMagFromHisto == True: textPrefix = ''; underline = ''
    else: textPrefix = 'TMP'; underline = '_' # 'TMP' stands for 'temporal' file.

    textfileMain = open(DirSaveOutput+'DistanceMu_Good_AfterCutoffs_Main_%s%s.txt'%(textPrefix, underline), 'w')
    textfile = open(DirSaveOutput+textPrefix+underline+'DistanceMu_All_BeforeCutoffs_.txt', 'w')
    textfil3 = open(DirSaveOutput+textPrefix+underline+'DistanceMu_Good_AfterCutoffs_z001_.txt', 'w')
    textfil4 = open(DirSaveOutput+textPrefix+underline+'RMS_tmp.txt', 'w')
    textfil5 = open(DirSaveOutput+textPrefix+underline+'DistanceMu_SNe_OutCutoffs_.txt', 'w')
    textfil6 = open(DirSaveOutput+textPrefix+underline+'DistanceMu_Good_AfterCutoffs_Main_NamesOnly_.txt', 'w')
    textfil7 = open(DirSaveOutput+textPrefix+underline+'SNe_all_before_any_cutoff_.txt', 'w')
    
    #------------------------------------------------------------------------------------

    #       Defining the text

    text_001 = '#    %s band  (template method) \n'%BandName
    text01b = '# MAIN DATA TABLE. Used to make further computations and plots \n'

    now = datetime.datetime.now() # Read the time and date right now
    text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
    text_Date   = '# On date: %s \n'%text_timenow
    text_Author = '# Data table created by: Arturo Avelino \n'
    text_script = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
    text_line = '#'+'-'*70 + '\n'

    text02 = '# This datatable is for information only. It is NOT used to compute or plot something else. \n'

    text01 = '# Apparent mag data used to construct the Hubble diagram: %s \n'%KindOfData4HD
    text1 = '# Template used: \n'
    text2 = '# %s \n'%DirTemplate[78:]
    text3 = '# Cosmology used to compute residuals: (Omega_M=%r, Omega_L=%r, w=%r, Ho=%r) \n'%(OmMFix, OmLFix, wFix, HoFix)

    text3_0_0 = '# 0.01 < z_cmb < %r \n'%zcmbUpperLim
    text3_0_1 = '# Cutoffs: z_cmb<%r | %r<dm15<%r | %r<EBVhost<%r | EBV_mw<%r | chi2_dof<%r | Residual<%r \n'%(zcmbUpperLim,
        dm15LowerLim,  dm15UpperLim, EBVhostMin, EBVhostMax, EBVMWLim, chi2_dof_Max, residualMax)
    text3_3 = '# Phase range used of the template: (%r, %r) days \n'%(PhaseMinTemp, PhaseMaxTemp)
    text3_4 = '# Minimal number of data in LC to be considered for the Hubble diagram: %r \n'%MinNumOfDataInLC
    text3_5 = '# Assumed uncertainty in the peculiar velocity: %r km/s \n'%vpecFix
    text4 = '# SN name                        z_CMB         \
error_zcmb        mu          error_mu       mu_residual    \
chi2_dof     (1=CfA,2=CSP)   appMag_TBmax  error_appMagTBmax  \
mu_LCDM  sigma_muLCDM_vPec AbsMagTBmax error_AbsMagTBmax  \
TBmax       Error_TBmax  PhaseB zhelio        err_zhelio   \
dm15         error_dm15      EBVhost      error_EBVhost   \
EBV_MW   error_EBV_MW  Alamb       err_Alamb      R_F          \
mu_Snoopy     error_muSnoopy   TBmax      Error_TBmax  \
appMag_TBmax  error_appMagTBmax  Notes   \n' 

    text4_1 = "# SN name                        z_CMB         \
error_zcmb          mu_LCDM  sigma_muLCDM_vPec     TBmax       \
Error_TBmax     zhelio        err_zhelio   dm15         \
error_dm15      EBVhost      error_EBVhost   EBV_MW   error_EBV_MW  \
Alamb       err_Alamb      R_F           mu_Snoopy     \
error_muSnoopy     TBmax           Error_TBmax #  Notes \n"
    
    # text4 = '# SN name                        z_CMB     error_zcmb        mu          error_mu      \
    # mu_residual     chi2_dof   (1=CfA,2=CSP)   appMag_TBmax   error_appMagTBmax   mu_LCDM   sigma_muLCDM_vPec  
    # dm15       error_dm15        EBVhost    error_EBVhost      \
    # EBV_MW     error_EBV_MW   mu_Snoopy   error_mu_Snoopy    TBmax     Error_TBmax    Notes \n' 

    #------------------------------------------------------------------ 

    textMethod = '# Method to determine the distance modulus?: %s \n'%Method 
    text_CorrMatrix = '# Use the correlation (exponential kernel) matrix in the template to compute the distance moduli?: (answer in the lines below) \n'
    text_CovMat_MeanMu  = '# CovMat_MeanMu = %s  # correlation matrix to compute the -mean-  distance moduli? \n'%CovMat_MeanMu
    text_CovMat_ErrorMu = '# CovMat_ErrorMu = %s # correlation matrix to compute the -error- distance moduli? \n'%CovMat_ErrorMu
    text_lkernel = '# l_kernel = %r. Lenght scale of the kernel in the covariance matrix for the residual. \n'%(l_kern)
    text_Template = "# Template = %s \n"%TemplateFile
    text_CovMatrix_MeanMu_1  = '# Total CovMat to compute Mean_mu:  CovMatrix_MeanMu = CovMatCorrData_4MeanMu  + Cov_appmag \n'
    text_CovMatrix_MeanMu_2  = '# Total CovMat to compute Mean_mu:  CovMatrix_MeanMu = CovMatCorrData_4MeanMu  + Cov_appmag + Cov_pecvel \n'
    text_CovMatrix_ErrorMu_1 = '# Total CovMat to compute Error_mu: CovMatrix_ErrorMu= CovMatCorrData_4ErrorMu + Cov_appmag + Cov_pecvel \n'
    text_CovMatrix_ErrorMu_2 = '# Total CovMat to compute Error_mu: CovMatrix_ErrorMu= CovMatCorrData_4ErrorMu + Cov_appmag \n'
    text_sigma2mu_1 = '# Variance mu: sigma2_mu = (1/Denominator_ErrorMu) \n'
    text_sigma2mu_2 = '# Variance mu: sigma2_mu = (1/Denominator_ErrorMu) +  sigma2_pec(zcmbInt, err_zcmb, vpecFix) \n'
    # text_IntrinsicDisp_1 = '# Intrinsic dispersion for the case (0 < z < %s) and used to obtain the total photometric \
    # distance modulus uncertainty:\n'%zcmbUpperLim
    # text_IntrinsicDisp_2 = '# Intrinsic dispersion = %s \n'%IntrinsicDisp

    #=====================================================================================

    #       Writting the text

    # 'DistanceMu_Good_AfterCutoffs_Main_.txt'

    textfileMain.write(text_001); textfileMain.write(text01b)

    textfileMain.write(text_line); 
    textfileMain.write(text_Author); textfileMain.write(text_Date); textfileMain.write(text_script);
    textfileMain.write(text_line); 

    textfileMain.write(text01); textfileMain.write(text1); textfileMain.write(text2); textfileMain.write(text3)
    textfileMain.write(text3_0_1)
    textfileMain.write(text3_3); textfileMain.write(text3_4); textfileMain.write(text3_5);
    # if CovMat_MeanMu == True: textfileMain.write(text_lkernel)
    textfileMain.write(text_line); 
    textfileMain.write(textMethod); 
    textfileMain.write(text_CorrMatrix); 
    textfileMain.write(text_CovMat_MeanMu); textfileMain.write(text_CovMat_ErrorMu); textfileMain.write(text_lkernel)
    textfileMain.write(text_Template)
    if Method==1:
        textfileMain.write(text_CovMatrix_MeanMu_1); textfileMain.write(text_CovMatrix_ErrorMu_1); 
        textfileMain.write(text_sigma2mu_1); 
    elif  Method==2:
        textfileMain.write(text_CovMatrix_MeanMu_1); textfileMain.write(text_CovMatrix_ErrorMu_2); 
        textfileMain.write(text_sigma2mu_2); 
    elif  Method==3:
        textfileMain.write(text_CovMatrix_MeanMu_2); textfileMain.write(text_CovMatrix_ErrorMu_1); 
        textfileMain.write(text_sigma2mu_1); 
    elif Method==4:
        textfileMain.write(text_CovMatrix_MeanMu_1); textfileMain.write(text_CovMatrix_ErrorMu_2); 
        textfileMain.write(text_sigma2mu_1);
    elif Method==5 or Method==6:
        textfileMain.write(text_CovMatrix_MeanMu_1); textfileMain.write(text_CovMatrix_ErrorMu_1); 
        textfileMain.write(text_sigma2mu_1); 
    elif Method==7:
        textfileMain.write('#       Distance modulus determination: \n')
        textfileMain.write('# CovMatrix_MeanMu = CovMatCorrData_4MeanMu + Cov_appmag \n')
        textfileMain.write('# appMag_TBmax = Numerator_MeanMu/Denominator_MeanMu \n')
        textfileMain.write('# mu_photo_Analytic = appMag_TBmax - <Average_NIRAbsMag_TBmax> \n')
        textfileMain.write('#       Uncertainty in the distance modulus determination: \n')
        textfileMain.write('# CovMatrix_ErrorMu = CovMatCorrData_4ErrorMu + Cov_appmag \n')
        textfileMain.write('# error_appMag_TBmax = np.sqrt(1/Denominator_ErrorMu)  \n')
        textfileMain.write('# sigma2_mu = (1/Denominator_ErrorMu) \n')
    if Method==6 or Method==7: 
        textfileMain.write('# Using a NORMALIZED template. Distance mu = m_TBmax - M_TBmax. \n')
        textfileMain.write('#       Assumed values for Average_NIRAbsMag_TBmax and error_Average_NIRAbsMag_TBmax: \n')
        textfileMain.write('# <Average_NIRAbsMag_TBmax> = %r. %s \n'%(Average_NIRAbsMag_TBmax, Notes1))
        textfileMain.write('# error_<Average_NIRAbsMag_TBmax> = %r. %s \n'%(error_Average_NIRAbsMag_TBmax, Notes1))
        # textfileMain.write(text_IntrinsicDisp_1); textfileMain.write(text_IntrinsicDisp_2);

    textfileMain.write(text_line); 
    textfileMain.write(text4)

    #----------------------------

    # 'DistanceMu_Good_AfterCutoffs_Main_NamesOnly_.txt'

    textfil6.write('# SN name \n')

    #----------------------------

    # 'DistanceMu_All_BeforeCutoffs_.txt'

    textfile.write(text_001); textfile.write(text02)
    textfile.write(text_line); 
    textfile.write(text_Author); textfile.write(text_Date); textfile.write(text_script);
    textfile.write(text_line); 
    textfile.write(text01); textfile.write(text1); textfile.write(text2); textfile.write(text3)
    textfile.write('# No restrictions on z_cmb, dm15, EBVhost, EBV_mw, chi2_dof, residuals \n')
    textfile.write(text3_3); textfile.write(text3_4); textfile.write(text3_5)
    # if CovMat_MeanMu == True: textfile.write(text_lkernel)
    textfile.write(text_line); 
    textfile.write(textMethod)
    textfile.write(text_CorrMatrix); 
    textfile.write(text_CovMat_MeanMu); textfile.write(text_CovMat_ErrorMu); textfile.write(text_lkernel)
    textfile.write(text_Template)
    if Method==1:
        textfile.write(text_CovMatrix_MeanMu_1); textfile.write(text_CovMatrix_ErrorMu_1); 
        textfile.write(text_sigma2mu_1); 
    elif  Method==2:
        textfile.write(text_CovMatrix_MeanMu_1); textfile.write(text_CovMatrix_ErrorMu_2); 
        textfile.write(text_sigma2mu_2); 
    elif  Method==3:
        textfile.write(text_CovMatrix_MeanMu_2); textfile.write(text_CovMatrix_ErrorMu_1); 
        textfile.write(text_sigma2mu_1); 
    elif Method==4:
        textfile.write(text_CovMatrix_MeanMu_1); textfile.write(text_CovMatrix_ErrorMu_2); 
        textfile.write(text_sigma2mu_1);
    elif Method==5 or Method==6:
        textfile.write(text_CovMatrix_MeanMu_1); textfile.write(text_CovMatrix_ErrorMu_1); 
        textfile.write(text_sigma2mu_1); 
    elif Method==7:
        textfile.write('#       Distance modulus determination: \n')
        textfile.write('# CovMatrix_MeanMu = CovMatCorrData_4MeanMu + Cov_appmag \n')
        textfile.write('# appMag_TBmax = Numerator_MeanMu/Denominator_MeanMu \n')
        textfile.write('# mu_photo_Analytic = appMag_TBmax - <Average_NIRAbsMag_TBmax> \n')
        textfile.write('#       Uncertainty in the distance modulus determination: \n')
        textfile.write('# CovMatrix_ErrorMu = CovMatCorrData_4ErrorMu + Cov_appmag \n')
        textfile.write('# error_appMag_TBmax = np.sqrt(1/Denominator_ErrorMu)  \n')
        textfile.write('# sigma2_mu = (1/Denominator_ErrorMu) \n')
    if Method==6 or Method==7: 
        textfile.write('# Using a NORMALIZED template. Distance mu = m_TBmax - M_TBmax. \n')
        textfile.write('#       Assumed values for Average_NIRAbsMag_TBmax and error_Average_NIRAbsMag_TBmax: \n')
        textfile.write('# <Average_NIRAbsMag_TBmax> = %r. %s \n'%(Average_NIRAbsMag_TBmax, Notes1))
        textfile.write('# error_<Average_NIRAbsMag_TBmax> = %r. %s \n'%(error_Average_NIRAbsMag_TBmax, Notes1))
        # textfile.write(text_IntrinsicDisp_1); textfile.write(text_IntrinsicDisp_2);

    textfile.write(text_line); 
    textfile.write(text4)

    #----------------------------

    textfil3.write(text02); textfil3.write(text01); textfil3.write(text1); textfil3.write(text2) 
    textfil3.write(text3)
    textfil3.write(text3_0_0); textfil3.write(text3_0_1)
    textfil3.write(text3_3); textfil3.write(text3_4); textfil3.write(text3_5)
    # if CovMat_MeanMu == True: textfil3.write(text_lkernel)
    textfil3.write(text_line); 
    textfil3.write(textMethod); 
    textfil3.write(text_CorrMatrix); 
    textfil3.write(text_CovMat_MeanMu); textfil3.write(text_CovMat_ErrorMu); textfil3.write(text_lkernel)
    textfil3.write(text_Template)
    if Method==1:
        textfil3.write(text_CovMatrix_MeanMu_1); textfil3.write(text_CovMatrix_ErrorMu_1); 
        textfil3.write(text_sigma2mu_1); 
    elif  Method==2:
        textfil3.write(text_CovMatrix_MeanMu_1); textfil3.write(text_CovMatrix_ErrorMu_2); 
        textfil3.write(text_sigma2mu_2); 
    elif  Method==3:
        textfil3.write(text_CovMatrix_MeanMu_2); textfil3.write(text_CovMatrix_ErrorMu_1); 
        textfil3.write(text_sigma2mu_1); 
    elif Method==4:
        textfil3.write(text_CovMatrix_MeanMu_1); textfil3.write(text_CovMatrix_ErrorMu_2); 
        textfil3.write(text_sigma2mu_1); 
    textfil3.write(text_line); 
    textfil3.write(text4)

    #----------------------------

    textfil4.write(text_line), textfil4.write('# <- Template range -> \n')
    textfil4.write('# Phase_min  Phase_max   RMS    # SNe   RMS(z>0.01)  # SNe  \n')

    #----------------------------

    textfil5.write('# SNe that are OUTSIDE the following cutoffs-criteria: \n')  
    textfil5.write(text3_0_0); textfil5.write(text3_0_1)
    textfil5.write(text02); textfil5.write(text01); textfil5.write(text1); textfil5.write(text2)
    textfil5.write(text3); 
    textfil5.write(text3_3); textfil5.write(text3_4); textfil5.write(text3_5)
    # if CovMat_MeanMu == True: textfil5.write(text_lkernel)
    textfil5.write(text_line); 
    textfil5.write(textMethod)
    textfil5.write(text_CorrMatrix); 
    textfil5.write(text_CovMat_MeanMu); textfil5.write(text_CovMat_ErrorMu); textfil5.write(text_lkernel)
    textfil5.write(text_Template)
    if Method==1:
        textfil5.write(text_CovMatrix_MeanMu_1); textfil5.write(text_CovMatrix_ErrorMu_1); 
        textfil5.write(text_sigma2mu_1); 
    elif  Method==2:
        textfil5.write(text_CovMatrix_MeanMu_1); textfil5.write(text_CovMatrix_ErrorMu_2); 
        textfil5.write(text_sigma2mu_2); 
    elif  Method==3:
        textfil5.write(text_CovMatrix_MeanMu_2); textfil5.write(text_CovMatrix_ErrorMu_1); 
        textfil5.write(text_sigma2mu_1); 
    elif Method==4:
        textfil5.write(text_CovMatrix_MeanMu_1); textfil5.write(text_CovMatrix_ErrorMu_2); 
        textfil5.write(text_sigma2mu_1); 
    textfil5.write(text_line); 
    textfil5.write(text4)
    
    #----------------------------
    
    # 'SNe_all_before_any_cutoff_.txt'
    
    # Copy/paste of section for 'DistanceMu_Good_AfterCutoffs_Main_.txt'
    textfil7.write(text_001); textfil7.write(text01b)

    textfil7.write(text_line); 
    textfil7.write(text_Author); textfil7.write(text_Date); textfil7.write(text_script);
    textfil7.write(text_line); 

    textfil7.write(text01); textfil7.write(text1); textfil7.write(text2); textfil7.write(text3)
    # textfil7.write(text3_0_1)
    textfil7.write(text3_3); 
    textfil7.write(text3_4); textfil7.write(text3_5);
    # if CovMat_MeanMu == True: textfil7.write(text_lkernel)
    textfil7.write(text_line); 
    textfil7.write(textMethod); 
    textfil7.write(text_CorrMatrix); 
    textfil7.write(text_CovMat_MeanMu); textfil7.write(text_CovMat_ErrorMu); textfil7.write(text_lkernel)
    textfil7.write(text_Template)
    if Method==1:
        textfil7.write(text_CovMatrix_MeanMu_1); textfil7.write(text_CovMatrix_ErrorMu_1); 
        textfil7.write(text_sigma2mu_1); 
    elif  Method==2:
        textfil7.write(text_CovMatrix_MeanMu_1); textfil7.write(text_CovMatrix_ErrorMu_2); 
        textfil7.write(text_sigma2mu_2); 
    elif  Method==3:
        textfil7.write(text_CovMatrix_MeanMu_2); textfil7.write(text_CovMatrix_ErrorMu_1); 
        textfil7.write(text_sigma2mu_1); 
    elif Method==4:
        textfil7.write(text_CovMatrix_MeanMu_1); textfil7.write(text_CovMatrix_ErrorMu_2); 
        textfil7.write(text_sigma2mu_1);
    elif Method==5 or Method==6:
        textfil7.write(text_CovMatrix_MeanMu_1); textfil7.write(text_CovMatrix_ErrorMu_1); 
        textfil7.write(text_sigma2mu_1); 
    elif Method==7:
        textfil7.write('#       Distance modulus determination: \n')
        textfil7.write('# CovMatrix_MeanMu = CovMatCorrData_4MeanMu + Cov_appmag \n')
        textfil7.write('# appMag_TBmax = Numerator_MeanMu/Denominator_MeanMu \n')
        textfil7.write('# mu_photo_Analytic = appMag_TBmax - <Average_NIRAbsMag_TBmax> \n')
        textfil7.write('#       Uncertainty in the distance modulus determination: \n')
        textfil7.write('# CovMatrix_ErrorMu = CovMatCorrData_4ErrorMu + Cov_appmag \n')
        textfil7.write('# error_appMag_TBmax = np.sqrt(1/Denominator_ErrorMu)  \n')
        textfil7.write('# sigma2_mu = (1/Denominator_ErrorMu) \n')
    if Method==6 or Method==7: 
        textfil7.write('# Using a NORMALIZED template. Distance mu = m_TBmax - M_TBmax. \n')
        textfil7.write('#       Assumed values for Average_NIRAbsMag_TBmax and error_Average_NIRAbsMag_TBmax: \n')
        textfil7.write('# <Average_NIRAbsMag_TBmax> = %r. %s \n'%(Average_NIRAbsMag_TBmax, Notes1))
        textfil7.write('# error_<Average_NIRAbsMag_TBmax> = %r. %s \n'%(error_Average_NIRAbsMag_TBmax, Notes1))
        # textfil7.write(text_IntrinsicDisp_1); textfil7.write(text_IntrinsicDisp_2);

    textfil7.write(text_line); 
    textfil7.write(text4_1) # Modified
    


# In[36]:


""" 
textfile.close()
textfileMain.close()
textfil3.close()
textfil4.close()
textfil5.close()
textfil6.close()
textfil7.close()
 """
0


# In[37]:


##############################################################################80


# # Main loop
# 
# ### Moving the template up or down to match the LC data. The amount of moving up or down corresponds to the photometric distance modulus.
# 
# ####  In this loop it is used ALL the apparent mag LCs located in "~/1_AllData_InitialFit/AppMag/" through the text file "ListSNe_AppMag_.txt" created and read above.

# In[38]:


#    PLOT OPTIONS
# Resolution in the Hubble diagram plots
ResolutionPlot_HD = 100  # dpi # ORIGINAL
fontsizePlot = 13

# Settings error bars:
MyPointSize = 6  
MyLineWidth = 1  
MyCapeSize = 2
fmt_int = '.'


# In[39]:


if NotebookToPlotOnly == False:
    
    # Initializing some quantities
    countGoodSNe = 0; countGoodSNe_z001 = 0; # Counters
    countCfA = 0;     countCSP = 0;     countOthers=0;  countSNe = 0 # Counters
    countCfA_z001= 0; countCSP_z001= 0; countOthers_z001=0 # Counters
    countNoCommented = 0; countCommented = 0;
    countBadSNe = 0; 

    # To be used to compute the RMS:
    mu_resid_array = []; mu_resid_z001_array = []; 
    chi2_list = []; chi2_z001_list = []; 
    SNeName_list = []; SNeName_z001_list = []; 

    # The phase range of the template to be used to determine the distance modulus
    lowerPhase = Template[indexLowerPhase][0]
    upperPhase = Template[indexUpperPhase-1][0]

    ############ THE LOOP ############################

    for j in range(len(list_SNe)): # loop over the SNe
    # for j in range(7):
    # for j in range(10):
    # for j in range(20, 40):

        # Upload the apparent magnitude data for a given SN
        magData = np.genfromtxt(DirAppMag+list_SNe[j][0])
        name = list_SNe[j][0]
        mask = list_SNe[j][1]

        zcmbInt = magData[0,0]
        dm15Int = magData[1,0]
        EBV_MW = magData[3,0]
        EBVhost = magData[4,0]

        err_zcmb = magData[0,1]
        err_dm15 = magData[1,1]
        err_EBVMW = magData[3,1]
        err_EBVhost = magData[4,1]

        # These quantities are not used to compute anything, they are just to be
        # written int the output file
        muSnoopy     = magData[2,0]
        err_muSnoopy = magData[2,1]
        TBmax      = magData[5,0]
        err_TBmax  = magData[5,1]
        zhelio     = magData[0,2]
        err_zhelio = magData[0,3]
        Alamb      = magData[3,2] 
        err_Alamb  = magData[3,3] 
        R_F        = magData[3,4]

        #-------- Peculiar velocity uncertainty -----------------        
        # Create the variable "snName" containing the first 8 
        # (or 7) letters of the SNe file name
        # I use "snName" to compare with the SNe names in 
        # 'SNeWithCepheidDistances.txt', so that
        # I will not compute a peculiar velocity uncertainty for those SNe.
        try:
            if   name[7] == '_': snName = name[:7] # To read correctly, e.g., "sn2011B_"
            elif name[7] != '_':
                # To read correctly, e.g., "snf20080514-002":
                if is_number(name[7]): snName = name[:15] 
                else: snName = name[:8]  # To read correctly, e.g., "sn1998bu" 
        except: snName = name[:6]  # To read correctly, e.g., "sn2011B"

        sigma_pecInt = 0 # Reset its value:
        # If this SNe has Cepheid distance, then use vpecFix=0 for it: 
        if snName in ListSNeCepheid['f0']: 
            sigma_pecInt = np.sqrt(sigma2_pec(zcmbInt, err_zcmb, 0))
        else: 
            sigma_pecInt = np.sqrt(sigma2_pec(zcmbInt, err_zcmb, vpecFix))

        #---------------------------------------
        # Write down the parameters that are independent of the band
        # needed to quantify the cutoffs

        # Theoretical distance mu from LCDM
        mu_LCDM_2 = DistanceMu(magData[0,0], OmMFix, wFix, HoFix)
                
        textfil7.write('%-30s   %.10f  %.10f        %.10f  %.10f  %13.5f  %13.5f   %.10f  %.8f  %.10f  %.10f  %13.10f  %.10f  %10.7f  %.7f  %.10f  %.10f  %.10f  %.10f  %.10f    %13.5f  %13.5f  # \n'%(
            list_SNe[j][0][:-4],
            zcmbInt, err_zcmb,
            mu_LCDM_2, sigma_pecInt,
            TBmax, err_TBmax,
            zhelio, err_zhelio, dm15Int, err_dm15,
            EBVhost, err_EBVhost, EBV_MW, err_EBVMW, Alamb, err_Alamb, R_F,
            muSnoopy, err_muSnoopy,
            TBmax, err_TBmax ))


        #---------------------------------------
        
        # If there are at least 3 photometric points in the 
        # app mag LC text file, then compute. 
        if len(magData) >= InfoSN_NumLinesSkip + MinNumOfDataInLC:

            #---------------------------
            # Find out if there are at least one datum inside the 
            # interpolated phase range of the template
            NumDataInRange = 0
            
            # Loop over all the phases for a given SNe
            for i in range(InfoSN_NumLinesSkip, len(magData)): 
                if magData[i,0] >= lowerPhase and magData[i,0] <= upperPhase:
                    NumDataInRange = NumDataInRange + 1

            #---------------------------
            # Fit the LCs with at least one data point within 
            # the phase range defined for the template.
            if NumDataInRange > 0:

                # Analytical minimization of chi2 to find the best 
                # estimate for distance modulus
                # and computing the uncertainty of distance 
                # modulus (analytic expression)

                # Resetting some quantities
                Numerator_MeanMu  = 0; Denominator_MeanMu= 0
                Numerator_ErrorMu = 0; Denominator_ErrorMu= 0
                mu_photo_Analytic=0; sigma2_mu=0; stdDev_mu=0
                appMag_TBmax=0; error_appMag_TBmax=0;
                AbsMag_s=0; error_AbsMag_s=0;
                mu_LCDM=0;

                CovMatrix_MeanMu = 0; InvCovMatrix_MeanMu = 0;
                CovMatrix_ErrorMu = 0; InvCovMatrix_ErrorMu =0;
                OffsetToMatchAppMagData = 0; appMag_TBmax = 0;

                #---------------------------    
                # Using the covariance matrix of the photometry to determine the distance modulus:

                # Create arrays for the data that is inside the interpolation 
                # range & within the range (lowerPhase, upperPhase).
                mag_array = [] # apparent magnitude data
                MagTemp_array = [] # template mean value (absolute magnitude)
                error2_appmag = [] # app magnitud errors
                error_temp = [] # population standard deviation of the template at different phases.
                phases = [] # phases inside the interpolation range & within the range (lowerPhase, upperPhase).
                phaseInt = 0
                countGoodPhases = 0; countLowerPhase = 0; 

                # Loop over all the photometry for a given SNe
                for i in range(InfoSN_NumLinesSkip, len(magData)): 
                    # Ignore data outside the phase interpolation range:
                    if magData[i,0] >= lowerPhase: 
                        if magData[i,0] <= upperPhase:
                            phaseInt = magData[i,0]
                            mag_array += [magData[i,3]]
                            error2_appmag += [magData[i,4]**2]
                            MagTemp_array += [MTemplInt(phaseInt)]
                            error_temp += [error_MTemplInt(phaseInt)]
                            phases += [phaseInt]
                            countGoodPhases = countGoodPhases + 1
                    else: countLowerPhase = countLowerPhase + 1

                error_temp = np.array(error_temp)

                #-------- Some covariance matrices --------

                # Diagonal variance matrix of the photometric data.
                Cov_appmag = np.diag(error2_appmag)

                # Square covariance matrix of the peculiar velocity uncertainty
                Cov_pecvel = np.ones((countGoodPhases,countGoodPhases))*(sigma_pecInt**2)


                #--------- MEAN photometric distance modulus -------------

                # Covariance matrix of the light-curve photometric data:
                CovMatCorrData_4MeanMu = np.zeros(shape=(countGoodPhases,countGoodPhases)) # Empty array
                for i in range(countGoodPhases):
                    if CovMat_MeanMu == True:
                        for k in range(countGoodPhases):
                            CovMatCorrData_4MeanMu[i,k] = error_temp[i]*error_temp[k]*np.exp( -(((phases[i] - phases[k])/l_kern)**2)/2 ) 
                    elif CovMat_MeanMu == False:
                        CovMatCorrData_4MeanMu[i,i] = error_temp[i]*error_temp[i]

                # Total covariance matrix
                if Method==1 or Method==2 or Method==4 or Method==5 or Method==6 or Method==7: 
                    CovMatrix_MeanMu = CovMatCorrData_4MeanMu + Cov_appmag
                elif Method==3: CovMatrix_MeanMu = CovMatCorrData_4MeanMu + Cov_appmag + Cov_pecvel

                # Inverse of the total covariance matrix
                InvCovMatrix_MeanMu = np.linalg.inv(CovMatrix_MeanMu)

                #---Computing the best estimated photometric distance modulus ---
                # Loop over all the photometry for a given SNe
                # Ignore data outside the interpolation range (lowerPhase, upperPhase):
                for i in range(countGoodPhases): 
                    # Ignore data outside the interpolation range:
                    Numerator_MeanMu = Numerator_MeanMu + (mag_array[i] - MagTemp_array[i])*np.sum(InvCovMatrix_MeanMu[i,:])
                    Denominator_MeanMu = Denominator_MeanMu + np.sum(InvCovMatrix_MeanMu[i,:])

                # Best estimated distance modulus:
                if Method==6 or Method==7: 
                    appMag_TBmax = Numerator_MeanMu/Denominator_MeanMu
                    mu_photo_Analytic = appMag_TBmax - Average_NIRAbsMag_TBmax

                else: 
                    appMag_TBmax = 0
                    mu_photo_Analytic = Numerator_MeanMu/Denominator_MeanMu   

                #--------- UNCERTAINTY in the photometric distance modulus -------------

                # Covariance matrix of the light-curve photometric data to compute ERROR_MU
                CovMatCorrData_4ErrorMu = np.zeros(shape=(countGoodPhases,countGoodPhases)) # Empty array
                for i in range(countGoodPhases):
                    if CovMat_ErrorMu == True:
                        for k in range(countGoodPhases):
                            CovMatCorrData_4ErrorMu[i,k] = error_temp[i]*error_temp[k]*np.exp( -(((phases[i] - phases[k])/l_kern)**2)/2 ) 
                    elif CovMat_ErrorMu == False:
                        CovMatCorrData_4ErrorMu[i,i] = error_temp[i]*error_temp[i]            

                # Total covariance matrix
                if Method==1 or Method==3 or Method==5 : 
                    CovMatrix_ErrorMu = CovMatCorrData_4ErrorMu + Cov_appmag + Cov_pecvel
                elif Method==2 or Method==4 or Method==6 or Method==7: 
                    CovMatrix_ErrorMu = CovMatCorrData_4ErrorMu + Cov_appmag                   

                # Inverse of the total covariance matrix
                InvCovMatrix_ErrorMu = np.linalg.inv(CovMatrix_ErrorMu)            

                #--- Uncertainty of the best estimated distance modulus---
                # Loop over all the photometry for a given SNe
                # Ignore data outside the interpolation range (lowerPhase, upperPhase):
                for i in range(countGoodPhases): 
                    # Ignore data outside the interpolation range:
                    Numerator_ErrorMu = Numerator_ErrorMu + (mag_array[i] - MagTemp_array[i])*np.sum(InvCovMatrix_ErrorMu[i,:])
                    Denominator_ErrorMu = Denominator_ErrorMu + np.sum(InvCovMatrix_ErrorMu[i,:])            

                if Method==7: # uncertainty in the estimated apparent magnitude at T_Bmax
                    error_appMag_TBmax = np.sqrt(1/Denominator_ErrorMu)
                else: error_appMag_TBmax = 0


                # Uncertainty in the estimated distance modulus
                if Method==1 or Method==3 or Method==4 or Method==5 or Method==6: 
                    sigma2_mu = (1/Denominator_ErrorMu)
                elif Method==2: 
                    sigma2_mu = (1/Denominator_ErrorMu) +  (sigma_pecInt**2)

                elif Method==7: 
                    # old1. sigma2_mu = (1/Denominator_ErrorMu) + error_Average_NIRAbsMag_TBmax**2
                    # old2. sigma2_mu = (1/Denominator_ErrorMu) + IntrinsicDisp**2
                    sigma2_mu = 1/Denominator_ErrorMu

                stdDev_mu = np.sqrt(sigma2_mu)               

                #------------------------------------------------------
                # The chi2 function: I define it just to compute the goodness-of-fit to data once 
                # I have computed the best estimated photometric distance modulus.
                # Note that the following definition need the value of "mu_photo_Analytic" computed above.
                chi2Int=0; mu_int=0;

                if Method==6 or Method==7: mu_int = appMag_TBmax
                else: mu_int = mu_photo_Analytic

                for i in range(countGoodPhases):
                    product1=0;
                    for k in range(countGoodPhases): 
                        product1 = product1 + InvCovMatrix_MeanMu[i,k]*(mag_array[k] - MagTemp_array[k] - mu_int)
                    chi2Int = chi2Int + (mag_array[i] - MagTemp_array[i] - mu_int)*product1

                # In case of considering the case where I have just one datum in the LC:
                if countGoodPhases==1: chi2_dof_Int = chi2Int/(countGoodPhases) 
                else: chi2_dof_Int = chi2Int/(countGoodPhases-1)

                #------------------------------------------------------

                # Theoretical distance mu from LCDM
                mu_LCDM = DistanceMu(magData[0,0], OmMFix, wFix, HoFix)
                
                # Compute the residual and absolute magnitude value
                mu_resid = mu_photo_Analytic - mu_LCDM

                AbsMag_s = appMag_TBmax - mu_LCDM 
                error_AbsMag_s = np.sqrt(error_appMag_TBmax**2 + sigma_pecInt**2)

                #------------------------------------------------------

                # Defining the subsample label

                if list_SNe[j][0][-9:-6] == 'CfA':
                    indexSample = 1
                elif list_SNe[j][0][-9:-6] == 'CSP':
                    indexSample = 2
                elif list_SNe[j][0][-12:-6] == 'Others':
                    indexSample = 3

                #######################################################

                # WRITTING THE RESULTS TO A TEXT FILE

                # Comment the masked SNe in ListSNe_AppMag_.txt:
                if mask==0: commentText = ''
                else: commentText = '##  '

                # Write the distance modulus for -all- the SNe, no matter any cutoffs.
                # This data will be used to write down the 
                # "DistanceMu_All_BeforeCutoffs_.txt" text file

                textfile.write('%s%-30s   %.10f  %.10f  %.10f  %.10f  %13.10f  %12.7f        %1.f          %12.9f  %.10f  %.10f  %.10f  %.10f  %.10f  %13.5f  %13.5f      0.00   %.10f  %.8f  %.10f  %.10f  %13.10f  %.10f  %10.7f  %.7f  %.10f  %.10f  %.10f  %.10f  %.10f    %13.5f  %13.5f  %.10f  %.10f # \n'%(
                    commentText, list_SNe[j][0][:-4],
                    zcmbInt, err_zcmb, mu_photo_Analytic, stdDev_mu, mu_resid, chi2_dof_Int, 
                    indexSample, appMag_TBmax, error_appMag_TBmax, 
                    mu_LCDM, sigma_pecInt,
                    AbsMag_s, error_AbsMag_s, TBmax, err_TBmax,
                    zhelio, err_zhelio, dm15Int, err_dm15, 
                    EBVhost, err_EBVhost, EBV_MW, err_EBVMW, Alamb, err_Alamb, R_F,
                    muSnoopy, err_muSnoopy,
                    TBmax, err_TBmax, appMag_TBmax, error_appMag_TBmax ))                        

                countSNe = countSNe + 1 # Counter

                # ------- QUALITY CUTOFFS:  -------
                # This data will be used to write down the 
                # "DistanceMu_Good_AfterCutoffs_Main_.txt" text file

                if (zcmbInt < zcmbUpperLim and dm15Int > dm15LowerLim and dm15Int < dm15UpperLim and 
                    EBVhost > EBVhostMin and EBVhost < EBVhostMax and EBV_MW < EBVMWLim and 
                    chi2_dof_Int < chi2_dof_Max and mu_resid < residualMax):
                    textfileMain.write('%s%-30s   %.10f  %.10f  %.10f  %.10f  %13.10f  %12.7f        %1.f          %12.9f  %.10f  %.10f  %.10f  %.10f  %.10f  %13.5f  %13.5f      0.00   %.10f  %.8f  %.10f  %.10f  %13.10f  %.10f  %10.7f  %.7f  %.10f  %.10f  %.10f  %.10f  %.10f    %13.5f  %13.5f  %.10f  %.10f # \n'%(
                        commentText, list_SNe[j][0][:-4],
                        zcmbInt, err_zcmb, mu_photo_Analytic, stdDev_mu, mu_resid, chi2_dof_Int, 
                        indexSample, appMag_TBmax, error_appMag_TBmax, 
                        mu_LCDM, sigma_pecInt,
                        AbsMag_s, error_AbsMag_s, TBmax, err_TBmax,
                        zhelio, err_zhelio, dm15Int, err_dm15, 
                        EBVhost, err_EBVhost, EBV_MW, err_EBVMW, Alamb, err_Alamb, R_F,
                        muSnoopy, err_muSnoopy,
                        TBmax, err_TBmax, appMag_TBmax, error_appMag_TBmax ))  

                    textfil6.write('%s \n'%list_SNe[j][0][:-6])

                    mu_resid_array += [mu_resid]
                    chi2_list += [chi2_dof_Int]
                    SNeName_list += [list_SNe[j][0][:-4]]

                    countGoodSNe = countGoodSNe + 1
                    if mask==0: countNoCommented = countNoCommented + 1;
                    else: countCommented = countCommented + 1; 

                    if indexSample == 1: countCfA =  countCfA + 1
                    elif indexSample == 2: countCSP =  countCSP + 1
                    elif indexSample == 3: countOthers =  countOthers + 1

                    #---- Cutoff redshift: z_CMB > 0.01 ---
                    # This data will be used to write down 
                    # the "DistanceMu_Good_AfterCutoffs_z001_.txt"

                    if zcmbInt > 0.01: 
                        textfil3.write('%s%-30s   %.10f  %.10f  %.10f  %.10f  %13.10f  %12.7f        %1.f          %12.9f  %.10f  %.10f  %.10f  %.10f  %.10f  %13.5f  %13.5f      0.00   %.10f  %.8f  %.10f  %.10f  %13.10f  %.10f  %10.7f  %.7f  %.10f  %.10f  %.10f  %.10f  %.10f    %13.5f  %13.5f  %.10f  %.10f # \n'%(
                            commentText, list_SNe[j][0][:-4],
                            zcmbInt, err_zcmb, mu_photo_Analytic, stdDev_mu, mu_resid, chi2_dof_Int, 
                            indexSample, appMag_TBmax, error_appMag_TBmax, 
                            mu_LCDM, sigma_pecInt,
                            AbsMag_s, error_AbsMag_s, TBmax, err_TBmax,
                            zhelio, err_zhelio, dm15Int, err_dm15, 
                            EBVhost, err_EBVhost, EBV_MW, err_EBVMW, Alamb, err_Alamb, R_F,
                            muSnoopy, err_muSnoopy,
                            TBmax, err_TBmax, appMag_TBmax, error_appMag_TBmax ))                                        

                        mu_resid_z001_array += [mu_resid]
                        chi2_z001_list += [chi2_dof_Int]
                        SNeName_z001_list += [list_SNe[j][0][:-4]]
                        countGoodSNe_z001 = countGoodSNe_z001 + 1

                        if indexSample == 1: countCfA_z001 =  countCfA_z001 + 1
                        elif indexSample == 2: countCSP_z001 =  countCSP_z001 + 1
                        elif indexSample == 3: countOthers_z001 =  countOthers_z001 + 1

                else:
                    # This data will be used to write down the "DistanceMu_SNe_OutCutoffs_.txt"
                    textfil5.write('%s%-30s   %.10f  %.10f  %.10f  %.10f  %13.10f  %12.7f        %1.f          %12.9f  %.10f  %.10f  %.10f  %.10f  %.10f  %13.5f  %13.5f      0.00   %.10f  %.8f  %.10f  %.10f  %13.10f  %.10f  %10.7f  %.7f  %.10f  %.10f  %.10f  %.10f  %.10f    %13.5f  %13.5f  %.10f  %.10f # \n'%(
                        commentText, list_SNe[j][0][:-4],
                        zcmbInt, err_zcmb, mu_photo_Analytic, stdDev_mu, mu_resid, chi2_dof_Int, 
                        indexSample, appMag_TBmax, error_appMag_TBmax, 
                        mu_LCDM, sigma_pecInt,
                        AbsMag_s, error_AbsMag_s, TBmax, err_TBmax,
                        zhelio, err_zhelio, dm15Int, err_dm15, 
                        EBVhost, err_EBVhost, EBV_MW, err_EBVMW, Alamb, err_Alamb, R_F,
                        muSnoopy, err_muSnoopy,
                        TBmax, err_TBmax, appMag_TBmax, error_appMag_TBmax ))                               

                    countBadSNe = countBadSNe + 1


                #######################################################
                #
                # Plotting the LC data and the fit with the template for a given SN

                # Offset to displace the template to match apparent magnitude data
                if Method==7: OffsetToMatchAppMagData = appMag_TBmax
                else: OffsetToMatchAppMagData = mu_photo_Analytic


                # Some definitions to plot the template in the range phase
                x = Template[indexLowerPhase:indexUpperPhase,0] 
                y_pred = Template[indexLowerPhase:indexUpperPhase,1] + OffsetToMatchAppMagData

                if TempType == 'MWA':
                    sigma = np.sqrt(Template[indexLowerPhase:indexUpperPhase,3])
                else: 
                    sigma = Template[indexLowerPhase:indexUpperPhase,2]

                # Creating the plot 
                plt.figure()

                # Standard deviation
                plt.fill(np.concatenate([x, x[::-1]]),
                        np.concatenate([y_pred - 1*sigma,
                                      (y_pred + 1*sigma)[::-1]]),
                        # np.concatenate([y_pred - 1 * sigma,
                        #                (y_pred + 1 * sigma)[::-1]]),
                        alpha=0.3, fc='g', ec='None')

                # PLOT GP TEMPLATE mean
                plt.plot(x, y_pred, 'k-', lw=2, alpha=0.4)

                #-------------------

                # Plot the photometric LC data (app mag) for the given SNe 

                plt.errorbar(magData[InfoSN_NumLinesSkip:,0], magData[InfoSN_NumLinesSkip:,3], 
                             magData[InfoSN_NumLinesSkip:,4], 
                             fmt=fmt_int, color='red', ms=MyPointSize, 
                                elinewidth=MyLineWidth, capsize=MyCapeSize
                             # old. ls='', fmt='r.', alpha=1,  markersize=8
                            )

                #-------------------

                plt.grid(alpha=0.6)

                #-- Defining the y-limits for plotting
                y_lowerlim = MinMagTempPlot + OffsetToMatchAppMagData
                y_upperlim = Template[indexMaxTempl,1]+ OffsetToMatchAppMagData

                plt.ylim(y_lowerlim+0.7, y_upperlim-0.5)
                plt.xlim(x_RangePlots)

                if Method==6 or Method==7:                     
                    plt.text(-7, y_lowerlim-0.7, 'Best fit: $m_{T_{Bmax}}$ = %r $\pm$ %r'%(
                        round(appMag_TBmax,3),round(np.sqrt(error_appMag_TBmax),3))) 

                plt.text(-7, y_lowerlim-0.5, '$\mu$ = %r $\pm$ %r'%(round(mu_photo_Analytic,3), 
                                                                    round(np.sqrt(sigma2_mu),3)))  
                plt.text(-7, y_lowerlim-0.3, '$z_{CMB}$ = %r'%round(magData[0,0],5))
                
                if deltamu_print:
                    plt.text(-7, y_lowerlim-0.1, '$\Delta \mu$ = %r'%round(mu_resid,3))
                if Chi2dofPrint:
                    plt.text(-7, y_lowerlim+0.1, '$\chi^2_{dof}$ = %r'%round(chi2_dof_Int,2))  

                plt.title(list_SNe[j][0][:-4])

                # plt.xlabel(r'phase = (MJD-$T_{Bmax}$)/(1+$z_{\rm hel}$)')
                plt.xlabel(r'(MJD-$T_{Bmax}$)/(1+$z_{\rm hel}$)', fontsize=fontsizePlot)
                plt.ylabel('apparent magnitude', fontsize=fontsizePlot)
                
                plt.tick_params(axis='x', labelsize=fontsizePlot-1)
                plt.tick_params(axis='y', labelsize=fontsizePlot-1)   
                
                plt.tight_layout()
                
                if (zcmbInt < zcmbUpperLim and dm15Int > dm15LowerLim and dm15Int < dm15UpperLim and 
                    EBVhost > EBVhostMin and EBVhost < EBVhostMax and EBV_MW < EBVMWLim and
                    chi2_dof_Int < chi2_dof_Max and mu_resid < residualMax):
                    if zcmbInt > 0.01:
                        plt.savefig(DirSaveOutput+'Plot_Fits/AfterCutoffs_z001/'+'%s_FitMu.png'%(
                            list_SNe[j][0][:-4]), dpi=ResolutionPlot_HD, format='png')
                    else:
                        plt.savefig(DirSaveOutput+'Plot_Fits/AfterCutoffs_z0/'+'%s_FitMu.png'%(
                            list_SNe[j][0][:-4]), dpi=ResolutionPlot_HD, format='png')
                else:
                    plt.savefig(DirSaveOutput+'Plot_Fits/OutCutoffs/'+'%s_FitMu.png'%(
                        list_SNe[j][0][:-4]), dpi=ResolutionPlot_HD, format='png')

                plt.close()

                # <--- End plotting data ----

                print 'Done:  %s'%list_SNe[j][0]
                # print ''

    # <--- End loop

    #######################################################

    #     RMS computation

    rms = np.sqrt(np.mean(np.square(mu_resid_array)))
    rms_z001 = np.sqrt(np.mean(np.square(mu_resid_z001_array)))
    # rms_z001_chi2_Res = np.sqrt(np.mean(np.square(mu_resid_z001_chi2_Res_array)))

    # Write to the RMS text file
    textfil4.write('%6.2f       %5.2f      %.4f  %3.f     %.4f      %3.f  \n'%(lowerPhase, upperPhase, rms, countGoodSNe, 
                     rms_z001, countGoodSNe_z001))

    #=======================================================

    # text7 = '# Total number of SNe in this list: %r (CfA=%r, CSP=%r, Others=%r) \n'%

    text7 = '# Total number of SNe in this list: %r (CfA=%r, CSP=%r, Others=%r) \n'%(
        countGoodSNe, countCfA, countCSP, countOthers)
    text8 = '# Total number of SNe with z >0.01: %r (CfA=%r, CSP=%r, Others=%r) \n'%(
        countGoodSNe_z001, countCfA_z001, countCSP_z001, countOthers_z001)
    text_08_01 = "# %s SNe no commented automatically (##). \n"%countNoCommented 
    text_08_02 = "# %s SNe were commented automatically (##). \n"%countCommented 
    text_22 = "# %s SNe were commented automatically (##). \n"%countCommented 

    text9  = '# RMS: %r (all the SNe in this file) \n'%round(rms,4)
    text10 = '# RMS: %r (SNe z > 0.01) \n'%round(rms_z001,4)

    text11 = '# Largest chi2: %r, SNe: %s \n'%(
        round(max(chi2_list),4), SNeName_list[chi2_list.index(max(chi2_list))]   )
    text12 = '# Largest chi2 (z>0.01): %r, SNe: %s \n'%(
        round(max(chi2_z001_list),4),
        SNeName_z001_list[chi2_z001_list.index(max(chi2_z001_list))] )
    # text12_1 = '# Largest chi2: %r, SNe: %s \n'%(round(max(chi2_z001_chi2_Res_list),4), 
    #   SNeName_z001_chi2_Res_list[chi2_z001_chi2_Res_list.index(max(chi2_z001_chi2_Res_list))] )

    text13 = '# Largest upper residual: %r, SNe: %s \n'%(
        round(max(mu_resid_array),4), SNeName_list[mu_resid_array.index(max(mu_resid_array))]   )
    text14 = '# Largest lower residual: %r, SNe: %s \n'%(
        round(min(mu_resid_array),4),SNeName_list[mu_resid_array.index(min(mu_resid_array))]  )

    text15 = '# Largest upper residual (z>0.01): %r, SNe: %s  \n'%(
        round(max(mu_resid_z001_array),4),
        SNeName_z001_list[mu_resid_z001_array.index(max(mu_resid_z001_array))]   )
    # text15_1 = '# Largest upper residual: %r, SNe: %s  \n'%(
    #    round(max(mu_resid_z001_chi2_Res_array),4), 
    #   SNeName_z001_chi2_Res_list[mu_resid_z001_chi2_Res_array.index(max(mu_resid_z001_chi2_Res_array))]   )

    text16 = '# Largest lower residual (z>0.01): %r, SNe: %s  \n'%(
        round(min(mu_resid_z001_array),4), 
        SNeName_z001_list[mu_resid_z001_array.index(min(mu_resid_z001_array))]   )
    # text16_1 = '# Largest lower residual: %r, SNe: %s  \n'%(
    #  round(min(mu_resid_z001_chi2_Res_array),4), 
    #  SNeName_z001_chi2_Res_list[mu_resid_z001_chi2_Res_array.index(min(mu_resid_z001_chi2_Res_array))]   )


    textfile.write('# Total number of SNe in this file: %r \n'%countSNe)

    #--------------------
    textfileMain.write(text_line);
    textfileMain.write(text7); textfileMain.write(text8); 
    textfileMain.write(text_08_01); textfileMain.write(text_08_02); 
    textfileMain.write(text9); textfileMain.write(text10)
    textfileMain.write(text11); textfileMain.write(text12); textfileMain.write(text13); 
    textfileMain.write(text14); textfileMain.write(text15); textfileMain.write(text16);
    #--------------------

    textfil3.write(text8); textfil3.write(text10); textfil3.write(text12); textfil3.write(text15)
    textfil3.write(text16)

    # textfil3_1.write(text8_1); textfil3_1.write(text10_1); textfil3_1.write(text12_1); textfil3_1.write(text15)
    # textfil3_1.write(text16)

    textfil5.write('# Total number of SNe in this file: %r \n'%countBadSNe)

    textfile.close()
    textfileMain.close()
    textfil3.close()
    textfil4.close()
    textfil5.close()
    textfil6.close()
    textfil7.close()

    print ' '
    print '--- All done smoothly, DistanceModuli_*.txt created ---'
    print 'Number of SNe before cutoffs = ', countSNe
    print 'Number of SNe after cutoffs = ', countGoodSNe
    print 'Number of SNe after cutoffs AND 0.01 < z < 0.06 = ', countGoodSNe_z001
    print text7, text_08_01, text_22


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[40]:


##############################################################################80


# # Reading the datatable I've just created above
# 
# #### This will be used to determine the cut based on $\chi^2$ and to determine the mean absolute value

# In[41]:


import os # To use command line like instructions
import glob # To read the files in my directory

#-- Array of minimum redshift to consider --
# zCMB_Min = np.array([0, 0.01])
# zCMB_Min_Label = ['z0', 'z001']

zCMB_Min = np.array([0])
zCMB_Min_Label = ['z0']

#-------------------------------------------------------

# Change the working directory where the data files are located
os.chdir(DirSaveOutput)

# Reset the main array container to avoid using the one from 
# the previous run.
DistMu_array = np.array([0])
if ScriptVersion == 'notebook':
    print "# 'DistMu_array' resetted."

# Reading the data files in that folder 
DistanceMuFiles = glob.glob('DistanceMu_Good_After*.txt')

# Check if 'DistanceModuli_Notes_.txt' is already there, otherwise read
# the 'DistanceModuli_.txt' file.
if 'DistanceMu_Good_AfterCutoffs_Main_Notes_.txt' in DistanceMuFiles:
    DistMu_array = np.genfromtxt(DirSaveOutput+
                        'DistanceMu_Good_AfterCutoffs_Main_Notes_.txt', 
                                 dtype=['S30',
               float,float,float,float,float,float,float,float,float,float,
               float,float,float,float,float,float,float,float,float,float,
               float,float,float,float,float,float,float,float,float,float,
               float,float,float])
    if ScriptVersion == 'notebook':
        print 'Reading the file:  < DistanceMu_Good_AfterCutoffs_Main_Notes_.txt >'
    
elif 'DistanceMu_Good_AfterCutoffs_Main_.txt' in DistanceMuFiles:
    DistMu_array = np.genfromtxt(DirSaveOutput+
                        'DistanceMu_Good_AfterCutoffs_Main_.txt',
                                 dtype=['S30',
               float,float,float,float,float,float,float,float,float,float,
               float,float,float,float,float,float,float,float,float,float,
               float,float,float,float,float,float,float,float,float,float,
               float,float,float])
    if ScriptVersion == 'notebook':
        print '# Reading the file: < DistanceMu_Good_AfterCutoffs_Main_.txt >'
    
elif 'DistanceMu_Good_AfterCutoffs_Main_TMP_Notes_.txt' in DistanceMuFiles:
    DistMu_array = np.genfromtxt(DirSaveOutput+
                        'DistanceMu_Good_AfterCutoffs_Main_TMP_Notes_.txt',
                                 dtype=['S30',
               float,float,float,float,float,float,float,float,float,float,
               float,float,float,float,float,float,float,float,float,float,
               float,float,float,float,float,float,float,float,float,float,
               float,float,float])
    if ScriptVersion == 'notebook':
        print '# Reading the file: < DistanceMu_Good_AfterCutoffs_Main_TMP_Notes_.txt >'
    
elif 'DistanceMu_Good_AfterCutoffs_Main_TMP_.txt' in DistanceMuFiles:
    DistMu_array = np.genfromtxt(DirSaveOutput+
                        'DistanceMu_Good_AfterCutoffs_Main_TMP_.txt',
                                 dtype=['S30',
               float,float,float,float,float,float,float,float,float,float,
               float,float,float,float,float,float,float,float,float,float,
               float,float,float,float,float,float,float,float,float,float,
               float,float,float])
    if ScriptVersion == 'notebook':
        print '# Reading the file: < DistanceMu_Good_AfterCutoffs_Main_TMP_.txt >'

else: print '# < DistanceMu_Good_AfterCutoffs_Main_.txt > file not found!!!'

print '# Number of SNe in the file = ', len(DistMu_array)


# In[42]:


# Reading the file: < DistanceMu_Good_AfterCutoffs_Main_.txt >
# Number of SNe in the file =  63


# In[43]:


#-----------------------------------------------------------------------------80


# ------

# # $\chi^2_{d.o.f.}$ histogram

# In[44]:


# Create 2 lists: for z<0 and z<0.01

if NotebookToPlotOnly == False:
    
    chi2dof_z0_list = []
    chi2dof_z001_list = []

    for i in range(len(DistMu_array)):

        name = DistMu_array[i][0]
        zcmbInt = DistMu_array[i][1]
        chi2dof_z0_list += [DistMu_array[i][6]]

        if zcmbInt > 0.01:
            chi2dof_z001_list += [DistMu_array[i][6]]


    if ScriptVersion == 'notebook':
        print 'len(chi2dof_z0_list) = ', len(chi2dof_z0_list)
        print 'len(chi2dof_z001_list) = ', len(chi2dof_z001_list)
        # len(chi2dof_z001_list) =  78


# In[45]:


# Plotting the histogram of chi2_dof

if NotebookToPlotOnly == False:
    
    # countChi2dofPlot = countChi2dofPlot + 1
    fontsizeText = 15

    BinSize = 45

    for zz in zCMB_Min:
        fig = plt.figure()

        if zz==0.0: plt.hist(chi2dof_z0_list, BinSize, histtype=u'barstacked')
        elif zz==0.01: plt.hist(chi2dof_z001_list, BinSize, histtype=u'barstacked')        

        plt.xlabel(r'$\chi^2_{\rm{d.o.f}}$', fontsize=fontsizeText+2)
        plt.ylabel('Frequency', fontsize=fontsizeText)
        plt.title(r'Histogram $\chi^2_{\rm{d.o.f}}$ (%r<z<%r)'%(zz, zcmbUpperLim), fontsize=fontsizeText)
        plt.text(15, 10, 'Bin size = %r'%BinSize)

        plt.grid(True)

        plt.savefig('Plot_histo_chi2dof_BinSize%r_z%r_.png'%(BinSize, zz))
        plt.close()


# In[46]:


if NotebookToPlotOnly == False: plt.close(); plt.close();


# In[47]:


# Plotting the values of chi2dof as lines

if NotebookToPlotOnly == False:
    
    # Array of zeros just for the y-axis values.
    zero_z0_np = np.zeros(len(chi2dof_z0_list))
    zero_z001_np = np.zeros(len(chi2dof_z001_list))

    for zz in zCMB_Min:
        fig = plt.figure(figsize=(12,4))

        if  zz==0.0:
            plt.plot(chi2dof_z0_list, zero_z0_np, ls='None', marker='|', ms=30, 
                     markeredgewidth=2, color='r',  alpha=1)
        elif zz==0.01:
            plt.plot(chi2dof_z001_list, zero_z001_np, ls='None', marker='|', ms=30, 
                     markeredgewidth=2, color='r',  alpha=1)

        plt.xlabel(r'$\chi^2_{\rm{d.o.f}}$', fontsize=fontsizeText+2)
        plt.ylabel('No meaning', fontsize=fontsizeText)
        plt.title(r'$\chi^2_{\rm{d.o.f}}$ (%r<z<%r)'%(zz,zcmbUpperLim), fontsize=fontsizeText)

        fig.tight_layout()

        plt.ylim(-1, 1)

        plt.grid(True)

        plt.savefig('Plot_histo_chi2dof_Points_z%r_.png'%(zz))
        plt.close()


# In[48]:


if NotebookToPlotOnly == False: plt.close(); plt.close();


# In[ ]:





# In[49]:


#-----------------------------------------------------------------------------80


# # Determining average Absolute magnitude  of the sample

# ### Plotting the histogram with the Gaussian estimation

# In[50]:


# Plotting the histogram with the Gaussian estimation

if NotebookToPlotOnly == False:
    
    from scipy.stats import norm
    import matplotlib.mlab as mlab
    
    # best fit of data.
    # 'norm.fit' simply compute the mean and standard devation of the sample.
    (Average_NIRAbsMag_TBmax, error_Average_NIRAbsMag_TBmax) = norm.fit(DistMu_array['f12'])

    #------ Plot -----------

    # the histogram of the data
    n, bins, patches = plt.hist(DistMu_array['f12'], 20, density=True, facecolor='green', alpha=0.5)
    # n, bins, patches = plt.hist(DistMu_array['f12'], 20, normed=False, facecolor='green', alpha=0.5)


    # add a 'best fit' Gaussian line
    y = mlab.normpdf(bins, Average_NIRAbsMag_TBmax, error_Average_NIRAbsMag_TBmax )
    l = plt.plot(bins, y, 'r--', linewidth=2)

    plt.xlabel(r'$\hat{m}_{\rm Bmax} - \mu_{\Lambda{\rm CDM}}$')
    plt.ylabel('Probability')
    plt.title(r'Histogram. (mean=%.2f, std dev = %.2f)' %(Average_NIRAbsMag_TBmax, error_Average_NIRAbsMag_TBmax))

    plt.grid(True)
    plt.savefig(DirSaveOutput+'Plot_histo_AbsMag_.png')
    plt.close()


# In[51]:


if NotebookToPlotOnly == False: plt.close(); plt.close();


# In[52]:


if NotebookToPlotOnly == False:
    
    # The results that I will use now.
    # NOTE: This value is independent of what I have assume about "Average_NIRAbsMag_TBmax" and
    # "error_Average_NIRAbsMag_TBmax" at the top of the notebook in the User section.
    # NOTE: The values of "Average_NIRAbsMag_TBmax" and "error_Average_NIRAbsMag_TBmax" do NOT 
    # depend on the value of "vpecFix".

    # ------------------------------
    now = datetime.datetime.now() # Read the time and date right now
    text_timenow = now.strftime("%Y-%m-%d; %H:%M hrs.")
    text_Date   = '# On date: %s \n'%text_timenow
    # ------------------------------

    print '#', '-'*30
    print '#  ', BandName, 'band | vpecFix =', vpecFix, 'km/s | 0 < z <', zcmbUpperLim
    # print '# (mean abs mag, std deviation)\n' 
    # print ""

    print "# Average_NIRAbsMag_TBmax = %.6f  # %s"%(Average_NIRAbsMag_TBmax, text_timenow)
    print "# error_Average_NIRAbsMag_TBmax = %.6f;"%(error_Average_NIRAbsMag_TBmax)


# In[ ]:





# #### Weighted averages

# In[53]:


# Inverse-variance weights array

# if NotebookToPlotOnly == False:
if 3<2:

    WeightsInvVar = 1/np.square(DistMu_array['f13'])

    # Inverse-variance weighted average
    WeightedAbsMag = np.average(DistMu_array['f12'], weights=WeightsInvVar )

    # Useful definitions to determine the weighted population standard deviation
    V1 = np.sum(WeightsInvVar)
    V2 = np.sum(np.square(WeightsInvVar))

    product2 = 0
    # Computing the unbiased weighted population standard deviation
    for i in range(len(DistMu_array['f12'])):
        # print i
        absMagInt = DistMu_array['f12'][i]
        product2 = product2 + WeightsInvVar[i]*(absMagInt-WeightedAbsMag)**2

    error_WeightedAbsMag = np.sqrt(product2/(V1 - (V2/V1)))


# In[ ]:





# In[54]:


# Compute some other values

# print '#'+'-'*30
# print '#   ', BandName, 'band   Pec velocity =', vpecFix, 'km/s'
# print '#', np.mean(DistMu_array['f12']), '= mean abs mag'
# print '#', np.median(DistMu_array['f12']), '= median' 
# print '#', WeightedAbsMag, '= Weighted Abs Mag'
# print '#', np.sqrt(1/np.sum(WeightsInvVar) ), '= uncertainty in the weighted average'
# print '#', np.std(DistMu_array['f12']), '= population standard deviation'
# print '#', error_WeightedAbsMag, '= unbiased weighted population standard deviation'


# ### RMS

# In[55]:



if ScriptVersion == 'notebook':
    # Notebook version

    import numpy as np
    import os # To use command line like instructions

    DirMyPyFuncs = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/SNIRfit/'

    # Change to this directory
    os.chdir(DirMyPyFuncs)

    # Auto -re-load  every given time (in case there are changes)
    # %load_ext autoreload
    # %autoreload 1

    # Load my python functions
    get_ipython().magic(u'aimport RMS')

    # Go back to the DirSaveOutput directory
    os.chdir(DirSaveOutput)

elif ScriptVersion == 'terminal':
    # Change to this directory
    print os.getcwd()
    os.chdir(DirSaveOutput)
    import RMS


# In[56]:


# Test that my python function works well

# xx = np.array([ 2.,  2.,  2.])
# sigma = np.array([ 2.,  2.,  2.])
# ww = np.array([ 3.,  3.,  3.])
# 
# print RMS.rms(xx), RMS.err_rms(xx,sigma)
# # 2.0 1.15470053838
# 
# print RMS.wrms(xx,ww), RMS.err_wrms(xx,ww,sigma)
# 2.0 0.816496580928


# In[57]:


# ?z_RMS.rms


# In[58]:


# ?z_RMS.sigma_rms


# In[ ]:







# In[59]:


##############################################################################80


# -------

# # HUBBLE DIAGRAM PLOTS
# 
# - Note that in all the rest of the the following cells, it is NOT computed or 
# used at all the values of (Average_NIRAbsMag_TBmax, error_Average_NIRAbsMag_TBmax).
# 
# 
# -  (SEE UPDATED INFORMATION BELOW) If I need to remove an outlier in any of 
# the "any YJHK", YJHK, YJH or JHK Gaussian-Process or template Hubble diagrams, then: 
# 
# 
# 1. first I need to "comment" that SN in 'DistanceMu_Good_AfterCutoffs_Main_.txt' in all the bands, 
# 2. then, rerun the "13_TotalDistanceMu_vXX.ipynb" iPython notebook to 
# recompute the covariance matrices and their inverse matrix from the remaining SNe.
# 3. plot the Hubble diagram using "11_DistanceMu_HubbleDiagram_vXX.ipynb" notebook (this section).
# 
# UPDATED INFORMATION ABOUT THE CASE OF REMOVING OUTLIERS IN THE "any YJHK", 
# YJHK, YJH or JHK Gaussian-Process or template Hubble diagrams:
# 
# The process described above is the correct to remove outliers, however, 
# for simplicity when making different experiments removing a lot of SNe 
# (specially when comparing the NIR vs SALT2 Hubble diagrams from SNe that 
# pass SALT2 cutoffs) I'm recalibrating the Hubble diagram residual by 
# simply determining the value that allows to minimize the weighted *average* 
# in the Hubble residual plot (exactly as done for the SALT2 Hubble diagram). 
# This is done by using the simple single instruction:
# 
# elif BandMax == 'NIRmax' and PlotTotalMu == True:
# 
#      delta_Mo = SimplexResult_1[0]
# 
# If I want to use the correct procedure to remove outliers in the "GP NIR 
# any band HD", I simply have to follow the instructions indicated and comment 
# the instruction above.
# The drawback of this method is that the uncertainties are less consistent 
# with scatter in the Hubble diagram, i.e., chi^2_dof is not close to one.
# 
# 
# 'DistanceModuli_Notes_.txt': The advantage of having a text file listing all 
# the SNe is that I can easily "comment" the outliers SNe in the Hubble 
# diagram and residual.

# ### USER 

# In[60]:


# Setting 

# Band to use as the reference peak maximum? Options: (NIRmax, Bmax).
# For the template method to determine the distance modulus I use BandMax = 'Bmax'
# For the GP method to determine the distance modulus I use BandMax = 'NIRmax'
# "Bmax_GP": Determine the distance moduli at B-band max but using the GP method.

# WATCH OUT!: Everytime I change between (NIRmax, Bmax, Snoopy, SALT2) 
# Hubble diagrams I must reset this notebook because in each case it 
# is used different values for (Ho, Omega_matter, w).

# ( Bmax , NIRmax , Bmax_GP , Snoopy , SALT2 ).
if ScriptVersion == 'notebook':
    BandMax = 'Bmax'  
elif ScriptVersion == 'terminal':
    BandMax = sys.argv[23]

#--------------------
#   LABELS IN THE PLOT'S TITLE for the OPTICAL

if ScriptVersion == 'notebook':
    if BandMax == 'SALT2': BandNameText = 'Optical' 
    if BandMax == 'Snoopy': BandNameText = 'Optical' 
    
elif ScriptVersion == 'terminal':
    BandNameText = sys.argv[24]

# if BandMax == 'Snoopy': BandName = 'Optical+NIR'  # TEMPORAL
# if BandMax == 'Snoopy': BandName = 'NIR'  # TEMPORAL
# if BandMax == 'Snoopy': BandName = 'Y'  # TEMPORAL

#    RAISINs
# if BandMax == 'Snoopy': BandName = 'Optical' # RAISINs
# if BandMax == 'Snoopy': BandName = 'Opt+NIR' # RAISINs

#------- "Total" distance modulus ---------

# Plot the "total" distance modulus derived from the three distance modulus computed 
# from each band?
# To put TRUE, first I need to have already computed the distance moduli for each SNe 
# in the YJHK bands!

if ScriptVersion == 'notebook':
    PlotTotalMu = False
elif ScriptVersion == 'terminal':
    PlotTotalMu = sys.argv[25] == 'True'

    
# "PlotTotalMu = True"  is only valid when BandMax = NIRmax, Bmax, Bmax_GP.
# "PlotTotalMu = False" if I'm plotting the SNoopy or SALT2 distance moduli.
if BandMax == 'Snoopy' or BandMax == 'SALT2':
    PlotTotalMu = False

# - "AllBands": If there is not distance mu estimation for a given SN, then 
# use the 3, 2 or 1 band that it was observed
# - "YJH": Consider the SNe that have distance mu estimation in these 
# 3 bands ONLY, and so on.

# ( AllBands , JH , YJH , JHK ,  YJHK )
if PlotTotalMu == True: 
    
    if ScriptVersion == 'notebook':
        BandsCombination = 'AllBands'  
        
    elif ScriptVersion == 'terminal':
        BandsCombination = sys.argv[26]
        
else: BandsCombination = ''
    
#####################################################

#    Plot RAISINs?
# Set the limits and labels in the plot if I'm plotting RAISINs.


if ScriptVersion == 'notebook':
    # False = plot the low-z sample.
    plot_raisins = False
    
elif ScriptVersion == 'terminal':
    # False = plot the low-z sample.
    plot_raisins =  sys.argv[27] == 'True'

#--------------------

if plot_raisins == True:
    zcmbUpperLim = 0.65  # (0.6, 0.65, anything else)
    xlimPlots = 0.2, zcmbUpperLim+0.003 # 
    ylimPlots = 40.0, 44.5
    
#####################################################

#    Weight definition to compute the wRMS

## Options: (wRMS_efit_ePec_eInt, wRMS_efit)
## 'wRMS_efit_ePec_eInt': w = 1/(sigma2_{fit, s}, sigma2_{mu_pec}, sigma2_{int})
## 'wRMS_efit': w = 1/(sigma2_{fit, s})

wRMS_type = 'wRMS_efit_ePec_eInt'
    
#####################################################

#    PLOT OPTIONS
# Resolution in the Hubble diagram plots
ResolutionPlot_HD = 150  # dpi

# Settings error bars:
MyPointSize = 6  
MyLineWidth = 1  
MyCapeSize = 2


# In[61]:


##############################################################################80


# ## Automatic

# In[62]:


# Apparent magnitude, determined from GP interpolation, template, snoopy.?
    
if BandMax == 'Bmax':
    AppMag_method = 'Template to compute the app mags'  
    #old HoFix = 73.24 # Valued used by default in the low-z paper
    OmMFix = 0.28 # Omega_Matter
    OmLFix = 0.72 # Omega_Lambda
    wFix = -1 # Dark energy EoS
    
elif BandMax == 'Bmax_GP':
    AppMag_method = 'Gaussian Process to compute the app mags'  
    #old HoFix = 73.24 # Valued used by default in the low-z paper
    OmMFix = 0.28 # Omega_Matter
    OmLFix = 0.72 # Omega_Lambda
    wFix = -1 # Dark energy EoS
    
elif BandMax == 'NIRmax':
    # AppMag_method = 'Gaussian Process to compute app mags'
    AppMag_method = 'Gaussian Process to compute the NIR app mags' # TEMPORAL
    #old HoFix = 73.24 # Valued used by default in the low-z paper
    OmMFix = 0.28 # Omega_Matter
    OmLFix = 0.72 # Omega_Lambda
    wFix = -1 # Dark energy EoS
    
elif BandMax == 'Snoopy':
    AppMag_method = 'SNooPy fit'
    # Redefing some variables: 
    # HoFix = 73.24 # km/s. Ho=72 is Snoopy's default value
    #old HoFix = 72.0 # km/s. Ho=72 is Snoopy's default value
    OmMFix = 0.28 # Omega_Matter
    OmLFix = 0.72 # Omega_Lambda
    wFix = -1 # Dark energy EoS
    
# ORIGINAL
elif BandMax == 'SALT2':
    AppMag_method = 'SALT2 fit'
    #old HoFix = 73.24 # Valued used by default in the low-z paper
    OmMFix = 0.28 # Omega_Matter
    OmLFix = 0.72 # Omega_Lambda
    wFix = -1 # Dark energy EoS  
    
# TEMPORAL FOR RAISIN DD proposal:
# elif BandMax == 'SALT2':
#     AppMag_method = '' # TEMPORAL
#     HoFix = 70.0 
#     OmMFix = 0.3 # Omega_Matter
#     OmLFix = 0.7 # Omega_Lambda
#     wFix = -1 # Dark energy EoS

#----------------------------------------------

print "# Hubble diagram type: %s"%BandMax
print "# PlotTotalMu? (when combined YJHK bands only) = %s \n"%PlotTotalMu
print "#     Cosmology to be used:"
print "# Ho = %s km/s/Mpc."%HoFix
print "# Omega_Matter = %s"%OmMFix
print "# Omega_Lambda = %s"%OmLFix


# In[63]:


# Reading the datafile with the TOTAL distance modulus inferred from
# the distance moduli coming from each band.

if PlotTotalMu and BandMax == 'NIRmax': Method_folder = 'GaussianProcess'
if PlotTotalMu and BandMax == 'Bmax': Method_folder = 'Template'
if PlotTotalMu and BandMax == 'Bmax_GP': Method_folder = 'GP_Bmax'

if PlotTotalMu: 
    # Change the working directory where the datafiles are located
    DirSaveOutput = DirMain_1+'AllBands/Plots/HubbleDiagram_vpec%s/%s/%s/'%(
        vpecFix,Method_folder, BandsCombination)
    if not os.path.exists(DirSaveOutput): os.makedirs(DirSaveOutput)
    
    DirDataTotalMu = DirMain_1+'AllBands/Plots/HubbleDiagram_vpec%s/'%vpecFix+Method_folder+'/'
    os.chdir(DirDataTotalMu)
   

    try:
        DistMu_array = np.genfromtxt('Table_TotalMu_%s_Notes_.txt'%(
                    BandsCombination),
                    dtype=['S30', float,float,float,
                            float,float,float,float])
        print "# File read: 'Table_TotalMu_%s_Notes_.txt'"%(
            BandsCombination)
    except: 
        DistMu_array = np.genfromtxt('Table_TotalMu_%s_.txt'%(
                    BandsCombination),
                    dtype=['S30', float,float,float,
                           float,float,float,float])
        print "# File read: 'Table_TotalMu_%s_.txt'"%BandsCombination
        
    BandNameText = 'All'
    
    
    print "# Location:\n", DirDataTotalMu
    print '# %s SNe in PlotTotalMu file.'%len(DistMu_array)
# Number of SNe in PlotTotalMu file: 30.


# #### Determine the value needed to have a weighted average of zero in the Hubble residual 
# 
# Used mainly (but not limited) for SALT2 Hubble diagram
# 

# In[64]:


# Define the average in the absolute value of the residual
# These functions are going to be minimized.

# AbsResidual values
mu_1_np = DistMu_array['f3'] 
z_1_np  = DistMu_array['f1']
res_mu_np_int = mu_1_np - DistanceMuVector(z_1_np, OmMFix, wFix, HoFix)

err_mu_1_np = DistMu_array['f4']

# Inverse-variance weighted average function to be minimized 
def WeightedAbsResidual_fun(deltaMo, UpLimit, LowLimit, AbsResidual_Roof):
    residuals_int = res_mu_np_int +  deltaMo
    WeightedAbsResidual_int = np.average(residuals_int, weights=err_mu_1_np )
    
    if deltaMo < UpLimit and deltaMo > LowLimit:  
        # For some unknown reason, simplex is search the maximum instead of
        # the minimimum, so I had to define the average absolute mag as 1/WeightedAbsResidual
        # WeightedAbsResidual_final = 1/WeightedAbsResidual_int
        WeightedAbsResidual_final = WeightedAbsResidual_int
        
    else: WeightedAbsResidual_final = AbsResidual_Roof
        
    return abs(WeightedAbsResidual_final)

#-----------------------------------------------------

# Simple average function to be minimized 
def AverageAbsResidual_fun(deltaMo, UpLimit, LowLimit, AbsResidual_Roof):
    residuals_int = res_mu_np_int + deltaMo
    # AverageAbsResidual_int = np.mean(residuals_int)
    AverageAbsResidual_int = np.sum(residuals_int)/len(residuals_int)
     
    if deltaMo < UpLimit and deltaMo > LowLimit: 
        # For some unknown reason, simplex is search the maximum instead of
        # the minimimum, so I had to define the average absolute mag as 1/AverageAbsResidual_int
        # AverageAbsResidual_final = 1/AverageAbsResidual_int
        AverageAbsResidual_final = AverageAbsResidual_int
        
    else: AverageAbsResidual_final = AbsResidual_Roof

    return abs(AverageAbsResidual_int)


# In[65]:


# Determine the value of the constant deltaMo in order to obtain a
# weighted average value of zero in the Hubble residual

import scipy.optimize

# Search limits:
# UpLimit = 0.2; LowLimit = -0.2
UpLimit = 3; LowLimit = -3

# Assume this value for deltaMo in case of the search is outside the range
# indicated above.
Residual_Roof = 10 

WeightedAbsResidual_Out = scipy.optimize.minimize_scalar(WeightedAbsResidual_fun, 
                                        args=(UpLimit, LowLimit, Residual_Roof))

AverageAbsResidual_Out = scipy.optimize.minimize_scalar(AverageAbsResidual_fun, 
                                        args=(UpLimit, LowLimit, Residual_Roof))

# Redefining the values:
SimplexResult_1 = [WeightedAbsResidual_Out['x'] ]
SimplexResult_2 = [AverageAbsResidual_Out['x'] ]

if ScriptVersion == 'notebook':
    print WeightedAbsResidual_Out
    print '# %.6f = value of delta_Mo that minimize the Hubble residual.'%SimplexResult_1[0]
    
# print AverageAbsResidual_Out
# print '#', SimplexResult_2, ' = value of delta_Mo that minimize the Hubble residual.'


# In[ ]:





# In[66]:


"""
     fun: 1.482459899373599e-12
    nfev: 33
     nit: 32
 success: True
       x: 5.0906660693874514e-07
# [5.0906660693874514e-07]  = value of delta_Mo that minimize the Hubble residual.
"""
0


# In[67]:


# Add a value to the theoretical distance modulus so that the weighted 
# average of the Hubble residual is zero. This is only for the SALT2 Hubble diagram.

if ScriptVersion == 'notebook':

    if BandMax == 'SALT2': # ORIGINAL
    # if BandMax == 'SALT2' or BandMax == 'Snoopy': # FOR RAISINS OPTICAL FITS
        delta_Mo = SimplexResult_1[0] 

    else: 
        delta_Mo_int1 = np.array([0.])
        delta_Mo = delta_Mo_int1[0]
    
    print '# delta_Mo = %s'%delta_Mo
    # print '# type of variable: %s'%type(delta_Mo)
    
elif ScriptVersion == 'terminal':
    minimize_residuals = sys.argv[28] == 'True'
    
    if minimize_residuals:
        delta_Mo = SimplexResult_1[0] 
    else: 
        delta_Mo_int1 = np.array([0.])
        delta_Mo = delta_Mo_int1[0]


# In[68]:


#-- Array of minimum redshift to consider --

# zCMB_Min = np.array([0, 0.01])
# zCMB_Min_Label = ['z0', 'z001']

zCMB_Min = np.array([0])
zCMB_Min_Label = ['z0']

#--------------------------
#   Color array
# ColorSamples = [True, False]
ColorSamples = [True]


# In[69]:


#-----------------------------------------------------------------------------80


# ### Plot Hubble diagram

# In[70]:


fontsizePlot = 10.5

# To plot the theoretical (spectroscopic) distance modulus
nbins1= 51
z1 = np.linspace(xlimPlots[0], xlimPlots[1], nbins1)  
mu1 = DistanceMuVector(z1, OmMFix, wFix, HoFix) + delta_Mo # ORIGINAL
# FOR RAISINS SNOOPY FIT:
# mu1 = DistanceMuVector(z1, OmMFix, wFix, HoFix) - delta_Mo 

# Count number of SNe for each subsample
count_CfA = 0; count_CSP = 0; count_all = 0

#---------------------------------
# Plot

for k in range(len(ColorSamples)): # Loop over the colors
    
    for j in range(len(zCMB_Min)):  # Loop over the 'zCMB_Min' array

        #------- Plotting the data -------
        fig = plt.figure()
        # fig = plt.figure(figsize=(8,6), dpi=80)
        # OLD. fig = plt.figure(figsize=(12, 8))
        
        # loop over 'DistanceMu_Good_AfterCutoffs_Main_.txt'
        for i in range(len(DistMu_array)): 

            zzInt = DistMu_array[i][1] # z_CMB for a given SNe
            # chi2dofInt = DistMu_array[i][6] # chi2 by d.o.f.
            # absolute value of the residual value:
            # residInt = abs(DistMu_array[i][5]) 
            sampleFlag = DistMu_array[i][7]  # Subsample           
                
            # Create the variable "snName" containing the first 8 
            # (or 7) letters of the SNe file name
            # I use "snName" to compare with the SNe names in 
            # 'SNeWithCepheidDistances.txt', so that
            # I will not compute a peculiar velocity uncertainty for those SNe.
            try:
                if   DistMu_array[i][0][7] == '_': 
                    # To read correctly, e.g., "sn2011B_"
                    snName = DistMu_array[i][0][:7] 
                elif DistMu_array[i][0][7] != '_':
                    # To read correctly, e.g., "snf20080514-002":
                    if is_number(DistMu_array[i][0][7]): snName = DistMu_array[i][0][:15] 
                    # To read correctly, e.g., "sn1998bu"
                    else: snName = DistMu_array[i][0][:8]     
            except: 
                # To read correctly, e.g., "sn2011B"
                snName = DistMu_array[i][0][:6]  
                
                
            # Make special symbol for SNe with Cepheid.
            if snName in ListSNeCepheid['f0']: fmt_int = '*'
            else: fmt_int = '.'
            
            if ColorSamples[k] and KindOfData4HD == 'AllSamples':
                # CSP data
                if sampleFlag == 2 and zzInt > zCMB_Min[j]:   
                    plt.errorbar(zzInt, DistMu_array[i][3], yerr=DistMu_array[i][4], 
                        fmt=fmt_int, color='blue', ms=MyPointSize, 
                        elinewidth=MyLineWidth, capsize=MyCapeSize)
                # CfA data
                elif sampleFlag == 1 and zzInt > zCMB_Min[j]: 
                    plt.errorbar(zzInt, DistMu_array[i][3], yerr=DistMu_array[i][4], 
                        fmt=fmt_int, color='red', ms=MyPointSize, 
                        elinewidth=MyLineWidth, capsize=MyCapeSize)
                # Others data
                elif sampleFlag == 3 and zzInt > zCMB_Min[j]: 
                    plt.errorbar(zzInt, DistMu_array[i][3], yerr=DistMu_array[i][4], 
                        fmt=fmt_int, color='green', ms=MyPointSize, 
                        elinewidth=MyLineWidth, capsize=MyCapeSize)
                # Anything else
                elif zzInt > zCMB_Min[j]: 
                    plt.errorbar(zzInt, DistMu_array[i][3], yerr=DistMu_array[i][4], 
                        fmt=fmt_int, color='black', ms=MyPointSize, 
                        elinewidth=MyLineWidth, capsize=MyCapeSize)
                    
            else:
                if zzInt > zCMB_Min[j]:  
                    if KindOfData4HD == 'CfA' and sampleFlag == 1:
                        plt.errorbar(zzInt, DistMu_array[i][3], yerr=DistMu_array[i][4], 
                            fmt=fmt_int, color='red', ms=MyPointSize, 
                            elinewidth=MyLineWidth, capsize=MyCapeSize)
                        # count_CfA = count_CfA + 1 # Counter
                    elif KindOfData4HD == 'CSP'  and sampleFlag == 2:
                        plt.errorbar(zzInt, DistMu_array[i][3], yerr=DistMu_array[i][4], 
                            fmt=fmt_int, color='blue', ms=MyPointSize, 
                            elinewidth=MyLineWidth, capsize=MyCapeSize)
                        # count_CSP = count_CSP + 1 # Counter
                    elif KindOfData4HD == 'Others'  and sampleFlag == 3:
                        plt.errorbar(zzInt, DistMu_array[i][3], yerr=DistMu_array[i][4], 
                            fmt=fmt_int, color='green', ms=MyPointSize, 
                            elinewidth=MyLineWidth, capsize=MyCapeSize)
                    elif KindOfData4HD == 'AllSamples':
                        plt.errorbar(zzInt, DistMu_array[i][3], yerr=DistMu_array[i][4], 
                            fmt=fmt_int, color='black', ms=MyPointSize, 
                            elinewidth=MyLineWidth, capsize=MyCapeSize)
                        # count_all = count_all + 1 # Counter    
             
        # Plotting the theory line
        plt.plot(z1, mu1, color='black')

        plt.xlim(xlimPlots)
        plt.ylim(ylimPlots)

        plt.xlabel('Redshift', fontsize=fontsizePlot)
        plt.ylabel(r'Distance modulus $\mu$', fontsize=fontsizePlot)
        
        #--- Plot Title ---
        if PlotTotalMu == False:
            plt.title('Hubble diagram (%s band)'%(BandNameText), 
                      fontsize=fontsizePlot+2)
        elif PlotTotalMu == True:
            if BandsCombination == 'AllBands':
                plt.title('Hubble diagram (all YJHK bands)', 
                          fontsize=fontsizePlot+2)
            else: plt.title('Hubble diagram (%s bands)'%(BandsCombination), 
                            fontsize=fontsizePlot+2)
        
        # plt.legend(loc='lower right')
        # plt.tick_params(labelsize=fontsizePlot)
        plt.tick_params(axis='x', labelsize=fontsizePlot)
        plt.tick_params(axis='y', labelsize=fontsizePlot+2)        

        plt.grid(True, ls='--', alpha=0.5)
        plt.tight_layout()
                
        # if   PlotTotalMu == True:  
        #     NIRbandsText_1 = BandsCombination+'_'
        #     NIRbandsText_2 = BandsCombination+'/'
        # elif PlotTotalMu == False: 
        #     NIRbandsText_1 = ''; NIRbandsText_2 = ''
            
        if   PlotTotalMu == True:  NIRbandsText = BandsCombination+'_'
        elif PlotTotalMu == False: NIRbandsText = ''

        if ColorSamples[k] and KindOfData4HD == 'AllSamples':
            plt.savefig(DirSaveOutput+'Plot_HD_%s_%s_%s_%s_%s_%s_Color_%s.png'%(
                KindOfData4HD, KindOfTemp,
                              KindOfTempSubgroup, 
                                   chi2_dof_Max_Label, residualMax_Label, 
                                      zCMB_Min_Label[j], NIRbandsText)
                             , dpi=ResolutionPlot_HD, format='png'
                       )
        else:
            plt.savefig(DirSaveOutput+'Plot_HD_%s_%s_%s_%s_%s_%s_%s.png'%(
                KindOfData4HD, KindOfTemp,
                                  KindOfTempSubgroup,
                                    chi2_dof_Max_Label, residualMax_Label, 
                                       zCMB_Min_Label[j], NIRbandsText) 
                                        , dpi=ResolutionPlot_HD, format='png'
                        )
        plt.close()

# print 'Number SNe from CfA: ', count_CfA
# print 'Number SNe from CSP: ', count_CSP
# print 'Number SNe from all: ', count_all


# In[71]:


plt.close();plt.close();plt.close();


# In[72]:


#-----------------------------------------------------------------------------80


# ## Intrinsic dispersion

# In[73]:


# Creation of arrays for mu_residuals.

mu_resid_z0_np = np.zeros(len(DistMu_array))
sigma2_appMagTBmax_z0_np = np.zeros(len(DistMu_array))
sigma2_pec_z0_np = np.zeros(len(DistMu_array))

mu_resid_z001 = []
sigma2_appMagTBmax_z001 = []
sigma2_pec_z001 = []

# loop over 'DistanceMu_Good_AfterCutoffs_Main_.txt'
for i in range(len(DistMu_array)): 
    
    zcmbInt     = DistMu_array[i][1]  # zcmb
    err_zcmbInt = DistMu_array[i][2]  # error_zcmb
    mu_resid_z0_np[i]   = DistMu_array[i][5] +delta_Mo  # Delta_mu
    
    # Create the variable "snName" containing the first 8 (or 7) 
    # letters of the SNe file name I use "snName" to compare with 
    # the SNe names in 'SNeWithCepheidDistances.txt', so that
    # I will not compute a peculiar velocity uncertainty for those SNe.
    try:
        if   DistMu_array[i][0][7] == '_': 
            # To read correctly, e.g., "sn2011B_"
            snName = DistMu_array[i][0][:7]  
        elif DistMu_array[i][0][7] != '_':
            # To read correctly, e.g., "snf20080514-002":
            if is_number(DistMu_array[i][0][7]): 
                snName = DistMu_array[i][0][:15] 
            # To read correctly, e.g., "sn1998bu"
            else: snName = DistMu_array[i][0][:8]    
    except: 
        # To read correctly, e.g., "sn2011B"
        snName = DistMu_array[i][0][:6]  
    
    
    #--- Determine the uncertainty in distance modulus 
    # due to the peculiar velocity uncertainty.---
    
    #-----------   For z > 0: --------------
    
    if PlotTotalMu == True: # For 'PlotTotalMu' *compute* "sigma_mu_pecVel"
        if snName in ListSNeCepheid['f0']: 
            sigma2_pec_z0_np[i] = sigma2_pec(zcmbInt, err_zcmbInt, 0)
        else: sigma2_pec_z0_np[i] = sigma2_pec(zcmbInt, err_zcmbInt, vpecFix)
    else: # For 'NIRmax', 'Bmax' 'Snoopy', 'SALT2' reads "sigma_mu_pecVel"
        sigma2_pec_z0_np[i] = (DistMu_array[i][11])**2   # (sigma_mu_pecVel)^2
    
    # Read the uncertainty on the distance mu due only to fitting the LCs:
    # old:
    # if PlotTotalMu==False and (BandMax == 'NIRmax' or 
    #                            BandMax == 'Bmax' or BandMax == 'Bmax_GP' or
    #                            BandMax == 'SALT2' or BandMax == 'Snoopy'):
    
    # uncertainty on mu due only to fitting the LC
    sigma2_appMagTBmax_z0_np[i] = (DistMu_array[i][4])**2.0  
    
    
    #-----------   For z > 0.01: --------------
    
    if zcmbInt > 0.01:
        mu_resid_z001 += [DistMu_array[i][5]]  # Delta_mu  
        
        if PlotTotalMu == True:
            if snName in ListSNeCepheid['f0']: 
                sigma2_pec_z001 += [sigma2_pec(zcmbInt, err_zcmbInt, 0)]
            else: sigma2_pec_z001 += [sigma2_pec(zcmbInt, err_zcmbInt, vpecFix)]
        else:  # For NIRmax', 'Bmax' 'Snoopy', 'SALT2' 
            sigma2_pec_z001 += [(DistMu_array[i][11])**2]  # (sigma_mu_pecVel)^2        
        
        # Read the uncertainty on the distance mu due only to fitting the LCs:
        # old:
        # if PlotTotalMu==False and (BandMax == 'NIRmax' or 
        #                            BandMax == 'Bmax'  or BandMax == 'Bmax_GP' or
        #                            BandMax == 'SALT2'  or BandMax == 'Snoopy'):
            
        # uncertainty on mu due only to fitting the LC
        sigma2_appMagTBmax_z001 += [(DistMu_array[i][4])**2.0]  
        
# Convert the list to np.arrays:
mu_resid_z001_np = np.array(mu_resid_z001)
sigma2_appMagTBmax_z001_np = np.array(sigma2_appMagTBmax_z001)
sigma2_pec_z001_np = np.array(sigma2_pec_z001)


# In[74]:


# sigma2_appMagTBmax_z0_np


# In[75]:


# sigma2_appMagTBmax_z001_np


# In[76]:


# print 'Number of SNe with z>0 :', len(mu_resid_z0_np)
# print 'Number of SNe with z>0.01 :', len(mu_resid_z001_np)


# In[77]:


# Checking the sized of the arrays: the left and right numbers have to be the same
# in a given row printed below.

if ScriptVersion == 'notebook':
    print ('#', len(mu_resid_z0_np),   len(sigma2_appMagTBmax_z0_np),   
           len(sigma2_pec_z0_np)   )
    print ('#', len(mu_resid_z001_np), len(sigma2_appMagTBmax_z001_np), 
           len(sigma2_pec_z001_np) )

# 119 119 119
# 95 95 95

# 63 63 63
# 49 49 49


# In[78]:


# Finding the best estimated value for sigma2Pred by minimizing the -2ln(Likelihood) function.

# Defining the function to compute the intrinsic dispersion (sigmaPred) instead
# of the square of the intrinsic dispersion (sigma2Pred): for the case of determining
# the instrisic dispersion.
# I obtain an error due to a very small value of sigmaPred, so that during minimizing 
# the likelihood function to determine sigmaPred, it is sampled some negative values.

# Likelihood to determine the intrinsic dispersion
# -2ln(Likelihood) function, Eq. (B.6) of Blondin et al 2011

# For the case of using a normalized template. This is the full Eq. (B.6) of Blondin et al 2011
def neg2lnLikeFull(sigmaPred, mu_resid_np, sigma2_pec_np, sigma2_appmagTBmax_np):
    sum1 = 0
    for i in range(len(mu_resid_np)):
        sum1 = (sum1 + np.log(sigma2_appmagTBmax_np[i] + sigmaPred**2 + sigma2_pec_np[i]) +  
                (mu_resid_np[i]**2)/(sigma2_appmagTBmax_np[i] + sigmaPred**2 + sigma2_pec_np[i]) )  
    return sum1

#----

# The truncated version.
def neg2lnLike(sigmaPred, mu_resid_np, sigma2_pec_np):
    sum1 = 0
    for i in range(len(mu_resid_np)):
        sum1 = (sum1 + np.log(sigmaPred**2 + sigma2_pec_np[i]) + 
                (mu_resid_np[i]**2)/(sigmaPred**2 + sigma2_pec_np[i]) ) 
    return sum1

#-----------------------------------------------------------------
from scipy.optimize import fmin as simplex

InitialGuess = 0.15


# old if Method==7 and PlotTotalMu==False and (BandMax == 'NIRmax' or BandMax == 'Bmax' or BandMax == 'Bmax_GP' #):
if Method==7 and (BandMax == 'NIRmax' or BandMax == 'Bmax' or BandMax == 'Bmax_GP' #):
                                         or BandMax == 'SALT2' or BandMax == 'Snoopy'): # ORIGINAL
                                         # or BandMax == 'SALT2'): # FOR RAISINS SNOOPY FIT
    SimplexResult_z0_2 = simplex(neg2lnLikeFull, InitialGuess, 
                                   args=(mu_resid_z0_np, sigma2_pec_z0_np, sigma2_appMagTBmax_z0_np))

    SimplexResult_z001_2 = simplex(neg2lnLikeFull, InitialGuess, 
                                   args=(mu_resid_z001_np, sigma2_pec_z001_np, sigma2_appMagTBmax_z001_np))

"""
# old 
elif PlotTotalMu==True: # ORIGINAL
# elif PlotTotalMu==True or BandMax == 'Snoopy' or BandMax == 'SALT2': # FOR RAISINS SNOOPY FIT
    SimplexResult_z0_2 = simplex(neg2lnLike, InitialGuess, 
                              args=(mu_resid_z0_np,   sigma2_pec_z0_np))
    SimplexResult_z001_2 = simplex(neg2lnLike, InitialGuess, 
                              args=(mu_resid_z001_np, sigma2_pec_z001_np)) """

SimplexResult_z0 = SimplexResult_z0_2**2
SimplexResult_z001 = SimplexResult_z001_2**2

# print 'Best estimated value of sigma_Pred (z>0) =', SimplexResult_z0_2
# print 'Best estimated value of sigma_Pred (z>0.01) =', SimplexResult_z001_2

if ScriptVersion == 'notebook': 
    print 'Best estimated value of sigma^2_Pred (z>0) =', SimplexResult_z0
    print 'Best estimated value of sigma_Pred (z>0) =', SimplexResult_z0_2

    print 
    print 'Best estimated value of sigma^2_Pred (z>0.01) =', SimplexResult_z001
    print 'Best estimated value of sigma_Pred (z>0.01) =', SimplexResult_z001_2


# In[ ]:





# In[79]:


# Computing the uncertainty on sigma_pred

# Define the Fisher information matrix, Eq. (B.7) of Blondin et al 2011

# For the case of using a normalized template. This is the full Eq. (B.7) of Blondin et al 2011
def FisherFuncFull(sigmaPred, mu_resid_np, sigma2_pec_np, sigma2_appmagTBmax_np):
    sum2 = 0
    for i in range(len(mu_resid_np)):
        sum2 = (sum2 + (mu_resid_np[i]**2)/(sigma2_appmagTBmax_np[i] + 
                sigmaPred**2 + sigma2_pec_np[i])**3 -
               1/(2*(sigma2_appmagTBmax_np[i] + sigmaPred**2 + sigma2_pec_np[i])**2)  )
    return sum2

#--------

def FisherFunc(sigmaPred, mu_resid_np, sigma2_pec_np):
    sum2 = 0
    for i in range(len(mu_resid_np)):
        sum2 = (sum2 + (mu_resid_np[i]**2)/(sigmaPred**2 + sigma2_pec_np[i])**3 -
               1/(2*(sigmaPred**2 + sigma2_pec_np[i])**2)  )
    return sum2

#-----------------------------------------------------------------

# The variance error of sigma^2_pred:
# old if Method==7 and PlotTotalMu==False and (BandMax == 'NIRmax' or BandMax == 'Bmax' or 
if Method==7 and (BandMax == 'NIRmax' or BandMax == 'Bmax' or 
                                         BandMax == 'Bmax_GP' or
                                        BandMax == 'SALT2'  or BandMax == 'Snoopy'):
    Var_sigma2_pred_z0 =   1/FisherFuncFull(SimplexResult_z0_2[0],   
                                    mu_resid_z0_np,   sigma2_pec_z0_np,
                                    sigma2_appMagTBmax_z0_np)
    Var_sigma2_pred_z001 = 1/FisherFuncFull(SimplexResult_z001_2[0], 
                                    mu_resid_z001_np, sigma2_pec_z001_np,
                                    sigma2_appMagTBmax_z001_np)
    
    
error_sigma_pred_z0   = np.sqrt(Var_sigma2_pred_z0  /(4*SimplexResult_z0[0]))
error_sigma_pred_z001 = np.sqrt(Var_sigma2_pred_z001/(4*SimplexResult_z001[0]))

if ScriptVersion == 'notebook':
    print 'Variance of sigma^2_pred (z>0) =', round(Var_sigma2_pred_z0,6)
    print 'Variance of sigma^2_pred (z>0.01) =', round(Var_sigma2_pred_z001,6)
    print '\nError_sigma_pred (z>0) =', round(error_sigma_pred_z0,3)
    print 'Error_sigma_pred (z>0.01) =', round(error_sigma_pred_z001,3)


# #### Report the intrinsic dispersion

# In[80]:


print '#'+'-'*50
# print "# %s band, vpecFix = %s km/s"%(Band, vpecFix)
# print "# Intrinsic dispersion for the case (0 < z < %s):"%(zcmb_Max)
print '# Intrinsic dispersion = %.6f +/- %.6f'%(np.sqrt(SimplexResult_z0[0]), 
                                                error_sigma_pred_z0)

# ------------------------------
now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d; %H:%M hrs.")
text_Date   = '# On date: %s \n'%text_timenow
# ------------------------------

print '# '
print '# %s | Date: %s '%(AppMag_method, text_timenow)
print "# IntrinsicDisp = %.5f # %s | vpecFix=%s km/s | 0<z<%s."%(
    np.sqrt(SimplexResult_z0[0]), BandNameText, vpecFix, zcmbUpperLim)


# In[81]:


#--------------------------------------------------
# Intrinsic dispersion = 0.07517578125 +/- 0.0629356479937
 
# Gaussian Process to compute app mags | Date: 2018-05-02; 17:04 hrs. 
# IntrinsicDisp = 0.07518 # K | vpecFix=150 km/s | 0<z<0.04.

#--------------------------------------------------
# Intrinsic dispersion = 0.12017578125 +/- 0.016725821387
 
# Gaussian Process to compute app mags | Date: 2018-02-22; 16:49 hrs. 
# IntrinsicDisp = 0.12018 # J | vpecFix=150 km/s | 0<z<0.04.


# In[82]:


##############################################################################80


# -------

# # $\chi^2_{\rm d.o.f.}$
# 
# #### Checking the consistency between the error bars of the residual distance modulus vs the scatter in the Hubble-diagram residual plot

# In[83]:


ratio_int = 0;  

for i in range(len(DistMu_array)):
    # print i
    
    mu_resid  = DistMu_array[i][5] +delta_Mo
    
    if   BandMax == 'Bmax':   IntrinsicDisp_int = np.sqrt(SimplexResult_z0[0])
    elif BandMax == 'NIRmax': IntrinsicDisp_int = np.sqrt(SimplexResult_z0[0])
    elif BandMax == 'Bmax_GP': IntrinsicDisp_int = np.sqrt(SimplexResult_z0[0])
    elif BandMax == 'SALT2':  IntrinsicDisp_int = np.sqrt(SimplexResult_z0[0])
    elif BandMax == 'Snoopy': IntrinsicDisp_int = np.sqrt(SimplexResult_z0[0])
     
    if  PlotTotalMu == False and (BandMax == 'NIRmax' or BandMax == 'Bmax' or
    # old. if  (BandMax == 'NIRmax' or BandMax == 'Bmax' or 
                                  BandMax == 'Bmax_GP' or
                                  BandMax == 'SALT2'  or BandMax == 'Snoopy'): 
        sigma_pecVel = DistMu_array[i][11]
        error_appMag = DistMu_array[i][4]
        ratio_int = ratio_int + ((mu_resid**2) / (error_appMag**2 + 
                                                  sigma_pecVel**2 + 
                                                  IntrinsicDisp_int**2) )

    elif PlotTotalMu == True:  
        sigma_pecVel = np.sqrt(sigma2_pec(zcmbInt, err_zcmbInt, vpecFix))
        error_appMag = DistMu_array[i][4]
        ratio_int = ratio_int + ((mu_resid**2) / (error_appMag**2 + 
                                                  sigma_pecVel**2 + 
                                                  IntrinsicDisp_int**2) )

    
chi2_dof_HD    = ratio_int / len(DistMu_array)


print '#'+'-'*40
print "# %s band | vpecFix = %s km/s | 0 < z < %s"%(BandNameText, vpecFix, zcmbUpperLim)
print '# %s.'%AppMag_method
print '# chi^2 = %.6f | Number of data = %s'%(ratio_int, len(DistMu_array)) 
print '# chi^2_dof = %.6f'%chi2_dof_HD


# In[84]:


#-----------------------------------------------------------------------------80


# ## Write the summary of values to a text file

# In[85]:


if (BandMax != 'NIRmax' or PlotTotalMu == True):
    textfile_1 = open(DirSaveOutput+"Summary_HDScatter_RMS_.txt", 'w')

    now = datetime.datetime.now() # Read the time and date right now
    text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
    text_Author = '# Data table created by: Arturo Avelino \n'
    text_Date   = '# On date: %s \n'%text_timenow
    text_script = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
    text_line = '#'+'-'*50 + '\n'

    textfile_1.write("# Summary of the scatter in the Hubble residual \n")
    textfile_1.write("# %s band \n"%BandNameText)

    textfile_1.write(text_line)
    textfile_1.write(text_Author); textfile_1.write(text_Date); 
    textfile_1.write(text_script); 
    textfile_1.write(text_line)

    if (BandMax == 'Bmax' and PlotTotalMu == False):
        
        try: 
            # This section is for the case that I've used 11.ipynb to compute
            # the distance moduli, so that the values above were previously defined
            # during the creation of the "DistanceMu_Good_AfterCutoffs_Main_.txt"
            # text file.
            textfile_1.write(text01); textfile_1.write(text1); textfile_1.write(text2);
            textfile_1.write(text3); textfile_1.write(text3_0_1); textfile_1.write(text3_3);
            textfile_1.write(text3_4); textfile_1.write(text3_5); 
            textfile_1.write(text_line)
            
        except:
            # This section is for the case that I'm using 11.ipynb just to plot 
            # in the template method.
            text01 = '# Apparent mag data used to construct the Hubble diagram: %s \n'%KindOfData4HD
            text1 = '# Template used: \n'
            text2 = '# %s \n'%DirTemplate[78:]
            text3 = '# Cosmology used to compute residuals: (Omega_M=%r, \
Omega_L=%r, w=%r, Ho=%r) \n'%(OmMFix, OmLFix, wFix, HoFix)

            text3_0_0 = '# 0.01 < z_cmb < %r \n'%zcmbUpperLim
            text3_0_1 = '# Cutoffs: z_cmb<%r | %r<dm15<%r | %r<EBVhost<%r \
| EBV_mw<%r | chi2_dof<%r | Residual<%r \n'%(zcmbUpperLim, dm15LowerLim,  
                dm15UpperLim, EBVhostMin, EBVhostMax, EBVMWLim, chi2_dof_Max, residualMax)
            text3_3 = '# Phase range used of the template: (%r, %r) days \n'%(PhaseMinTemp, PhaseMaxTemp)
            text3_4 = '# Minimal number of data in LC to be considered \
for the Hubble diagram: %r \n'%MinNumOfDataInLC
            text3_5 = '# Assumed uncertainty in the peculiar velocity: %r km/s \n'%vpecFix

            textfile_1.write(text01); textfile_1.write(text1); textfile_1.write(text2);
            textfile_1.write(text3); textfile_1.write(text3_0_1); textfile_1.write(text3_3);
            textfile_1.write(text3_4); textfile_1.write(text3_5); 
            textfile_1.write(text_line)
    
    textfile_1.write("# %s km/s = peculiar-velocity uncertainty. \n"%vpecFix)
        
    textfile_1.write('# %.5f # = delta_Mo. Valued ADDED to the theoretical distance \n\
# modulus so that the weighted average of the Hubble residual is zero. \n'%delta_Mo)
        
    textfile_1.write('# Intrinsic dispersion for the case (0 < z < %s) and used to obtain the \n\
# total photometric distance modulus uncertainty: \n'%zcmbUpperLim)
    
    textfile_1.write('%14.6f  %.6f  0.0  0.0 # Intrinsic dispersion and its uncertainty for the case \
0 < z_cmb < %s \n'%(np.sqrt(SimplexResult_z0[0]), error_sigma_pred_z0, zcmbUpperLim))
    textfile_1.write('%14.6f  %.6f  0.0  0.0 # Intrinsic dispersion and its uncertainty for the case \
0.01 < z_cmb < %s \n'%(np.sqrt(SimplexResult_z001[0]), error_sigma_pred_z001, zcmbUpperLim))

    textfile_1.write('%14.6f  %.6f  0.0  0.0 # (AverageAbsMag_atMax, err_AverageAbsMag_atMax) \
    \n'%(Average_NIRAbsMag_TBmax, error_Average_NIRAbsMag_TBmax))

    textfile_1.write('%14.6f  %-8.0f  0.0  0.0 # (chi^2, Number of SNe) for the case (0 < z < %s) \
    \n'%(ratio_int, len(DistMu_array), zcmbUpperLim))
    textfile_1.write('%14.6f  0         0.0  0.0 # chi^2_dof \n'%chi2_dof_HD)

    textfile_1.write(text_line)
    textfile_1.close()
    
print "# File written."


# In[86]:


##############################################################################80


# # Plot residual

# In[87]:




# Compute the standard errors of RMS from bootstrap:
errRMS_bootstrap = True
loopsize = 2000

# Plot error bars of the Hubble residuals?
# This option sometimes is useful to better visualize any trend
# in the Hubble residual plot.
# NOTE: This only apply to Y-band data, I need to implement the
# instruction in the other NIR bands.
# plotErrorBars_list = [True, False]
plotErrorBars_list = [True]

fontsizePlot = 10.5

# Label the outliers.
# In one plot, print the name of those SNe with one value in "label_array"
# and in another plot, print those SNe with the other value.
# Best option: "np.array([0.25, 3])"
label_array = np.array([0.01, 3])

#---------------------------------
#    Text files

if   PlotTotalMu == False: textfile7 = open(DirSaveOutput+'Verbose_%s.txt'%BandNameText, 'w')
elif PlotTotalMu == True:  textfile7 = open(DirSaveOutput+'Verbose_%s.txt'%BandsCombination, 'w')

# Append data to the already existing file 'Summary_HDScatter_.txt', but
# if the file doesn't exist then create it:
textfile_8 = open(DirSaveOutput+'Summary_HDScatter_RMS_.txt', 'a')

textfile7.write('# RMS SUMMARY \n')
textfile7.write(' \n')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_line = '#'+'-'*50 + '\n'

# Add a value to the theoretical distance modulus so that the weighted
# average of the Hubble residual is zero. This is only for the SALT2 Hubble diagram.
textfile7.write('%s # = delta_Mo. Valued ADDED to the theoretical distance modulus so that \
the weighted average of the Hubble residual is zero: \n'%delta_Mo)

textfile7.write('%s band | vpecFix = %s km/s \n'%(BandNameText, vpecFix))
textfile7.write(' \n')

textfile7.write('intrinsic dispersion sigma_Pred (z>0): %r +/- %r \n'%
                (np.sqrt(SimplexResult_z0[0]), error_sigma_pred_z0))
textfile7.write('intrinsic dispersion sigma_Pred (z>0.01): %r +/- %r \n'%
                (np.sqrt(SimplexResult_z001[0]), error_sigma_pred_z001))

textfile7.write('\n Checking the consistency between the error bars of the residual distance modulus vs the scatter in the Hubble-diagram residual plot (case 0 < z < %s) \n'%zcmbUpperLim)
textfile7.write('chi^2 = %s | Number of data = %s \n'%(ratio_int, len(DistMu_array)))
textfile7.write('chi^2_dof = %s \n'%chi2_dof_HD)
textfile7.write(text_line)

textfile7.write('# Data table created by: Arturo Avelino \n')
textfile7.write('# On date: %s \n'%text_timenow)
textfile7.write('# Script used: %s (version %s | last update: %s)\n'%(
    code_name, version_code, last_update))
textfile7.write(text_line)

#---------------------------------

#    Theoretical values

# Array plot the -theoretical- (spectroscopic) distance modulus
nbins1= 151
# old. z1 = np.linspace(xlimPlots[0], xlimPlots[1], nbins1)
z1 = np.linspace(xlimPlots[0], xlimPlots[1]+0.05, nbins1)
mu0 = np.zeros(len(z1))

# Array to plot the -theoretical- peculiar velocity uncertanty
# 0.001 is the average error_zcmb in Andy's compilation
# old. error_z1 = 0.001 * np.ones(len(z1))
error_z1 = err_zCMB_fix * np.ones(len(z1))
sigma_pec_np = np.sqrt(sigma2_pec(z1, error_z1, vpecFix))


#--------------------------------------------------------60
# Plot

# Loop over plotting or not the error bars.
for i3 in range(len(plotErrorBars_list)):

    plotErrorBars = plotErrorBars_list[i3]

    for ll in label_array: # Loop to label the name of the outlier SNe
        labelOutlierLim = ll # Label the outliers with residual >= this quantity

        for k in range(len(ColorSamples)): # Loop over the colors

            for j in range(len(zCMB_Min)):  # Loop over the 'zCMB_Min' array

                # fig = plt.figure(figsize=(12, 8))
                fig = plt.figure()

                # Plotting theory: Flat Universe
                plt.plot(z1, mu0, color='black', lw=2, alpha=0.5)

                # Plot the peculiar velocity uncertainty

                plt.plot(z1, sigma_pec_np, ls='--', color='black')
                plt.plot(z1, -sigma_pec_np, ls='--', color='black')


                # Reset all the values
                # To be used to compute the RMS:
                mu_resid_all = []; error_mu_all = []; sigma2_pec_all=[];
                mu_resid_CSP = []; error_mu_CSP = []; sigma2_pec_CSP=[]; countCSP = 0;
                mu_resid_CfA = []; error_mu_CfA = []; sigma2_pec_CfA=[]; countCfA = 0;
                mu_resid_Others = []; error_mu_Others = []; sigma2_pec_Others=[];
                countOthers = 0
                mu_resid_AnyOther = []; error_mu_AnyOther = []; sigma2_pec_AnyOther=[];
                countAnyOther = 0

                mu_resid_all_np = 0;  error_mu_all_np = 0; sigma2_pec_all_np=0;
                mu_resid_CSP_np = 0;  error_mu_CSP_np = 0; sigma2_pec_CSP_np=0;
                mu_resid_CfA_np = 0;  error_mu_CfA_np = 0; sigma2_pec_CfA_np=0;
                mu_resid_Others_np=0;error_mu_Others_np=0;sigma2_pec_Others_np=0;
                mu_resid_AnyOther_np=0;error_mu_AnyOther_np=0;sigma2_pec_AnyOther_np=0;

                weight_int_np = np.array([]); sigma2_total_np = np.array([]);

                rms = 0;     ratio = 0;     rms_weight =0;
                rms_CfA = 0; ratio_CfA = 0; rms_CfA_weight =0;
                rms_CSP = 0; ratio_CSP = 0; rms_CSP_weight =0;
                rms_Others = 0; ratio_Others = 0; rms_Others_weight =0;
                rms_AnyOther = 0; ratio_AnyOther = 0; rms_AnyOther_weight =0;

                ratio_1 = 0; var_wrms2 = 0; err_wrms = 0;
                ratio_CSP_1 = 0; var_CSP_wrms2 = 0; err_CSP_wrms = 0; rms_err_CSP=0;
                ratio_CfA_1 = 0; var_CfA_wrms2 = 0; err_CfA_wrms = 0; rms_err_CfA=0;
                ratio_Others_1 = 0; var_Others_wrms2 = 0;
                err_Others_wrms = 0; rms_err_Others=0;
                ratio_AnyOther_1 = 0; var_AnyOther_wrms2 = 0;
                err_AnyOther_wrms = 0; rms_err_AnyOther=0;

                #----  Plotting the data ----

                # loop over 'DistanceMu_Good_AfterCutoffs_Main_.txt'
                for i in range(len(DistMu_array)):

                    zzInt = DistMu_array[i][1] # z_CMB for a given SNe
                    err_zz= DistMu_array[i][2] # z_CMB for a given SNe
                    residInt = abs(DistMu_array[i][5]) # residual (absolute value)
                    error_mu = DistMu_array[i][4] # uncertainty in the distance modulus, mu.
                    mu_resid_int = DistMu_array[i][5]+delta_Mo # Delta_mu

                    # Create the variable "snName" containing the first 8
                    # (or 7) letters of the SNe file name
                    # I use "snName" to compare with the SNe names in
                    # 'SNeWithCepheidDistances.txt', so that
                    # I will not compute a peculiar velocity uncertainty for those SNe.

                    if plot_raisins == 'False':
                        try:
                            if   DistMu_array[i][0][7] == '_':
                                # To read correctly, e.g., "sn2011B_":
                                snName = DistMu_array[i][0][:7]
                            elif DistMu_array[i][0][7] != '_':
                                # To read correctly, e.g., "snf20080514-002":
                                if is_number(DistMu_array[i][0][7]):
                                    snName = DistMu_array[i][0][:15]
                                # To read correctly, e.g., "sn1998bu":
                                else: snName = DistMu_array[i][0][:8]
                        # To read correctly, e.g., "sn2011B"
                        except: snName = DistMu_array[i][0][:6]

                    else:
                        # FOR RAISINS
                        snName = DistMu_array[i][0]

                    #-----------------------------
                    # Define the SN name used to label the outliers

                    snName_2 = '' # reset

                    if plot_raisins == 'False':
                        snName_2 = DistMu_array[i][0][2:8]
                    else:
                        # FOR RAISINS SNOOPY
                        snName_2 = DistMu_array[i][0]

                    #-----------------------------

                    # Define the characters to use to print outliers:
                    ini_char_outlier = 3
                    end_char_outlier = 10

                    #-----------------------------

                    if PlotTotalMu == True:
                        if snName in ListSNeCepheid['f0']:
                            sigma2_pec_int = sigma2_pec(zzInt, err_zz, 0)
                        else: sigma2_pec_int = sigma2_pec(zzInt, err_zz, vpecFix)
                    else:
                        sigma2_pec_int = (DistMu_array[i][11])**2  # (sigma_mu_pecVel)^2

                    sampleFlag = DistMu_array[i][7] # Distinguish between CfA or CSP

                    # Make special symbol for SNe with Cepheid.
                    if snName in ListSNeCepheid['f0']: fmt_int = '*'
                    else: fmt_int = '.'

                    #====================================================

                    if ColorSamples[k] and KindOfData4HD == 'AllSamples':

                        #--- CSP data:
                        if sampleFlag == 2 and zzInt > zCMB_Min[j]:
                            if plotErrorBars:
                                plt.errorbar(zzInt, mu_resid_int, yerr=error_mu,
                                    fmt=fmt_int, color='blue', ms=MyPointSize,
                                    elinewidth=MyLineWidth, capsize=MyCapeSize)
                            else:
                                plt.plot(zzInt, mu_resid_int,
                                    marker=fmt_int, color='blue', ms=MyPointSize)

                            mu_resid_CSP += [mu_resid_int]
                            error_mu_CSP += [error_mu]
                            sigma2_pec_CSP += [sigma2_pec_int]
                            countCSP = countCSP + 1

                            # Label the SNe with residual larger than a given value.
                            if residInt > labelOutlierLim:
                                plt.text(zzInt+0.002, mu_resid_int-0.015,
                                DistMu_array[i][0][ini_char_outlier:end_char_outlier],
                                         fontsize=(fontsizePlot-4), color='blue')

                        #--- CfA data
                        elif sampleFlag == 1 and zzInt > zCMB_Min[j]:
                            plt.errorbar(zzInt, mu_resid_int, yerr=error_mu,
                                fmt=fmt_int, color='red',  ms=MyPointSize,
                                elinewidth=MyLineWidth, capsize=MyCapeSize)
                            mu_resid_CfA += [mu_resid_int]
                            error_mu_CfA += [error_mu]
                            sigma2_pec_CfA += [sigma2_pec_int]
                            countCfA = countCfA + 1

                            # Label the SNe with residual larger than a given value.
                            if residInt > labelOutlierLim:
                                plt.text(zzInt+0.002, mu_resid_int-0.015,
                                DistMu_array[i][0][ini_char_outlier:end_char_outlier],
                                         fontsize=(fontsizePlot-4), color='red')

                        #--- Others data
                        elif sampleFlag == 3 and zzInt > zCMB_Min[j]:
                            plt.errorbar(zzInt, mu_resid_int, yerr=error_mu,
                                fmt=fmt_int, color='green', ms=MyPointSize,
                                elinewidth=MyLineWidth, capsize=MyCapeSize)
                            mu_resid_Others += [mu_resid_int]
                            error_mu_Others += [error_mu]
                            sigma2_pec_Others += [sigma2_pec_int]
                            countOthers = countOthers + 1

                            # Label the SNe with residual larger than a given value.
                            if residInt > labelOutlierLim:
                                plt.text(zzInt+0.002, mu_resid_int-0.015,
                                DistMu_array[i][0][ini_char_outlier:end_char_outlier],
                                         fontsize=(fontsizePlot-4), color='green')

                        #--- Any other kind of data
                        elif zzInt > zCMB_Min[j]:
                            plt.errorbar(zzInt, mu_resid_int, yerr=error_mu,
                                fmt=fmt_int, color='black', ms=MyPointSize,
                                elinewidth=MyLineWidth, capsize=MyCapeSize)
                            mu_resid_AnyOther += [mu_resid_int]
                            error_mu_AnyOther += [error_mu]
                            sigma2_pec_AnyOther += [sigma2_pec_int]
                            countAnyOther = countAnyOther + 1

                            # Label the SNe with residual larger than a given value.
                            if residInt > labelOutlierLim:
                                plt.text(zzInt+0.002, mu_resid_int-0.015,
                                DistMu_array[i][0][ini_char_outlier:end_char_outlier],
                                         fontsize=(fontsizePlot-4), color='black')


                    else:
                        if  zzInt > zCMB_Min[j]:

                            if KindOfData4HD == 'CfA' and sampleFlag == 1:
                                plt.errorbar(zzInt, mu_resid_int, yerr=error_mu,
                                    fmt=fmt_int, color='red', ms=MyPointSize,
                                elinewidth=MyLineWidth, capsize=MyCapeSize)
                                mu_resid_all += [mu_resid_int]
                                error_mu_all += [error_mu]
                                sigma2_pec_all += [sigma2_pec_int]
                                # Label the SNe with residual larger than a given value.
                                if residInt > labelOutlierLim:
                                    plt.text(zzInt+0.002, mu_resid_int-0.015,
                                    DistMu_array[i][0][ini_char_outlier:end_char_outlier],
                                             fontsize=(fontsizePlot-4), color='red')

                            elif KindOfData4HD == 'CSP'  and sampleFlag == 2:
                                plt.errorbar(zzInt, mu_resid_int, yerr=error_mu,
                                    fmt=fmt_int, color='blue', ms=MyPointSize,
                                    elinewidth=MyLineWidth, capsize=MyCapeSize)
                                mu_resid_all += [mu_resid_int]
                                error_mu_all += [error_mu]
                                sigma2_pec_all += [sigma2_pec_int]
                                # Label the SNe with residual larger than a given value.
                                if residInt > labelOutlierLim:
                                    plt.text(zzInt+0.002, mu_resid_int-0.015,
                                    DistMu_array[i][0][ini_char_outlier:end_char_outlier],
                                             fontsize=(fontsizePlot-4), color='blue')

                            elif KindOfData4HD == 'Others'  and sampleFlag == 3:
                                plt.errorbar(zzInt, mu_resid_int, yerr=error_mu,
                                    fmt=fmt_int, color='green', ms=MyPointSize,
                                    elinewidth=MyLineWidth, capsize=MyCapeSize)
                                mu_resid_all += [mu_resid_int]
                                error_mu_all += [error_mu]
                                sigma2_pec_all += [sigma2_pec_int]
                                # Label the SNe with residual larger than a given value.
                                if residInt > labelOutlierLim:
                                    plt.text(zzInt+0.002, mu_resid_int-0.015,
                                    DistMu_array[i][0][ini_char_outlier:end_char_outlier],
                                             fontsize=(fontsizePlot-4), color='green')

                            elif KindOfData4HD == 'AllSamples':
                                plt.errorbar(zzInt, mu_resid_int, yerr=error_mu,
                                    fmt=fmt_int, color='black', ms=MyPointSize,
                                    elinewidth=MyLineWidth, capsize=MyCapeSize)
                                mu_resid_all += [mu_resid_int]
                                error_mu_all += [error_mu]
                                sigma2_pec_all += [sigma2_pec_int]
                                # Label the SNe with residual larger than a given value.
                                if residInt > labelOutlierLim:
                                    plt.text(zzInt+0.002, mu_resid_int-0.015,
                                    DistMu_array[i][0][ini_char_outlier:end_char_outlier],
                                             fontsize=(fontsizePlot-4), color='black')

                if ColorSamples[k] and KindOfData4HD == 'AllSamples':
                    mu_resid_all = mu_resid_CSP + mu_resid_CfA + mu_resid_Others + mu_resid_AnyOther
                    error_mu_all = error_mu_CSP + error_mu_CfA + error_mu_Others + error_mu_AnyOther
                    sigma2_pec_all = sigma2_pec_CSP + sigma2_pec_CfA + sigma2_pec_Others + sigma2_pec_AnyOther

                ############ Computing the (RMS, WRMS) residuals ##########

                mu_resid_all_np = np.array(mu_resid_all)
                error_mu_all_np = np.array(error_mu_all)
                sigma2_pec_all_np = np.array(sigma2_pec_all)

                if len(mu_resid_all) > 0:

                    # reset
                    weight_int_np = np.array([]); sigma2_total_np = np.array([]);

                    if wRMS_type == 'wRMS_efit':
                        weight_int_np = 1.0/(error_mu_all_np**2.0)
                        sigma2_total_np = error_mu_all_np**2.0

                    elif wRMS_type == 'wRMS_efit_ePec_eInt':
                        weight_int_np = 1.0/( (error_mu_all_np**2.0) +
                                                            sigma2_pec_all_np +
                                                            SimplexResult_z0 )
                        sigma2_total_np = ((error_mu_all_np**2.0) +
                                                            sigma2_pec_all_np +
                                                            SimplexResult_z0)

                    rms_weight = RMS.wrms(mu_resid_all_np, weight_int_np);
                    rms = RMS.rms(mu_resid_all_np)


                    if errRMS_bootstrap:
                        err_wrms = RMS.err_wrms_boot(mu_resid_all_np, weight_int_np,
                            loopsize);
                        rms_err = RMS.err_rms_boot(mu_resid_all_np,loopsize)
                    else:
                        try:
                            err_wrms = RMS.err_wrms(mu_resid_all_np, weight_int_np,
                                sigma2_total_np);
                        except: err_wrms = -1.0

                        try:rms_err=RMS.err_rms(mu_resid_all_np,np.sqrt(sigma2_total_np))
                        except:rms_err = -1.0

                    # try:
                        # From Eq. (A.1) of Blondin, Mandel, Kirshner:
                        # var_wrms2 = (2.0/((np.sum(ratio_1))**2.0 ) ) * np.sum(np.square(ratio))
                        # err_wrms = np.sqrt(var_wrms2)/(2.0*rms_weight)
                    # except: err_wrms = -1.0

                    textfile7.write('RMS all (z > %r) = %r (%r WRMS) \n'%(zCMB_Min[j],
                                            round(rms,5), round(rms_weight,5)))
                    textfile7.write('   \n')

                    if ll == label_array[0] and k == 0:
                        textfile_8.write("%.6f  %.6f  %.6f  %.6f # (wRMS, err_wRMS, RMS, RMS_err) \
for all. Case: z > %s.\n"%(rms_weight, err_wrms, rms, rms_err, zCMB_Min[j]))
                    print 'RMS (z > %r, color=%s): %r (%r WRMS)'%(zCMB_Min[j],
                            ColorSamples[k], round(rms,5), round(rms_weight,5))
                    print ' '

                #----------------------------

                mu_resid_CSP_np = np.array(mu_resid_CSP)
                error_mu_CSP_np = np.array(error_mu_CSP)
                sigma2_pec_CSP_np = np.array(sigma2_pec_CSP)

                if len(mu_resid_CSP) > 0:

                    # reset
                    weight_int_np = np.array([]); sigma2_total_np = np.array([]);

                    if wRMS_type == 'wRMS_efit':
                        weight_int_np = 1.0/(error_mu_CSP_np**2.0)
                        sigma2_total_np = error_mu_CSP_np**2.0

                    elif wRMS_type == 'wRMS_efit_ePec_eInt':
                        weight_int_np = 1.0/( (error_mu_CSP_np**2.0) +
                                                            sigma2_pec_CSP_np +
                                                            SimplexResult_z0 )
                        sigma2_total_np = ((error_mu_CSP_np**2.0) +
                                                            sigma2_pec_CSP_np +
                                                            SimplexResult_z0)

                    rms_CSP_weight = RMS.wrms(mu_resid_CSP_np, weight_int_np);
                    rms_CSP = RMS.rms(mu_resid_CSP_np)

                    if errRMS_bootstrap:
                        err_CSP_wrms = RMS.err_wrms_boot(mu_resid_CSP_np, weight_int_np,
                            loopsize);
                        rms_err_CSP = RMS.err_rms_boot(mu_resid_CSP_np,loopsize)
                    else:
                        try:
                            err_CSP_wrms = RMS.err_wrms(mu_resid_CSP_np, weight_int_np,
                                sigma2_total_np);
                        except: err_CSP_wrms = -1.0

                        try:rms_err_CSP=RMS.err_rms(mu_resid_CSP_np,np.sqrt(sigma2_total_np))
                        except: rms_err_CSP = -1.0

#                     try:
#                         # From Eq. (A.1) of Blondin, Mandel, Kirshner:
#                         var_CSP_wrms2 = (2.0/((np.sum(ratio_CSP_1))**2.0 ) )*
#                                         np.sum(np.square(ratio_CSP))
#                         err_CSP_wrms = np.sqrt(var_CSP_wrms2)/(2.0*rms_CSP_weight)
#                     except: err_CSP_wrms = -1.0

                    textfile7.write('RMS CSP (z > %r) = %r (%r WRMS) \n'%(zCMB_Min[j],
                                            round(rms_CSP,5), round(rms_CSP_weight, 5)))
                    textfile7.write('   \n')
                    if ll == label_array[0] and k == 0:
                        textfile_8.write("%.6f  %.6f  %.6f  %.6f # \
(wRMS, err_wRMS, RMS, RMS_err) for CSP subsample. Case: z > %s.\n"%(
                            rms_CSP_weight, err_CSP_wrms, rms_CSP, rms_err_CSP, zCMB_Min[j]))
                    print 'RMS CSP (z > %r, color=%s): %r (%r WRMS)'%(zCMB_Min[j],
                            ColorSamples[k], round(rms_CSP,5), round(rms_CSP_weight, 5))
                    print ' '

                #----------------------------

                mu_resid_CfA_np = np.array(mu_resid_CfA)
                error_mu_CfA_np = np.array(error_mu_CfA)
                sigma2_pec_CfA_np = np.array(sigma2_pec_CfA)

                if len(mu_resid_CfA) > 0:

                    # reset
                    weight_int_np = np.array([]); sigma2_total_np = np.array([]);

                    if wRMS_type == 'wRMS_efit':
                        weight_int_np = 1.0/(error_mu_CfA_np**2.0)
                        sigma2_total_np = error_mu_CfA_np**2.0

                    elif wRMS_type == 'wRMS_efit_ePec_eInt':
                        weight_int_np = 1.0/( (error_mu_CfA_np**2.0) +
                                                            sigma2_pec_CfA_np +
                                                            SimplexResult_z0 )
                        sigma2_total_np = ((error_mu_CfA_np**2.0) +
                                                            sigma2_pec_CfA_np +
                                                            SimplexResult_z0)

                    rms_CfA_weight = RMS.wrms(mu_resid_CfA_np, weight_int_np);
                    rms_CfA = RMS.rms(mu_resid_CfA_np)

                    if errRMS_bootstrap:
                        err_CfA_wrms = RMS.err_wrms_boot(mu_resid_CfA_np, weight_int_np,
                            loopsize);
                        rms_err_CfA = RMS.err_rms_boot(mu_resid_CfA_np,loopsize)
                    else:
                        try:
                            err_CfA_wrms = RMS.err_wrms(mu_resid_CfA_np, weight_int_np,
                                sigma2_total_np);
                        except: err_CfA_wrms = -1.0

                        try:rms_err_CfA=RMS.err_rms(mu_resid_CfA_np,np.sqrt(sigma2_total_np))
                        except:rms_err_CfA = -1.0

#                     try:
#                         # From Eq. (A.1) of Blondin, Mandel, Kirshner:
#                         var_CfA_wrms2 = (2.0/((np.sum(ratio_CfA_1))**2.0 ) )
#                                            * np.sum(np.square(ratio_CfA))
#                         err_CfA_wrms = np.sqrt(var_CfA_wrms2)/(2.0*rms_CfA_weight)
#                     except: err_CfA_wrms = -1.0

                    textfile7.write('RMS CfA (z > %r) = %r (%r WRMS) \n'%(zCMB_Min[j],
                                            round(rms_CfA,5), round(rms_CfA_weight,5)))
                    textfile7.write('   \n')
                    if ll == label_array[0] and k == 0:
                        textfile_8.write("%.6f  %.6f  %.6f  %.6f # (wRMS, err_wRMS, RMS, RMS_err) \
for CfA subsample. Case: z > %s.\n"%(rms_CfA_weight, err_CfA_wrms, rms_CfA, rms_err_CfA, zCMB_Min[j]))
                    print 'RMS CfA (z > %r, color=%s): %r (%r WRMS)'%(zCMB_Min[j],
                            ColorSamples[k], round(rms_CfA,5), round(rms_CfA_weight,5))

                #----------------------------

                mu_resid_Others_np = np.array(mu_resid_Others)
                error_mu_Others_np = np.array(error_mu_Others)
                sigma2_pec_Others_np = np.array(sigma2_pec_Others)

                if len(mu_resid_Others) > 0:

                    # reset
                    weight_int_np = np.array([]); sigma2_total_np = np.array([]);

                    if wRMS_type == 'wRMS_efit':
                        weight_int_np = 1.0/(error_mu_Others_np**2.0)
                        sigma2_total_np = error_mu_Others_np**2.0

                    elif wRMS_type == 'wRMS_efit_ePec_eInt':
                        weight_int_np = 1.0/( (error_mu_Others_np**2.0) +
                                                            sigma2_pec_Others_np +
                                                            SimplexResult_z0 )
                        sigma2_total_np = ((error_mu_Others_np**2.0) +
                                                            sigma2_pec_Others_np +
                                                            SimplexResult_z0)

                    rms_Others_weight = RMS.wrms(mu_resid_Others_np, weight_int_np);
                    rms_Others = RMS.rms(mu_resid_Others_np)


                    if errRMS_bootstrap:
                        err_Others_wrms = RMS.err_wrms_boot(mu_resid_Others_np, weight_int_np,
                            loopsize);
                        rms_err_Others = RMS.err_rms_boot(mu_resid_Others_np,loopsize)
                    else:
                        try:
                            err_Others_wrms = RMS.err_wrms(mu_resid_Others_np, weight_int_np,
                                sigma2_total_np);
                        except: err_Others_wrms = -1.0

                        try:
                            rms_err_Others = RMS.err_rms(mu_resid_Others_np,
                                                 np.sqrt(sigma2_total_np))
                        except: rms_err_Others = -1.0

#                     try:
#                         # From Eq. (A.1) of Blondin, Mandel, Kirshner:
#                         var_Others_wrms2 = (2.0/((np.sum(ratio_Others_1))**2.0 ) ) *
#                                            np.sum(np.square(ratio_Others))
#                         err_Others_wrms = np.sqrt(var_Others_wrms2)/(2.0*rms_Others_weight)
#                     except: err_Others_wrms = -1.0

                    textfile7.write('RMS Others (z > %r) = %r (%r WRMS)  \n'%(zCMB_Min[j],
                                            round(rms_Others,5), round(rms_Others_weight,5)))
                    textfile7.write('   \n')
                    if ll == label_array[0] and k == 0:
                        textfile_8.write("%.6f  %.6f  %.6f  %.6f # (wRMS, err_wRMS, RMS, RMS_err) \
for Others subsamples. Case: z > %s.\n"%(rms_Others_weight, err_Others_wrms, rms_Others, rms_err_Others,
                                         zCMB_Min[j]))
                    print 'RMS Others (z > %r, color=%s): %r (%r WRMS)'%(zCMB_Min[j],
                            ColorSamples[k], round(rms_Others,5), round(rms_Others_weight,5))

                ####################################################

                if plot_raisins == True: plt.xlim(0.2, 0.65)
                else: plt.xlim(xlimPlots)
                plt.ylim(ylimPlots_residual)

                #--- Labeling ----------------------------------------

                # PRINTING THE WEIGHTED RMS AND INTRINSIC DISPERSION IN THE PLOTS.

                RightPlotLimitText_1 = 0.0135 # old: 0.05, 0.057
                if WRMS_label == True and j==0:

                    # ORIGINAL
                    plt.text(xlimPlots[0]+RightPlotLimitText_1-0.005, ylimPlots_residual[1]-0.15,
                             r'%r SN, RMS=%.2f$\pm$%.2f (wRMS = %.2f$\pm$%.2f), $\sigma_{\rm int}$=%.3f$\pm$%.3f'%(
                                 len(mu_resid_all),
                            rms, rms_err, rms_weight, err_wrms, np.sqrt(SimplexResult_z0),
                            error_sigma_pred_z0), fontsize=fontsizePlot-2)

                    # TEMPORAL: for RAISIN2-DD proposal
#                     plt.text(xlimPlots[0]+RightPlotLimitText_1-0.005, ylimPlots_residual[1]-0.15,
#                              r'wRMS = %.2f$\pm$%.2f ($\sigma_{\rm int}$=%.3f$\pm$%.3f)'%(
#                                  rms_weight, err_wrms, np.sqrt(SimplexResult_z0), error_sigma_pred_z0),
#                              fontsize=fontsizePlot-2)

                elif WRMS_label == True and j==1:
                    plt.text(xlimPlots[0]+RightPlotLimitText_1-0.005, ylimPlots_residual[1]-0.15,
                             r'%r SN, RMS=%.2f$\pm$%.2f (wRMS = %.2f$\pm$%.2f), $\sigma_{\rm int}$=%.3f$\pm$%.3f'%(
                                 len(mu_resid_all),
                            rms, rms_err, rms_weight, err_wrms, np.sqrt(SimplexResult_z001),
                            error_sigma_pred_z001), fontsize=fontsizePlot-1)


                RightPlotLimitText_2 = 0.0235 #
                if (WRMS_label == True and WRMS_subsamples==True and
                    ColorSamples[k] and KindOfData4HD == 'AllSamples'):

                    if plot_raisins == True:

                        #  ORIGINAL:
                        plt.text(xlimPlots[0]+RightPlotLimitText_2, ylimPlots_residual[1]-0.25,
                                 '%r DES, RMS=%.2f$\pm$%.2f (%.2f$\pm$%.2f WRMS)'%(countCfA,
                                            rms_CfA, rms_err_CfA, rms_CfA_weight, err_CfA_wrms),
                                 color='red', fontsize=fontsizePlot-4)
                        plt.text(xlimPlots[0]+RightPlotLimitText_2, ylimPlots_residual[1]-0.35,
                                 '%r PS1, RMS=%.2f$\pm$%.2f (%.2f$\pm$%.2f WRMS)'%(countCSP,
                                            rms_CSP, rms_err_CSP, rms_CSP_weight, err_CSP_wrms),
                                 color='blue', fontsize=fontsizePlot-4)

                        # TEMPORAL 2
#                         plt.text(xlimPlots[0]+RightPlotLimitText_2, ylimPlots_residual[1]-0.25,
#                                  'DES, wRMS = %.2f$\pm$%.2f'%(
#                                             rms_CfA_weight, err_CfA_wrms),
#                                  color='red', fontsize=fontsizePlot-2)
#                         plt.text(xlimPlots[0]+RightPlotLimitText_2, ylimPlots_residual[1]-0.35,
#                                  'PS1, wRMS = %.2f$\pm$%.2f'%(
#                                             rms_CSP_weight, err_CSP_wrms),
#                                  color='blue', fontsize=fontsizePlot-2)


                    else: # plot the low-z:
                        plt.text(xlimPlots[0]+RightPlotLimitText_2, ylimPlots_residual[1]-0.25,
                                 '%r CfA, RMS=%.2f$\pm$%.2f (%.2f$\pm$%.2f WRMS)'%(countCfA,
                                rms_CfA, rms_err_CfA, rms_CfA_weight, err_CfA_wrms),
                                 color='red', fontsize=fontsizePlot-4)
                        plt.text(xlimPlots[0]+RightPlotLimitText_2, ylimPlots_residual[1]-0.35,
                                 '%r CSP, RMS=%.2f$\pm$%.2f (%.2f$\pm$%.2f WRMS)'%(countCSP,
                                rms_CSP, rms_err_CSP, rms_CSP_weight, err_CSP_wrms),
                                 color='blue', fontsize=fontsizePlot-4)
                        plt.text(xlimPlots[0]+(RightPlotLimitText_2-0.001), ylimPlots_residual[1]-0.45,
                                 '%r Others, RMS=%.2f$\pm$%.2f (%.2f$\pm$%.2f WRMS)'%(countOthers,
                                rms_Others, rms_err_Others, rms_Others_weight,err_Others_wrms),
                                 color='green', fontsize=fontsizePlot-4)


                RightPlotLimitText_3 = 0.0235  # old: 0.032
                # -NO- PRINTING THE WEIGHTED RMSs IN THE PLOTS.
                if WRMS_label == False and RMS_simple==True:

                    plt.text(xlimPlots[0]+RightPlotLimitText_3, ylimPlots_residual[1]- 0.15, # 0.12,
                             '%r SN, RMS=%r$\pm$%.2f'%(len(mu_resid_all),
                            rms, rms_err), fontsize=fontsizePlot-4)

                    if ColorSamples[k] and KindOfData4HD == 'AllSamples':
                        plt.text(xlimPlots[0]+RightPlotLimitText_3, ylimPlots_residual[1]- 0.25, # 0.20,
                                 '%r CfA, RMS=%.2f$\pm$%.2f'%(countCfA,
                                rms_CfA, rms_err_CfA),
                                 color='red', fontsize=fontsizePlot-4)
                        plt.text(xlimPlots[0]+RightPlotLimitText_3, ylimPlots_residual[1]- 0.35, # 0.28,
                                 '%r CSP, RMS=%.2f$\pm$%.2f'%(countCSP,
                                rms_CSP, rms_err_CSP),
                                 color='blue', fontsize=fontsizePlot-4)
                        plt.text(xlimPlots[0]+RightPlotLimitText_3, ylimPlots_residual[1]- 0.45, # 0.36,
                                 '%r Others, RMS=%.2f$\pm$%.2f'%(countOthers,
                                rms_Others, rms_err_Others),
                                 color='green', fontsize=fontsizePlot-4)

                #-----------------------------------------------

                #--- Print the peculiar velocity, method, chi2_dof: -------

                # ORIGINAL
                # Peculiar velocity uncertainty
                plt.text(xlimPlots[0]+RightPlotLimitText_2, ylimPlots_residual[0]+0.35,
                         r'$\sigma_{\rm pec}$ = %r km/s'%(vpecFix), color='black',
                         fontsize=fontsizePlot-3)

                # chi2_dof: Consistency between error bars and scatter in the data
                # TEMPORAL: commented temporally
                # if zCMB_Min[j] == 0:
                #     plt.text(xlimPlots[0]+RightPlotLimitText_2, ylimPlots_residual[0]+0.25,
                #          r'$\chi^2_{\rm d.o.f.}$ = %.2f'%chi2_dof_HD, color='black',
                #              fontsize=fontsizePlot-3)

                # Method to determine distance moduli
                plt.text(xlimPlots[0]+RightPlotLimitText_2-0.01, ylimPlots_residual[0]+0.15,
                         r'%s'%AppMag_method, color='black', fontsize=fontsizePlot-3)

                #---- Label title and axes ------

                if   BandMax == 'NIRmax': BandMaxText = '$T_{%s}$ max'%BandNameText
                elif BandMax == 'Bmax' or BandMax == 'Bmax_GP':   BandMaxText = '$T_B$ max'
                elif BandMax == 'Snoopy': BandMaxText = 'Snoopy'
                elif BandMax == 'SALT2':  BandMaxText = 'SALT2'


                plt.xlabel('Redshift', fontsize=fontsizePlot)
                plt.ylabel(r'$\mu - \mu_{\rm \Lambda CDM}$', fontsize=fontsizePlot) # ORIGINAL
                # plt.ylabel(r'(B - %s)'%BandNameText, fontsize=fontsizePlot) # TEMPORAL

                #--- Plot Title ---

                if PlotTotalMu == False:

                    if plot_raisins == True:
                        # ORIGINAL:
                        plt.title('%s Hubble Residuals from %s'%(BandNameText, BandMaxText),
                               fontsize=fontsizePlot+2)
                        # TEMPORAL:
                        # plt.title('Preliminary RAISIN Hubble residuals',
                        #       fontsize=fontsizePlot+2)
                    else: # Low-z
                        plt.title('%s-band Hubble Residuals from %s'%(BandNameText, BandMaxText), # ORIGINAL
                        # plt.title('(B-%s) colors, from B and %s max'%(BandNameText,BandNameText), # TEMPORAL
                              fontsize=fontsizePlot+2)

                elif PlotTotalMu == True:
                    if BandsCombination == 'AllBands':
                        if BandMax == 'Bmax' or BandMax == 'Bmax_GP':
                            plt.title('Any YJHK bands Hubble Residuals from %s'%(BandMaxText),
                                      fontsize=fontsizePlot+2)
                        elif BandMax == 'NIRmax':
                            plt.title(r'Any YJHK bands Hubble Residuals from $T_{\rm NIR}$ max',
                                      fontsize=fontsizePlot+2)
                    else:
                        if BandMax == 'Bmax' or BandMax == 'Bmax_GP':
                            plt.title('%s bands Hubble Residuals from %s'%(
                                BandsCombination, BandMaxText),
                                      fontsize=fontsizePlot+2)
                        elif BandMax == 'NIRmax':
                            plt.title(r'%s bands Hubble Residuals from $T_{\rm NIR}$ max'%(
                                BandsCombination), fontsize=fontsizePlot+2)

                #--- Axes ---
                # plt.tick_params(axis='x', labelsize=fontsizePlot)
                # plt.tick_params(axis='y', labelsize=fontsizePlot+2)

                plt.grid(True, ls='--', alpha=0.5)
                # plt.tight_layout()

                #----- Save the figures --------

                # In the file name, append 'text' if the labeling the Hubble diagram outliers:
                if labelOutlierLim < 0.7: textLabel = 'label_'
                else: textLabel = ''

                if   PlotTotalMu == True:  NIRbandsText = BandsCombination+'_'
                elif PlotTotalMu == False: NIRbandsText = ''

                Append_savetext_1 = '' # reset
                if plotErrorBars: Append_savetext_1 = ''
                else: Append_savetext_1 = 'NoBars_'

                if ColorSamples[k] and KindOfData4HD == 'AllSamples':
                    plt.savefig(DirSaveOutput+'Plot_Residual_%s_%s_%s_%s_%s_%s_Color_%s%s%s.png'%(
                        KindOfData4HD, KindOfTemp,
                                KindOfTempSubgroup, chi2_dof_Max_Label,
                        residualMax_Label,
                        zCMB_Min_Label[j],
                                NIRbandsText, textLabel,Append_savetext_1),
                                dpi=ResolutionPlot_HD)

                    # EPS format
#                     plt.savefig(DirSaveOutput+'Plot_Residual_%s_%s_%s_%s_%s_%s_Color_%s%s%s.eps'%(
#                         KindOfData4HD, KindOfTemp,
#                                 KindOfTempSubgroup, chi2_dof_Max_Label,
#                         residualMax_Label,
#                         zCMB_Min_Label[j],
#                                 NIRbandsText, textLabel,Append_savetext_1),format='eps')

                else:
                    plt.savefig(DirSaveOutput+'Plot_Residual_%s_%s_%s_%s_%s_%s_%s%s%s.png'%(
                        KindOfData4HD, KindOfTemp,
                                KindOfTempSubgroup, chi2_dof_Max_Label,
                        residualMax_Label,
                        zCMB_Min_Label[j],
                                NIRbandsText, textLabel,Append_savetext_1),
                                dpi=ResolutionPlot_HD)
                plt.close()

textfile7.close()
textfile_8.close()

#----------------------

if ScriptVersion == 'notebook': print "\n # All done."



# In[88]:


plt.close();plt.close();plt.close();

textfile7.close();textfile7.close();textfile7.close();
textfile_8.close();textfile_8.close();textfile_8.close();


# In[89]:


print "\n#       All done smoothly\n"
print "############################################################\n\n"


# In[ ]:




