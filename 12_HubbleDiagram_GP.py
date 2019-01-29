#!/usr/bin/env python
# coding: utf-8

# # Hubble diagram from Gaussian Processes fitted apparent magnitude
#
# I have to run this notebook 3 times:
#
# 1) The first time to determine:
# - the apparent magnitudes (and their uncertainties, sigma_m) at t_NIRmax infered from the GP fit
# - uncertainty in the photometric distance moduli of the SNe sample defined as:
#     error_mu_photometric^2 = sigma_m^2
# where sigma_m is the uncertainty in the apparent magnitude at t_Bmax computed in this first run.
# - the absolute magnitude for each SN from AbsMag = appMag - mu_LCDM(z)
# - the mean absolute magnitude at t_Bmax, mean_AbsMag, with its standard deviation of the SNe sample. Write down these values in (AverageAbsMag_atMax, err_AverageAbsMag_atMax).
#
# Run this notebook until the end of section "Determining average Absolute magnitude of the sample".
# Set: AbsMagFromHisto = False
#
# 2) The second time to determine:
# - the photometric distance moduli of the SNe sample defined as mu_photometric = appMag - mean_AbsMag
# - the intrinsic dispersion of the Hubble-diagram residual, sigma_intrinsic.
#
# Run the entire notebook
# Set: AbsMagFromHisto = True
#
#     NOTES
# - When plotting the Hubble diagram and the residuals plots and computing the intrisic dispersion using 11_DistanceMu_HubbleDiagram_*.ipynb, it is NOT necessary to set the values of (AverageNIRMax_AbsMag, err_AverageNIRMax_AbsMag ) there, because in that section it is not necessary these values.
#
#

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import quad as intquad
from shutil import copyfile # To copy/paste data files.

import os
import glob # To read the name of the files in a given directory

#--------------------------------------------------------60
code_created_by = 'Arturo_Avelino'
# On date: 2017.01.10 (yyyy.mm.dd)
code_name = '12_HubbleDiagram_GP.ipynb'
version_code = '0.1.16'
last_update = '2019.01.28'
#--------------------------------------------------------60

##############################################################################80

# # User interface

HoFix = 73.24 # 72  # Hubble constant (km/(s Mpc))
# HoFix = 72.78 # # TEMPORAL: value reported by Dhawan et al 2017.

# Peculiar velocity. km/s. Options (150, 200, 300, 400)
# Must be integer number.
vpecFix = 250

Band = 'K'    # What band to fit:(Y, J, H ,K)

# Band to use as the reference peak maximum? Options: (NIRmax, Bmax)
BandMax = 'Bmax'

# Mean Absolute magnitude determined from histogram of 'appMagTmax_s - mu_s'?:
# First run the notebook with this option with setting the value to "False"
# in order to determine the mean abs mag, then once I get that number, run a
# second time the notebook with this option as "True" to compute the final
# tables and results.
# If "True" it is generated the final 'DistanceMu_Good_AfterCutoffs_Main_.txt'
# text file, otherwise if "False" generate temporal text files.

AbsMagFromHisto = True

#-----------------------------------------------------------------------------80
# zcmb upper cutoff. Options: (0.04, 0.06, 0.09)
zcmb_Max = 0.04

#-- EBV cutoff
EBVhostMin = -0.4 # -0.4 # host galaxy
EBVhostMax = 0.4 # 0.4 # host galaxy.
EBVMWLim = 1.0 # Milky-Way galaxy

#-- dm15 cutoff
dm15LowerLim = 0.8 # I assume 0.79.
dm15UpperLim = 1.6

#-----------------------------------------------------------------------------80
#--- Fixed values ---
OmMFix = 0.28 # Omega_Matter
OmLFix = 0.72 # Omega_Lambda
wFix = -1 # Dark energy EoS
c = 299792.458  # Speed of light (km/s)
cc = 299792.458  # Speed of light (km/s)

# Characteristics of the fitting for each band

# Directory where the "Settings_GPFit_.txt" file is located.
# From it I read the value of the kernel 'length' hyperpavrameter.
# DirSettings = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/\
# 10Compute/TheTemplates/%s_band/Std_filters/2_Selection_FlatPrior/\
# AllSamples_appMag/Goods'%Band

#-------------------
# Light-curve types to use to construct the Hubble diagram:
# KindOfData4HD = 'CfA'
# KindOfData4HD = 'CSP'
# KindOfData4HD = 'Others'
KindOfData4HD = 'AllSamples'

#================================================================

#    J band

if Band == 'J':

    # NOTE: The values of "AverageAbsMag_atMax" and
    # "err_AverageAbsMag_atMax" do NOT
    # depend on the value of "vpecFix".

    if BandMax == 'NIRmax':

        if zcmb_Max == 0.04:

            # ------------------------------
            #   J band | Band Max = NIRmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.5169458627;  # 2018-07-12; 23:47 hrs.
            # err_AverageAbsMag_atMax = 0.14124319086;

            # ------------------------------
            #   J band | Band Max = NIRmax | 0 < z < 0.04
            AverageAbsMag_atMax = -18.515116;  # 2018-09-29; 13:02 hrs.
            err_AverageAbsMag_atMax = 0.139428;

        elif zcmb_Max == 0.06:

            # Using Foley+Cepheid+Special cases  z_CMB values
            AverageAbsMag_atMax = -18.5527768929;
            err_AverageAbsMag_atMax = 0.143477501392;

        # Search for the NIR maximum app mag of the light curve
        # in the following rest-frame days range:
        MinPhase = -8.5   # restframe days after T_Bmax
        MaxPhase = 1   # restframe days after T_Bmax

        # Discard if the maximum is located at phase >= phaseB_atMax_Cutoff:
        phaseB_atMax_Cutoff = -2.5

        # Ok if the light curve is truncated (instead of a maximum)
        # at phaseB =< phaseB_truncatedLC days.
        phaseB_truncatedLC = phaseB_atMax_Cutoff

    elif BandMax == 'Bmax':

        if zcmb_Max == 0.04:

            # ------------------------------
            #   J band | Band Max = Bmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.3728286471;  # 2018-06-21; 14:46 hrs.
            # err_AverageAbsMag_atMax = 0.158568172895;

            # ------------------------------
            #   J band | Band Max = Bmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.342148902;  # 2018-07-13; 01:54 hrs.
            # err_AverageAbsMag_atMax = 0.152401896151;

            # ------------------------------
            #   J band | Band Max = Bmax | 0 < z < 0.04
            AverageAbsMag_atMax = -18.341034;  # 2018-09-30; 18:42 hrs.
            err_AverageAbsMag_atMax = 0.151497;

        elif zcmb_Max == 0.06:

            # Using Foley+Cepheid+Special cases  z_CMB values
            AverageAbsMag_atMax = -18.4100281935;
            err_AverageAbsMag_atMax = 0.192207626319;

        # Find the NIR maximum app mag of the light curve
        # in a given rest-frame days range. For B max there is not
        # a search for the maximum because I've already determine
        # the phase = 0 from the optical bands.
        MinPhase = -1.   # restframe days after T_Bmax
        MaxPhase = 1.  # restframe days after T_Bmax

        # Discard if the maximum is located at phase >= phaseB_atMax_Cutoff:
        phaseB_atMax_Cutoff = 0.5

        # Ok if the light curve is truncated (instead of a maximum)
        # at phaseB =< phaseB_truncatedLC days.
        phaseB_truncatedLC = phaseB_atMax_Cutoff

#================================================================

#    Y band

elif Band == 'Y':

    # NOTE: The values of "AverageAbsMag_atMax" and
    # "err_AverageAbsMag_atMax" do NOT
    # depend on the value of "vpecFix".

    if BandMax == 'NIRmax':

        if zcmb_Max == 0.04:

            # ------------------------------
            #   Y band | Band Max = NIRmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.392380;  # 2018-08-09; 18:54 hrs.
            # err_AverageAbsMag_atMax = 0.111230;

            # ------------------------------
            #   Y band | Band Max = NIRmax | 0 < z < 0.04
            AverageAbsMag_atMax = -18.392629;  # 2018-09-30; 17:34 hrs.
            err_AverageAbsMag_atMax = 0.111422;

        elif zcmb_Max == 0.06:

            AverageAbsMag_atMax =  nan;
            err_AverageAbsMag_atMax = nan;

        # Search for the NIR maximum app mag of the light curve
        # in the following rest-frame days range:
        MinPhase = -8.5   # restframe days after T_Bmax
        MaxPhase = 1   # restframe days after T_Bmax

        # Discard if the maximum is located at phase >= phaseB_atMax_Cutoff:
        phaseB_atMax_Cutoff = -2.5

        # Ok if the light curve is truncated (instead of a maximum)
        # at phaseB =< phaseB_truncatedLC days.
        phaseB_truncatedLC = phaseB_atMax_Cutoff

    elif BandMax == 'Bmax':

        if zcmb_Max == 0.04:

            # ------------------------------
            #   Y band | Band Max = Bmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.2358845862;  # 2018-06-21; 14:48 hrs.
            # err_AverageAbsMag_atMax = 0.141701608988;

            # ------------------------------
            #   Y band | Band Max = Bmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.1630574138;  # 2018-07-13; 01:57 hrs.
            # err_AverageAbsMag_atMax = 0.11783879715;

            # ------------------------------
            #   Y band | Band Max = Bmax | 0 < z < 0.04
            AverageAbsMag_atMax = -18.162957;  # 2018-09-30; 18:49 hrs.
            err_AverageAbsMag_atMax = 0.117983;

        elif zcmb_Max == 0.06:

            AverageAbsMag_atMax = nan;
            err_AverageAbsMag_atMax = nan;

        # Find the NIR maximum app mag of the light curve
        # in a given rest-frame days range. For B max there is not
        # a search for the maximum because I've already determine
        # the phase = 0 from the optical bands.
        MinPhase = -1.   # restframe days after T_Bmax
        MaxPhase = 1.  # restframe days after T_Bmax

        # Discard if the maximum is located at phase >= phaseB_atMax_Cutoff:
        phaseB_atMax_Cutoff = 0.5

        # Ok if the light curve is truncated (instead of a maximum)
        # at phaseB =< phaseB_truncatedLC days.
        phaseB_truncatedLC = phaseB_atMax_Cutoff

#================================================================

#    H band

elif Band == 'H':

    # NOTE: The values of "AverageAbsMag_atMax" and
    # "err_AverageAbsMag_atMax" do NOT
    # depend on the value of "vpecFix".

    if BandMax == 'NIRmax':

        if zcmb_Max == 0.04:

            # ------------------------------
            #   H band | Band Max = NIRmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.3023154545;  # 2018-07-13; 00:19 hrs.
            # err_AverageAbsMag_atMax = 0.113290021133;

            # ------------------------------
            #   H band | Band Max = NIRmax | 0 < z < 0.04
            AverageAbsMag_atMax = -18.302525;  # 2018-09-30; 17:49 hrs.
            err_AverageAbsMag_atMax = 0.113689;

        elif zcmb_Max == 0.06:

            # Using Foley+Cepheid+Special cases  z_CMB values.
            AverageAbsMag_atMax = nan;
            err_AverageAbsMag_atMax = nan;

        # Search for the NIR maximum app mag of the light curve
        # in the following rest-frame days range:
        MinPhase = -8.5   # restframe days after T_Bmax
        MaxPhase = 1   # restframe days after T_Bmax

        # Discard if the maximum is located at phase >= phaseB_atMax_Cutoff:
        phaseB_atMax_Cutoff = -2.5

        # Ok if the light curve is truncated (instead of a maximum)
        # at phaseB =< phaseB_truncatedLC days.
        phaseB_truncatedLC = phaseB_atMax_Cutoff

    elif BandMax == 'Bmax':

        if zcmb_Max == 0.04:

            # ------------------------------
            #   H band | Band Max = Bmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.1691600889;  # 2018-06-21; 14:50 hrs.
            # err_AverageAbsMag_atMax = 0.1197098627;

            # ------------------------------
            #   H band | Band Max = Bmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.1869583409;  # 2018-07-13; 03:56 hrs.
            # err_AverageAbsMag_atMax = 0.120716918869;

            # ------------------------------
            #   H band | Band Max = Bmax | 0 < z < 0.04
            AverageAbsMag_atMax = -18.186980;  # 2018-09-30; 18:56 hrs.
            err_AverageAbsMag_atMax = 0.120082;

        elif zcmb_Max == 0.06:

            # Using Foley+Cepheid+Special cases  z_CMB values.
            AverageAbsMag_atMax = nan;
            err_AverageAbsMag_atMax = nan;

        # Find the NIR maximum app mag of the light curve
        # in a given rest-frame days range. For B max there is not
        # a search for the maximum because I've already determine
        # the phase = 0 from the optical bands.
        MinPhase = -1.   # restframe days after T_Bmax
        MaxPhase = 1.  # restframe days after T_Bmax

        # Discard if the maximum is located at phase >= phaseB_atMax_Cutoff:
        phaseB_atMax_Cutoff = 0.5

        # Ok if the light curve is truncated (instead of a maximum)
        # at phaseB =< phaseB_truncatedLC days.
        phaseB_truncatedLC = phaseB_atMax_Cutoff

#================================================================

#    K band

elif Band == 'K':

    # NOTE: The values of "AverageAbsMag_atMax" and
    # "err_AverageAbsMag_atMax" do NOT
    # depend on the value of "vpecFix".

    if BandMax == 'NIRmax':

        if zcmb_Max == 0.04:

            # ------------------------------
            #   K band | Band Max = NIRmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.4352868462;  # 2018-05-02; 17:02 hrs.
            # err_AverageAbsMag_atMax = 0.17068523975;

            # ------------------------------
            #   K band | Band Max = NIRmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.3746893846;  # 2018-07-13; 00:32 hrs.
            # err_AverageAbsMag_atMax = 0.184788076388;

            # ------------------------------
            #   K band | Band Max = NIRmax | 0 < z < 0.04
            AverageAbsMag_atMax = -18.368747;  # 2018-09-30; 17:57 hrs.
            err_AverageAbsMag_atMax = 0.179273;

        elif zcmb_Max == 0.06:

            # Using Foley+Cepheid+Special cases  z_CMB values.
            AverageAbsMag_atMax = nan ;
            err_AverageAbsMag_atMax = nan;

        # Search for the NIR maximum app mag of the light curve
        # in the following rest-frame days range:
        MinPhase = -9   # restframe days after T_Bmax
        MaxPhase = 1   # restframe days after T_Bmax

        # Discard if the maximum is located at phase >= phaseB_atMax_Cutoff:
        phaseB_atMax_Cutoff = 0

        # Ok if the light curve is truncated (instead of a maximum)
        # at phaseB =< phaseB_truncatedLC days.
        phaseB_truncatedLC = phaseB_atMax_Cutoff

    elif BandMax == 'Bmax':

        if zcmb_Max == 0.04:

            # ------------------------------
            #   K band | Band Max = Bmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.3489594615;  # 2018-06-21; 14:53 hrs.
            # err_AverageAbsMag_atMax = 0.170252048275;

            # ------------------------------
            #   K band | Band Max = Bmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.2878237647;  # 2018-07-13; 02:01 hrs.
            # err_AverageAbsMag_atMax = 0.202175569;

            # ------------------------------
            #   K band | Band Max = Bmax | 0 < z < 0.04
            # AverageAbsMag_atMax = -18.283445;  # 2018-07-13; 04:00 hrs.
            # err_AverageAbsMag_atMax = 0.177039421209;

            # ------------------------------
            #   K band | Band Max = Bmax | 0 < z < 0.04
            AverageAbsMag_atMax = -18.283624;  # 2018-09-30; 19:04 hrs.
            err_AverageAbsMag_atMax = 0.169889;

        elif zcmb_Max == 0.06:

            # Using Foley+Cepheid+Special cases  z_CMB values.
            AverageAbsMag_atMax = nan ;
            err_AverageAbsMag_atMax = nan;

        # Find the NIR maximum app mag of the light curve
        # in a given rest-frame days range. For B max there is not
        # a search for the maximum because I've already determine
        # the phase = 0 from the optical bands.
        MinPhase = -1.   # restframe days after T_Bmax
        MaxPhase = 1.  # restframe days after T_Bmax

        # Discard if the maximum is located at phase >= phaseB_atMax_Cutoff:
        phaseB_atMax_Cutoff = 0.5

        # Ok if the light curve is truncated (instead of a maximum)
        # at phaseB =< phaseB_truncatedLC days.
        phaseB_truncatedLC = phaseB_atMax_Cutoff

##############################################################################80

# # Automatic

# #### Get the name of this ipython notebook
# To print it in the output text files as reference

# %%javascript
# var kernel = IPython.notebook.kernel;
# var thename = window.document.getElementById("notebook_name").innerHTML;
# var command = "NotebookName = " + "'"+thename+".ipynb"+"'";
# kernel.execute(command);

# print '#', (NotebookName)
# 12_HubbleDiagram_GP_v1_10.ipynb

# Get the current date and time
import datetime

# Read the time and date now
now = datetime.datetime.now()

shiftNum = 71

#- Function to convert from index to days (phase).
def index2day(index):
    day = (index-shiftNum)/2.
    return day

#- Function to convert from days (phase) to index.
def day2index(day):
    index = 2*day + shiftNum
    return index

# print 'Testing the functions:', index2day(105), ',', day2index(-11)
# Testing the functions: 17.0 , 49

#-------------------
# Defining the minimal/maximal rows based on the days:

MinRow = int(day2index(MinPhase)) # 54 # (53 = -9 days) (54 = -8.5 days)
MaxRow = int(day2index(MaxPhase)) # 66 # (75 = 2 days), (66 = -2.5 days),
# (65 = -3 days), (64 = -3.5 days), (63 = -4 days)

print MinRow, MaxRow
# 53 73

# ### $\mu_{\Lambda{\rm CDM}}$

# Inverse of the dimensionless Hubble parameter
def InvEHubblePar(z, OmM, wde):
    "Dimensionless Hubble parameter"
    InvEHubbleParInt = 1.0/(np.sqrt(OmM*(1.0+z)**3.0 + (1.0-OmM)*(1.+z)**(3.*(1.+wde))))
    return InvEHubbleParInt

# ---- The luminosity distance ----
def LumDistanceVec(z, OmM, wde, Ho):
    "Luminosity distance"
    LumDistanceVecInt = 0.
    LumDistanceVecInt = cc*(1.+z)*intquad(InvEHubblePar, 0., z, args=(OmM, wde))[0]/Ho
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
ztest1 = 0.01

print 'Checking that the functions work well:', DistanceMu(ztest1, OmMFix, wFix, HoFix)
# Checking that the functions work well: 33.1141460988 # Ho=72
# Checking that the functions work well: 33.0773926577 # Ho=73.24

# #### Function to identify string or number

# Function to identify if a string is an integer number or a letter.
# This will be used in the dictionary construction to properly read some SN names.

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

# Tests
print is_number('5'), is_number('e')
# True False

# ### For the instrinsic dispersion computation¶
#
# ##### Sigma peculiar

# sigma^2_mu from the peculiar velocity uncertainty
# This function is used to determine in the sections "Intrinsic dispersion" and "Optical RMS", to
# determine the intrinsic dispersion.

def sigma2_pec(zcmb, err_zcmb, vpec):
    sigma2_pecInt = ((5/(zcmb*np.log(10)))*np.sqrt((vpec/cc)**2 + err_zcmb**2))**2
    return sigma2_pecInt

# Test
sigma2_pec(0.0109942726, 0.0010420420, 150)
# 0.052125037354329877

# Likelihood to determine the intrinsic dispersion

# -2ln(Likelihood) function
# Eq. (B.6) of Blondin et al 2011

# 'sigma2Pred' is the square of the intrinsic dispersion.

def neg2lnLike(sigma2Pred, mu_resid_np, sigma2_pec_np):
    sum1 = 0
    for i in range(len(mu_resid_np)):
        sum1 = (sum1 + np.log(sigma2Pred + sigma2_pec_np[i]) +
                (mu_resid_np[i]**2)/(sigma2Pred + sigma2_pec_np[i]) )
    return sum1

# For the case of using a normalized template. This is the full Eq. (B.6) of Blondin et al 2011
def neg2lnLikeFull(sigma2Pred, mu_resid_np, sigma2_pec_np, sigma2_appmagTBmax_np):
    sum1 = 0
    for i in range(len(mu_resid_np)):
        sum1 = (sum1 + np.log(sigma2_appmagTBmax_np[i] + sigma2Pred + sigma2_pec_np[i]) +
                (mu_resid_np[i]**2)/(sigma2_appmagTBmax_np[i] + sigma2Pred + sigma2_pec_np[i]) )
    return sum1

# Test
# neg2lnLike(0.15, mu_resid_z0_np, sigma2_pec_z0_np)
# -98.297121896077712

# Finding the uncertainty in the determination of sigma_int

# Define the Fisher information matrix
# Eq. (B.7) of Blondin et al 2011

def FisherFunc(sigma2Pred, mu_resid_np, sigma2_pec_np):
    sum2 = 0
    for i in range(len(mu_resid_np)):
        sum2 = (sum2 + (mu_resid_np[i]**2)/(sigma2Pred + sigma2_pec_np[i])**3 -
               1/(2*(sigma2Pred + sigma2_pec_np[i])**2)  )
    return sum2

# For the case of using a normalized template. This is the full Eq. (B.7) of Blondin et al 2011
def FisherFuncFull(sigma2Pred, mu_resid_np, sigma2_pec_np, sigma2_appmagTBmax_np):
    sum2 = 0
    for i in range(len(mu_resid_np)):
        sum2 = (sum2 + (mu_resid_np[i]**2)/(sigma2_appmagTBmax_np[i] + sigma2Pred + sigma2_pec_np[i])**3 -
               1/(2*(sigma2_appmagTBmax_np[i] + sigma2Pred + sigma2_pec_np[i])**2)  )
    return sum2

# Test
# FisherFunc(0.049875, mu_resid_z0_np, sigma2_pec_z0_np)

#-----------------------------------------------------------------------------80

# ### Cepheids

DirSNeWithCepheid = '/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/AndyLCComp_2018_02/MyNotes/'

# From the Cepheid SNe list the only part that I use in the entire
# notebook is the first column, i.e., the SN Name column.
ListSNeCepheid = np.genfromtxt(DirSNeWithCepheid+
                               'SNeWithCepheidDistances.txt', dtype=['S10',
                                float,float,float,float,float,float])

print "# %s SNe with Cepheid distances in Andy's compilation."%len(ListSNeCepheid['f0'])
ListSNeCepheid['f0']

# Defining directories

MainDir = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/TheTemplates/'

DirGPFits = MainDir+Band+'_band/Std_filters/2_Selection_FlatPrior/%s_appMag_vpec_0/Goods/'%(KindOfData4HD)

DirSaveOutput = MainDir+Band+'_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/'+KindOfData4HD+'/vpec%s_%s/'%(vpecFix,BandMax)

DirSavePlots = DirSaveOutput+'/Plots_GPFits/'

#- Force the creation of the directory to save the outputs.
#- "If the subdirectory does not exist then create it"
import os # To use command line like instructions
if not os.path.exists(DirSaveOutput): os.makedirs(DirSaveOutput)
if not os.path.exists(DirSavePlots): os.makedirs(DirSavePlots)

import os # To use command line like instructions
import glob # To read the files in my directory

# Change the working directory where the data files are located
os.chdir(DirGPFits)

#--- Reading the data files in the app mag folder
if KindOfData4HD == 'CfA':
    list_SNe2 = glob.glob('*'+'CfA_'+Band+'.txt')
elif KindOfData4HD == 'CSP':
    list_SNe2 = glob.glob('*'+'CSP_'+Band+'.txt')
elif KindOfData4HD == 'Others':
    list_SNe2 = glob.glob('*'+'Others_'+Band+'.txt')
elif KindOfData4HD == 'AllSamples':
    list_SNe2 = glob.glob('*'+Band+'.txt')

"""
#--- Reading the data files in the app mag folder
if KindOfData4HD == 'CfA':
    list_SNe2 = glob.glob('*'+'CfA_'+Band+'_GP_mean_sigma_Filled.dat')
elif KindOfData4HD == 'CSP':
    list_SNe2 = glob.glob('*'+'CSP_'+Band+'_GP_mean_sigma_Filled.dat')
elif KindOfData4HD == 'Others':
    list_SNe2 = glob.glob('*'+'Others_'+Band+'_GP_mean_sigma_Filled.dat')
elif KindOfData4HD == 'AllSamples':
    list_SNe2 = glob.glob('*'+Band+'_GP_mean_sigma_Filled.dat')
"""

print 'Number of -%s- SNe with GP fit: %s'%(KindOfData4HD, len(list_SNe2))

#-----------------------------------------------------------------------------80

#     CREATE &/or READ THE NAME OF THE SNe FROM THE TEXT FILE

# Change the working directory where the data files are located
os.chdir(DirSaveOutput)

# Reading the data files in that folder
Textfiles = glob.glob('ListSNe_AppMag_*.txt')

if 'ListSNe_AppMag_Notes_.txt' in Textfiles:
    list_SNe = np.genfromtxt(DirSaveOutput+'ListSNe_AppMag_Notes_.txt', dtype=['S35', int])
    print '# Reading: <ListSNe_AppMag_Notes_.txt>'

elif 'ListSNe_AppMag_.txt' in Textfiles:
    list_SNe = np.genfromtxt(DirSaveOutput+'ListSNe_AppMag_.txt', dtype=['S35', int])
    print '# Reading: <ListSNe_AppMag_.txt>'

else: # if ListSNe_AppMag_.txt doesn't exist, then create it and read it.
    print '# ListSNe_AppMag_.txt not found. Creating it.'
    # Create a list text file with the name of the SNe of the apparent mag LCs.
    list_SNe_file = open(DirSaveOutput+'ListSNe_AppMag_.txt', 'w')
    list_SNe_file.write('#        SN list with GP fit \n')
    list_SNe_file.write('# Note that these SNe have passed -a broader version- of my cutoffs. \n')
    list_SNe_file.write('# These SNe are ALL those located in: \n')
    list_SNe_file.write('# %s \n'%DirGPFits)

    #------
    now = datetime.datetime.now() # Read the time and date right now
    text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd).")
    text_Date   = '# On date: %s \n'%text_timenow
    text_Author = '# Data table created by: Arturo Avelino \n'
    text_script = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
    text_line = '#'+'-'*50 + '\n'

    list_SNe_file.write(text_line);
    list_SNe_file.write(text_Author); list_SNe_file.write(text_Date); list_SNe_file.write(text_script);
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
    """
    for name in list_SNe2:
        list_SNe_file.write('%-32s  0 \n'%name)
        # list_SNe_file.write(name)
        # print name


    list_SNe_file.write('# Number of SNe in this list = %s \n'%(len(list_SNe2)))
    list_SNe_file.close()
    """
    #------------------------

    # Read the list I've just created:

    list_SNe = np.genfromtxt(DirSaveOutput+'ListSNe_AppMag_.txt', dtype=['S35', int])
    print '# Reading: < ListSNe_AppMag_.txt >'

print '# Number of SNe in the list =', len(list_SNe)
print '# Number of masked SNe in the list: ', sum(list_SNe['f1'])

# ListSNe_AppMag_.txt not found. Creating it.
# Number of SNe in the list = 132
# Number of masked SNe in the list:  0

##############################################################################80

# ### Main loop
#
# #### Compute apparent magnitude and distance modulus

# Using the GP fitting to the apparent magnitude light curves

# Create a list of (MJD_NIRmax, appMag_NIRmax, error_appMag_NIRmax), in
# addition to some other relevant quantities.
os.chdir(DirGPFits)

# Mean Absolute magnitude determined from histogram of 'appMagTmax_s - mu_s'?:
# If so it is generated the final 'DistanceMu_Good_AfterCutoffs_Main_.txt' text file, otherwise
# generate temporal text files.
if  AbsMagFromHisto == True: textPrefix = ''; underline = ''
else: textPrefix = 'TMP'; underline = '_' # 'TMP' stands for 'temporal' file.

file_results = open(DirSaveOutput+'DistanceMu_Good_AfterCutoffs_Main_%s%s.txt'%(
    textPrefix, underline), 'w')
file_all = open(DirSaveOutput+'DistanceMu_All_BeforeCutoffs_%s%s.txt'%(
    textPrefix, underline), 'w')

#--------------------------------------------------
file_results.write('# Apparent magnitudes and distance moduli computed from the GP fit directly. \n')
file_results.write('# <AverageAbsMag_atMax> = %s +/- %s used. \n'%(AverageAbsMag_atMax, err_AverageAbsMag_atMax))
file_results.write('# Peculiar velocity uncertainty assumed: %s km/s  \n'%vpecFix)
file_results.write('# Band used as reference of maximum: %s \n'%BandMax)
file_results.write('# Phase_Bmax range to look for the maximum: Minimum = %s days (%s rows), \
Maximum = %s days (%s rows) \n'%(MinPhase, MinRow, MaxPhase, MaxRow))
file_results.write('# Discard if the maximum is located at phase > phaseB_atMax_Cutoff = %s days. \n'%phaseB_atMax_Cutoff)
file_results.write('# Ok if the light curve is truncated (instead of a maximum) at phaseB =< \
phaseB_truncatedLC = %s days. \n'%phaseB_truncatedLC)
file_results.write('# Upper limit zcmb cutoff = %s \n'%zcmb_Max)
# file_results.write('# Intrinsic dispersion for the case (0 < z < %s) and used to obtain the \
# total photometric distance modulus uncertainty: \n'%zcmb_Max)
# file_results.write('# Intrinsic dispersion = %s \n'%IntrinsicDisp)

file_results.write('# phaseB_truncatedLC = %s days \n'%phaseB_truncatedLC)
file_results.write("# \n")
file_results.write("# Cutoffs: \n")
file_results.write("# z_cmb<%s | %s<dm15<%s | %s<EBV_host<%s | EBV_MW<%s \n"%(
    zcmb_Max,dm15LowerLim, dm15UpperLim, EBVhostMin, EBVhostMax, EBVMWLim))

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd).")
text_Author = '# Data table created by: Arturo Avelino \n'
text_Date   = '# On date: %s \n'%text_timenow
text_script = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
text_line = '#'+'-'*80 + '\n'

file_results.write(text_line)
file_results.write(text_Author); file_results.write(text_Date);
file_results.write(text_script);
file_results.write(text_line)

file_results.write('#  SN name                      zcmb          err_zcmb      mu_GP      err_mu_GP \
mu_GP_residual  chi2_dof FlatCode  NIRMax_appMag  ErrNIRMax_appMag   mu_LCDM  sigma_pecInt  NIRMax_AbsMag  \
ErrNIRMax_AbsMag  MJD_NIRmax   err_MJD_Bmax   phaseB_NIRMax  zhelio        err_zhelio       dm15      \
err_dm15   EBVhost  err_EBVhost  EBV_MW    err_EBV_MW  Alamb    err_Alamb   R_F     mu_snoopy  \
err_mu_snoopy    MJD_Bmax     err_MJD_Bmax   BMax_appMag  ErrBMax_appMag \n')

#-----------------------------------------------------------------------------80
# Copy/paste of the text above.

file_all.write('# Apparent magnitudes and distance moduli computed from the GP fit directly. \n')
file_all.write('# <AverageAbsMag_atMax> = %s +/- %s used. \n'%(AverageAbsMag_atMax, err_AverageAbsMag_atMax))
file_all.write('# Peculiar velocity uncertainty assumed: %s km/s  \n'%vpecFix)
file_all.write('# Band used as reference of maximum: %s \n'%BandMax)
file_all.write('# Phase_Bmax range to look for the maximum: Minimum = %s days (%s rows), \
Maximum = %s days (%s rows) \n'%(MinPhase, MinRow, MaxPhase, MaxRow))
file_all.write('# Discard if the maximum is located at phase > phaseB_atMax_Cutoff = %s days. \n'%phaseB_atMax_Cutoff)
file_all.write('# Ok if the light curve is truncated (instead of a maximum) at phaseB =< \
phaseB_truncatedLC = %s days. \n'%phaseB_truncatedLC)
file_all.write('# Upper limit zcmb cutoff = %s \n'%zcmb_Max)
# file_all.write('# Intrinsic dispersion for the case (0 < z < %s) and used to obtain the \
# total photometric distance modulus uncertainty: \n'%zcmb_Max)
# file_all.write('# Intrinsic dispersion = %s \n'%IntrinsicDisp)

file_all.write('# phaseB_truncatedLC = %s days \n'%phaseB_truncatedLC)
file_all.write("# \n")
file_all.write("# Cutoffs: \n")
file_all.write("# z_cmb<%s | %s<dm15<%s | %s<EBV_host<%s | EBV_MW<%s \n"%(
    zcmb_Max,dm15LowerLim, dm15UpperLim, EBVhostMin, EBVhostMax, EBVMWLim))

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd).")
text_Author = '# Data table created by: Arturo Avelino \n'
text_Date   = '# On date: %s \n'%text_timenow
text_script = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
text_line = '#'+'-'*80 + '\n'

file_all.write(text_line)
file_all.write(text_Author); file_all.write(text_Date);
file_all.write(text_script);
file_all.write(text_line)

file_all.write('#  SN name                      zcmb          err_zcmb      mu_GP      err_mu_GP \
mu_GP_residual  chi2_dof FlatCode  NIRMax_appMag  ErrNIRMax_appMag   mu_LCDM  sigma_pecInt  NIRMax_AbsMag  \
ErrNIRMax_AbsMag  MJD_NIRmax   err_MJD_Bmax   phaseB_NIRMax  zhelio        err_zhelio       dm15      \
err_dm15   EBVhost  err_EBVhost  EBV_MW    err_EBV_MW  Alamb    err_Alamb   R_F     mu_snoopy  \
err_mu_snoopy    MJD_Bmax     err_MJD_Bmax   BMax_appMag  ErrBMax_appMag \n')

#-----------------------------------------------------------------------------80
# The loop

# Reset variables and arrays
countTotal = 0;
countCommented = 0; countNoCommented = 0;
countCommented_cutoffs = 0;
mu_resid_array = []; SNeName_list = [];

for name, mask in list_SNe: # Loop over SNe.
    # print name
    SN_GPfit  = np.genfromtxt(name[:-4]+'_GP_mean_sigma_Filled.dat')
    SN_LCdata = np.genfromtxt(name)

    zcmb      = SN_LCdata[0,0] # CMB redshift
    err_zcmb  = SN_LCdata[0,1] # CMB redshift
    zhelio    = SN_LCdata[0,2] # helio redshift
    mu_LCDM   = SN_LCdata[2,2] # distance mu, LCDM
    MJD_Bmax  = SN_LCdata[5,0] # MJD at t_Bmax

    # These values are just to be printed in the text file, i.e., they are not used
    # for computing anything here.
    mu_snoopy     = SN_LCdata[2,0] # distance mu, snoopy
    err_mu_snoopy = SN_LCdata[2,1] # error distance mu, snoopy
    err_zhelio    = SN_LCdata[0,3] # helio redshift
    dm15          = SN_LCdata[1,0] #
    err_dm15      = SN_LCdata[1,1] #
    EBVhost       = SN_LCdata[4,0] #
    err_EBVhost   = SN_LCdata[4,1] #
    EBV_MW        = SN_LCdata[3,0] #
    err_EBV_MW    = SN_LCdata[3,1] #
    err_MJD_Bmax = SN_LCdata[5,1] # MJD at t_Bmax
    Alamb      = SN_LCdata[3,2]
    err_Alamb  = SN_LCdata[3,3]
    R_F        = SN_LCdata[3,4]

    # The sample flag:
    sampleFlag = name[-9:-6]
    if sampleFlag == "CfA": FlatCode = 1
    elif sampleFlag == "CSP": FlatCode = 2
    elif sampleFlag == "ers": FlatCode = 3

    #-------- Peculiar velocity uncertainty -----------------
    # Create the variable "snName" containing the first 8 (or 7) letters of the SNe file name
    # I use "snName" to compare with the SNe names in 'SNeWithCepheidDistances.txt', so that
    # I will not compute a peculiar velocity uncertainty for those SNe.
    if   name[7] == '_': snName = name[:7] # To read correctly, e.g., "sn2011B"
    elif name[7] != '_':
        if is_number(name[7]): snName = name[:15] # To read correctly, e.g., "snf20080514-002"
        else: snName = name[:8]  # To read correctly, e.g., "sn1998bu"


    sigma_pecInt = 0 # Reset its value:
    # If this SNe has Cepheid distance, then use vpecFix=0 for it:
    if snName in ListSNeCepheid['f0']: sigma_pecInt = np.sqrt(sigma2_pec(zcmb, err_zcmb, 0))
    else: sigma_pecInt = np.sqrt(sigma2_pec(zcmb, err_zcmb, vpecFix))

    #-------- At B maximum ------------------------------

    BMax_appMag    = SN_GPfit[70,1]  # B-band app mag
    ErrBMax_appMag = SN_GPfit[70,2]  # B-band app mag error

    if BandMax == 'Bmax':

        appMag_atMax    = BMax_appMag  # B-band app mag
        err_appMag_atMax = ErrBMax_appMag  # B-band app mag error

        phaseB_atMax = 0
        MJD_atMax = MJD_Bmax

        flag_kindOfMax = 1

    #-------- At NIR maximum ------------------------------

    elif BandMax == 'NIRmax':

        # Find the NIR maximum app mag of the light curve in a given rest-frame days range.
        appMag_atMax = min(SN_GPfit[MinRow:MaxRow,1])

        # Find the array index of the maximum
        IndexNIRMax_appMag = np.where(appMag_atMax == SN_GPfit[MinRow:MaxRow,1])

        # phase (relative to B band) of the NIR maximum
        phaseB_atMax    = SN_GPfit[MinRow:MaxRow,0][IndexNIRMax_appMag[0][0]]
        # NIR maximum apparent mag
        # ok, old. NIRMax_appMag2  = SN_GPfit[MinRow:MaxRow,1][IndexNIRMax_appMag[0][0]]
        # error_app mag at NIR max
        err_appMag_atMax = SN_GPfit[MinRow:MaxRow,2][IndexNIRMax_appMag[0][0]]

        # MJD of the NIR maximum
        MJD_atMax = phaseB_atMax*(1+zhelio) + MJD_Bmax

        #--- Distinguishing if the value found is a maximum or is truncated light curves ---
        # NIR apparent mag at 1/2 day before t_NIR_max
        appMag_BeforeMax = SN_GPfit[MinRow:MaxRow,1][IndexNIRMax_appMag[0][0]-1]

        # If appMag_BeforeMax > 39 then it means that the found maximum is
        # actually a truncated light-curve then flag it with the number "2",
        # otherwise flag it with "1" (meaning that it is an actual maximum).
        if appMag_BeforeMax < 39 : flag_kindOfMax = 1   # = a maximum.
        else: flag_kindOfMax = 2  # = a truncated light curve.

    #-------- Absolute magnitudes at max ------------

    # NIR Absolute magnitude at band max
    AbsMag_atMax = appMag_atMax - mu_LCDM

    # Uncertainty in NIR Absolute magnitude at band max.
    # This value is not used for computation, I determine its value
    # just to write it down in the output text file.
    err_AbsMag_atMax = np.sqrt(err_appMag_atMax**2 + sigma_pecInt**2)

    #-------- Photometric distance mu from GP fit ----------------

    # Distance modulus from the Gaussian Processes fit at maximum
    mu_GP = appMag_atMax - AverageAbsMag_atMax

    # Uncertainty in the distance modulus from the Gaussian Processes fit at maximum
    err_mu_GP = err_appMag_atMax

    # Residual
    mu_GP_residual = mu_GP - mu_LCDM

    #--- write the data BEFORE cutoffs --->
    file_all.write('%-30s  %.10f  %.10f  %.6f  %.6f    %9.6f        1.0       %.0f       %.6f      %.6f         %.6f  %.6f    %.6f      %.6f         %.6f  %.6f      %10.6f      %.9f   %.9f      %.6f  %.6f  %9.6f  %.6f    %.6f  %.6f    %.6f  %.6f  %.6f    %.6f   %.6f     %.6f  %.6f      %.6f    %.6f        # \n'%(name[:-4],
            zcmb, err_zcmb, mu_GP, err_mu_GP, mu_GP_residual, FlatCode,
            appMag_atMax, err_appMag_atMax, mu_LCDM, sigma_pecInt,
            AbsMag_atMax, err_AbsMag_atMax,
            MJD_atMax, err_MJD_Bmax, phaseB_atMax,
            zhelio, err_zhelio,
            dm15, err_dm15, EBVhost, err_EBVhost, EBV_MW, err_EBV_MW, Alamb, err_Alamb, R_F,
            mu_snoopy, err_mu_snoopy,
            MJD_Bmax, err_MJD_Bmax, BMax_appMag, ErrBMax_appMag ))
        # <--- write the data ---

    #---------- Cutoff and write the results to a data text file ------------------

    # Filter the SNe where there is not LC data at NIR-max or B-max. In the '_GP_mean_sigma_Filled.dat' text
    # file I put the value of 40 at those phases.
    # Comment the masked SNe in ListSNe_AppMag_.txt:
    if mask==0: commentText = ''
    else: commentText = '##  '

    if appMag_atMax<39 and zcmb < zcmb_Max and phaseB_atMax <= phaseB_atMax_Cutoff:

        # Comment the SNe that fail the cutoffs.
        if (dm15 >= dm15LowerLim and dm15 <= dm15UpperLim and
            EBVhost >= EBVhostMin and EBVhost <= EBVhostMax and
            EBV_MW < EBVMWLim):
            comment_text_1 = ''
            flag_cutoff = 0
        else:
            comment_text_1 = '#_ '
            flag_cutoff = 1
            countCommented_cutoffs += 1

        # Text to be written at the end of the line (in the "Notes"
        # column) about the kind of maximum:
        if flag_kindOfMax == 1 and phaseB_atMax <= phaseB_atMax_Cutoff:
            comment_text_2 = 'NIR maximum at phase = %s days.'%phaseB_atMax
        elif flag_kindOfMax == 2 and phaseB_atMax <= phaseB_truncatedLC:
            comment_text_2 = 'LC truncated at phase = %s days.'%phaseB_atMax

        #--- write the data -AFTER- cutoffs --->
        file_results.write('%s%s%-30s  %.10f  %.10f  %.6f  %.6f    %9.6f        1.0       %.0f       %.6f      %.6f         %.6f  %.6f    %.6f      %.6f         %.6f  %.6f      %10.6f      %.9f   %.9f      %.6f  %.6f  %9.6f  %.6f    %.6f  %.6f    %.6f  %.6f  %.6f    %.6f   %.6f     %.6f  %.6f      %.6f    %.6f        # %s \n'%(commentText, comment_text_1, name[:-4],
            zcmb, err_zcmb, mu_GP, err_mu_GP, mu_GP_residual, FlatCode,
            appMag_atMax, err_appMag_atMax, mu_LCDM, sigma_pecInt,
            AbsMag_atMax, err_AbsMag_atMax,
            MJD_atMax, err_MJD_Bmax, phaseB_atMax,

            zhelio, err_zhelio,
            dm15, err_dm15, EBVhost, err_EBVhost, EBV_MW, err_EBV_MW, Alamb, err_Alamb, R_F,
            mu_snoopy, err_mu_snoopy,
            MJD_Bmax, err_MJD_Bmax, BMax_appMag, ErrBMax_appMag,
            comment_text_2
            ))
        # <--- write the data ---

        # Add this good SNe to an array to later determine the SN with the
        # largest values of distance-modulus residuals, i.e., the outliers.
        mu_resid_array += [mu_GP_residual];
        SNeName_list += [snName];

        copyfile(DirGPFits+name[:-4]+'_GP_plot.png', DirSavePlots+name[:-4]+'_GP_plot.png')
        print name[:-4]

        countTotal = countTotal + 1
        if mask==0 and flag_cutoff==0:
            countNoCommented = countNoCommented + 1;
        else: countCommented = countCommented + 1;

file_results.write(text_line)
file_all.write(text_line)

text_10 = '# %s SNe in this list. \n'%countTotal
text_11 = "# %s SNe no commented automatically (## or #_). \n"%countNoCommented
text_12 = "# %s SNe were commented automatically (## or #_). \n"%countCommented
text_12_2 = "# %s SNe were commented automatically because didn't pass the \
cutoffs (#_). \n"%countCommented_cutoffs

text_13 = '# Largest upper residual: %r, SNe: %s \n'%(
    round(max(mu_resid_array),4), SNeName_list[mu_resid_array.index(max(mu_resid_array))]   )
text_14 = '# Largest lower residual: %r, SNe: %s \n'%(
    round(min(mu_resid_array),4),SNeName_list[mu_resid_array.index(min(mu_resid_array))]  )

file_results.write(text_10);
file_results.write(text_11); file_results.write(text_12); file_results.write(text_12_2);
file_results.write(text_13); file_results.write(text_14);

file_all.write(text_10);
file_all.write(text_11); file_all.write(text_12); file_all.write(text_12_2);
file_all.write(text_13); file_all.write(text_14);

print ' '
print text_10, text_11, text_12, text_13,text_14

file_results.close();
file_all.close();

file_results.close();file_results.close();file_results.close();
file_all.close();file_all.close();file_all.close();

##############################################################################80

# ## Reading the datatable I've just created above

import os # To use command line like instructions
import glob # To read the files in my directory

# Change the working directory where the data files are located
os.chdir(DirSaveOutput)

# Reading the data files in that folder
DistanceMuFiles = glob.glob('DistanceMu_Good_AfterCutoffs_Main*.txt')

#--------------------------------------------------
# Check if 'DistanceModuli_Notes_.txt' is already there, otherwise read
# the 'DistanceModuli_.txt' file.
if 'DistanceMu_Good_AfterCutoffs_Main_Notes_.txt' in DistanceMuFiles:
    DistMu_array = np.genfromtxt(DirSaveOutput+'DistanceMu_Good_AfterCutoffs_Main_Notes_.txt',
                                 dtype=['S30',
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float])
    print '# Reading the file:  < DistanceMu_Good_AfterCutoffs_Main_Notes_.txt >'

elif 'DistanceMu_Good_AfterCutoffs_Main_.txt' in DistanceMuFiles:
    DistMu_array = np.genfromtxt(DirSaveOutput+'DistanceMu_Good_AfterCutoffs_Main_.txt',
                                 dtype=['S30',
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float])
    print '# Reading the file: < DistanceMu_Good_AfterCutoffs_Main_.txt >'

elif 'DistanceMu_Good_AfterCutoffs_Main_TMP_Notes_.txt' in DistanceMuFiles:
    DistMu_array = np.genfromtxt(DirSaveOutput+'DistanceMu_Good_AfterCutoffs_Main_TMP_Notes_.txt',
                                 dtype=['S30',
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float])
    print '# Reading the file: < DistanceMu_Good_AfterCutoffs_Main_TMP_Notes_.txt >'

elif 'DistanceMu_Good_AfterCutoffs_Main_TMP_.txt' in DistanceMuFiles:
    DistMu_array = np.genfromtxt(DirSaveOutput+'DistanceMu_Good_AfterCutoffs_Main_TMP_.txt',
                                 dtype=['S30',
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float])
    print '# Reading the file: < DistanceMu_Good_AfterCutoffs_Main_TMP_.txt >'
else: print '# < DistanceMu_Good_AfterCutoffs_Main_.txt > file not found!!!'

print "# %s SNe in this file."%len(DistMu_array)

# Reading the file: < DistanceMu_Good_AfterCutoffs_Main_TMP_.txt >
# 63 SNe in this file.

#-----------------------------------------------------------------------------80

# ## Determining average Absolute magnitude of the sample

# Plotting the histogram with the Gaussian estimation

from scipy.stats import norm
import matplotlib.mlab as mlab

# best fit of data.
# 'norm.fit' simply compute the mean and standard devation of the sample.
(AverageAbsMag_atMax, err_AverageAbsMag_atMax) = norm.fit(DistMu_array['f12'])

#------ Plot -----------

# the histogram of the data
n, bins, patches = plt.hist(DistMu_array['f12'], 20, normed=True, facecolor='green', alpha=0.5)
# n, bins, patches = plt.hist(DistMu_array['f12'], 20, normed=False, facecolor='green', alpha=0.5)

# add a 'best fit' Gaussian line
y = mlab.normpdf(bins, AverageAbsMag_atMax, err_AverageAbsMag_atMax )
l = plt.plot(bins, y, 'r--', linewidth=2)

plt.xlabel(r'$\hat{m}_{\rm Bmax} - \mu_{\Lambda{\rm CDM}}$')
plt.ylabel('Probability')
plt.title(r'Histogram. (mean=%.2f, std dev = %.2f)' %(AverageAbsMag_atMax, err_AverageAbsMag_atMax))

plt.grid(True)
plt.savefig(DirSaveOutput+'Plot_histo_AbsMag_GP_%s_.png'%vpecFix)
plt.close()

# The results that I will use now.
# NOTE: This value is independent of what I have assume about "AverageAbsMag_atMax" and
# "err_AverageAbsMag_atMax" at the top of the notebook in the User section because
# I've just recomputed it in the cell above.
# NOTE: The values of "AverageAbsMag_atMax" and "err_AverageAbsMag_atMax" do NOT
# depend on the value of "vpecFix".

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d; %H:%M hrs.")
text_Date   = '# On date: %s \n'%text_timenow
print '#', '-'*30
print '#  ', Band, 'band | Band Max =', BandMax, '| 0 < z <', zcmb_Max
# print '# (mean abs mag, std deviation)\n'
# print ""

print "# AverageAbsMag_atMax = %.6f;  # %s"%(AverageAbsMag_atMax, text_timenow)
print "# err_AverageAbsMag_atMax = %.6f;"%(err_AverageAbsMag_atMax)

#--------------------------------------------------------60

# #### Weighted averages

if 2<1:
    # Inverse-variance weights array
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

    # print 'Weighted Abs Mag =', WeightedAbsMag
    # print 'Weighted Std dev =', error_WeightedAbsMag

# Compute some other values

# print '-'*30
# print '   ', Band, 'band   Pec velocity =', vpecFix, 'km/s. Band Max =', BandMax
# print np.mean(DistMu_array['f12']), ' = mean abs mag'
# print np.median(DistMu_array['f12']), '= median'
# print WeightedAbsMag, '= Weighted Abs Mag'
# print np.sqrt(1/np.sum(WeightsInvVar) ), '= uncertainty in the weighted average'
# print np.std(DistMu_array['f12']), '= population standard deviation'
# print error_WeightedAbsMag, '= unbiased weighted population standard deviation'

# 0 < z < 0.04
"""

------------------------------
    J band   Pec velocity = 150 km/s. Band Max = NIRmax
-18.566632254  = mean abs mag
-18.587673 = median
-18.586175105 = Weighted Abs Mag
0.0107267650756 = uncertainty in the weighted average
0.170036484735 = population standard deviation
0.145823909562 = unbiased weighted population standard deviation

"""
0

##############################################################################80

# # Intrinsic dispersion

# Creation of arrays for mu_residuals.

mu_resid_z0_np = np.zeros(len(DistMu_array))
sigma2_appMagTBmax_z0_np = np.zeros(len(DistMu_array))
sigma2_pec_z0_np = np.zeros(len(DistMu_array))

mu_resid_z001 = []
sigma2_appMagTBmax_z001 = []
sigma2_pec_z001 = []

PlotTotalMu = False # I silence this variable for now.

for i in range(len(DistMu_array)): # loop over 'DistanceMu_Good_AfterCutoffs_Main_.txt'

    zcmbInt     = DistMu_array[i][1]  # zcmb
    err_zcmbInt = DistMu_array[i][2]  # error_zcmb
    mu_resid_z0_np[i]   = DistMu_array[i][5]  # Delta_mu


    # Create the variable "snName" containing the first 8 (or 7) letters of the SNe file name
    # I use "snName" to compare with the SNe names in 'SNeWithCepheidDistances.txt', so that
    # I will not compute a peculiar velocity uncertainty for those SNe.
    try:
        if   DistMu_array[i][0][7] == '_':
            snName = DistMu_array[i][0][:7]  # To read correctly, e.g., "sn2011B_"
        elif DistMu_array[i][0][7] != '_':
            # To read correctly, e.g., "snf20080514-002"
            if is_number(DistMu_array[i][0][7]): snName = DistMu_array[i][0][:15]
            else: snName = DistMu_array[i][0][:8]  # To read correctly, e.g., "sn1998bu"
    except: snName = DistMu_array[i][0][:6]  # To read correctly, e.g., "sn2011B"

    if PlotTotalMu == True:
        if snName in ListSNeCepheid['f0']:
            sigma2_pec_z0_np[i] = sigma2_pec(zcmbInt, err_zcmbInt, 0)
        else: sigma2_pec_z0_np[i] = sigma2_pec(zcmbInt, err_zcmbInt, vpecFix)
    else:
        sigma2_pec_z0_np[i] = (DistMu_array[i][11])**2   # (sigma_mu_pecVel)^2

    if PlotTotalMu==False:
        # error variance of the app mag at TBmax
        sigma2_appMagTBmax_z0_np[i] = (DistMu_array[i][4])**2

    if zcmbInt > 0.01:
        mu_resid_z001 += [DistMu_array[i][5]]  # Delta_mu

        if PlotTotalMu == True:
            if snName in ListSNeCepheid['f0']:
                sigma2_pec_z001 += [sigma2_pec(zcmbInt, err_zcmbInt, 0)]
            else: sigma2_pec_z001 += [sigma2_pec(zcmbInt, err_zcmbInt, vpecFix)]
        else:
            sigma2_pec_z001 += [(DistMu_array[i][11])**2]  # (sigma_mu_pecVel)^2

        if PlotTotalMu==False:
            # error variance of the app mag at TBmax
            sigma2_appMagTBmax_z001 += [(DistMu_array[i][4])**2]

# Convert the list to np.arrays:
mu_resid_z001_np = np.array(mu_resid_z001)
sigma2_appMagTBmax_z001_np = np.array(sigma2_appMagTBmax_z001)
sigma2_pec_z001_np = np.array(sigma2_pec_z001)

# print 'Number of SNe with z>0 :', len(mu_resid_z0_np)
# print 'Number of SNe with z>0.01 :', len(mu_resid_z001_np)

# Checking the sized of the arrays: the left and right numbers have to be the same
# in a given row printed below.

print '#', len(mu_resid_z0_np),   len(sigma2_appMagTBmax_z0_np),   len(sigma2_pec_z0_np)
print '#', len(mu_resid_z001_np), len(sigma2_appMagTBmax_z001_np), len(sigma2_pec_z001_np)

# 63 63 63
# 49 49 49

# Finding the best estimated value for sigma2Pred by
# minimizing the -2ln(Likelihood) function

from scipy.optimize import fmin as simplex

Method=7 # I silence this variable for now.

InitialGuess = 0.15**2

if Method==7 and PlotTotalMu==False:
    SimplexResult_z0 = simplex(neg2lnLikeFull, InitialGuess,
                               args=(mu_resid_z0_np, sigma2_pec_z0_np, sigma2_appMagTBmax_z0_np))
    SimplexResult_z001 = simplex(neg2lnLikeFull, InitialGuess,
                                 args=(mu_resid_z001_np, sigma2_pec_z001_np, sigma2_appMagTBmax_z001_np))

else:
    SimplexResult_z0 = simplex(neg2lnLike, InitialGuess,
                               args=(mu_resid_z0_np, sigma2_pec_z0_np))
    SimplexResult_z001 = simplex(neg2lnLike, InitialGuess,
                                 args=(mu_resid_z001_np, sigma2_pec_z001_np))

print ' '
print 'Best estimated value of sigma^2_Pred (z>0) =', SimplexResult_z0
print 'Best estimated value of sigma_Pred (z>0) =', np.sqrt(SimplexResult_z0)
print ' '
print 'Best estimated value of sigma^2_Pred (z>0.01) =', SimplexResult_z001
print 'Best estimated value of sigma_Pred (z>0.01) =', np.sqrt(SimplexResult_z001)

# SPECIAL CASE:

"""
# When the intrinsic dispersion is very close to zero.

# Finding the best estimated value for sigma2Pred by
# minimizing the -2ln(Likelihood) function

# Defining the function to compute the intrinsic dispersion (sigmaPred) instead
# of the square of the intrinsic dispersion (sigma2Pred): for the case of determining
# the instrisic dispersion from total distance modulus of YJHK band only and 0.01<z<0.06,
# I obtain an error due to a very small value of sigmaPred, so that during minimizing
# the likelihood function to determine sigmaPred, it is sampled some negative values.

# For the case of using a normalized template. This is the full Eq. (B.6) of Blondin et al 2011
def neg2lnLikeFull_2(sigmaPred, mu_resid_np, sigma2_pec_np, sigma2_appmagTBmax_np):
    sum1 = 0
    for i in range(len(mu_resid_np)):
        sum1 = (sum1 + np.log(sigma2_appmagTBmax_np[i] + sigmaPred**2 + sigma2_pec_np[i]) +
                (mu_resid_np[i]**2)/(sigma2_appmagTBmax_np[i] + sigmaPred**2 + sigma2_pec_np[i]) )
    return sum1

#-------------------
from scipy.optimize import fmin as simplex

SimplexResult_z0_a = simplex(neg2lnLikeFull_2, 0.15,
                               args=(mu_resid_z0_np, sigma2_pec_z0_np, sigma2_appMagTBmax_z0_np))

SimplexResult_z001_a = simplex(neg2lnLikeFull_2, 0.15,
                               args=(mu_resid_z001_np, sigma2_pec_z001_np, sigma2_appMagTBmax_z001_np))

SimplexResult_z001 = SimplexResult_z0_a**2
SimplexResult_z001 = SimplexResult_z001_a**2

print 'Best estimated value of sigma_Pred (z>0) =', SimplexResult_z0_a
print 'Best estimated value of sigma_Pred (z>0.01) =', SimplexResult_z001_a

"""
0

# Computing the uncertainty on sigma_pred

# The variance error of sigma^2_pred:

if Method==7 and PlotTotalMu==False:
    Var_sigma2_pred_z0 =   1/FisherFuncFull(SimplexResult_z0[0],
                                            mu_resid_z0_np,   sigma2_pec_z0_np, sigma2_appMagTBmax_z0_np)
    Var_sigma2_pred_z001 = 1/FisherFuncFull(SimplexResult_z001[0],
                                            mu_resid_z001_np, sigma2_pec_z001_np, sigma2_appMagTBmax_z001_np)
else:
    Var_sigma2_pred_z0 =   1/FisherFunc(SimplexResult_z0[0],   mu_resid_z0_np,   sigma2_pec_z0_np)
    Var_sigma2_pred_z001 = 1/FisherFunc(SimplexResult_z001[0], mu_resid_z001_np, sigma2_pec_z001_np)

error_sigma_pred_z0   = np.sqrt(Var_sigma2_pred_z0  /(4*SimplexResult_z0[0]))
error_sigma_pred_z001 = np.sqrt(Var_sigma2_pred_z001/(4*SimplexResult_z001[0]))

print 'Variance of sigma^2_pred (z>0) =', round(Var_sigma2_pred_z0,6)
print 'Variance of sigma^2_pred (z>0.01) =', round(Var_sigma2_pred_z001,6)
print '\nError_sigma_pred (z>0) =', round(error_sigma_pred_z0,3)
print 'Error_sigma_pred (z>0.01) =', round(error_sigma_pred_z001,3)

# SPECIAL CASE:

"""
# Computing the uncertainty on sigma_pred

# When the intrinsic dispersion is very close to zero.

# For the case of using a normalized template. This is the full Eq. (B.7) of Blondin et al 2011
def FisherFuncFull_2(sigmaPred, mu_resid_np, sigma2_pec_np, sigma2_appmagTBmax_np):
    sum2 = 0
    for i in range(len(mu_resid_np)):
        sum2 = (sum2 + (mu_resid_np[i]**2)/(sigma2_appmagTBmax_np[i] + sigmaPred**2 + sigma2_pec_np[i])**3 -
               1/(2*(sigma2_appmagTBmax_np[i] + sigmaPred**2 + sigma2_pec_np[i])**2)  )
    return sum2

#--------------------------------------------------
if Method==7 and PlotTotalMu==False:
    Var_sigma2_pred_z0 =   1/FisherFuncFull_2(SimplexResult_z0_a[0],
                                            mu_resid_z0_np,   sigma2_pec_z0_np, sigma2_appMagTBmax_z0_np)
    Var_sigma2_pred_z001 = 1/FisherFuncFull_2(SimplexResult_z001_a[0],
                                            mu_resid_z001_np, sigma2_pec_z001_np, sigma2_appMagTBmax_z001_np)

error_sigma_pred_z0   = np.sqrt(Var_sigma2_pred_z0  /(4*SimplexResult_z0_a[0]**2))
error_sigma_pred_z001 = np.sqrt(Var_sigma2_pred_z001/(4*SimplexResult_z001_a[0]**2))

print 'Variance of sigma^2_pred (z>0) =', round(Var_sigma2_pred_z0,6)
print 'Variance of sigma^2_pred (z>0.01) =', round(Var_sigma2_pred_z001,6)
print '\nError_sigma_pred (z>0) =', round(error_sigma_pred_z0,3)
print 'Error_sigma_pred (z>0.01) =', round(error_sigma_pred_z001,3)

"""
0

# SPECIAL CASE:

"""
# Setting by hand the intrinsic dispersion to zero when it is negative and very close to zero

SimplexResult_z0 = [0]

SimplexResult_z0[0] = 0
error_sigma_pred_z0 = 0

"""
0

# #### Report the intrinsic dispersion

print '#'+'-'*50
# print "# %s band, vpecFix = %s km/s"%(Band, vpecFix)
# print "# Intrinsic dispersion for the case (0 < z < %s):"%(zcmb_Max)
print '# Intrinsic dispersion = %s +/- %s'%(np.sqrt(SimplexResult_z0[0]),
                                            error_sigma_pred_z0)

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d; %H:%M hrs.")
print ' '
print '# Date: %s '%(text_timenow)
print "# IntrinsicDisp = %.5f # %s | vpecFix=%s km/s | 0<z<%s."%(
    np.sqrt(SimplexResult_z0[0]), Band, vpecFix, zcmb_Max)

#--------------------------------------------------
# Intrinsic dispersion = 0.103380365641 +/- 0.0181455651003

# Date: 2018-02-22; 17:53 hrs.
# IntrinsicDisp = 0.10338 # J | vpecFix=150 km/s | 0<z<0.04.

#--------------------------------------------------
# Intrinsic dispersion = 0.170408516513 +/- 0.0363568720916

# Date: 2018-02-22; 22:10 hrs.
# IntrinsicDisp = 0.17041 # Y | vpecFix=150 km/s | 0<z<0.04.

##############################################################################80

# ## Check the consistency between the error bars of the residual distance modulus vs the scatter in the Hubble-diagram residual plot

ratio_int = 0;

for i in range(len(DistMu_array)):
    # print i

    mu_resid     = DistMu_array[i][5]
    error_appMag = DistMu_array[i][9]
    sigma_pecVel = DistMu_array[i][11]

    ratio_int = ratio_int + ((mu_resid**2) / (error_appMag**2 +
                            sigma_pecVel**2 + SimplexResult_z0[0]) )

chi2_dof_HD    = ratio_int / len(DistMu_array)

print '#'+'-'*40
print "# %s band | vpecFix = %s km/s | 0 < z < %s"%(Band, vpecFix, zcmb_Max)
print '# chi^2 =', ratio_int, '| Number of data =', len(DistMu_array)
print '# chi^2_dof =', chi2_dof_HD

#--------------------------------------------------
# J band | vpecFix = 150 km/s | 0 < z < 0.04
# chi^2 = 57.9577533241 | Number of data = 63
# chi^2_dof = 0.919964338478

#--------------------------------------------------------60

# # Write the summary of values to a text file

textfile_1 = open(DirSaveOutput+"Summary_HDScatter_RMS_.txt", 'w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd).")
text_Author = '# Data table created by: Arturo Avelino \n'
text_Date   = '# On date: %s \n'%text_timenow
text_script = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
text_line = '#'+'-'*50 + '\n'

textfile_1.write("# Summary of the scatter in the Hubble residual \n")
textfile_1.write("# %s band \n"%Band)

textfile_1.write(text_line)
textfile_1.write(text_Author); textfile_1.write(text_Date);
textfile_1.write(text_script);
textfile_1.write(text_line)

textfile_1.write('# Apparent magnitudes and distance moduli computed from the GP fit directly. \n')
textfile_1.write('# Peculiar velocity uncertainty assumed: %s km/s  \n'%vpecFix)
textfile_1.write('# Band used as reference of maximum: %s \n'%BandMax)
textfile_1.write('# Phase_Bmax range to look for the maximum: Minimum = %s days (%s rows), \
Maximum = %s days (%s rows) \n'%(MinPhase, MinRow, MaxPhase, MaxRow))
textfile_1.write('# Discard if the maximum is located at phase > \
phaseB_atMax_Cutoff = %s days. \n'%phaseB_atMax_Cutoff)
textfile_1.write('# Ok if the light curve is truncated (instead of a maximum) at phaseB =< \
phaseB_truncatedLC = %s days. \n'%phaseB_truncatedLC)
textfile_1.write('# phaseB_truncatedLC = %s days \n'%phaseB_truncatedLC)
textfile_1.write('# Upper limit zcmb cutoff = %s \n'%zcmb_Max)

textfile_1.write(text_line)
textfile_1.write('# Intrinsic dispersion for the case (0 < z < %s) and used to obtain the \n\
# total photometric distance modulus uncertainty: \n'%zcmb_Max)

textfile_1.write('%14.6f  %.6f  0.0 # Intrinsic dispersion and its uncertainty for the case \
0 < z_cmb < %s \n'%(np.sqrt(SimplexResult_z0[0]), error_sigma_pred_z0, zcmb_Max))

# Sometimes I get error when computing SimplexResult_z001, these lines help
# to alleviate the writting of the file.
# ---->>
try:
    intDisp001 = np.sqrt(SimplexResult_z001[0])
    err_intDisp001 = error_sigma_pred_z001
except:
    print "Hola"

if intDisp001>0.:
    textfile_1.write('%14.6f  %.6f  0.0  0.0 # Intrinsic dispersion and its uncertainty for the case \
# 0.01 < z_cmb < %s \n'%(intDisp001, err_intDisp001, zcmb_Max))
else:
    textfile_1.write('-1         -1  0.0  0.0 # Intrinsic dispersion and its uncertainty for the case \
# 0.01 < z_cmb < %s \n')
# <<-----

textfile_1.write('%14.6f  %.6f  0.0  0.0 # (AverageAbsMag_atMax, err_AverageAbsMag_atMax) \n'%(
       AverageAbsMag_atMax, err_AverageAbsMag_atMax))

textfile_1.write('%14.6f  %-8.0f  0.0  0.0 # (chi^2, Number of SNe) for the case (0 < z < %s) \
\n'%(ratio_int, len(DistMu_array), zcmb_Max))
textfile_1.write('%14.6f  0         0.0  0.0  # chi^2_dof \n'%chi2_dof_HD)

textfile_1.write(text_line)
textfile_1.close()

#--------------------------------------------------------60

# ## Plot the app-mag Hubble diagram from (zcmb, $m_{\rm NIR max}$)

fontsizePlot = 16

nbins1= 51
z1 = np.linspace(0.01, 0.1, nbins1)
mu1 = DistanceMuVector(z1, OmMFix, wFix, HoFix)

#-------------------
# Plot

fig = plt.figure()

# Plotting the data
plt.errorbar(DistMu_array['f1'], DistMu_array['f8'], yerr=DistMu_array['f9'],
             fmt='.', color=[1,0,0],  ms=10, elinewidth=1.5, capthick=1.5)

# Plotting the theory
# plt.plot(z1, mu1, color='black')

# plt.xlim(xlimPlots)
plt.xlim(0,0.06)

# Labeling
plt.xlabel('Redshift', fontsize=fontsizePlot)
plt.ylabel('Apparent magnitude', fontsize=fontsizePlot)
plt.title('Hubble diagram, %s band'%Band)

plt.grid(True)
plt.tight_layout()

plt.savefig(DirSaveOutput+'Plot_HubbleDiagram_appMag_%s_.png'%vpecFix, format='png')
plt.close()

#-------------------
print "# All done."

