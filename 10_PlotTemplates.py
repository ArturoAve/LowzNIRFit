#!/usr/bin/env python
# coding: utf-8

# # Plotting NIR template
#

# # User inputs
#
# ### Here is the only part where I have to setup the things for a specific band
#
# The rest of the notebook can be run without any additional manipulation.

import numpy as np
from matplotlib import pyplot as plt
#---- Setting directories -----
import glob # To read the files in my directory
import os # To use command line like instructions

#--------------------------------------------------------60
code_created_by = 'Arturo_Avelino'
# On date: 2017.01.10 (yyyy.mm.dd)
code_name = '10_PlotTemplates.ipynb'
version_code = '0.7.2'
last_update = '2019.01.28'
#--------------------------------------------------------60

##############################################################################80

# # User

#--- Kind of data to be fitted ---

Band = 'K'   # (Y, J, H, K)

# NOTE: To plot the normalized template I have to GP fit the light curves
# assuming a peculiar velocity uncertainty of vpec = 0, then to plot this
# case select "KindOfData = 'AllSamples_vpec_0'" and "Normalized = True".

# KindOfData = 'CfA'
# KindOfData = 'CSP'
# KindOfData = 'Others'
# KindOfData = 'AllSamples_vpec_0'
KindOfData = 'AllSamples_vpec_150'

# Plot the normalized template?:
Normalized = True

# Just for plotting purposes to plot the error bars of
# the photometric data.
# When "Normalized = False" then the value of "velPecuFix" has to be the
# same than the one used in the GP fit of the ABS-mag LCs.
velPecuFix = 150  # 0, 150 , 300 # km/sec

# Error in the estimation of redshift to plot the error bars of the photometric data.
zerrorFix = 0.001

# Minimal redshift
zcmb_min = 0.0
# zcmb_min = 0.01

#-------------------
# Indicate the technique used to compute the GP fitting in the '1_AllData_InitialFit' step.
# This information is just to label the '3_Template' folder where I will save the output.
# TempPrior_Hyper = Using a template prior (computed from the Moving Windows
# Average template) for the Gaussian Process fitting, then determining the
# hyperparameters using all the LCs simultaneously.
# FlatPrior= Assuming a flat prior at ~ -17 Abs mag, then computing the hyperparameters
# for each LC independently.
# For the low-z paper I use the option 'FlatPrior'.

Technique = 'FlatPrior'    # ('FlatPrior', 'TempPrior_Hyper')

#-------------------
#-- If there are -already existing- plots overlapping the GP fit and template
# in '2_Selection/[Sample]/Good/', then copy those plots to "/3_Template/[Sample]/etc".
# Check in '2_Selection/[Sample]/Good/' the extension used and write down it below.
TemplateOverplotted = 'TempAllz001'
# TemplateOverplotted = 'Temp'
# TemplateOverplotted = 'TempAllz0'
# TemplateOverplotted = 'TempCSPz001'
# TemplateOverplotted = 'TempCSPz0'

#=========================================

# (FIX) For the confidence-interval plots
CL_size = 1 # For the confidence-interval plots (1.=68.3% CL or 1.9600=95% CL )
CL_label = 68.3

LegendFontSize = 11
PlotResolution = 110 # dpi

#--- Filter system
FilterSyst = 'Std_filters/'
# FilterSyst = 'CSP_filters/'

##############################################################################80

# # AUTOMATIC
#
# # Loading data

#-------------------
#    J band:

if Band == 'J':
    # x_IntervalPlots = -10, 60 #
    x_IntervalPlots = -10, 45 #
    yIntervalPlotsUnnorma = -14, -20
    yIntervalPlotsNorma = 4, -1
    y_Interval_Residual = -0.8, 0.8

#-------------------
#    Y band:

if Band == 'Y':
    # x_IntervalPlots = -10, 60 #
    x_IntervalPlots = -10, 43 #
    yIntervalPlotsUnnorma = -15.5, -19.5
    yIntervalPlotsNorma = 4, -1
    y_Interval_Residual = -0.8, 0.8

#-------------------
# H band:

if Band == 'H':
    # x_IntervalPlots = -10, 60
    x_IntervalPlots = -10, 40 #
    yIntervalPlotsUnnorma = -15.5, -20
    yIntervalPlotsNorma = 4, -1
    y_Interval_Residual = -0.8, 0.8

#-------------------
# K band:

if Band == 'K':
    # x_IntervalPlots = -10, 55
    x_IntervalPlots = -7.5, 40   #
    yIntervalPlotsUnnorma = -17, -20.5
    yIntervalPlotsNorma = 2, -1
    y_Interval_Residual = -0.8, 0.8

#--------------------------------------------------------60

# Definitions based on zcmb_min

if zcmb_min==0.0:
    zCMBFolder = 'z_gr_0/'
elif zcmb_min==0.01:
    zCMBFolder = 'z_gr_001/'
else: zCMBFolder = 'z_gr_/'

# The text for the zcmb_min
if zcmb_min==0.0:
    zCMBText = '(z > 0)'
elif zcmb_min==0.01:
    zCMBText = '(z > 0.01)'
else: zCMBText = ' '

#--------------------------------------------------
TitleForPlots = Band+' band template' + zCMBText

#--------------------------------------------------------60

import os # To use command line like instructions

MainDirectory =  '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/TheTemplates/'+Band+'_band/'

if Normalized == True:
    NormaFolder = 'Normalized'
elif Normalized == False:
    NormaFolder = 'Unnormalized'

DirTempFiles = (MainDirectory+FilterSyst+'3_Template_'+Technique+'/'+KindOfData+
                '/'+ zCMBFolder)
DirSaveData  = (MainDirectory+FilterSyst+'3_Template_'+Technique+'/'+KindOfData+
                '/'+ zCMBFolder+NormaFolder+'/')

if not os.path.exists(DirSaveData): os.makedirs(DirSaveData)

print 'Directory to save the outputs: \n'
print DirSaveData

# Some useful definitions based on if I'm plotting the normalized or unnormalized template
if Normalized == True:
    NormaSuffix = '_Norma'
    yIntervalPlots = yIntervalPlotsNorma
elif Normalized == False:
    NormaSuffix = ''
    yIntervalPlots = yIntervalPlotsUnnorma

#--------------------------------------------------------60

# ### Loading the list of SNe to be plotted

import glob # To read the files in my directory

# Change the working directory where the data files are located
os.chdir(DirTempFiles)
# Reading the data files in that folder
ListFiles = glob.glob('SN_list_template_*.txt')

# Check if 'SN_list_template_Notes_.txt' is already there, otherwise read
# the 'SN_list_template_.txt' file.
if 'SN_list_template%s_Notes_.txt'%(NormaSuffix) in ListFiles:
    list_SNe = np.genfromtxt(DirTempFiles+'SN_list_template%s_Notes_.txt'%(NormaSuffix), dtype='S100')
    print 'Reading the file:  < SN_list_template%s_Notes_.txt >'%(NormaSuffix)
elif 'SN_list_template%s_.txt'%(NormaSuffix) in ListFiles:
    list_SNe = np.genfromtxt(DirTempFiles+'SN_list_template%s_.txt'%(NormaSuffix), dtype='S100')
    print 'Reading the file: < SN_list_template%s_.txt >'%(NormaSuffix)
else: print '< SN_list_template%s_.txt > file not found'%(NormaSuffix)

print 'Number of SNe in the list:', len(list_SNe)

list_SNe[1]

# ### Template: Confidence interval

Template_GP_Hier = np.genfromtxt(DirTempFiles+'Template_phase_mu_stdError_FromR%s.dat'%(NormaSuffix))

x = np.atleast_2d(Template_GP_Hier[:,0]).T
y_pred = Template_GP_Hier[:,1]
sigma = Template_GP_Hier[:,2]

# ### Template: Standard deviation

Template_GP_Hier_2 = np.genfromtxt(DirTempFiles+'Template_phase_mu_tau_FromR%s.dat'%(NormaSuffix))

x_2 = np.atleast_2d(Template_GP_Hier_2[:,0]).T
y_pred_2 = Template_GP_Hier_2[:,1]
sigma_2 = Template_GP_Hier_2[:,2]

# ### Individual Gaussian Process mean an variance for each SN

# ### Loading the light-curve data

cc = 299792.458  # Speed of light

#---- Setting directories -----
import glob # To read the files in my directory
import os # To use command line like instructions

#-------------------
#- Setting the where the LC data is located:
DirLCs = MainDirectory+FilterSyst+'2_Selection_'+Technique+'/'+KindOfData+'/Goods/'
os.chdir(DirLCs)

#- Reading the LC data file names
# list_SNe = glob.glob(ExtensionDataFiles) # OK. OLD

#-------------------
#- "If the subdirectory does not exist then create it"
# if not os.path.exists(DirSaveData): os.makedirs(DirSaveData)

print 'Directory to read the individual LC data: \n'
print DirLCs
print '\n Number of LCs to be plotted:', len(list_SNe)

# Number of LCs to be plotted: 50

# Peculiar velocity uncertainty: $\sigma_{\mu}(v_{pec}, z, \sigma_z)$

# Peculiar velocity uncertainty

# Wood-Vasey+08 definition
def sigma_muPecu(vpec, z, zerror):
    "Sigma_mu from the peculiar velocity"
    sigma_mu = (5/(z*np.log(10))) * np.sqrt((vpec/cc)**2 + zerror**2)
    return sigma_mu

sigma_muPecu(150, 0.01, zerrorFix)
# 0.2428116195960309

#--------------------------------------------------------60

# #### Combininig all data in a single array

# Combininig all the data in a single array

# This lines have information about the SN
skipFirstLines = 10

# Number of characters to remove the suffix "_GP_mean_sigma_Filled_norma.dat"
# or "_GP_mean_sigma_Filled.dat"
ncharsuffix = 0 # reset
if Normalized: ncharsuffix = 31
else: ncharsuffix = 25

#-------------------
#  Defining the first element to append the rest of the data

snfilename = '' # reset
snfilename = '%s.txt'%list_SNe[0][0:(len(list_SNe[0]) - ncharsuffix)]

# First_list_SNe = list_SNe+'.txt'
DataLC = 0
DataLC = np.genfromtxt(snfilename, skip_header=skipFirstLines)

if Normalized:
    AbsMag_TBmax = np.genfromtxt(list_SNe[0])[70][1]
else: AbsMag_TBmax = 0

print 'Test. First SN:', list_SNe[0]
print 'Test. Number of data in the first LC:', len(DataLC)

# print DataLC

zz = np.genfromtxt(list_SNe[0])[0,0]

#- Summing the errors in quadratures:
error_M_PecVel = 0
error_M_PecVel = np.sqrt(DataLC[:,2]**2 + sigma_muPecu(velPecuFix, zz, zerrorFix)**2)
print 'Test. Number of data errors in the first LC:', len(error_M_PecVel)

#-------------------
# Appending the rest of the SNe

DataLC_temp = 0
for i in range(1, len(list_SNe)):

    snfilename = '' # reset
    snfilename = '%s.txt'%list_SNe[i][0:(len(list_SNe[i]) - ncharsuffix)]

    DataLC_temp = np.genfromtxt(snfilename, skip_header=skipFirstLines)
    zz_temp = np.genfromtxt(snfilename)[0,0]
    error_M_PecVel_temp = np.sqrt(DataLC_temp[:,2]**2 + sigma_muPecu(velPecuFix, zz_temp, zerrorFix)**2)

    DataLC = np.vstack((DataLC, DataLC_temp))
    error_M_PecVel = np.hstack((error_M_PecVel, error_M_PecVel_temp))

print 'Length of final light-curve array: (', len(DataLC), ',', len(error_M_PecVel), ',', len(DataLC.T), ')'

# Test. First SN: sn1998bu__U_69_B_9_J_band
# Test. Number of data in the first LC: 18
# Test. Number of data errors in the first LC: 18
# Length of final light-curve array: ( 1174 , 1174 , 3 )

# ### OPTIONAL: Loading the moving window average template

# Loading the data of the weighted average template
"""
Template_WeightMeanSmooth = np.genfromtxt(MainDirectory+'4_WoodVaseyProcedure/Plots/'+'Template_WeightedMean_StdErrorMean_SmoothW_6_Box_3_Step_05.dat')

len(Template_WeightMeanSmooth), len(Template_WeightMeanSmooth.T)
# (125, 3)
"""
0

##############################################################################80

# # PLOTTING
#
# # Overlaying all the individual GP plots and the template

# ### Plotting all the GP functions together computed for each SN

# Plotting all the GP functions together computed for each SN

ii =0
for file in list_SNe:
    # print file
    ii =ii+1

    snname = '' # reset
    snname = '%s'%file[0:(len(file) - ncharsuffix)]

    DataGP=np.genfromtxt(snname+'_GP_mean_sigma.dat')

    if Normalized:
        AbsMag_TBmax = np.genfromtxt(snname+
                        '_GP_mean_sigma_Filled.dat')[70][1]
    else: AbsMag_TBmax = 0

    plt.plot(DataGP[:,0], DataGP[:,1]-AbsMag_TBmax, label=''.format(i=ii),
            linewidth = 0.6)

plt.xlim(x_IntervalPlots)

plt.ylim(yIntervalPlots)

plt.title(TitleForPlots)

#-------------------
plt.grid(True)

plt.xlabel(r'Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
plt.ylabel('Absolute magnitude')

plt.savefig(DirSaveData+'Plot_GP_All_Individual_.png', dpi=(PlotResolution+20))
# plt.show()
plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ### Plotting all the GP functions together computed for each SN + data

# Plotting all the GP functions together computed for each SN

if Normalized == False:
    fig = plt.figure()

    ii =0
    for file in list_SNe:
        # print file
        ii =ii+1

        snname = '' # reset
        snname = '%s'%file[0:(len(file) - ncharsuffix)]

        DataGP=np.genfromtxt(snname+'_GP_mean_sigma.dat')

        if Normalized: AbsMag_TBmax = np.genfromtxt(snname+'_GP_mean_sigma_Filled.dat')[70][1]
        else: AbsMag_TBmax = 0

        plt.plot(DataGP[:,0], DataGP[:,1]-AbsMag_TBmax, label=''.format(i=ii))

    #---------------------
    # THE LIGHT-CURVE DATA. Including peculiar velocity

    plt.errorbar(DataLC[:,0], DataLC[:,1], error_M_PecVel, fmt='grey', ls='', markersize=5, alpha=0.3)

    #---------------------

    plt.xlim(x_IntervalPlots)

    plt.ylim(yIntervalPlots)
    plt.title(TitleForPlots)

    plt.grid(True)

    plt.xlabel(r'Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
    plt.ylabel('Absolute magnitude')

    plt.savefig(DirSaveData+'Plot_GP_All_Individual_Data_.png', dpi=PlotResolution)
    # plt.show()
    plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ## Plot the abs mag data only

# Plot the function, the prediction and the confidence interval based on
# the MSE

if Normalized == False:

    fig = plt.figure()

    #-----------------------------------------------
    # THE LIGHT-CURVE DATA. Including peculiar velocity

    plt.errorbar(DataLC[:,0], DataLC[:,1], error_M_PecVel, fmt='r.', markersize=5, alpha=0.3)

    #---------------------

    plt.xlim(x_IntervalPlots)
    plt.ylim(yIntervalPlots)
    plt.title(TitleForPlots)

    #---------------------

    plt.grid(True)

    # plt.legend(loc='lower left', fontsize=LegendFontSize)
    plt.xlabel('Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
    plt.ylabel('Absolute magnitude')

    # plt.savefig(DirSaveData+'Plot_GP_Template_.png')
    plt.savefig(DirSaveData+'Plot_Data_%s_.png'%Band, dpi=PlotResolution)
    # plt.show()
    plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ## Plot the abs mag data only, with colors

# Plotting all the GP functions together computed for each SN

if Normalized == False:
    skipFirstLines = 10 # This lines have information about the SN

    fig = plt.figure()

    ii = 0
    DataLC_temp = 0
    for i in range(1, len(list_SNe)):
        ii =ii+1

        DataLCInt = np.genfromtxt(list_SNe[i], skip_header=skipFirstLines)
        errorM_PecVel_Int = np.sqrt(DataLCInt[:,2]**2 + sigma_muPecu(velPecuFix, zz_temp, zerrorFix)**2)

        plt.errorbar(DataLCInt[:,0], DataLCInt[:,1], errorM_PecVel_Int,
                     ls='None', markersize=5, alpha=0.7, label=''.format(i=ii))

    plt.xlim(x_IntervalPlots)

    plt.ylim(yIntervalPlots)
    plt.title(TitleForPlots)

    #---------------------
    plt.grid(True)

    plt.xlabel(r'Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
    plt.ylabel('Absolute magnitude')

    plt.savefig(DirSaveData+'Plot_Data_Color_%s_.png'%Band, dpi=PlotResolution)
    # plt.show()
    plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ## Plot the abs mag data and GP strings, with colors

# Plotting all the GP functions together computed for each SN

if Normalized == False:
    skipFirstLines = 10 # This lines have information about the SN

    fig = plt.figure()

    ii = 0
    DataLC_temp = 0
    for i in range(1, len(list_SNe)):
        ii =ii+1

        labelInt = ''.format(i=ii)

        # The GP string
        DataGP=np.genfromtxt(list_SNe[i][:-4]+'_GP_mean_sigma.dat')
        # plt.plot(DataGP[:,0], DataGP[:,1], label=''.format(i=ii))
        plt.plot(DataGP[:,0], DataGP[:,1], label=labelInt)

        # The data
        DataLCInt = np.genfromtxt(list_SNe[i], skip_header=skipFirstLines)
        errorM_PecVel_Int = np.sqrt(DataLCInt[:,2]**2 + sigma_muPecu(velPecuFix, zz_temp, zerrorFix)**2)
        plt.errorbar(DataLCInt[:,0], DataLCInt[:,1], errorM_PecVel_Int,
                     ls='None', markersize=5, alpha=0.7,
                     # label=''.format(i=ii))
                     label=labelInt)

    plt.xlim(x_IntervalPlots)

    plt.ylim(yIntervalPlots)
    plt.title(TitleForPlots)

    #---------------------
    plt.grid(True)

    plt.xlabel(r'Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
    plt.ylabel('Absolute magnitude')

    plt.savefig(DirSaveData+'Plot_Data_Color2_%s_.png'%Band, dpi=PlotResolution)
    # plt.show()
    plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ### Plotting all the GP function VARIANCES together computed for each SN

# Plotting all the GP functions together computed for each SN

# Defining the alpha transparency based on if I use all the data
# or just a subsample.
if KindOfData == '2_AllSubsamples/':
    alphaTransparency = 0.04
else:
    alphaTransparency = 0.07

for file in list_SNe:
    # print file

    snname = '' # reset
    snname = '%s'%file[0:(len(file) - ncharsuffix)]

    DataGP=np.genfromtxt(snname+'_GP_mean_sigma.dat')

    if Normalized: AbsMag_TBmax = np.genfromtxt(snname+
                                    '_GP_mean_sigma_Filled.dat')[70][1]
    else: AbsMag_TBmax = 0

    xx = np.atleast_2d(DataGP[:,0]).T
    yy_pred = DataGP[:,1]-AbsMag_TBmax
    sigmas = DataGP[:,2]

    plt.fill(np.concatenate([xx, xx[::-1]]),
        np.concatenate([yy_pred - CL_size * sigmas,
                      (yy_pred + CL_size * sigmas)[::-1]]),
        # np.concatenate([y_pred - 1 * sigma,
        #                (y_pred + 1 * sigma)[::-1]]),
        alpha=alphaTransparency, fc='b', ec='None')

plt.xlim(x_IntervalPlots)

plt.ylim(yIntervalPlots)
plt.title(TitleForPlots+' (GP std dev)')

#-------------------
plt.grid(True)

plt.xlabel(r'Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
plt.ylabel('Absolute magnitude')

plt.savefig(DirSaveData+'Plot_GP_All_Individual_Var_.png', dpi=PlotResolution)
# plt.show()
plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ### TO EDIT
# ### $\Delta_{\rm m15}$. Plotting all the GP functions together computed for each SN with a gradient of color for $\Delta_{\rm m15}$ parameter

# Plotting all the GP functions together computed for each SN

if 2 < 1:

    for file in list_SNe:
        # print file

        snname = '' # reset
        snname = '%s'%file[0:(len(file) - ncharsuffix)]

        DataGP=np.genfromtxt(snname+'_GP_mean_sigma.dat')

        if Normalized: AbsMag_TBmax = np.genfromtxt(snname+
                                '_GP_mean_sigma_Filled.dat')[70][1]
        else: AbsMag_TBmax = 0

        DataLC_Int = np.genfromtxt(file)
        dm15_value = DataLC_Int[1,0]
        plt.plot(DataGP[:,0], DataGP[:,1]-AbsMag_TBmax, lw=2
                 # , c=(1-dm15_value/2, 1-dm15_value/2, 1-dm15_value/2)
                 , c=(1-dm15_value/2, 0, 1-dm15_value/2)
                 # , c='blue'
                 # , alpha=1.9/dm15_value
                 , alpha=dm15_value/2
                )

    # Creating the legend of colors by hand
    dm15Artificial = np.arange(0.8, 1.7, 0.2)
    for i in dm15Artificial:
        plt.plot((x_IntervalPlots[1]-13,x_IntervalPlots[1]-7),
                 (yIntervalPlots[1]+i, yIntervalPlots[1]+i),
                 lw=3, c=(1-i/2, 0, 1-i/2), alpha=i/2  )
        plt.text(x_IntervalPlots[1]-6, (yIntervalPlots[1]+i+0.1), '%r'%round(i,2))

    plt.text(x_IntervalPlots[1]-10, (yIntervalPlots[1]+0.4), '$\Delta m15 $', fontsize=14)

    plt.xlim(x_IntervalPlots)

    plt.ylim(yIntervalPlots)
    plt.title(TitleForPlots+r' ($\Delta {\rm m15}$)')

    # plt.text(40,-19.4, '----- 1', color=(0, 0., 1))
    # plt.text(40,-19.2, '----- 0.8', color=(0, 0., 0.8))

    #---------------------
    plt.grid(True)

    plt.xlabel(r'Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
    plt.ylabel('Absolute magnitude')

    # plt.savefig(DirSaveData+'Plot_GP_All_Individual_dm15_.png', dpi=PlotResolution)
    plt.savefig(DirSaveData+'Plot_GP_All_Individual_dm15_.png')

    # plt.show()
    plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ### $z_{\rm CMB}$. Plotting all the GP functions together computed for each SN with a gradient of color for $z_{\rm CMB}$ parameter

# Plotting all the GP functions together computed for each SN

factorInt = 12

for file in list_SNe:
    # print file

    snname = '' # reset
    snname = '%s'%file[0:(len(file) - ncharsuffix)]

    DataGP=np.genfromtxt(snname+'_GP_mean_sigma.dat')

    if Normalized: AbsMag_TBmax = np.genfromtxt(snname+
                                '_GP_mean_sigma_Filled.dat')[70][1]
    else: AbsMag_TBmax = 0

    DataLC_Int = np.genfromtxt(file)
    zcmb_value = DataLC_Int[0,0]
    plt.plot(DataGP[:,0], DataGP[:,1]-AbsMag_TBmax, lw=2
             # , c=(0, 0, 1-zcmb_value*factorInt)
             , c='blue'
             , alpha=zcmb_value*factorInt
            )

# Creating the legend of colors by hand
zCMBArtificial = np.arange(0.002, 0.08, 0.01)
for i in zCMBArtificial:
    plt.plot((x_IntervalPlots[1]-13,x_IntervalPlots[1]-7),
             (yIntervalPlots[1]+i*(factorInt+6)+0.5, yIntervalPlots[1]+i*(factorInt+6)+0.5),
             lw=2.5, color='blue', alpha = i*factorInt )
    plt.text(x_IntervalPlots[1]-6, (yIntervalPlots[1]+i*(factorInt+6)+0.58), '%r'%round(i,3))


plt.text(x_IntervalPlots[1]-10, (yIntervalPlots[1]+0.36), 'z_CMB')

plt.xlim(x_IntervalPlots)

plt.ylim(yIntervalPlots)
plt.title(TitleForPlots+r' ($z_{\rm CMB}$)')

# plt.text(40,-19.4, '----- 1', color=(0, 0., 1))
# plt.text(40,-19.2, '----- 0.8', color=(0, 0., 0.8))

#-------------------
plt.grid(True)

plt.xlabel(r'Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
plt.ylabel('Absolute magnitude')

plt.savefig(DirSaveData+'Plot_GP_All_Individual_zCMB_.png', dpi=PlotResolution)
# plt.show()
plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ### TO EDIT
# ### Plotting all the GP functions together computed for each SN and overlaying both snakes.

if 2< 1:

    # Plotting all the GP functions together computed for each SN

    ii =0
    for file in list_SNe:
        # print file
        ii =ii +1
        DataGP=np.genfromtxt(file[:-4]+'_GP_mean_sigma.dat')

        if Normalized: AbsMag_TBmax = np.genfromtxt(file[:-4]+'_GP_mean_sigma_Filled.dat')[70][1]
        else: AbsMag_TBmax = 0

        plt.plot(DataGP[:,0], DataGP[:,1]-AbsMag_TBmax, label=''.format(i=ii))

    #-----------------------------------------------
    #      PLOT GP TEMPLATE

    # Standard deviation

    plt.fill(np.concatenate([x_2, x_2[::-1]]),
            # np.concatenate([y_pred_2 - 1.9600 * sigma,
            #               (y_pred_2 + 1.9600 * sigma)[::-1]]),
            np.concatenate([y_pred_2 - 1. * sigma_2,
                           (y_pred_2 + 1. * sigma_2)[::-1]]),
            alpha=0.5, fc='g', ec='None', label='Standard deviation')
    #-----------------------------------------------

    # Confidence interval
    plt.fill(np.concatenate([x, x[::-1]]),
            np.concatenate([y_pred - CL_size * sigma,
                          (y_pred + CL_size * sigma)[::-1]]),
            # np.concatenate([y_pred - 1 * sigma,
            #                (y_pred + 1 * sigma)[::-1]]),
            alpha=0.7, fc='b', ec='None', label='{0}% confidence interval GP'.format(CL_label))

    # PLOT GP TEMPLATE mean
    plt.plot(x, y_pred, 'k-', lw=2, label=u'Gaussian Process mean (GP)')

    #-----------------------------------------------

    plt.xlim(x_IntervalPlots)

    if Normalized: plt.ylim(yIntervalPlotsNorma)
    else: plt.ylim(yIntervalPlots)

    plt.title(TitleForPlots)

    #---------------------

    plt.grid(True)

    plt.xlabel(r'Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
    plt.ylabel('Absolute magnitude')
    plt.legend(loc='lower left', fontsize=LegendFontSize)

    plt.savefig(DirSaveData+'Plot_GP_All_Template_CL_.png', dpi=PlotResolution)
    # plt.show()
    plt.close()

plt.close();plt.close();plt.close();

##############################################################################80

# # Plotting the template computed with GP and hierarchical

# ### Plotting the Confidence Interval  and the template. OK

# Plot the function, the prediction and the confidence interval based on
# the MSE

fig = plt.figure()
# plt.errorbar(X.ravel(), y, dy, fmt='r.', markersize=10, label=u'Observations')
plt.plot(x, y_pred, 'k-', lw=2, label=u'Gaussian Process mean (GP)')
plt.fill(np.concatenate([x, x[::-1]]),
        # np.concatenate([y_pred - 1.9600 * sigma,
        #               (y_pred + 1.9600 * sigma)[::-1]]),
        np.concatenate([y_pred - CL_size * sigma,
                       (y_pred + CL_size * sigma)[::-1]]),
        alpha=.5, fc='b', ec='None', label='{0}% confidence interval GP'.format(CL_label))

plt.xlim(x_IntervalPlots)
plt.ylim(yIntervalPlots)
plt.title(TitleForPlots)

#-------------------
plt.grid(True)

plt.legend(loc='lower left', fontsize=LegendFontSize)
plt.xlabel('Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
plt.ylabel('Absolute magnitude')

plt.savefig(DirSaveData+'Plot_GP_Template_CL_.png', dpi=PlotResolution)
# plt.show()
plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ### Plotting the Standard Deviation and the template. OK

# Plot the function, the prediction and the confidence interval based on
# the MSE

fig = plt.figure()
# plt.errorbar(X.ravel(), y, dy, fmt='r.', markersize=10, label=u'Observations')
plt.plot(x_2, y_pred_2, 'k-', lw=2, label=u'Gaussian Process mean (GP)')
plt.fill(np.concatenate([x_2, x_2[::-1]]),
        # np.concatenate([y_pred_2 - 1.9600 * sigma,
        #               (y_pred_2 + 1.9600 * sigma)[::-1]]),
        np.concatenate([y_pred_2 - 1. * sigma_2,
                       (y_pred_2 + 1. * sigma_2)[::-1]]),
        alpha=.5, fc='g', ec='None', label='Standard deviation')

plt.xlim(x_IntervalPlots)
plt.ylim(yIntervalPlots)
plt.title(TitleForPlots)

#-------------------
plt.grid(True)

plt.legend(loc='lower left', fontsize=LegendFontSize)
plt.xlabel('Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
plt.ylabel('Absolute magnitude')

# plt.savefig(DirSaveData+'Plot_GP_Template_.png')
plt.savefig(DirSaveData+'Plot_GP_Template_StdDev_.png', dpi=PlotResolution)
# plt.show()
plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ### Overlay snakes: the confidence interval and standard deviation snakes. OK

# Plot the function, the prediction and the confidence interval based on
# the MSE

fig = plt.figure()
# Standard deviation
plt.plot(x_2, y_pred_2, 'k-', lw=2, label=u'Gaussian Process mean (GP)')
plt.fill(np.concatenate([x_2, x_2[::-1]]),
        # np.concatenate([y_pred_2 - 1.9600 * sigma,
        #               (y_pred_2 + 1.9600 * sigma)[::-1]]),
        np.concatenate([y_pred_2 - 1. * sigma_2,
                       (y_pred_2 + 1. * sigma_2)[::-1]]),
        alpha=0.5, fc='g', ec='None', label='Standard deviation')

# Confidence interval
plt.plot(x, y_pred, 'k-', lw=2)
plt.fill(np.concatenate([x, x[::-1]]),
        # np.concatenate([y_pred - 1.9600 * sigma,
        #               (y_pred + 1.9600 * sigma)[::-1]]),
        np.concatenate([y_pred - CL_size * sigma,
                       (y_pred + CL_size * sigma)[::-1]]),
        alpha=0.5, fc='b', ec='None', label='{0}% confidence interval GP'.format(CL_label))

plt.xlim(x_IntervalPlots)
plt.ylim(yIntervalPlots)
plt.title(TitleForPlots)

#-------------------
plt.grid(True)

plt.legend(loc='lower left', fontsize=LegendFontSize)
plt.xlabel('Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
plt.ylabel('Absolute magnitude')

# plt.savefig(DirSaveData+'Plot_GP_Template_.png')
plt.savefig(DirSaveData+'Plot_GP_Template_CL_StdDev_.png', dpi=PlotResolution)
# plt.show()
plt.close()

plt.close();plt.close();plt.close();

##############################################################################80

# # Overplot: GP template and the light curve data
#
# Optional: Overplot also the weighted-mean template

# ## Plot: LC data  + confidence interval snake

# Plot the function, the prediction and the confidence interval based on
# the MSE

if Normalized == False:
    fig = plt.figure()

    #-----------------------------------------------
    #      PLOT GP TEMPLATE

    plt.fill(np.concatenate([x, x[::-1]]),
            # np.concatenate([y_pred - 1.9600 * sigma,
            #               (y_pred + 1.9600 * sigma)[::-1]]),
            np.concatenate([y_pred - CL_size * sigma,
                           (y_pred + CL_size * sigma)[::-1]]),
            alpha=.6, fc='b', ec='None', label='{0}% confidence interval for GP'.format(CL_label))

    #-----------------------------------------------
    #      PLOT SMOOTHED WEIGHTED MEAN TEMPLATE
    """
    # Mean template: Wood-Vasey algorithm
    plt.plot(Template_WeightMeanSmooth[:,0], Template_WeightMeanSmooth[:,1],
             color='red', lw=2, ls='-', alpha=1, label='Smooth weighted average (SWA)')

    plt.fill(np.concatenate([Template_WeightMeanSmooth[:,0], Template_WeightMeanSmooth[:,0][::-1]]),
            np.concatenate([Template_WeightMeanSmooth[:,1] - 1.960 * np.array(Template_WeightMeanSmooth[:,2]),
                           (Template_WeightMeanSmooth[:,1] + 1.960 * np.array(Template_WeightMeanSmooth[:,2]) )[::-1]]),
            alpha=1, fc='orange', ec='None', label='95% confidence interval for SWA')
    """
    #-----------------------------------------------
    # THE LIGHT-CURVE DATA. Including peculiar velocity

    plt.errorbar(DataLC[:,0], DataLC[:,1], error_M_PecVel, ls='', fmt='grey', alpha=0.3)

    #-----------------------------------------------

    # PLOT GP TEMPLATE mean
    plt.plot(x, y_pred, 'k-', lw=2, label=u'Gaussian Process mean (GP)')

    plt.xlim(x_IntervalPlots)
    plt.ylim(yIntervalPlots)
    plt.title(TitleForPlots)

    #---------------------

    plt.grid(True)

    plt.legend(loc='lower left', fontsize=LegendFontSize)
    plt.xlabel('Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
    plt.ylabel('Absolute magnitude')

    plt.savefig(DirSaveData+'Plot_GP_Template_CL_Data_.png', dpi=PlotResolution)
    # plt.show()
    plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ## Plot: LC data  + standard deviation snake

# Plot the function, the prediction and the confidence interval based on
# the MSE

if Normalized == False:
    fig = plt.figure()

    #-----------------------------------------------
    #      PLOT GP TEMPLATE

    plt.fill(np.concatenate([x_2, x_2[::-1]]),
            # np.concatenate([y_pred_2 - 1.9600 * sigma,
            #               (y_pred_2 + 1.9600 * sigma)[::-1]]),
            np.concatenate([y_pred_2 - 1. * sigma_2,
                           (y_pred_2 + 1. * sigma_2)[::-1]]),
            alpha=.5, fc='g', ec='None', label='Standard deviation')

    #-----------------------------------------------
    #      PLOT SMOOTHED WEIGHTED MEAN TEMPLATE
    """
    # Mean template: Wood-Vasey algorithm
    plt.plot(Template_WeightMeanSmooth[:,0], Template_WeightMeanSmooth[:,1],
             color='red', lw=2, ls='-', alpha=1, label='Smooth weighted average (SWA)')

    plt.fill(np.concatenate([Template_WeightMeanSmooth[:,0], Template_WeightMeanSmooth[:,0][::-1]]),
            np.concatenate([Template_WeightMeanSmooth[:,1] - 1.960 * np.array(Template_WeightMeanSmooth[:,2]),
                           (Template_WeightMeanSmooth[:,1] + 1.960 * np.array(Template_WeightMeanSmooth[:,2]) )[::-1]]),
            alpha=1, fc='orange', ec='None', label='95% confidence interval for SWA')
    """
    #-----------------------------------------------
    # THE LIGHT-CURVE DATA. Including peculiar velocity

    plt.errorbar(DataLC[:,0], DataLC[:,1], error_M_PecVel, ls='', fmt='grey', alpha=0.3)

    #-----------------------------------------------

    # PLOT GP TEMPLATE mean
    plt.plot(x_2, y_pred_2, 'k-', lw=2, label=u'Gaussian Process mean (GP)')

    plt.xlim(x_IntervalPlots)
    plt.ylim(yIntervalPlots)
    plt.title(TitleForPlots)

    #---------------------

    plt.grid(True)

    plt.legend(loc='lower left', fontsize=LegendFontSize)
    plt.xlabel('Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
    plt.ylabel('Absolute magnitude')

    plt.savefig(DirSaveData+'Plot_GP_Template_StdDev_Data_.png', dpi=PlotResolution)
    # plt.show()
    plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ## Overlay plots of the LC data and both confidence interval and standard deviation snakes

# Plot the function, the prediction and the confidence interval based on
# the MSE

if Normalized == False:
    fig = plt.figure()
    # Standard deviation
    plt.plot(x_2, y_pred_2, 'k-', lw=2, label=u'Gaussian Process mean (GP)')
    plt.fill(np.concatenate([x_2, x_2[::-1]]),
            # np.concatenate([y_pred_2 - 1.9600 * sigma,
            #               (y_pred_2 + 1.9600 * sigma)[::-1]]),
            np.concatenate([y_pred_2 - 1. * sigma_2,
                           (y_pred_2 + 1. * sigma_2)[::-1]]),
            alpha=0.5, fc='g', ec='None', label='Standard deviation')

    # Confidence interval
    plt.plot(x, y_pred, 'k-', lw=2)
    plt.fill(np.concatenate([x, x[::-1]]),
            # np.concatenate([y_pred - 1.9600 * sigma,
            #               (y_pred + 1.9600 * sigma)[::-1]]),
            np.concatenate([y_pred - CL_size * sigma,
                           (y_pred + CL_size * sigma)[::-1]]),
            alpha=0.5, fc='b', ec='None', label='{0}% confidence interval GP'.format(CL_label))

    #-----------------------------------------------
    # THE LIGHT-CURVE DATA. Including peculiar velocity

    plt.errorbar(DataLC[:,0], DataLC[:,1], error_M_PecVel, ls='', fmt='grey', alpha=0.3)

    #---------------------

    plt.xlim(x_IntervalPlots)
    plt.ylim(yIntervalPlots)
    plt.title(TitleForPlots)

    #---------------------

    plt.grid(True)

    plt.legend(loc='lower left', fontsize=LegendFontSize)
    plt.xlabel('Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
    plt.ylabel('Absolute magnitude')

    # plt.savefig(DirSaveData+'Plot_GP_Template_.png')
    plt.savefig(DirSaveData+'Plot_GP_Template_CL_StdDev_Data_.png', dpi=PlotResolution)
    # plt.show()
    plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ### Residual plot. OK

# Plot the function, the prediction and the confidence interval based on
# the MSE

phaseInt = np.linspace(-10, 60, 10)
zeroInt = np.zeros(len(phaseInt))

fig = plt.figure()

# Standard deviation
plt.fill(np.concatenate([x_2, x_2[::-1]]),
        np.concatenate([0 - 1. * sigma_2,
                       (0 + 1. * sigma_2)[::-1]]),
        alpha=0.5, fc='g', ec='None', label='Standard deviation')

# Confidence interval
plt.fill(np.concatenate([x, x[::-1]]),
        np.concatenate([0 - CL_size * sigma,
                       (0 + CL_size * sigma)[::-1]]),
        alpha=0.5, fc='b', ec='None', label='{0}% confidence interval GP'.format(CL_label))

#--------------------------------------------------
# THE LIGHT-CURVE DATA. Including peculiar velocity

# plt.errorbar(DataLC[:,0], DataLC[:,1], error_M_PecVel, ls='', fmt='grey', alpha=0.3)

#-------------------
# A black line in the center of the residual plot

plt.plot(phaseInt, zeroInt, color='black', lw=2)

#-------------------
plt.xlim(x_IntervalPlots)
plt.ylim(y_Interval_Residual)
plt.title(TitleForPlots)

#-------------------
plt.grid(True)

plt.legend(loc='lower left', fontsize=LegendFontSize)
plt.xlabel('Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
# plt.ylabel('Residual (mag)') # old
plt.ylabel('Sample Std Dev (mag)')

# plt.savefig(DirSaveData+'Plot_GP_Template_.png')
plt.savefig(DirSaveData+'Plot_GP_Template_CL_StdDev_Residual_.png', dpi=PlotResolution)
# plt.show()
plt.close()

plt.close();plt.close();plt.close();

#-----------------------------------------------------------------------------80

# ## Template + Residual. OK

TitleForPlots2 = Band+' band template'

# Array to create the central line:
phaseInt = np.linspace(-10, 60, 10)
zeroInt = np.zeros(len(phaseInt))

plt.clf()
# fig, axes = plt.subplots(2, 1, sharex=True,  gridspec_kw = {'width_ratios':[3, 1]})
fig, axes = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[3, 1]}  , figsize=(9, 9),)

# Close the space between the subplots
plt.subplots_adjust(hspace = .005)

#--------------------------------------------------
#     TEMPLATE

#- Standard deviation
axes[0].plot(x_2, y_pred_2, 'k-', lw=2, label='Gaussian Process mean (GP)')
axes[0].fill(np.concatenate([x_2, x_2[::-1]]),
        np.concatenate([y_pred_2 - 1. * sigma_2,
                       (y_pred_2 + 1. * sigma_2)[::-1]]),
        alpha=0.5, fc='g', ec='None', label='Standard deviation')

#- Confidence interval
axes[0].plot(x, y_pred, 'k-', lw=2)
axes[0].fill(np.concatenate([x, x[::-1]]),
        np.concatenate([y_pred - CL_size * sigma,
                       (y_pred + CL_size * sigma)[::-1]]),
        alpha=0.5, fc='b', ec='None', label='{0}% confidence interval GP'.format(CL_label))

axes[0].legend(loc='lower left', fontsize=LegendFontSize) # no working! :(
# axes[0].legend() # no working! :(
# axes[0].get_legend() # no working! :(

#-- THE LIGHT-CURVE DATA. Including peculiar velocity

if Normalized == False:
    axes[0].errorbar(DataLC[:,0], DataLC[:,1], error_M_PecVel, ls='', fmt='grey', alpha=0.3)

#--------

axes[0].grid(True)

axes[0].set_ylim(yIntervalPlots)
axes[0].set_ylabel('Absolute magnitude')
axes[0].set_title(TitleForPlots2)

#--------------------------------------------------
#     RESIDUAL

# The central black line:
axes[1].plot(phaseInt, zeroInt, color='black', lw=2)

#- Standard deviation
axes[1].fill(np.concatenate([x_2, x_2[::-1]]),
        np.concatenate([0 - 1. * sigma_2,
                       (0 + 1. * sigma_2)[::-1]]),
        alpha=0.5, fc='g', ec='None', label='Standard deviation')

#- Confidence interval
axes[1].fill(np.concatenate([x, x[::-1]]),
        np.concatenate([0 - CL_size * sigma,
                       (0 + CL_size * sigma)[::-1]]),
        alpha=0.5, fc='b', ec='None', label='{0}% confidence interval GP'.format(CL_label))

axes[1].grid(True)
axes[1].set_xlim(x_IntervalPlots)
# axes[1].set_ylabel('Residual (mag)')
axes[1].set_ylabel('Sample Std Dev (mag)')
axes[1].set_xlabel('Phase = (MJD - $T_{Bmax}$)/(1+$z_{hel}$)')
# plt.legend(loc='upper left', fontsize=LegendFontSize)

plt.savefig(DirSaveData+'Plot_GP_Template_CL_StdDev_Both_.png', dpi=PlotResolution)
plt.close()

#-------------------
plt.close();plt.close();plt.close();

print "# All done :)"

