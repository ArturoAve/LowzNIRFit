#!/usr/bin/env python
# coding: utf-8

# ## Create the final figures and latex table for the paper

import numpy as np
import os # To use command line like instructions

# To read arguments in command line
# Used in the ".py" version of this notebook.
import sys

#--------------------------------------------------------60
code_created_by = 'Arturo_Avelino'
# On date: 2018-06-30 (yyyy.mm.dd)
code_name = '15_Latex_scatter.ipynb'
version_code = '0.1.7'
last_update = '2019.01.28'
#--------------------------------------------------------60

##############################################################################80

# # User

## Terminal or notebook version of this script?
ScriptVersion = 'terminal' # ( terminal , notebook )

#--------------------------------------------------------60
if ScriptVersion == 'notebook':

    #   1st case

    ## Peculiar velocity
    vpecFix = 150 # km/s # Notebook version.
    # print 'vpecFix:', vpecFix, type(vpecFix)

    ## Write the simple RMS also?
    WriteSimple_RMS = True # Notebook version.
    # print 'WriteSimple_RMS:', WriteSimple_RMS, type(WriteSimple_RMS)

    ## Write also the chi2_dof values:
    chi2dof_print = False # Notebook version.
    # print 'chi2dof_print:', chi2dof_print, type(chi2dof_print)

    #----------------------------------------

    #    2nd case

    ## Write a second column with the intrinsic dispersion computed
    ## assuming another value for the peculiar velocity uncertainty?
    Write_2nd_vpec = True # Notebook version.
    # print 'Write_2nd_vpec:', Write_2nd_vpec, type(Write_2nd_vpec)

    ## What is the other value of the peculiar velocity uncertainty?:
    vpecFix_2 = 250  # km/s # Notebook version.
    # print 'vpecFix_2:', vpecFix_2, type(vpecFix_2)

    ## Write the wRMS and its uncertainty also?
    Write_wRMS_2nd = False # Notebook version.
    # print 'Write_wRMS_2nd:', Write_wRMS_2nd, type(Write_wRMS_2nd)

    ## Write the simple RMS also?
    WriteSimple_RMS_2nd = False # Notebook version.
    # print 'WriteSimple_RMS_2nd:', WriteSimple_RMS_2nd, type(WriteSimple_RMS_2nd)

    DirSaveOutput = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Latex/Plots_LatexTables_ForThePaper/'

#-----------------------------------------------------------------------------80

elif ScriptVersion == 'terminal':

    #   1st case

    ## Peculiar velocity
    vpecFix = int(sys.argv[1])
    # print 'vpecFix:', vpecFix, type(vpecFix)

    ## Write the simple RMS also?
    WriteSimple_RMS = sys.argv[2] == 'True'
    # print 'WriteSimple_RMS:', WriteSimple_RMS, type(WriteSimple_RMS)

    ## Write also the chi2_dof values:
    chi2dof_print = sys.argv[3] == 'True'
    # print 'chi2dof_print:', chi2dof_print, type(chi2dof_print)

    #----------------------------------------

    #    2nd case

    ## Write a second column with the intrinsic dispersion computed
    ## assuming another value for the peculiar velocity uncertainty?
    Write_2nd_vpec = sys.argv[4] == 'True'
    # print 'Write_2nd_vpec:', Write_2nd_vpec, type(Write_2nd_vpec)

    ## What is the other value of the peculiar velocity uncertainty?:
    vpecFix_2 = int(sys.argv[5])
    # print 'vpecFix_2:', vpecFix_2, type(vpecFix_2)

    ## Write the wRMS and its uncertainty also?
    Write_wRMS_2nd = sys.argv[6] == 'True'
    # print 'Write_wRMS_2nd:', Write_wRMS_2nd, type(Write_wRMS_2nd)

    ## Write the simple RMS also?
    WriteSimple_RMS_2nd = sys.argv[7] == 'True'
    # print 'WriteSimple_RMS_2nd:', WriteSimple_RMS_2nd, type(WriteSimple_RMS_2nd)

    DirSaveOutput = sys.argv[8]

###################################################################

MainDir = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/TheTemplates/'

# Dir to save the  python output code generated to create the main body
# of this section
DirSavePythonOutput = DirSaveOutput

###################################################################

"""
vpecFix: 150 <type 'int'>
WriteSimple_RMS: True <type 'bool'>
chi2dof_print: False <type 'bool'>
Write_2nd_vpec: True <type 'bool'>
vpecFix_2: 250 <type 'int'>
Write_wRMS_2nd: False <type 'bool'>
WriteSimple_RMS_2nd: False <type 'bool'>
"""
0

##############################################################################80

# # AUTOMATIC

# Create the directory to save the output if it doesn't exist yet.
if not os.path.exists(DirSaveOutput): os.makedirs(DirSaveOutput)

# Get the current date and time
import datetime

# Read the time and date now
now = datetime.datetime.now()

# Information about the author and date to write down as header
# in each latex table.

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd).")
text_Date   = '%% On date: %s \n'%text_timenow
text_Author = '% Data table created by: Arturo Avelino \n'
text_script = '%% Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
text_line = '%'+'-'*60 + '\n'

# ### Get the name of this ipython notebook
# To print it in the output text files as reference

# %%javascript
# var kernel = IPython.notebook.kernel;
# var thename = window.document.getElementById("notebook_name").innerHTML;
# var command = "NotebookName = " + "'"+thename+".ipynb"+"'";
# kernel.execute(command);

# print '#', (NotebookName)
# Update_zcmb_in_SNANA_datafiles_v1_0.ipynb

##############################################################################80

# # Latex table: Intrinsic dispersion $\sigma_{\rm int}$ and wRMS

# ### Load the data

#################################################

# Y band

#    vpec = 150 km/s

Y_HubbleDir_Tp =  MainDir+'Y_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi3_EBVh0.4_Method7_MinData3_vpec%s_ok/plots_HD/'%vpecFix

Y_HubbleDir_GP =  MainDir+'Y_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_NIRmax/plots_HD/'%vpecFix

Y_HubbleDir_GP_tBmax =  MainDir+'Y_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_Bmax/plots_HD/'%vpecFix

Y_wRMS_Tp = np.genfromtxt(Y_HubbleDir_Tp+'Summary_HDScatter_RMS_.txt')
Y_wRMS_GP = np.genfromtxt(Y_HubbleDir_GP+'Summary_HDScatter_RMS_.txt')
Y_wRMS_GP_tBmax = np.genfromtxt(Y_HubbleDir_GP_tBmax+'Summary_HDScatter_RMS_.txt')

#-------------------
#    vpec = 250 km/s

if Write_2nd_vpec:

    Y_HubbleDir_Tp_2 =  MainDir+'Y_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi3_EBVh0.4_Method7_MinData3_vpec%s_ok/plots_HD/'%vpecFix_2

    Y_HubbleDir_GP_2 =  MainDir+'Y_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_NIRmax/plots_HD/'%vpecFix_2

    Y_HubbleDir_GP_tBmax_2 =  MainDir+'Y_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_Bmax/plots_HD/'%vpecFix_2

    Y_wRMS_Tp_2 = np.genfromtxt(Y_HubbleDir_Tp_2+'Summary_HDScatter_RMS_.txt')
    Y_wRMS_GP_2 = np.genfromtxt(Y_HubbleDir_GP_2+'Summary_HDScatter_RMS_.txt')
    Y_wRMS_GP_tBmax_2 = np.genfromtxt(Y_HubbleDir_GP_tBmax_2+'Summary_HDScatter_RMS_.txt')

#################################################

# J band

#    vpec = 150 km/s

J_HubbleDir_Tp = MainDir+'J_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi_1e6_EBVh0.4_Method7_MinData3_vpec%s_ok/plots_HD/'%vpecFix

J_HubbleDir_GP =  MainDir+'J_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_NIRmax/plots_HD/'%vpecFix

_HubbleDir_GP_tBmax =  MainDir+'J_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_Bmax/plots_HD/'%vpecFix

J_wRMS_Tp = np.genfromtxt(J_HubbleDir_Tp+'Summary_HDScatter_RMS_.txt')
J_wRMS_GP = np.genfromtxt(J_HubbleDir_GP+'Summary_HDScatter_RMS_.txt')
J_wRMS_GP_tBmax = np.genfromtxt(_HubbleDir_GP_tBmax+'Summary_HDScatter_RMS_.txt')

#-------------------
#    vpec = 250 km/s

if Write_2nd_vpec:

    J_HubbleDir_Tp_2 = MainDir+'J_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi_1e6_EBVh0.4_Method7_MinData3_vpec%s_ok/plots_HD/'%vpecFix_2

    J_HubbleDir_GP_2 =  MainDir+'J_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_NIRmax/plots_HD/'%vpecFix_2

    _HubbleDir_GP_tBmax_2 =  MainDir+'J_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_Bmax/plots_HD/'%vpecFix_2

    J_wRMS_Tp_2 = np.genfromtxt(J_HubbleDir_Tp_2+'Summary_HDScatter_RMS_.txt')
    J_wRMS_GP_2 = np.genfromtxt(J_HubbleDir_GP_2+'Summary_HDScatter_RMS_.txt')
    J_wRMS_GP_tBmax_2 = np.genfromtxt(_HubbleDir_GP_tBmax_2+'Summary_HDScatter_RMS_.txt')

#################################################

# H band

H_HubbleDir_Tp = MainDir+'H_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi_1e6_EBVh0.4_Method7_MinData3_vpec%s_ok/plots_HD/'%vpecFix

H_HubbleDir_GP =  MainDir+'H_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_NIRmax/plots_HD/'%vpecFix

H_HubbleDir_GP_tBmax =  MainDir+'H_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_Bmax/plots_HD/'%vpecFix

H_wRMS_Tp = np.genfromtxt(H_HubbleDir_Tp+'Summary_HDScatter_RMS_.txt')
H_wRMS_GP = np.genfromtxt(H_HubbleDir_GP+'Summary_HDScatter_RMS_.txt')
H_wRMS_GP_tBmax = np.genfromtxt(H_HubbleDir_GP_tBmax+'Summary_HDScatter_RMS_.txt')

#-------------------
#    vpec = 250 km/s

if Write_2nd_vpec:

    H_HubbleDir_Tp_2 = MainDir+'H_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi_1e6_EBVh0.4_Method7_MinData3_vpec%s_ok/plots_HD/'%vpecFix_2

    H_HubbleDir_GP_2 =  MainDir+'H_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_NIRmax/plots_HD/'%vpecFix_2

    H_HubbleDir_GP_tBmax_2 =  MainDir+'H_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_Bmax/plots_HD/'%vpecFix_2

    H_wRMS_Tp_2 = np.genfromtxt(H_HubbleDir_Tp_2+'Summary_HDScatter_RMS_.txt')
    H_wRMS_GP_2 = np.genfromtxt(H_HubbleDir_GP_2+'Summary_HDScatter_RMS_.txt')
    H_wRMS_GP_tBmax_2 = np.genfromtxt(H_HubbleDir_GP_tBmax_2+'Summary_HDScatter_RMS_.txt')

#################################################

# K band

K_HubbleDir_Tp = MainDir+'K_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi4_EBVh0.4_Method7_MinData3_vpec%s_ok/plots_HD/'%vpecFix

K_HubbleDir_GP =  MainDir+'K_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_NIRmax/plots_HD/'%vpecFix

K_HubbleDir_GP_tBmax =  MainDir+'K_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_Bmax/plots_HD/'%vpecFix

K_wRMS_Tp = np.genfromtxt(K_HubbleDir_Tp+'Summary_HDScatter_RMS_.txt')
K_wRMS_GP = np.genfromtxt(K_HubbleDir_GP+'Summary_HDScatter_RMS_.txt')
K_wRMS_GP_tBmax = np.genfromtxt(K_HubbleDir_GP_tBmax+'Summary_HDScatter_RMS_.txt')

#-------------------
#    vpec = 250 km/s

if Write_2nd_vpec:

    K_HubbleDir_Tp_2 = MainDir+'K_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi4_EBVh0.4_Method7_MinData3_vpec%s_ok/plots_HD/'%vpecFix_2

    K_HubbleDir_GP_2 =  MainDir+'K_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_NIRmax/plots_HD/'%vpecFix_2

    K_HubbleDir_GP_tBmax_2 =  MainDir+'K_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_Bmax/plots_HD/'%vpecFix_2

    K_wRMS_Tp_2 = np.genfromtxt(K_HubbleDir_Tp_2+'Summary_HDScatter_RMS_.txt')
    K_wRMS_GP_2 = np.genfromtxt(K_HubbleDir_GP_2+'Summary_HDScatter_RMS_.txt')
    K_wRMS_GP_tBmax_2 = np.genfromtxt(K_HubbleDir_GP_tBmax_2+'Summary_HDScatter_RMS_.txt')

#################################################

# Any_YJHK

AnyNIR_HubbleDir_Tp = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/Template/AllBands/'%vpecFix
AnyNIR_HubbleDir_Tp_GPsubsample = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/Template__GP_subsample/AllBands/'%vpecFix
AnyNIR_HubbleDir_GP = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GaussianProcess/AllBands/'%vpecFix
AnyNIR_HubbleDir_GP_tBmax = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GP_Bmax/AllBands/'%vpecFix

AnyNIR_wRMS_Tp = np.genfromtxt(AnyNIR_HubbleDir_Tp+'Summary_HDScatter_RMS_.txt')
AnyNIR_wRMS_Tp_GPsubsample = np.genfromtxt(AnyNIR_HubbleDir_Tp_GPsubsample+'Summary_HDScatter_RMS_.txt')
AnyNIR_wRMS_GP = np.genfromtxt(AnyNIR_HubbleDir_GP+'Summary_HDScatter_RMS_.txt')
AnyNIR_wRMS_GP_tBmax = np.genfromtxt(AnyNIR_HubbleDir_GP_tBmax+'Summary_HDScatter_RMS_.txt')

#-------------------
#    vpec = 250 km/s

if Write_2nd_vpec:

    AnyNIR_HubbleDir_Tp_2 = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/Template/AllBands/'%vpecFix_2
    AnyNIR_HubbleDir_Tp_GPsubsample_2 = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/Template__GP_subsample/AllBands/'%vpecFix_2
    AnyNIR_HubbleDir_GP_2 = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GaussianProcess/AllBands/'%vpecFix_2
    AnyNIR_HubbleDir_GP_tBmax_2 = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GP_Bmax/AllBands/'%vpecFix_2

    AnyNIR_wRMS_Tp_2 = np.genfromtxt(AnyNIR_HubbleDir_Tp_2+'Summary_HDScatter_RMS_.txt')
    AnyNIR_wRMS_Tp_GPsubsample_2 = np.genfromtxt(AnyNIR_HubbleDir_Tp_GPsubsample_2+'Summary_HDScatter_RMS_.txt')
    AnyNIR_wRMS_GP_2 = np.genfromtxt(AnyNIR_HubbleDir_GP_2+'Summary_HDScatter_RMS_.txt')
    AnyNIR_wRMS_GP_tBmax_2 = np.genfromtxt(AnyNIR_HubbleDir_GP_tBmax_2+'Summary_HDScatter_RMS_.txt')

#################################################

# JH bands

JH_HubbleDir_Tp = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/Template/JH/'%vpecFix
JH_HubbleDir_GP = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GaussianProcess/JH/'%vpecFix
JH_HubbleDir_GP_tBmax = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GP_Bmax/JH/'%vpecFix

JH_wRMS_Tp = np.genfromtxt(JH_HubbleDir_Tp+'Summary_HDScatter_RMS_.txt')
JH_wRMS_GP = np.genfromtxt(JH_HubbleDir_GP+'Summary_HDScatter_RMS_.txt')
JH_wRMS_GP_tBmax = np.genfromtxt(JH_HubbleDir_GP_tBmax+'Summary_HDScatter_RMS_.txt')

#-------------------
#    vpec = 250 km/s

if Write_2nd_vpec:

    JH_HubbleDir_Tp_2 = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/Template/JH/'%vpecFix_2
    JH_HubbleDir_GP_2 = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GaussianProcess/JH/'%vpecFix_2
    JH_HubbleDir_GP_tBmax_2 = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GP_Bmax/JH/'%vpecFix_2

    JH_wRMS_Tp_2 = np.genfromtxt(JH_HubbleDir_Tp_2+'Summary_HDScatter_RMS_.txt')
    JH_wRMS_GP_2 = np.genfromtxt(JH_HubbleDir_GP_2+'Summary_HDScatter_RMS_.txt')
    JH_wRMS_GP_tBmax_2 = np.genfromtxt(JH_HubbleDir_GP_tBmax_2+'Summary_HDScatter_RMS_.txt')

#################################################

# YJH bands

YJH_HubbleDir_Tp = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/Template/YJH/'%vpecFix
YJH_HubbleDir_GP = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GaussianProcess/YJH/'%vpecFix
YJH_HubbleDir_GP_tBmax = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GP_Bmax/YJH/'%vpecFix

YJH_wRMS_Tp = np.genfromtxt(YJH_HubbleDir_Tp+'Summary_HDScatter_RMS_.txt')
YJH_wRMS_GP = np.genfromtxt(YJH_HubbleDir_GP+'Summary_HDScatter_RMS_.txt')
YJH_wRMS_GP_tBmax = np.genfromtxt(YJH_HubbleDir_GP_tBmax+'Summary_HDScatter_RMS_.txt')

#-------------------
#    vpec = 250 km/s

if Write_2nd_vpec:

    YJH_HubbleDir_Tp_2 = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/Template/YJH/'%vpecFix_2
    YJH_HubbleDir_GP_2 = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GaussianProcess/YJH/'%vpecFix_2
    YJH_HubbleDir_GP_tBmax_2 = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GP_Bmax/YJH/'%vpecFix_2

    YJH_wRMS_Tp_2 = np.genfromtxt(YJH_HubbleDir_Tp_2+'Summary_HDScatter_RMS_.txt')
    YJH_wRMS_GP_2 = np.genfromtxt(YJH_HubbleDir_GP_2+'Summary_HDScatter_RMS_.txt')
    YJH_wRMS_GP_tBmax_2 = np.genfromtxt(YJH_HubbleDir_GP_tBmax_2+'Summary_HDScatter_RMS_.txt')

#################################################

# JHK bands

JHK_HubbleDir_Tp = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/Template/JHK/'%vpecFix
JHK_HubbleDir_GP = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GaussianProcess/JHK/'%vpecFix
JHK_HubbleDir_GP_tBmax = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GP_Bmax/JHK/'%vpecFix

JHK_wRMS_Tp = np.genfromtxt(JHK_HubbleDir_Tp+'Summary_HDScatter_RMS_.txt')
JHK_wRMS_GP = np.genfromtxt(JHK_HubbleDir_GP+'Summary_HDScatter_RMS_.txt')
JHK_wRMS_GP_tBmax = np.genfromtxt(JHK_HubbleDir_GP_tBmax+'Summary_HDScatter_RMS_.txt')

#-------------------
#    vpec = 250 km/s

if Write_2nd_vpec:

    JHK_HubbleDir_Tp_2 = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/Template/JHK/'%vpecFix_2
    JHK_HubbleDir_GP_2 = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GaussianProcess/JHK/'%vpecFix_2
    JHK_HubbleDir_GP_tBmax_2 = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/GP_Bmax/JHK/'%vpecFix_2

    JHK_wRMS_Tp_2 = np.genfromtxt(JHK_HubbleDir_Tp_2+'Summary_HDScatter_RMS_.txt')
    JHK_wRMS_GP_2 = np.genfromtxt(JHK_HubbleDir_GP_2+'Summary_HDScatter_RMS_.txt')
    JHK_wRMS_GP_tBmax_2 = np.genfromtxt(JHK_HubbleDir_GP_tBmax_2+'Summary_HDScatter_RMS_.txt')

#################################################

# SALT2

# DirSalt2Mu = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/\
# 10Compute/TheTemplates/AllBands/Plots/HubbleDiagram/SALT2/plots_HD/'

DirSalt2Mu = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/TheTemplates/AllBands/Plots/HubbleDiagram_vpec%s/SALT2/plots_HD/'%vpecFix
SALT2_RMSData = np.genfromtxt(DirSalt2Mu+"Summary_HDScatter_RMS_.txt")

#-------------------
#    vpec = 250 km/s

DirSalt2Mu_2 = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/TheTemplates/AllBands/Plots/HubbleDiagram_vpec%s/SALT2/plots_HD/'%vpecFix_2
SALT2_RMSData_2 = np.genfromtxt(DirSalt2Mu_2+"Summary_HDScatter_RMS_.txt")

#################################################

# Snoopy

# DirSnoopyMu = MainDir+'AllBands/Plots/HubbleDiagram/Snoopy_opt/plots_HD/'

DirSnoopyMu = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/TheTemplates/AllBands/Plots/HubbleDiagram_vpec%s/Snoopy_opt/plots_HD/'%vpecFix
Snoopy_RMSData = np.genfromtxt(DirSnoopyMu+"Summary_HDScatter_RMS_.txt")

#-------------------
#    vpec = 250 km/s

DirSnoopyMu_2 = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/TheTemplates/AllBands/Plots/HubbleDiagram_vpec%s/Snoopy_opt/plots_HD/'%vpecFix_2
Snoopy_RMSData_2 = np.genfromtxt(DirSnoopyMu_2+"Summary_HDScatter_RMS_.txt")

#-------------------
print "# All the 'Summary_HDScatter_RMS_.txt' files read with no issues."

#-----------------------------------------------------------------------------80

# ###  Main cell: Write the latex table on the intrinsic scatter.
#
# Main cell of this section

# Write the latex table on the intrinsic scatter.

textfile_1 = open(DirSaveOutput+'Table_Latex_scatter_.tex', 'w')
textfile_2 = open(DirSaveOutput+'Table_Latex_scatter_.dat', 'w')
textfile_3 = open(DirSaveOutput+'Table_Latex_scatter_GPNIRmax.dat', 'w')
textfile_4 = open(DirSaveOutput+'Table_Latex_scatter_GPBmax.dat', 'w')
textfile_5 = open(DirSaveOutput+'Table_Latex_scatter_Template.dat', 'w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_Date   = '%% On date: %s \n'%text_timenow

textfile_1.write(text_line)
textfile_1.write(text_Author); textfile_1.write(text_Date);
textfile_1.write(text_script); textfile_1.write(text_line);
textfile_1.write(' \n')

textfile_2.write('#'+text_line)
textfile_2.write('#'+text_Author); textfile_2.write('#'+text_Date);
textfile_2.write('#'+text_script); textfile_2.write('#'+text_line);

textfile_3.write('#'+text_line)
textfile_3.write('#'+text_Author); textfile_3.write('#'+text_Date);
textfile_3.write('#'+text_script); textfile_3.write('#'+text_line);

textfile_4.write('#'+text_line)
textfile_4.write('#'+text_Author); textfile_4.write('#'+text_Date);
textfile_4.write('#'+text_script); textfile_4.write('#'+text_line);

textfile_5.write('#'+text_line)
textfile_5.write('#'+text_Author); textfile_5.write('#'+text_Date);
textfile_5.write('#'+text_script); textfile_5.write('#'+text_line);

if chi2dof_print:

    Y_chi2dof_Tp_1 = "(%.2f)"%Y_wRMS_Tp[4][0]
    Y_chi2dof_GP_1 = "(%.2f)"%Y_wRMS_GP[4][0]
    if Write_2nd_vpec:
        Y_chi2dof_Tp_2 = " "
        Y_chi2dof_GP_2 = " "

    J_chi2dof_Tp_1 = "(%.2f)"%J_wRMS_Tp[4][0]
    J_chi2dof_GP_1 = "(%.2f)"%J_wRMS_GP[4][0]
    if Write_2nd_vpec:
        J_chi2dof_Tp_2 = " "
        J_chi2dof_GP_2 = " "

    H_chi2dof_Tp_1 = "(%.2f)"%H_wRMS_Tp[4][0]
    H_chi2dof_GP_1 = "(%.2f)"%H_wRMS_GP[4][0]
    if Write_2nd_vpec:
        H_chi2dof_Tp_2 = " "
        H_chi2dof_GP_2 = " "

    K_chi2dof_Tp_1 = "(%.2f)"%K_wRMS_Tp[4][0]
    K_chi2dof_GP_1 = "(%.2f)"%K_wRMS_GP[4][0]
    if Write_2nd_vpec:
        K_chi2dof_Tp_2 = " "
        K_chi2dof_GP_2 = " "

    AnyNIR_chi2dof_Tp_1 = "(%.2f)"%AnyNIR_wRMS_Tp[4][0]
    AnyNIR_chi2dof_GP_1 = "(%.2f)"%AnyNIR_wRMS_GP[4][0]
    if Write_2nd_vpec:
        AnyNIR_chi2dof_Tp_2 = " "
        AnyNIR_chi2dof_GP_2 = " "

    JH_chi2dof_Tp_1 = "(%.2f)"%JH_wRMS_Tp[4][0]
    JH_chi2dof_GP_1 = "(%.2f)"%JH_wRMS_GP[4][0]
    if Write_2nd_vpec:
        JH_chi2dof_Tp_2 = " "
        JH_chi2dof_GP_2 = " "

    YJH_chi2dof_Tp_1 = "(%.2f)"%YJH_wRMS_Tp[4][0]
    YJH_chi2dof_GP_1 = "(%.2f)"%YJH_wRMS_GP[4][0]
    if Write_2nd_vpec:
        YJH_chi2dof_Tp_2 = " "
        YJH_chi2dof_GP_2 = " "

    JHK_chi2dof_Tp_1 = "(%.2f)"%JHK_wRMS_Tp[4][0]
    JHK_chi2dof_GP_1 = "(%.2f)"%JHK_wRMS_GP[4][0]
    if Write_2nd_vpec:
        JHK_chi2dof_Tp_2 = " "
        JHK_chi2dof_GP_2 = " "

    SALT2_chi2dof_1  = "(%.2f)"%SALT2_RMSData[4][0]
    Snoopy_chi2dof_1 = "(%.2f)"%Snoopy_RMSData[4][0]
    if Write_2nd_vpec:
        SALT2_chi2dof_2  = " "
        Snoopy_chi2dof_2 = " "

else:

    Y_chi2dof_Tp_1 = " "; Y_chi2dof_GP_1 = " "; Y_chi2dof_GP_tBmax_1 = " ";
    Y_chi2dof_Tp_2 = " "; Y_chi2dof_GP_2 = " "; Y_chi2dof_GP_tBmax_2 = " ";

    J_chi2dof_Tp_1 = " "; J_chi2dof_GP_1 = " "; J_chi2dof_GP_tBmax_1 = " ";
    J_chi2dof_Tp_2 = " "; J_chi2dof_GP_2 = " "; J_chi2dof_GP_tBmax_2 = " ";

    H_chi2dof_Tp_1 = " "; H_chi2dof_GP_1 = " "; H_chi2dof_GP_tBmax_1 = " ";
    H_chi2dof_Tp_2 = " "; H_chi2dof_GP_2 = " "; H_chi2dof_GP_tBmax_2 = " ";

    K_chi2dof_Tp_1 = " "; K_chi2dof_GP_1 = " "; K_chi2dof_GP_tBmax_1 = " ";
    K_chi2dof_Tp_2 = " "; K_chi2dof_GP_2 = " "; K_chi2dof_GP_tBmax_2 = " ";

    AnyNIR_chi2dof_Tp_1 = " "; AnyNIR_chi2dof_GP_1 = " "; AnyNIR_chi2dof_GP_tBmax_1 = " ";
    AnyNIR_chi2dof_Tp_2 = " "; AnyNIR_chi2dof_GP_2 = " "; AnyNIR_chi2dof_GP_tBmax_2 = " ";

    AnyNIR_chi2dof_Tp_GPsubsample_1 = " ";
    AnyNIR_chi2dof_Tp_GPsubsample_2 = " ";

    JH_chi2dof_Tp_1 = " "; JH_chi2dof_GP_1 = " "; JH_chi2dof_GP_tBmax_1 = " ";
    JH_chi2dof_Tp_2 = " "; JH_chi2dof_GP_2 = " "; JH_chi2dof_GP_tBmax_2 = " ";

    YJH_chi2dof_Tp_1 = " "; YJH_chi2dof_GP_1 = " "; YJH_chi2dof_GP_tBmax_1 = " ";
    YJH_chi2dof_Tp_2 = " "; YJH_chi2dof_GP_2 = " "; YJH_chi2dof_GP_tBmax_2 = " ";

    JHK_chi2dof_Tp_1 = " "; JHK_chi2dof_GP_1 = " "; JHK_chi2dof_GP_tBmax_1 = " ";
    JHK_chi2dof_Tp_2 = " "; JHK_chi2dof_GP_2 = " "; JHK_chi2dof_GP_tBmax_2 = " ";

    SALT2_chi2dof_1 = " "; SALT2_chi2dof_2 = " ";
    Snoopy_chi2dof_1 = " "; Snoopy_chi2dof_2 = " ";

##############################################################################80

textfile_1.write("% Intrinsic dispersion as a function of the uncertainty in the peculiar velocity, and the wRMS \n")
textfile_1.write(text_line)

text_cols = '# Band   Method     N_sn  s_int   e_s_int   s_int   \
e_s_int   wRMS    e_wRMS   RMS     e_RMS\n'

textfile_2.write(text_cols);textfile_3.write(text_cols);
textfile_4.write(text_cols);textfile_5.write(text_cols);

#-----------------------------------------------------------------------------80

# reset
salt2_RMS_1_txt=''; salt2_wRMS_2_txt=''; salt2_RMS_2_txt=''
if WriteSimple_RMS:
    salt2_RMS_1_txt = '& $%.3f \\pm %.3f$'%(SALT2_RMSData[5][2],SALT2_RMSData[5][3])
    salt2_RMS_1_dat = '%7.3f %7.3f'%(SALT2_RMSData[5][2],SALT2_RMSData[5][3])
if Write_wRMS_2nd:
    salt2_wRMS_2_txt = '& $%.3f \\pm %.3f$'%(SALT2_RMSData_2[5][0], SALT2_RMSData_2[5][1])
if WriteSimple_RMS_2nd:
    salt2_RMS_2_txt = '& $%.3f \\pm %.3f$'%(SALT2_RMSData_2[5][2],SALT2_RMSData_2[5][3])

textfile_1.write("Optical $BVR$ & SALT2    & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n"%(
    SALT2_RMSData[3][1],
    SALT2_RMSData[0][0], SALT2_RMSData[0][1], SALT2_chi2dof_1,
    SALT2_RMSData_2[0][0], SALT2_RMSData_2[0][1], SALT2_chi2dof_2,
    SALT2_RMSData[5][0], SALT2_RMSData[5][1], salt2_wRMS_2_txt,
    salt2_RMS_1_txt, salt2_RMS_2_txt ))
                           #          #
textfile_2.write('BVR      SALT2      %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    SALT2_RMSData[3][1],
    SALT2_RMSData[0][0], SALT2_RMSData[0][1], SALT2_chi2dof_1,
    SALT2_RMSData_2[0][0], SALT2_RMSData_2[0][1], SALT2_chi2dof_2,
    SALT2_RMSData[5][0], SALT2_RMSData[5][1], salt2_wRMS_2_txt,
    salt2_RMS_1_dat, salt2_RMS_2_txt ))

#--------------------------------------------------
# reset
snoopy_RMS_1_txt=''; snoopy_wRMS_2_txt=''; snoopy_RMS_2_txt=''
if WriteSimple_RMS:
    snoopy_RMS_1_txt = '& $%.3f \\pm %.3f$'%(Snoopy_RMSData[5][2],Snoopy_RMSData[5][3])
    snoopy_RMS_1_dat = '%7.3f %7.3f'%(Snoopy_RMSData[5][2],Snoopy_RMSData[5][3])
if Write_wRMS_2nd:
    snoopy_wRMS_2_txt = '& $%.3f \\pm %.3f$'%(Snoopy_RMSData_2[5][0],Snoopy_RMSData_2[5][1])
if WriteSimple_RMS_2nd:
    snoopy_RMS_2_txt = '& $%.3f \\pm %.3f$'%(Snoopy_RMSData_2[5][2],Snoopy_RMSData_2[5][3])

textfile_1.write("Optical $BVR$ & SNooPy   & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n"%(
    Snoopy_RMSData[3][1],
    Snoopy_RMSData[0][0], Snoopy_RMSData[0][1], Snoopy_chi2dof_1,
    Snoopy_RMSData_2[0][0], Snoopy_RMSData_2[0][1], Snoopy_chi2dof_2,
    Snoopy_RMSData[5][0], Snoopy_RMSData[5][1], snoopy_wRMS_2_txt,
    snoopy_RMS_1_txt, snoopy_RMS_2_txt ))
                                      #
textfile_2.write('BVR      SNooPy     %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    Snoopy_RMSData[3][1],
    Snoopy_RMSData[0][0], Snoopy_RMSData[0][1], Snoopy_chi2dof_1,
    Snoopy_RMSData_2[0][0], Snoopy_RMSData_2[0][1], Snoopy_chi2dof_2,
    Snoopy_RMSData[5][0], Snoopy_RMSData[5][1], snoopy_wRMS_2_txt,
    snoopy_RMS_1_dat, snoopy_RMS_2_txt ))

textfile_1.write("\\hline \n")
textfile_1.write('%------------------------------------- \n')

#-----------------------------------------------------------------------------80
# CASE "any YJHK_s & Template (GP subsample)"

# reset
AnyNIR_RMS_Tp_GPsubsample_1_txt='';
AnyNIR_wRMS_Tp_GPsubsample_2_txt='';
AnyNIR_RMS_Tp_GPsubsample_2_txt='';
if WriteSimple_RMS:
    AnyNIR_RMS_Tp_GPsubsample_1_txt = '& $%.3f \\pm %.3f$'%(
        AnyNIR_wRMS_Tp_GPsubsample[5][2],AnyNIR_wRMS_Tp_GPsubsample[5][3])
    AnyNIR_RMS_Tp_GPsubsample_1_dat = '%7.3f %7.3f'%(
        AnyNIR_wRMS_Tp_GPsubsample[5][2],AnyNIR_wRMS_Tp_GPsubsample[5][3])
if Write_wRMS_2nd:
    AnyNIR_wRMS_Tp_GPsubsample_2_txt = '& $%.3f \\pm %.3f$'%(
        AnyNIR_wRMS_Tp_GPsubsample_2[5][0],
        AnyNIR_wRMS_Tp_GPsubsample_2[5][1])
if WriteSimple_RMS_2nd:
    AnyNIR_RMS_Tp_GPsubsample_2_txt = '& $%.3f \\pm %.3f$'%(
        AnyNIR_wRMS_Tp_GPsubsample_2[5][2],AnyNIR_wRMS_Tp_GPsubsample_2[5][3])

textfile_1.write("any $YJHK_s$ & Template  & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n"%(
    AnyNIR_wRMS_Tp_GPsubsample[3][1],
    AnyNIR_wRMS_Tp_GPsubsample[0][0],
    AnyNIR_wRMS_Tp_GPsubsample[0][1], AnyNIR_chi2dof_Tp_GPsubsample_1,
    AnyNIR_wRMS_Tp_GPsubsample_2[0][0],
    AnyNIR_wRMS_Tp_GPsubsample_2[0][1], AnyNIR_chi2dof_Tp_GPsubsample_2,
    AnyNIR_wRMS_Tp_GPsubsample[5][0], AnyNIR_wRMS_Tp_GPsubsample[5][1],
    AnyNIR_wRMS_Tp_GPsubsample_2_txt,
    AnyNIR_RMS_Tp_GPsubsample_1_txt,
    AnyNIR_RMS_Tp_GPsubsample_2_txt ))

textfile_2.write('AnyNIR   Tp_GPsampl %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    AnyNIR_wRMS_Tp_GPsubsample[3][1],
    AnyNIR_wRMS_Tp_GPsubsample[0][0],
    AnyNIR_wRMS_Tp_GPsubsample[0][1], AnyNIR_chi2dof_Tp_GPsubsample_1,
    AnyNIR_wRMS_Tp_GPsubsample_2[0][0],
    AnyNIR_wRMS_Tp_GPsubsample_2[0][1], AnyNIR_chi2dof_Tp_GPsubsample_2,
    AnyNIR_wRMS_Tp_GPsubsample[5][0], AnyNIR_wRMS_Tp_GPsubsample[5][1],
    AnyNIR_wRMS_Tp_GPsubsample_2_txt,
    AnyNIR_RMS_Tp_GPsubsample_1_dat,
    AnyNIR_RMS_Tp_GPsubsample_2_txt ))

textfile_1.write("\\hline \n")
textfile_1.write('%------------------------------------- \n')

###############################################
# ---- PART GENERATED WITH THE CODE BELOW ---------->>

#-----------------------------------------------------------------------------80
# Data table created by: Arturo Avelino
# On date: 2018-12-20 (yyyy-mm-dd); 12:22 hrs.
# Script used: 15_Latex_scatter.ipynb
#-----------------------------------------------------------------------------80
#
# reset
Y_RMS_Tp_1_txt=''; Y_wRMS_Tp_2_txt=''; Y_RMS_Tp_2_txt='';
if WriteSimple_RMS:
    Y_RMS_Tp_1_txt = '& $%.3f \\pm %.3f$'%(Y_wRMS_Tp[5][2], Y_wRMS_Tp[5][3] )
    Y_RMS_Tp_1_dat = '%7.3f %7.3f'%(Y_wRMS_Tp[5][2], Y_wRMS_Tp[5][3] )
if Write_wRMS_2nd:
    Y_wRMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(Y_wRMS_Tp_2[5][0], Y_wRMS_Tp_2[5][1])
if WriteSimple_RMS_2nd:
    Y_RMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(Y_wRMS_Tp_2[5][2], Y_wRMS_Tp_2[5][3] )

textfile_1.write('$Y$     & Template     & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    Y_wRMS_Tp[3][1],
    Y_wRMS_Tp[0][0], Y_wRMS_Tp[0][1], Y_chi2dof_Tp_1,
    Y_wRMS_Tp_2[0][0], Y_wRMS_Tp_2[0][1], Y_chi2dof_Tp_2,
    Y_wRMS_Tp[5][0], Y_wRMS_Tp[5][1], Y_wRMS_Tp_2_txt,
    Y_RMS_Tp_1_txt, Y_RMS_Tp_2_txt ))

textfile_2.write('Y        Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    Y_wRMS_Tp[3][1],
    Y_wRMS_Tp[0][0], Y_wRMS_Tp[0][1], Y_chi2dof_Tp_1,
    Y_wRMS_Tp_2[0][0], Y_wRMS_Tp_2[0][1], Y_chi2dof_Tp_2,
    Y_wRMS_Tp[5][0], Y_wRMS_Tp[5][1], Y_wRMS_Tp_2_txt,
    Y_RMS_Tp_1_dat, Y_RMS_Tp_2_txt ))

textfile_5.write('Y        Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    Y_wRMS_Tp[3][1],
    Y_wRMS_Tp[0][0], Y_wRMS_Tp[0][1], Y_chi2dof_Tp_1,
    Y_wRMS_Tp_2[0][0], Y_wRMS_Tp_2[0][1], Y_chi2dof_Tp_2,
    Y_wRMS_Tp[5][0], Y_wRMS_Tp[5][1], Y_wRMS_Tp_2_txt,
    Y_RMS_Tp_1_dat, Y_RMS_Tp_2_txt ))

#-------------------
# reset
Y_RMS_GP_1_txt=''; Y_wRMS_GP_2_txt=''; Y_RMS_GP_2_txt='';
if WriteSimple_RMS:
    Y_RMS_GP_1_txt = '& $%.3f \\pm %.3f$'%(Y_wRMS_GP[5][2], Y_wRMS_GP[5][3] )
    Y_RMS_GP_1_dat = '%7.3f %7.3f'%(Y_wRMS_GP[5][2], Y_wRMS_GP[5][3] )
if Write_wRMS_2nd:
    Y_wRMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(Y_wRMS_GP_2[5][0], Y_wRMS_GP_2[5][1])
if WriteSimple_RMS_2nd:
    Y_RMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(Y_wRMS_GP_2[5][2], Y_wRMS_GP_2[5][3] )

textfile_1.write('$Y$     & GP (NIR max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    Y_wRMS_GP[3][1],
    Y_wRMS_GP[0][0], Y_wRMS_GP[0][1], Y_chi2dof_GP_1,
    Y_wRMS_GP_2[0][0], Y_wRMS_GP_2[0][1], Y_chi2dof_GP_2,
    Y_wRMS_GP[5][0], Y_wRMS_GP[5][1], Y_wRMS_GP_2_txt,
    Y_RMS_GP_1_txt, Y_RMS_GP_2_txt ))

textfile_2.write('Y        GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    Y_wRMS_GP[3][1],
    Y_wRMS_GP[0][0], Y_wRMS_GP[0][1], Y_chi2dof_GP_1,
    Y_wRMS_GP_2[0][0], Y_wRMS_GP_2[0][1], Y_chi2dof_GP_2,
    Y_wRMS_GP[5][0], Y_wRMS_GP[5][1], Y_wRMS_GP_2_txt,
    Y_RMS_GP_1_dat, Y_RMS_GP_2_txt ))

textfile_3.write('Y        GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    Y_wRMS_GP[3][1],
    Y_wRMS_GP[0][0], Y_wRMS_GP[0][1], Y_chi2dof_GP_1,
    Y_wRMS_GP_2[0][0], Y_wRMS_GP_2[0][1], Y_chi2dof_GP_2,
    Y_wRMS_GP[5][0], Y_wRMS_GP[5][1], Y_wRMS_GP_2_txt,
    Y_RMS_GP_1_dat, Y_RMS_GP_2_txt ))

#-------------------
# reset
Y_RMS_GP_tBmax_1_txt=''; Y_wRMS_GP_tBmax_2_txt=''; Y_RMS_GP_tBmax_2_txt='';
if WriteSimple_RMS:
    Y_RMS_GP_tBmax_1_txt = '& $%.3f \\pm %.3f$'%(Y_wRMS_GP_tBmax[5][2], Y_wRMS_GP_tBmax[5][3] )
    Y_RMS_GP_tBmax_1_dat = '%7.3f %7.3f'%(Y_wRMS_GP_tBmax[5][2], Y_wRMS_GP_tBmax[5][3] )
if Write_wRMS_2nd:
    Y_wRMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(Y_wRMS_GP_tBmax_2[5][0], Y_wRMS_GP_tBmax_2[5][1])
if WriteSimple_RMS_2nd:
    Y_RMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(Y_wRMS_GP_tBmax_2[5][2], Y_wRMS_GP_tBmax_2[5][3] )

textfile_1.write('$Y$     & GP ($B$ max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    Y_wRMS_GP_tBmax[3][1],
    Y_wRMS_GP_tBmax[0][0], Y_wRMS_GP_tBmax[0][1], Y_chi2dof_GP_tBmax_1,
    Y_wRMS_GP_tBmax_2[0][0], Y_wRMS_GP_tBmax_2[0][1], Y_chi2dof_GP_tBmax_2,
    Y_wRMS_GP_tBmax[5][0], Y_wRMS_GP_tBmax[5][1], Y_wRMS_GP_tBmax_2_txt,
    Y_RMS_GP_tBmax_1_txt, Y_RMS_GP_tBmax_2_txt ))

textfile_2.write('Y        GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    Y_wRMS_GP_tBmax[3][1],
    Y_wRMS_GP_tBmax[0][0], Y_wRMS_GP_tBmax[0][1], Y_chi2dof_GP_tBmax_1,
    Y_wRMS_GP_tBmax_2[0][0], Y_wRMS_GP_tBmax_2[0][1], Y_chi2dof_GP_tBmax_2,
    Y_wRMS_GP_tBmax[5][0], Y_wRMS_GP_tBmax[5][1], Y_wRMS_GP_tBmax_2_txt,
    Y_RMS_GP_tBmax_1_dat, Y_RMS_GP_tBmax_2_txt ))

textfile_4.write('Y        GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    Y_wRMS_GP_tBmax[3][1],
    Y_wRMS_GP_tBmax[0][0], Y_wRMS_GP_tBmax[0][1], Y_chi2dof_GP_tBmax_1,
    Y_wRMS_GP_tBmax_2[0][0], Y_wRMS_GP_tBmax_2[0][1], Y_chi2dof_GP_tBmax_2,
    Y_wRMS_GP_tBmax[5][0], Y_wRMS_GP_tBmax[5][1], Y_wRMS_GP_tBmax_2_txt,
    Y_RMS_GP_tBmax_1_dat, Y_RMS_GP_tBmax_2_txt ))

#-------------------
textfile_1.write('%------------------------------------- \n')
# reset
J_RMS_Tp_1_txt=''; J_wRMS_Tp_2_txt=''; J_RMS_Tp_2_txt='';
if WriteSimple_RMS:
    J_RMS_Tp_1_txt = '& $%.3f \\pm %.3f$'%(J_wRMS_Tp[5][2], J_wRMS_Tp[5][3] )
    J_RMS_Tp_1_dat = '%7.3f %7.3f'%(J_wRMS_Tp[5][2], J_wRMS_Tp[5][3] )
if Write_wRMS_2nd:
    J_wRMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(J_wRMS_Tp_2[5][0], J_wRMS_Tp_2[5][1])
if WriteSimple_RMS_2nd:
    J_RMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(J_wRMS_Tp_2[5][2], J_wRMS_Tp_2[5][3] )

textfile_1.write('$J$     & Template     & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    J_wRMS_Tp[3][1],
    J_wRMS_Tp[0][0], J_wRMS_Tp[0][1], J_chi2dof_Tp_1,
    J_wRMS_Tp_2[0][0], J_wRMS_Tp_2[0][1], J_chi2dof_Tp_2,
    J_wRMS_Tp[5][0], J_wRMS_Tp[5][1], J_wRMS_Tp_2_txt,
    J_RMS_Tp_1_txt, J_RMS_Tp_2_txt ))

textfile_2.write('J        Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    J_wRMS_Tp[3][1],
    J_wRMS_Tp[0][0], J_wRMS_Tp[0][1], J_chi2dof_Tp_1,
    J_wRMS_Tp_2[0][0], J_wRMS_Tp_2[0][1], J_chi2dof_Tp_2,
    J_wRMS_Tp[5][0], J_wRMS_Tp[5][1], J_wRMS_Tp_2_txt,
    J_RMS_Tp_1_dat, J_RMS_Tp_2_txt ))

textfile_5.write('J        Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    J_wRMS_Tp[3][1],
    J_wRMS_Tp[0][0], J_wRMS_Tp[0][1], J_chi2dof_Tp_1,
    J_wRMS_Tp_2[0][0], J_wRMS_Tp_2[0][1], J_chi2dof_Tp_2,
    J_wRMS_Tp[5][0], J_wRMS_Tp[5][1], J_wRMS_Tp_2_txt,
    J_RMS_Tp_1_dat, J_RMS_Tp_2_txt ))

#-------------------
# reset
J_RMS_GP_1_txt=''; J_wRMS_GP_2_txt=''; J_RMS_GP_2_txt='';
if WriteSimple_RMS:
    J_RMS_GP_1_txt = '& $%.3f \\pm %.3f$'%(J_wRMS_GP[5][2], J_wRMS_GP[5][3] )
    J_RMS_GP_1_dat = '%7.3f %7.3f'%(J_wRMS_GP[5][2], J_wRMS_GP[5][3] )
if Write_wRMS_2nd:
    J_wRMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(J_wRMS_GP_2[5][0], J_wRMS_GP_2[5][1])
if WriteSimple_RMS_2nd:
    J_RMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(J_wRMS_GP_2[5][2], J_wRMS_GP_2[5][3] )

textfile_1.write('$J$     & GP (NIR max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    J_wRMS_GP[3][1],
    J_wRMS_GP[0][0], J_wRMS_GP[0][1], J_chi2dof_GP_1,
    J_wRMS_GP_2[0][0], J_wRMS_GP_2[0][1], J_chi2dof_GP_2,
    J_wRMS_GP[5][0], J_wRMS_GP[5][1], J_wRMS_GP_2_txt,
    J_RMS_GP_1_txt, J_RMS_GP_2_txt ))

textfile_2.write('J        GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    J_wRMS_GP[3][1],
    J_wRMS_GP[0][0], J_wRMS_GP[0][1], J_chi2dof_GP_1,
    J_wRMS_GP_2[0][0], J_wRMS_GP_2[0][1], J_chi2dof_GP_2,
    J_wRMS_GP[5][0], J_wRMS_GP[5][1], J_wRMS_GP_2_txt,
    J_RMS_GP_1_dat, J_RMS_GP_2_txt ))

textfile_3.write('J        GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    J_wRMS_GP[3][1],
    J_wRMS_GP[0][0], J_wRMS_GP[0][1], J_chi2dof_GP_1,
    J_wRMS_GP_2[0][0], J_wRMS_GP_2[0][1], J_chi2dof_GP_2,
    J_wRMS_GP[5][0], J_wRMS_GP[5][1], J_wRMS_GP_2_txt,
    J_RMS_GP_1_dat, J_RMS_GP_2_txt ))

#-------------------
# reset
J_RMS_GP_tBmax_1_txt=''; J_wRMS_GP_tBmax_2_txt=''; J_RMS_GP_tBmax_2_txt='';
if WriteSimple_RMS:
    J_RMS_GP_tBmax_1_txt = '& $%.3f \\pm %.3f$'%(J_wRMS_GP_tBmax[5][2], J_wRMS_GP_tBmax[5][3] )
    J_RMS_GP_tBmax_1_dat = '%7.3f %7.3f'%(J_wRMS_GP_tBmax[5][2], J_wRMS_GP_tBmax[5][3] )
if Write_wRMS_2nd:
    J_wRMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(J_wRMS_GP_tBmax_2[5][0], J_wRMS_GP_tBmax_2[5][1])
if WriteSimple_RMS_2nd:
    J_RMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(J_wRMS_GP_tBmax_2[5][2], J_wRMS_GP_tBmax_2[5][3] )

textfile_1.write('$J$     & GP ($B$ max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    J_wRMS_GP_tBmax[3][1],
    J_wRMS_GP_tBmax[0][0], J_wRMS_GP_tBmax[0][1], J_chi2dof_GP_tBmax_1,
    J_wRMS_GP_tBmax_2[0][0], J_wRMS_GP_tBmax_2[0][1], J_chi2dof_GP_tBmax_2,
    J_wRMS_GP_tBmax[5][0], J_wRMS_GP_tBmax[5][1], J_wRMS_GP_tBmax_2_txt,
    J_RMS_GP_tBmax_1_txt, J_RMS_GP_tBmax_2_txt ))

textfile_2.write('J        GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    J_wRMS_GP_tBmax[3][1],
    J_wRMS_GP_tBmax[0][0], J_wRMS_GP_tBmax[0][1], J_chi2dof_GP_tBmax_1,
    J_wRMS_GP_tBmax_2[0][0], J_wRMS_GP_tBmax_2[0][1], J_chi2dof_GP_tBmax_2,
    J_wRMS_GP_tBmax[5][0], J_wRMS_GP_tBmax[5][1], J_wRMS_GP_tBmax_2_txt,
    J_RMS_GP_tBmax_1_dat, J_RMS_GP_tBmax_2_txt ))

textfile_4.write('J        GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    J_wRMS_GP_tBmax[3][1],
    J_wRMS_GP_tBmax[0][0], J_wRMS_GP_tBmax[0][1], J_chi2dof_GP_tBmax_1,
    J_wRMS_GP_tBmax_2[0][0], J_wRMS_GP_tBmax_2[0][1], J_chi2dof_GP_tBmax_2,
    J_wRMS_GP_tBmax[5][0], J_wRMS_GP_tBmax[5][1], J_wRMS_GP_tBmax_2_txt,
    J_RMS_GP_tBmax_1_dat, J_RMS_GP_tBmax_2_txt ))

#-------------------
textfile_1.write('%------------------------------------- \n')
# reset
H_RMS_Tp_1_txt=''; H_wRMS_Tp_2_txt=''; H_RMS_Tp_2_txt='';
if WriteSimple_RMS:
    H_RMS_Tp_1_txt = '& $%.3f \\pm %.3f$'%(H_wRMS_Tp[5][2], H_wRMS_Tp[5][3] )
    H_RMS_Tp_1_dat = '%7.3f %7.3f'%(H_wRMS_Tp[5][2], H_wRMS_Tp[5][3] )
if Write_wRMS_2nd:
    H_wRMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(H_wRMS_Tp_2[5][0], H_wRMS_Tp_2[5][1])
if WriteSimple_RMS_2nd:
    H_RMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(H_wRMS_Tp_2[5][2], H_wRMS_Tp_2[5][3] )

textfile_1.write('$H$     & Template     & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    H_wRMS_Tp[3][1],
    H_wRMS_Tp[0][0], H_wRMS_Tp[0][1], H_chi2dof_Tp_1,
    H_wRMS_Tp_2[0][0], H_wRMS_Tp_2[0][1], H_chi2dof_Tp_2,
    H_wRMS_Tp[5][0], H_wRMS_Tp[5][1], H_wRMS_Tp_2_txt,
    H_RMS_Tp_1_txt, H_RMS_Tp_2_txt ))

textfile_2.write('H        Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    H_wRMS_Tp[3][1],
    H_wRMS_Tp[0][0], H_wRMS_Tp[0][1], H_chi2dof_Tp_1,
    H_wRMS_Tp_2[0][0], H_wRMS_Tp_2[0][1], H_chi2dof_Tp_2,
    H_wRMS_Tp[5][0], H_wRMS_Tp[5][1], H_wRMS_Tp_2_txt,
    H_RMS_Tp_1_dat, H_RMS_Tp_2_txt ))

textfile_5.write('H        Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    H_wRMS_Tp[3][1],
    H_wRMS_Tp[0][0], H_wRMS_Tp[0][1], H_chi2dof_Tp_1,
    H_wRMS_Tp_2[0][0], H_wRMS_Tp_2[0][1], H_chi2dof_Tp_2,
    H_wRMS_Tp[5][0], H_wRMS_Tp[5][1], H_wRMS_Tp_2_txt,
    H_RMS_Tp_1_dat, H_RMS_Tp_2_txt ))

#-------------------
# reset
H_RMS_GP_1_txt=''; H_wRMS_GP_2_txt=''; H_RMS_GP_2_txt='';
if WriteSimple_RMS:
    H_RMS_GP_1_txt = '& $%.3f \\pm %.3f$'%(H_wRMS_GP[5][2], H_wRMS_GP[5][3] )
    H_RMS_GP_1_dat = '%7.3f %7.3f'%(H_wRMS_GP[5][2], H_wRMS_GP[5][3] )
if Write_wRMS_2nd:
    H_wRMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(H_wRMS_GP_2[5][0], H_wRMS_GP_2[5][1])
if WriteSimple_RMS_2nd:
    H_RMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(H_wRMS_GP_2[5][2], H_wRMS_GP_2[5][3] )

textfile_1.write('$H$     & GP (NIR max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    H_wRMS_GP[3][1],
    H_wRMS_GP[0][0], H_wRMS_GP[0][1], H_chi2dof_GP_1,
    H_wRMS_GP_2[0][0], H_wRMS_GP_2[0][1], H_chi2dof_GP_2,
    H_wRMS_GP[5][0], H_wRMS_GP[5][1], H_wRMS_GP_2_txt,
    H_RMS_GP_1_txt, H_RMS_GP_2_txt ))

textfile_2.write('H        GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    H_wRMS_GP[3][1],
    H_wRMS_GP[0][0], H_wRMS_GP[0][1], H_chi2dof_GP_1,
    H_wRMS_GP_2[0][0], H_wRMS_GP_2[0][1], H_chi2dof_GP_2,
    H_wRMS_GP[5][0], H_wRMS_GP[5][1], H_wRMS_GP_2_txt,
    H_RMS_GP_1_dat, H_RMS_GP_2_txt ))

textfile_3.write('H        GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    H_wRMS_GP[3][1],
    H_wRMS_GP[0][0], H_wRMS_GP[0][1], H_chi2dof_GP_1,
    H_wRMS_GP_2[0][0], H_wRMS_GP_2[0][1], H_chi2dof_GP_2,
    H_wRMS_GP[5][0], H_wRMS_GP[5][1], H_wRMS_GP_2_txt,
    H_RMS_GP_1_dat, H_RMS_GP_2_txt ))

#-------------------
# reset
H_RMS_GP_tBmax_1_txt=''; H_wRMS_GP_tBmax_2_txt=''; H_RMS_GP_tBmax_2_txt='';
if WriteSimple_RMS:
    H_RMS_GP_tBmax_1_txt = '& $%.3f \\pm %.3f$'%(H_wRMS_GP_tBmax[5][2], H_wRMS_GP_tBmax[5][3] )
    H_RMS_GP_tBmax_1_dat = '%7.3f %7.3f'%(H_wRMS_GP_tBmax[5][2], H_wRMS_GP_tBmax[5][3] )
if Write_wRMS_2nd:
    H_wRMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(H_wRMS_GP_tBmax_2[5][0], H_wRMS_GP_tBmax_2[5][1])
if WriteSimple_RMS_2nd:
    H_RMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(H_wRMS_GP_tBmax_2[5][2], H_wRMS_GP_tBmax_2[5][3] )

textfile_1.write('$H$     & GP ($B$ max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    H_wRMS_GP_tBmax[3][1],
    H_wRMS_GP_tBmax[0][0], H_wRMS_GP_tBmax[0][1], H_chi2dof_GP_tBmax_1,
    H_wRMS_GP_tBmax_2[0][0], H_wRMS_GP_tBmax_2[0][1], H_chi2dof_GP_tBmax_2,
    H_wRMS_GP_tBmax[5][0], H_wRMS_GP_tBmax[5][1], H_wRMS_GP_tBmax_2_txt,
    H_RMS_GP_tBmax_1_txt, H_RMS_GP_tBmax_2_txt ))

textfile_2.write('H        GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    H_wRMS_GP_tBmax[3][1],
    H_wRMS_GP_tBmax[0][0], H_wRMS_GP_tBmax[0][1], H_chi2dof_GP_tBmax_1,
    H_wRMS_GP_tBmax_2[0][0], H_wRMS_GP_tBmax_2[0][1], H_chi2dof_GP_tBmax_2,
    H_wRMS_GP_tBmax[5][0], H_wRMS_GP_tBmax[5][1], H_wRMS_GP_tBmax_2_txt,
    H_RMS_GP_tBmax_1_dat, H_RMS_GP_tBmax_2_txt ))

textfile_4.write('H        GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    H_wRMS_GP_tBmax[3][1],
    H_wRMS_GP_tBmax[0][0], H_wRMS_GP_tBmax[0][1], H_chi2dof_GP_tBmax_1,
    H_wRMS_GP_tBmax_2[0][0], H_wRMS_GP_tBmax_2[0][1], H_chi2dof_GP_tBmax_2,
    H_wRMS_GP_tBmax[5][0], H_wRMS_GP_tBmax[5][1], H_wRMS_GP_tBmax_2_txt,
    H_RMS_GP_tBmax_1_dat, H_RMS_GP_tBmax_2_txt ))

#-------------------
textfile_1.write('%------------------------------------- \n')
# reset
K_RMS_Tp_1_txt=''; K_wRMS_Tp_2_txt=''; K_RMS_Tp_2_txt='';
if WriteSimple_RMS:
    K_RMS_Tp_1_txt = '& $%.3f \\pm %.3f$'%(K_wRMS_Tp[5][2], K_wRMS_Tp[5][3] )
    K_RMS_Tp_1_dat = '%7.3f %7.3f'%(K_wRMS_Tp[5][2], K_wRMS_Tp[5][3] )
if Write_wRMS_2nd:
    K_wRMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(K_wRMS_Tp_2[5][0], K_wRMS_Tp_2[5][1])
if WriteSimple_RMS_2nd:
    K_RMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(K_wRMS_Tp_2[5][2], K_wRMS_Tp_2[5][3] )

textfile_1.write('$K$     & Template     & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    K_wRMS_Tp[3][1],
    K_wRMS_Tp[0][0], K_wRMS_Tp[0][1], K_chi2dof_Tp_1,
    K_wRMS_Tp_2[0][0], K_wRMS_Tp_2[0][1], K_chi2dof_Tp_2,
    K_wRMS_Tp[5][0], K_wRMS_Tp[5][1], K_wRMS_Tp_2_txt,
    K_RMS_Tp_1_txt, K_RMS_Tp_2_txt ))

textfile_2.write('K        Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    K_wRMS_Tp[3][1],
    K_wRMS_Tp[0][0], K_wRMS_Tp[0][1], K_chi2dof_Tp_1,
    K_wRMS_Tp_2[0][0], K_wRMS_Tp_2[0][1], K_chi2dof_Tp_2,
    K_wRMS_Tp[5][0], K_wRMS_Tp[5][1], K_wRMS_Tp_2_txt,
    K_RMS_Tp_1_dat, K_RMS_Tp_2_txt ))

textfile_5.write('K        Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    K_wRMS_Tp[3][1],
    K_wRMS_Tp[0][0], K_wRMS_Tp[0][1], K_chi2dof_Tp_1,
    K_wRMS_Tp_2[0][0], K_wRMS_Tp_2[0][1], K_chi2dof_Tp_2,
    K_wRMS_Tp[5][0], K_wRMS_Tp[5][1], K_wRMS_Tp_2_txt,
    K_RMS_Tp_1_dat, K_RMS_Tp_2_txt ))

#-------------------
# reset
K_RMS_GP_1_txt=''; K_wRMS_GP_2_txt=''; K_RMS_GP_2_txt='';
if WriteSimple_RMS:
    K_RMS_GP_1_txt = '& $%.3f \\pm %.3f$'%(K_wRMS_GP[5][2], K_wRMS_GP[5][3] )
    K_RMS_GP_1_dat = '%7.3f %7.3f'%(K_wRMS_GP[5][2], K_wRMS_GP[5][3] )
if Write_wRMS_2nd:
    K_wRMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(K_wRMS_GP_2[5][0], K_wRMS_GP_2[5][1])
if WriteSimple_RMS_2nd:
    K_RMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(K_wRMS_GP_2[5][2], K_wRMS_GP_2[5][3] )

textfile_1.write('$K$     & GP (NIR max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    K_wRMS_GP[3][1],
    K_wRMS_GP[0][0], K_wRMS_GP[0][1], K_chi2dof_GP_1,
    K_wRMS_GP_2[0][0], K_wRMS_GP_2[0][1], K_chi2dof_GP_2,
    K_wRMS_GP[5][0], K_wRMS_GP[5][1], K_wRMS_GP_2_txt,
    K_RMS_GP_1_txt, K_RMS_GP_2_txt ))

textfile_2.write('K        GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    K_wRMS_GP[3][1],
    K_wRMS_GP[0][0], K_wRMS_GP[0][1], K_chi2dof_GP_1,
    K_wRMS_GP_2[0][0], K_wRMS_GP_2[0][1], K_chi2dof_GP_2,
    K_wRMS_GP[5][0], K_wRMS_GP[5][1], K_wRMS_GP_2_txt,
    K_RMS_GP_1_dat, K_RMS_GP_2_txt ))

textfile_3.write('K        GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    K_wRMS_GP[3][1],
    K_wRMS_GP[0][0], K_wRMS_GP[0][1], K_chi2dof_GP_1,
    K_wRMS_GP_2[0][0], K_wRMS_GP_2[0][1], K_chi2dof_GP_2,
    K_wRMS_GP[5][0], K_wRMS_GP[5][1], K_wRMS_GP_2_txt,
    K_RMS_GP_1_dat, K_RMS_GP_2_txt ))

#-------------------
# reset
K_RMS_GP_tBmax_1_txt=''; K_wRMS_GP_tBmax_2_txt=''; K_RMS_GP_tBmax_2_txt='';
if WriteSimple_RMS:
    K_RMS_GP_tBmax_1_txt = '& $%.3f \\pm %.3f$'%(K_wRMS_GP_tBmax[5][2], K_wRMS_GP_tBmax[5][3] )
    K_RMS_GP_tBmax_1_dat = '%7.3f %7.3f'%(K_wRMS_GP_tBmax[5][2], K_wRMS_GP_tBmax[5][3] )
if Write_wRMS_2nd:
    K_wRMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(K_wRMS_GP_tBmax_2[5][0], K_wRMS_GP_tBmax_2[5][1])
if WriteSimple_RMS_2nd:
    K_RMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(K_wRMS_GP_tBmax_2[5][2], K_wRMS_GP_tBmax_2[5][3] )

textfile_1.write('$K$     & GP ($B$ max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    K_wRMS_GP_tBmax[3][1],
    K_wRMS_GP_tBmax[0][0], K_wRMS_GP_tBmax[0][1], K_chi2dof_GP_tBmax_1,
    K_wRMS_GP_tBmax_2[0][0], K_wRMS_GP_tBmax_2[0][1], K_chi2dof_GP_tBmax_2,
    K_wRMS_GP_tBmax[5][0], K_wRMS_GP_tBmax[5][1], K_wRMS_GP_tBmax_2_txt,
    K_RMS_GP_tBmax_1_txt, K_RMS_GP_tBmax_2_txt ))

textfile_2.write('K        GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    K_wRMS_GP_tBmax[3][1],
    K_wRMS_GP_tBmax[0][0], K_wRMS_GP_tBmax[0][1], K_chi2dof_GP_tBmax_1,
    K_wRMS_GP_tBmax_2[0][0], K_wRMS_GP_tBmax_2[0][1], K_chi2dof_GP_tBmax_2,
    K_wRMS_GP_tBmax[5][0], K_wRMS_GP_tBmax[5][1], K_wRMS_GP_tBmax_2_txt,
    K_RMS_GP_tBmax_1_dat, K_RMS_GP_tBmax_2_txt ))

textfile_4.write('K        GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    K_wRMS_GP_tBmax[3][1],
    K_wRMS_GP_tBmax[0][0], K_wRMS_GP_tBmax[0][1], K_chi2dof_GP_tBmax_1,
    K_wRMS_GP_tBmax_2[0][0], K_wRMS_GP_tBmax_2[0][1], K_chi2dof_GP_tBmax_2,
    K_wRMS_GP_tBmax[5][0], K_wRMS_GP_tBmax[5][1], K_wRMS_GP_tBmax_2_txt,
    K_RMS_GP_tBmax_1_dat, K_RMS_GP_tBmax_2_txt ))

#-------------------
textfile_1.write('%------------------------------------- \n')
# reset
AnyNIR_RMS_Tp_1_txt=''; AnyNIR_wRMS_Tp_2_txt=''; AnyNIR_RMS_Tp_2_txt='';
if WriteSimple_RMS:
    AnyNIR_RMS_Tp_1_txt = '& $%.3f \\pm %.3f$'%(AnyNIR_wRMS_Tp[5][2], AnyNIR_wRMS_Tp[5][3] )
    AnyNIR_RMS_Tp_1_dat = '%7.3f %7.3f'%(AnyNIR_wRMS_Tp[5][2], AnyNIR_wRMS_Tp[5][3] )
if Write_wRMS_2nd:
    AnyNIR_wRMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(AnyNIR_wRMS_Tp_2[5][0], AnyNIR_wRMS_Tp_2[5][1])
if WriteSimple_RMS_2nd:
    AnyNIR_RMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(AnyNIR_wRMS_Tp_2[5][2], AnyNIR_wRMS_Tp_2[5][3] )

textfile_1.write('any $YJHK_s$     & Template     & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    AnyNIR_wRMS_Tp[3][1],
    AnyNIR_wRMS_Tp[0][0], AnyNIR_wRMS_Tp[0][1], AnyNIR_chi2dof_Tp_1,
    AnyNIR_wRMS_Tp_2[0][0], AnyNIR_wRMS_Tp_2[0][1], AnyNIR_chi2dof_Tp_2,
    AnyNIR_wRMS_Tp[5][0], AnyNIR_wRMS_Tp[5][1], AnyNIR_wRMS_Tp_2_txt,
    AnyNIR_RMS_Tp_1_txt, AnyNIR_RMS_Tp_2_txt ))

textfile_2.write('AnyNIR   Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    AnyNIR_wRMS_Tp[3][1],
    AnyNIR_wRMS_Tp[0][0], AnyNIR_wRMS_Tp[0][1], AnyNIR_chi2dof_Tp_1,
    AnyNIR_wRMS_Tp_2[0][0], AnyNIR_wRMS_Tp_2[0][1], AnyNIR_chi2dof_Tp_2,
    AnyNIR_wRMS_Tp[5][0], AnyNIR_wRMS_Tp[5][1], AnyNIR_wRMS_Tp_2_txt,
    AnyNIR_RMS_Tp_1_dat, AnyNIR_RMS_Tp_2_txt ))

textfile_5.write('AnyNIR   Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    AnyNIR_wRMS_Tp[3][1],
    AnyNIR_wRMS_Tp[0][0], AnyNIR_wRMS_Tp[0][1], AnyNIR_chi2dof_Tp_1,
    AnyNIR_wRMS_Tp_2[0][0], AnyNIR_wRMS_Tp_2[0][1], AnyNIR_chi2dof_Tp_2,
    AnyNIR_wRMS_Tp[5][0], AnyNIR_wRMS_Tp[5][1], AnyNIR_wRMS_Tp_2_txt,
    AnyNIR_RMS_Tp_1_dat, AnyNIR_RMS_Tp_2_txt ))

#-------------------
# reset
AnyNIR_RMS_GP_1_txt=''; AnyNIR_wRMS_GP_2_txt=''; AnyNIR_RMS_GP_2_txt='';
if WriteSimple_RMS:
    AnyNIR_RMS_GP_1_txt = '& $%.3f \\pm %.3f$'%(AnyNIR_wRMS_GP[5][2], AnyNIR_wRMS_GP[5][3] )
    AnyNIR_RMS_GP_1_dat = '%7.3f %7.3f'%(AnyNIR_wRMS_GP[5][2], AnyNIR_wRMS_GP[5][3] )
if Write_wRMS_2nd:
    AnyNIR_wRMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(AnyNIR_wRMS_GP_2[5][0], AnyNIR_wRMS_GP_2[5][1])
if WriteSimple_RMS_2nd:
    AnyNIR_RMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(AnyNIR_wRMS_GP_2[5][2], AnyNIR_wRMS_GP_2[5][3] )

textfile_1.write('any $YJHK_s$     & GP (NIR max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    AnyNIR_wRMS_GP[3][1],
    AnyNIR_wRMS_GP[0][0], AnyNIR_wRMS_GP[0][1], AnyNIR_chi2dof_GP_1,
    AnyNIR_wRMS_GP_2[0][0], AnyNIR_wRMS_GP_2[0][1], AnyNIR_chi2dof_GP_2,
    AnyNIR_wRMS_GP[5][0], AnyNIR_wRMS_GP[5][1], AnyNIR_wRMS_GP_2_txt,
    AnyNIR_RMS_GP_1_txt, AnyNIR_RMS_GP_2_txt ))

textfile_2.write('AnyNIR   GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    AnyNIR_wRMS_GP[3][1],
    AnyNIR_wRMS_GP[0][0], AnyNIR_wRMS_GP[0][1], AnyNIR_chi2dof_GP_1,
    AnyNIR_wRMS_GP_2[0][0], AnyNIR_wRMS_GP_2[0][1], AnyNIR_chi2dof_GP_2,
    AnyNIR_wRMS_GP[5][0], AnyNIR_wRMS_GP[5][1], AnyNIR_wRMS_GP_2_txt,
    AnyNIR_RMS_GP_1_dat, AnyNIR_RMS_GP_2_txt ))

textfile_3.write('AnyNIR   GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    AnyNIR_wRMS_GP[3][1],
    AnyNIR_wRMS_GP[0][0], AnyNIR_wRMS_GP[0][1], AnyNIR_chi2dof_GP_1,
    AnyNIR_wRMS_GP_2[0][0], AnyNIR_wRMS_GP_2[0][1], AnyNIR_chi2dof_GP_2,
    AnyNIR_wRMS_GP[5][0], AnyNIR_wRMS_GP[5][1], AnyNIR_wRMS_GP_2_txt,
    AnyNIR_RMS_GP_1_dat, AnyNIR_RMS_GP_2_txt ))

#-------------------
# reset
AnyNIR_RMS_GP_tBmax_1_txt=''; AnyNIR_wRMS_GP_tBmax_2_txt=''; AnyNIR_RMS_GP_tBmax_2_txt='';
if WriteSimple_RMS:
    AnyNIR_RMS_GP_tBmax_1_txt = '& $%.3f \\pm %.3f$'%(AnyNIR_wRMS_GP_tBmax[5][2], AnyNIR_wRMS_GP_tBmax[5][3] )
    AnyNIR_RMS_GP_tBmax_1_dat = '%7.3f %7.3f'%(AnyNIR_wRMS_GP_tBmax[5][2], AnyNIR_wRMS_GP_tBmax[5][3] )
if Write_wRMS_2nd:
    AnyNIR_wRMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(AnyNIR_wRMS_GP_tBmax_2[5][0], AnyNIR_wRMS_GP_tBmax_2[5][1])
if WriteSimple_RMS_2nd:
    AnyNIR_RMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(AnyNIR_wRMS_GP_tBmax_2[5][2], AnyNIR_wRMS_GP_tBmax_2[5][3] )

textfile_1.write('any $YJHK_s$     & GP ($B$ max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    AnyNIR_wRMS_GP_tBmax[3][1],
    AnyNIR_wRMS_GP_tBmax[0][0], AnyNIR_wRMS_GP_tBmax[0][1], AnyNIR_chi2dof_GP_tBmax_1,
    AnyNIR_wRMS_GP_tBmax_2[0][0], AnyNIR_wRMS_GP_tBmax_2[0][1], AnyNIR_chi2dof_GP_tBmax_2,
    AnyNIR_wRMS_GP_tBmax[5][0], AnyNIR_wRMS_GP_tBmax[5][1], AnyNIR_wRMS_GP_tBmax_2_txt,
    AnyNIR_RMS_GP_tBmax_1_txt, AnyNIR_RMS_GP_tBmax_2_txt ))

textfile_2.write('AnyNIR   GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    AnyNIR_wRMS_GP_tBmax[3][1],
    AnyNIR_wRMS_GP_tBmax[0][0], AnyNIR_wRMS_GP_tBmax[0][1], AnyNIR_chi2dof_GP_tBmax_1,
    AnyNIR_wRMS_GP_tBmax_2[0][0], AnyNIR_wRMS_GP_tBmax_2[0][1], AnyNIR_chi2dof_GP_tBmax_2,
    AnyNIR_wRMS_GP_tBmax[5][0], AnyNIR_wRMS_GP_tBmax[5][1], AnyNIR_wRMS_GP_tBmax_2_txt,
    AnyNIR_RMS_GP_tBmax_1_dat, AnyNIR_RMS_GP_tBmax_2_txt ))

textfile_4.write('AnyNIR   GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    AnyNIR_wRMS_GP_tBmax[3][1],
    AnyNIR_wRMS_GP_tBmax[0][0], AnyNIR_wRMS_GP_tBmax[0][1], AnyNIR_chi2dof_GP_tBmax_1,
    AnyNIR_wRMS_GP_tBmax_2[0][0], AnyNIR_wRMS_GP_tBmax_2[0][1], AnyNIR_chi2dof_GP_tBmax_2,
    AnyNIR_wRMS_GP_tBmax[5][0], AnyNIR_wRMS_GP_tBmax[5][1], AnyNIR_wRMS_GP_tBmax_2_txt,
    AnyNIR_RMS_GP_tBmax_1_dat, AnyNIR_RMS_GP_tBmax_2_txt ))

#-------------------
textfile_1.write('%------------------------------------- \n')
# reset
JH_RMS_Tp_1_txt=''; JH_wRMS_Tp_2_txt=''; JH_RMS_Tp_2_txt='';
if WriteSimple_RMS:
    JH_RMS_Tp_1_txt = '& $%.3f \\pm %.3f$'%(JH_wRMS_Tp[5][2], JH_wRMS_Tp[5][3] )
    JH_RMS_Tp_1_dat = '%7.3f %7.3f'%(JH_wRMS_Tp[5][2], JH_wRMS_Tp[5][3] )
if Write_wRMS_2nd:
    JH_wRMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(JH_wRMS_Tp_2[5][0], JH_wRMS_Tp_2[5][1])
if WriteSimple_RMS_2nd:
    JH_RMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(JH_wRMS_Tp_2[5][2], JH_wRMS_Tp_2[5][3] )

textfile_1.write('$JH$      & Template     & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    JH_wRMS_Tp[3][1],
    JH_wRMS_Tp[0][0], JH_wRMS_Tp[0][1], JH_chi2dof_Tp_1,
    JH_wRMS_Tp_2[0][0], JH_wRMS_Tp_2[0][1], JH_chi2dof_Tp_2,
    JH_wRMS_Tp[5][0], JH_wRMS_Tp[5][1], JH_wRMS_Tp_2_txt,
    JH_RMS_Tp_1_txt, JH_RMS_Tp_2_txt ))

textfile_2.write('JH       Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    JH_wRMS_Tp[3][1],
    JH_wRMS_Tp[0][0], JH_wRMS_Tp[0][1], JH_chi2dof_Tp_1,
    JH_wRMS_Tp_2[0][0], JH_wRMS_Tp_2[0][1], JH_chi2dof_Tp_2,
    JH_wRMS_Tp[5][0], JH_wRMS_Tp[5][1], JH_wRMS_Tp_2_txt,
    JH_RMS_Tp_1_dat, JH_RMS_Tp_2_txt ))

textfile_5.write('JH       Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    JH_wRMS_Tp[3][1],
    JH_wRMS_Tp[0][0], JH_wRMS_Tp[0][1], JH_chi2dof_Tp_1,
    JH_wRMS_Tp_2[0][0], JH_wRMS_Tp_2[0][1], JH_chi2dof_Tp_2,
    JH_wRMS_Tp[5][0], JH_wRMS_Tp[5][1], JH_wRMS_Tp_2_txt,
    JH_RMS_Tp_1_dat, JH_RMS_Tp_2_txt ))

#-------------------
# reset
JH_RMS_GP_1_txt=''; JH_wRMS_GP_2_txt=''; JH_RMS_GP_2_txt='';
if WriteSimple_RMS:
    JH_RMS_GP_1_txt = '& $%.3f \\pm %.3f$'%(JH_wRMS_GP[5][2], JH_wRMS_GP[5][3] )
    JH_RMS_GP_1_dat = '%7.3f %7.3f'%(JH_wRMS_GP[5][2], JH_wRMS_GP[5][3] )
if Write_wRMS_2nd:
    JH_wRMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(JH_wRMS_GP_2[5][0], JH_wRMS_GP_2[5][1])
if WriteSimple_RMS_2nd:
    JH_RMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(JH_wRMS_GP_2[5][2], JH_wRMS_GP_2[5][3] )

textfile_1.write('$JH$      & GP (NIR max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    JH_wRMS_GP[3][1],
    JH_wRMS_GP[0][0], JH_wRMS_GP[0][1], JH_chi2dof_GP_1,
    JH_wRMS_GP_2[0][0], JH_wRMS_GP_2[0][1], JH_chi2dof_GP_2,
    JH_wRMS_GP[5][0], JH_wRMS_GP[5][1], JH_wRMS_GP_2_txt,
    JH_RMS_GP_1_txt, JH_RMS_GP_2_txt ))

textfile_2.write('JH       GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    JH_wRMS_GP[3][1],
    JH_wRMS_GP[0][0], JH_wRMS_GP[0][1], JH_chi2dof_GP_1,
    JH_wRMS_GP_2[0][0], JH_wRMS_GP_2[0][1], JH_chi2dof_GP_2,
    JH_wRMS_GP[5][0], JH_wRMS_GP[5][1], JH_wRMS_GP_2_txt,
    JH_RMS_GP_1_dat, JH_RMS_GP_2_txt ))

textfile_3.write('JH       GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    JH_wRMS_GP[3][1],
    JH_wRMS_GP[0][0], JH_wRMS_GP[0][1], JH_chi2dof_GP_1,
    JH_wRMS_GP_2[0][0], JH_wRMS_GP_2[0][1], JH_chi2dof_GP_2,
    JH_wRMS_GP[5][0], JH_wRMS_GP[5][1], JH_wRMS_GP_2_txt,
    JH_RMS_GP_1_dat, JH_RMS_GP_2_txt ))

#-------------------
# reset
JH_RMS_GP_tBmax_1_txt=''; JH_wRMS_GP_tBmax_2_txt=''; JH_RMS_GP_tBmax_2_txt='';
if WriteSimple_RMS:
    JH_RMS_GP_tBmax_1_txt = '& $%.3f \\pm %.3f$'%(JH_wRMS_GP_tBmax[5][2], JH_wRMS_GP_tBmax[5][3] )
    JH_RMS_GP_tBmax_1_dat = '%7.3f %7.3f'%(JH_wRMS_GP_tBmax[5][2], JH_wRMS_GP_tBmax[5][3] )
if Write_wRMS_2nd:
    JH_wRMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(JH_wRMS_GP_tBmax_2[5][0], JH_wRMS_GP_tBmax_2[5][1])
if WriteSimple_RMS_2nd:
    JH_RMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(JH_wRMS_GP_tBmax_2[5][2], JH_wRMS_GP_tBmax_2[5][3] )

textfile_1.write('$JH$      & GP ($B$ max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    JH_wRMS_GP_tBmax[3][1],
    JH_wRMS_GP_tBmax[0][0], JH_wRMS_GP_tBmax[0][1], JH_chi2dof_GP_tBmax_1,
    JH_wRMS_GP_tBmax_2[0][0], JH_wRMS_GP_tBmax_2[0][1], JH_chi2dof_GP_tBmax_2,
    JH_wRMS_GP_tBmax[5][0], JH_wRMS_GP_tBmax[5][1], JH_wRMS_GP_tBmax_2_txt,
    JH_RMS_GP_tBmax_1_txt, JH_RMS_GP_tBmax_2_txt ))

textfile_2.write('JH       GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    JH_wRMS_GP_tBmax[3][1],
    JH_wRMS_GP_tBmax[0][0], JH_wRMS_GP_tBmax[0][1], JH_chi2dof_GP_tBmax_1,
    JH_wRMS_GP_tBmax_2[0][0], JH_wRMS_GP_tBmax_2[0][1], JH_chi2dof_GP_tBmax_2,
    JH_wRMS_GP_tBmax[5][0], JH_wRMS_GP_tBmax[5][1], JH_wRMS_GP_tBmax_2_txt,
    JH_RMS_GP_tBmax_1_dat, JH_RMS_GP_tBmax_2_txt ))

textfile_4.write('JH       GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    JH_wRMS_GP_tBmax[3][1],
    JH_wRMS_GP_tBmax[0][0], JH_wRMS_GP_tBmax[0][1], JH_chi2dof_GP_tBmax_1,
    JH_wRMS_GP_tBmax_2[0][0], JH_wRMS_GP_tBmax_2[0][1], JH_chi2dof_GP_tBmax_2,
    JH_wRMS_GP_tBmax[5][0], JH_wRMS_GP_tBmax[5][1], JH_wRMS_GP_tBmax_2_txt,
    JH_RMS_GP_tBmax_1_dat, JH_RMS_GP_tBmax_2_txt ))

#-------------------
textfile_1.write('%------------------------------------- \n')
# reset
YJH_RMS_Tp_1_txt=''; YJH_wRMS_Tp_2_txt=''; YJH_RMS_Tp_2_txt='';
if WriteSimple_RMS:
    YJH_RMS_Tp_1_txt = '& $%.3f \\pm %.3f$'%(YJH_wRMS_Tp[5][2], YJH_wRMS_Tp[5][3] )
    YJH_RMS_Tp_1_dat = '%7.3f %7.3f'%(YJH_wRMS_Tp[5][2], YJH_wRMS_Tp[5][3] )
if Write_wRMS_2nd:
    YJH_wRMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(YJH_wRMS_Tp_2[5][0], YJH_wRMS_Tp_2[5][1])
if WriteSimple_RMS_2nd:
    YJH_RMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(YJH_wRMS_Tp_2[5][2], YJH_wRMS_Tp_2[5][3] )

textfile_1.write('$YJH$     & Template     & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    YJH_wRMS_Tp[3][1],
    YJH_wRMS_Tp[0][0], YJH_wRMS_Tp[0][1], YJH_chi2dof_Tp_1,
    YJH_wRMS_Tp_2[0][0], YJH_wRMS_Tp_2[0][1], YJH_chi2dof_Tp_2,
    YJH_wRMS_Tp[5][0], YJH_wRMS_Tp[5][1], YJH_wRMS_Tp_2_txt,
    YJH_RMS_Tp_1_txt, YJH_RMS_Tp_2_txt ))

textfile_2.write('YJH      Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    YJH_wRMS_Tp[3][1],
    YJH_wRMS_Tp[0][0], YJH_wRMS_Tp[0][1], YJH_chi2dof_Tp_1,
    YJH_wRMS_Tp_2[0][0], YJH_wRMS_Tp_2[0][1], YJH_chi2dof_Tp_2,
    YJH_wRMS_Tp[5][0], YJH_wRMS_Tp[5][1], YJH_wRMS_Tp_2_txt,
    YJH_RMS_Tp_1_dat, YJH_RMS_Tp_2_txt ))

textfile_5.write('YJH      Template   %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    YJH_wRMS_Tp[3][1],
    YJH_wRMS_Tp[0][0], YJH_wRMS_Tp[0][1], YJH_chi2dof_Tp_1,
    YJH_wRMS_Tp_2[0][0], YJH_wRMS_Tp_2[0][1], YJH_chi2dof_Tp_2,
    YJH_wRMS_Tp[5][0], YJH_wRMS_Tp[5][1], YJH_wRMS_Tp_2_txt,
    YJH_RMS_Tp_1_dat, YJH_RMS_Tp_2_txt ))

#-------------------
# reset
YJH_RMS_GP_1_txt=''; YJH_wRMS_GP_2_txt=''; YJH_RMS_GP_2_txt='';
if WriteSimple_RMS:
    YJH_RMS_GP_1_txt = '& $%.3f \\pm %.3f$'%(YJH_wRMS_GP[5][2], YJH_wRMS_GP[5][3] )
    YJH_RMS_GP_1_dat = '%7.3f %7.3f'%(YJH_wRMS_GP[5][2], YJH_wRMS_GP[5][3] )
if Write_wRMS_2nd:
    YJH_wRMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(YJH_wRMS_GP_2[5][0], YJH_wRMS_GP_2[5][1])
if WriteSimple_RMS_2nd:
    YJH_RMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(YJH_wRMS_GP_2[5][2], YJH_wRMS_GP_2[5][3] )

textfile_1.write('$YJH$     & GP (NIR max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    YJH_wRMS_GP[3][1],
    YJH_wRMS_GP[0][0], YJH_wRMS_GP[0][1], YJH_chi2dof_GP_1,
    YJH_wRMS_GP_2[0][0], YJH_wRMS_GP_2[0][1], YJH_chi2dof_GP_2,
    YJH_wRMS_GP[5][0], YJH_wRMS_GP[5][1], YJH_wRMS_GP_2_txt,
    YJH_RMS_GP_1_txt, YJH_RMS_GP_2_txt ))

textfile_2.write('YJH      GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    YJH_wRMS_GP[3][1],
    YJH_wRMS_GP[0][0], YJH_wRMS_GP[0][1], YJH_chi2dof_GP_1,
    YJH_wRMS_GP_2[0][0], YJH_wRMS_GP_2[0][1], YJH_chi2dof_GP_2,
    YJH_wRMS_GP[5][0], YJH_wRMS_GP[5][1], YJH_wRMS_GP_2_txt,
    YJH_RMS_GP_1_dat, YJH_RMS_GP_2_txt ))

textfile_3.write('YJH      GP_NIRmax  %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    YJH_wRMS_GP[3][1],
    YJH_wRMS_GP[0][0], YJH_wRMS_GP[0][1], YJH_chi2dof_GP_1,
    YJH_wRMS_GP_2[0][0], YJH_wRMS_GP_2[0][1], YJH_chi2dof_GP_2,
    YJH_wRMS_GP[5][0], YJH_wRMS_GP[5][1], YJH_wRMS_GP_2_txt,
    YJH_RMS_GP_1_dat, YJH_RMS_GP_2_txt ))

#-------------------
# reset
YJH_RMS_GP_tBmax_1_txt=''; YJH_wRMS_GP_tBmax_2_txt=''; YJH_RMS_GP_tBmax_2_txt='';
if WriteSimple_RMS:
    YJH_RMS_GP_tBmax_1_txt = '& $%.3f \\pm %.3f$'%(YJH_wRMS_GP_tBmax[5][2], YJH_wRMS_GP_tBmax[5][3] )
    YJH_RMS_GP_tBmax_1_dat = '%7.3f %7.3f'%(YJH_wRMS_GP_tBmax[5][2], YJH_wRMS_GP_tBmax[5][3] )
if Write_wRMS_2nd:
    YJH_wRMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(YJH_wRMS_GP_tBmax_2[5][0], YJH_wRMS_GP_tBmax_2[5][1])
if WriteSimple_RMS_2nd:
    YJH_RMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(YJH_wRMS_GP_tBmax_2[5][2], YJH_wRMS_GP_tBmax_2[5][3] )

textfile_1.write('$YJH$     & GP ($B$ max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    YJH_wRMS_GP_tBmax[3][1],
    YJH_wRMS_GP_tBmax[0][0], YJH_wRMS_GP_tBmax[0][1], YJH_chi2dof_GP_tBmax_1,
    YJH_wRMS_GP_tBmax_2[0][0], YJH_wRMS_GP_tBmax_2[0][1], YJH_chi2dof_GP_tBmax_2,
    YJH_wRMS_GP_tBmax[5][0], YJH_wRMS_GP_tBmax[5][1], YJH_wRMS_GP_tBmax_2_txt,
    YJH_RMS_GP_tBmax_1_txt, YJH_RMS_GP_tBmax_2_txt ))

textfile_2.write('YJH      GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    YJH_wRMS_GP_tBmax[3][1],
    YJH_wRMS_GP_tBmax[0][0], YJH_wRMS_GP_tBmax[0][1], YJH_chi2dof_GP_tBmax_1,
    YJH_wRMS_GP_tBmax_2[0][0], YJH_wRMS_GP_tBmax_2[0][1], YJH_chi2dof_GP_tBmax_2,
    YJH_wRMS_GP_tBmax[5][0], YJH_wRMS_GP_tBmax[5][1], YJH_wRMS_GP_tBmax_2_txt,
    YJH_RMS_GP_tBmax_1_dat, YJH_RMS_GP_tBmax_2_txt ))

textfile_4.write('YJH      GP_Bmax    %3.0f %7.3f %7.3f %s %7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    YJH_wRMS_GP_tBmax[3][1],
    YJH_wRMS_GP_tBmax[0][0], YJH_wRMS_GP_tBmax[0][1], YJH_chi2dof_GP_tBmax_1,
    YJH_wRMS_GP_tBmax_2[0][0], YJH_wRMS_GP_tBmax_2[0][1], YJH_chi2dof_GP_tBmax_2,
    YJH_wRMS_GP_tBmax[5][0], YJH_wRMS_GP_tBmax[5][1], YJH_wRMS_GP_tBmax_2_txt,
    YJH_RMS_GP_tBmax_1_dat, YJH_RMS_GP_tBmax_2_txt ))

#-------------------
textfile_1.write('%------------------------------------- \n')
# reset
JHK_RMS_Tp_1_txt=''; JHK_wRMS_Tp_2_txt=''; JHK_RMS_Tp_2_txt='';
if WriteSimple_RMS:
    JHK_RMS_Tp_1_txt = '& $%.3f \\pm %.3f$'%(JHK_wRMS_Tp[5][2], JHK_wRMS_Tp[5][3] )
    JHK_RMS_Tp_1_dat = '%7.3f %7.3f'%(JHK_wRMS_Tp[5][2], JHK_wRMS_Tp[5][3] )
if Write_wRMS_2nd:
    JHK_wRMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(JHK_wRMS_Tp_2[5][0], JHK_wRMS_Tp_2[5][1])
if WriteSimple_RMS_2nd:
    JHK_RMS_Tp_2_txt = '& $%.3f \\pm %.3f$'%(JHK_wRMS_Tp_2[5][2], JHK_wRMS_Tp_2[5][3] )

textfile_1.write('%% $JHK_s$     & Template     & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    JHK_wRMS_Tp[3][1],
    JHK_wRMS_Tp[0][0], JHK_wRMS_Tp[0][1], JHK_chi2dof_Tp_1,
    JHK_wRMS_Tp_2[0][0], JHK_wRMS_Tp_2[0][1], JHK_chi2dof_Tp_2,
    JHK_wRMS_Tp[5][0], JHK_wRMS_Tp[5][1], JHK_wRMS_Tp_2_txt,
    JHK_RMS_Tp_1_txt, JHK_RMS_Tp_2_txt ))

textfile_2.write('# JHK      Template   %3.0f %7.3f %7.3f %s \
%7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    JHK_wRMS_Tp[3][1],
    JHK_wRMS_Tp[0][0], JHK_wRMS_Tp[0][1], JHK_chi2dof_Tp_1,
    JHK_wRMS_Tp_2[0][0], JHK_wRMS_Tp_2[0][1], JHK_chi2dof_Tp_2,
    JHK_wRMS_Tp[5][0], JHK_wRMS_Tp[5][1], JHK_wRMS_Tp_2_txt,
    JHK_RMS_Tp_1_dat, JHK_RMS_Tp_2_txt ))

textfile_5.write('# JHK      Template   %3.0f %7.3f %7.3f %s \
%7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    JHK_wRMS_Tp[3][1],
    JHK_wRMS_Tp[0][0], JHK_wRMS_Tp[0][1], JHK_chi2dof_Tp_1,
    JHK_wRMS_Tp_2[0][0], JHK_wRMS_Tp_2[0][1], JHK_chi2dof_Tp_2,
    JHK_wRMS_Tp[5][0], JHK_wRMS_Tp[5][1], JHK_wRMS_Tp_2_txt,
    JHK_RMS_Tp_1_dat, JHK_RMS_Tp_2_txt ))

#-------------------
# reset
JHK_RMS_GP_1_txt=''; JHK_wRMS_GP_2_txt=''; JHK_RMS_GP_2_txt='';
if WriteSimple_RMS:
    JHK_RMS_GP_1_txt = '& $%.3f \\pm %.3f$'%(JHK_wRMS_GP[5][2], JHK_wRMS_GP[5][3] )
    JHK_RMS_GP_1_dat = '%7.3f %7.3f'%(JHK_wRMS_GP[5][2], JHK_wRMS_GP[5][3] )
if Write_wRMS_2nd:
    JHK_wRMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(JHK_wRMS_GP_2[5][0], JHK_wRMS_GP_2[5][1])
if WriteSimple_RMS_2nd:
    JHK_RMS_GP_2_txt = '& $%.3f \\pm %.3f$'%(JHK_wRMS_GP_2[5][2], JHK_wRMS_GP_2[5][3] )

textfile_1.write('%% $JHK_s$     & GP (NIR max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    JHK_wRMS_GP[3][1],
    JHK_wRMS_GP[0][0], JHK_wRMS_GP[0][1], JHK_chi2dof_GP_1,
    JHK_wRMS_GP_2[0][0], JHK_wRMS_GP_2[0][1], JHK_chi2dof_GP_2,
    JHK_wRMS_GP[5][0], JHK_wRMS_GP[5][1], JHK_wRMS_GP_2_txt,
    JHK_RMS_GP_1_txt, JHK_RMS_GP_2_txt ))

textfile_2.write('# JHK      GP_NIRmax  %3.0f %7.3f %7.3f %s \
%7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    JHK_wRMS_GP[3][1],
    JHK_wRMS_GP[0][0], JHK_wRMS_GP[0][1], JHK_chi2dof_GP_1,
    JHK_wRMS_GP_2[0][0], JHK_wRMS_GP_2[0][1], JHK_chi2dof_GP_2,
    JHK_wRMS_GP[5][0], JHK_wRMS_GP[5][1], JHK_wRMS_GP_2_txt,
    JHK_RMS_GP_1_dat, JHK_RMS_GP_2_txt ))

textfile_3.write('# JHK      GP_NIRmax  %3.0f %7.3f %7.3f %s \
%7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    JHK_wRMS_GP[3][1],
    JHK_wRMS_GP[0][0], JHK_wRMS_GP[0][1], JHK_chi2dof_GP_1,
    JHK_wRMS_GP_2[0][0], JHK_wRMS_GP_2[0][1], JHK_chi2dof_GP_2,
    JHK_wRMS_GP[5][0], JHK_wRMS_GP[5][1], JHK_wRMS_GP_2_txt,
    JHK_RMS_GP_1_dat, JHK_RMS_GP_2_txt ))

#-------------------
# reset
JHK_RMS_GP_tBmax_1_txt=''; JHK_wRMS_GP_tBmax_2_txt=''; JHK_RMS_GP_tBmax_2_txt='';
if WriteSimple_RMS:
    JHK_RMS_GP_tBmax_1_txt = '& $%.3f \\pm %.3f$'%(JHK_wRMS_GP_tBmax[5][2], JHK_wRMS_GP_tBmax[5][3] )
    JHK_RMS_GP_tBmax_1_dat = '%7.3f %7.3f'%(JHK_wRMS_GP_tBmax[5][2], JHK_wRMS_GP_tBmax[5][3] )
if Write_wRMS_2nd:
    JHK_wRMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(JHK_wRMS_GP_tBmax_2[5][0], JHK_wRMS_GP_tBmax_2[5][1])
if WriteSimple_RMS_2nd:
    JHK_RMS_GP_tBmax_2_txt = '& $%.3f \\pm %.3f$'%(JHK_wRMS_GP_tBmax_2[5][2], JHK_wRMS_GP_tBmax_2[5][3] )

textfile_1.write('%% $JHK_s$     & GP ($B$ max) & %.0f & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s & $%.3f \\pm %.3f$ %s %s %s \\\\ \n'%(
    JHK_wRMS_GP_tBmax[3][1],
    JHK_wRMS_GP_tBmax[0][0], JHK_wRMS_GP_tBmax[0][1], JHK_chi2dof_GP_tBmax_1,
    JHK_wRMS_GP_tBmax_2[0][0], JHK_wRMS_GP_tBmax_2[0][1], JHK_chi2dof_GP_tBmax_2,
    JHK_wRMS_GP_tBmax[5][0], JHK_wRMS_GP_tBmax[5][1], JHK_wRMS_GP_tBmax_2_txt,
    JHK_RMS_GP_tBmax_1_txt, JHK_RMS_GP_tBmax_2_txt ))

textfile_2.write('# JHK      GP_Bmax    %3.0f %7.3f %7.3f %s \
%7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    JHK_wRMS_GP_tBmax[3][1],
    JHK_wRMS_GP_tBmax[0][0], JHK_wRMS_GP_tBmax[0][1], JHK_chi2dof_GP_tBmax_1,
    JHK_wRMS_GP_tBmax_2[0][0], JHK_wRMS_GP_tBmax_2[0][1], JHK_chi2dof_GP_tBmax_2,
    JHK_wRMS_GP_tBmax[5][0], JHK_wRMS_GP_tBmax[5][1], JHK_wRMS_GP_tBmax_2_txt,
    JHK_RMS_GP_tBmax_1_dat, JHK_RMS_GP_tBmax_2_txt ))

textfile_4.write('# JHK      GP_Bmax    %3.0f %7.3f %7.3f %s \
%7.3f %7.3f %s %7.3f %7.3f %s %s %s \n'%(
    JHK_wRMS_GP_tBmax[3][1],
    JHK_wRMS_GP_tBmax[0][0], JHK_wRMS_GP_tBmax[0][1], JHK_chi2dof_GP_tBmax_1,
    JHK_wRMS_GP_tBmax_2[0][0], JHK_wRMS_GP_tBmax_2[0][1], JHK_chi2dof_GP_tBmax_2,
    JHK_wRMS_GP_tBmax[5][0], JHK_wRMS_GP_tBmax[5][1], JHK_wRMS_GP_tBmax_2_txt,
    JHK_RMS_GP_tBmax_1_dat, JHK_RMS_GP_tBmax_2_txt ))

#-------------------
textfile_1.write('%------------------------------------- \n')

# << --------- PART GENERATED WITH THE CODE BELOW ---
###################################################

textfile_1.close(); textfile_2.close(); textfile_3.close();
textfile_4.close(); textfile_5.close();

if ScriptVersion == 'notebook':
    print "All done with no issues :)"

textfile_1.close();textfile_1.close();textfile_1.close();
textfile_1.close();textfile_1.close();textfile_1.close();

textfile_2.close();textfile_2.close();textfile_2.close();
textfile_2.close();textfile_2.close();textfile_2.close();

textfile_3.close();textfile_3.close();textfile_3.close();
textfile_4.close();textfile_4.close();textfile_4.close();
textfile_5.close();textfile_5.close();textfile_5.close();

#-----------------------------------------------------------------------------80

# ### Code to generate the NIR part of the python script written above.
#
# #### This is very useful!!

# Comment this cell if I'm NOT working on generating the python code to
# be used in the cell above.

band_list = ['Y', 'J','H','K', 'AnyNIR', 'JH', 'YJH', 'JHK']
band_text_list = ['$Y$', '$J$','$H$','$K$',
                  'any $YJHK_s$',
                  '$JH$ ', '$YJH$', '$JHK_s$']

method_list = ['Tp', 'GP', 'GP_tBmax']
method_text_nir_list = ['Template    ', 'GP (NIR max)', 'GP ($B$ max)']
method_text_nir_list2 = ['Template', 'GP_NIRmax', 'GP_Bmax']

######################################################
#    Create the text file to save the generated python script

textfile_1 = open(DirSavePythonOutput+'Code_Table_Latex_scatter_.py','w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_Date_2   = '# On date: %s \n'%text_timenow
text_Author_2 = '# Data table created by: Arturo Avelino \n'
text_script_2 = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
text_line_2= '#'+'-'*60 + '\n'

textfile_1.write(text_line_2)
textfile_1.write(text_Author_2); textfile_1.write(text_Date_2);
textfile_1.write(text_script_2);
textfile_1.write(text_line_2)
textfile_1.write('#\n')

######################################################

for jj in range(len(band_list)):
    bb = band_list[jj]
    bb_text = band_text_list[jj]

    for ii in range(len(method_list)):
        i1 = method_list[ii]
        i1_text = method_text_nir_list[ii]
        i2_text = method_text_nir_list2[ii]

        ## Comment the results from 'JHK' because they produce no sense values.
        if bb == 'JHK':
            comment1 = '%% '; comment2 = '# ';
        else:
            comment1 = ''; comment2 = '';

        textfile_1.write("# reset \n")
        textfile_1.write("%s_RMS_%s_1_txt=''; %s_wRMS_%s_2_txt=''; %s_RMS_%s_2_txt='';\n"%(
        bb,i1,bb,i1,bb,i1))
        textfile_1.write("if WriteSimple_RMS:\n")
        textfile_1.write("    %s_RMS_%s_1_txt = '& $%%.3f \\\\pm %%.3f$'%%(%s_wRMS_%s[5][2], %s_wRMS_%s[5][3] )\n"%(bb,i1,bb,i1, bb,i1))
        textfile_1.write("    %s_RMS_%s_1_dat = '%%7.3f %%7.3f'%%(%s_wRMS_%s[5][2], %s_wRMS_%s[5][3] )\n"%(bb,i1,bb,i1, bb,i1))
        textfile_1.write("if Write_wRMS_2nd:\n")
        textfile_1.write("    %s_wRMS_%s_2_txt = '& $%%.3f \\\\pm %%.3f$'%%(%s_wRMS_%s_2[5][0], %s_wRMS_%s_2[5][1])\n"%(bb,i1,bb,i1,bb,i1))
        textfile_1.write("if WriteSimple_RMS_2nd:\n")
        textfile_1.write("    %s_RMS_%s_2_txt = '& $%%.3f \\\\pm %%.3f$'%%(%s_wRMS_%s_2[5][2], %s_wRMS_%s_2[5][3] )\n"%(bb,i1,bb,i1,bb,i1))
        textfile_1.write("\n")

        textfile_1.write("textfile_1.write('%s%s     & %s & %%.0f & $%%.3f \\\\pm %%.3f$ %%s & \\\n"%(
        comment1, bb_text, i1_text))
        textfile_1.write("$%.3f \\\\pm %.3f$ %s & $%.3f \\\\pm %.3f$ %s %s %s \\\\\\\\ \\n'%(\n")
        textfile_1.write("    %s_wRMS_%s[3][1], \n"%(bb,i1))
        textfile_1.write("    %s_wRMS_%s[0][0], %s_wRMS_%s[0][1], %s_chi2dof_%s_1, \n"%(
        bb,i1,bb,i1,bb,i1))
        textfile_1.write("    %s_wRMS_%s_2[0][0], %s_wRMS_%s_2[0][1], %s_chi2dof_%s_2, \n"%(
        bb,i1,bb,i1,bb,i1))
        textfile_1.write("    %s_wRMS_%s[5][0], %s_wRMS_%s[5][1], %s_wRMS_%s_2_txt, \n"%(
        bb,i1,bb,i1,bb,i1))
        textfile_1.write("    %s_RMS_%s_1_txt, %s_RMS_%s_2_txt )) \n"%(bb,i1,bb,i1))
        textfile_1.write("\n")

        #----------------------------------------------------------
        #    WRITE THE TABLES IN DATA FORMAT INSTEAD OF LATEX FORMAT

        textfile_1.write("textfile_2.write('%s%-8s %-10s %%3.0f %%7.3f %%7.3f %%s \\\n"%(
        comment2, bb, i2_text))
        textfile_1.write("%7.3f %7.3f %s %7.3f %7.3f %s %s %s \\n'%(\n")

        # Copy/paste from the paragraph above:
        textfile_1.write("    %s_wRMS_%s[3][1], \n"%(bb,i1))
        textfile_1.write("    %s_wRMS_%s[0][0], %s_wRMS_%s[0][1], %s_chi2dof_%s_1, \n"%(
        bb,i1,bb,i1,bb,i1))
        textfile_1.write("    %s_wRMS_%s_2[0][0], %s_wRMS_%s_2[0][1], %s_chi2dof_%s_2, \n"%(
        bb,i1,bb,i1,bb,i1))
        textfile_1.write("    %s_wRMS_%s[5][0], %s_wRMS_%s[5][1], %s_wRMS_%s_2_txt, \n"%(
        bb,i1,bb,i1,bb,i1))
        textfile_1.write("    %s_RMS_%s_1_dat, %s_RMS_%s_2_txt )) \n"%(bb,i1,bb,i1))
        textfile_1.write("\n")

        if i1 == 'GP':
            textfile_1.write("textfile_3.write('%s%-8s %-10s %%3.0f %%7.3f %%7.3f %%s \\\n"%(
            comment2, bb, i2_text))
            textfile_1.write("%7.3f %7.3f %s %7.3f %7.3f %s %s %s \\n'%(\n")

        elif i1 == 'GP_tBmax':
            textfile_1.write("textfile_4.write('%s%-8s %-10s %%3.0f %%7.3f %%7.3f %%s \\\n"%(
            comment2, bb, i2_text))
            textfile_1.write("%7.3f %7.3f %s %7.3f %7.3f %s %s %s \\n'%(\n")

        elif i1 == 'Tp':
            textfile_1.write("textfile_5.write('%s%-8s %-10s %%3.0f %%7.3f %%7.3f %%s \\\n"%(
            comment2, bb, i2_text))
            textfile_1.write("%7.3f %7.3f %s %7.3f %7.3f %s %s %s \\n'%(\n")

        # Copy/paste from the paragraph above:
        textfile_1.write("    %s_wRMS_%s[3][1], \n"%(bb,i1))
        textfile_1.write("    %s_wRMS_%s[0][0], %s_wRMS_%s[0][1], %s_chi2dof_%s_1, \n"%(
        bb,i1,bb,i1,bb,i1))
        textfile_1.write("    %s_wRMS_%s_2[0][0], %s_wRMS_%s_2[0][1], %s_chi2dof_%s_2, \n"%(
        bb,i1,bb,i1,bb,i1))
        textfile_1.write("    %s_wRMS_%s[5][0], %s_wRMS_%s[5][1], %s_wRMS_%s_2_txt, \n"%(
        bb,i1,bb,i1,bb,i1))
        textfile_1.write("    %s_RMS_%s_1_dat, %s_RMS_%s_2_txt )) \n"%(bb,i1,bb,i1))


        #----------------------------------------------------------

        textfile_1.write("\n#%s\n\n"%('-'*20))
    textfile_1.write("textfile_1.write('%------------------------------------- \\n') \n")

textfile_1.close()
if ScriptVersion == 'notebook':
    print "All done with no issues :)"

textfile_1.close();textfile_1.close();textfile_1.close();
textfile_1.close();textfile_1.close();textfile_1.close();

##############################################################################80

# ## Latex Table:  Optical vs NIR intrinsic scatter (Andy's table)

band_list = ['Y', 'J','H','K', 'AnyNIR', 'JH', 'YJH', 'JHK']

textfile_1 = open(DirSaveOutput+'Table_Latex_scatter_Opt_vs_NIR_.tex', 'w')
textfile_2 = open(DirSaveOutput+'Table_Latex_scatter_Opt_vs_NIR_.dat', 'w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_Date   = '%% On date: %s \n'%text_timenow

textfile_1.write(text_line)
textfile_1.write(text_Author); textfile_1.write(text_Date);
textfile_1.write(text_script);
textfile_1.write(text_line)
textfile_1.write(' \n')

textfile_2.write('#'+text_line)
textfile_2.write('#'+text_Author); textfile_2.write('#'+text_Date);
textfile_2.write('#'+text_script);
textfile_2.write('#'+text_line)

txt_hline = "\\hline \n"

# textfile_1.write("Optical $BVR$ Method - NIR band(s) & $\\Delta\\sigmaint$ & $n$-$\sigma$ \\\\ \n")
# textfile_1.write(txt_hline)

###################################################
# ---- PART GENERATED WITH THE CODE BELOW ---------->>

#-----------------------------------------------------------------------------80
# Data table created by: Arturo Avelino
# On date: 2018-12-19 (yyyy-mm-dd); 13:21 hrs.
# Script used: 15_Latex_scatter.ipynb
#-----------------------------------------------------------------------------80
#
salt2_y     = SALT2_RMSData[0][0] - Y_wRMS_GP[0][0]
err_salt2_y = np.sqrt(SALT2_RMSData[0][1]**2 + Y_wRMS_GP[0][1]**2)
ratio_y = salt2_y/err_salt2_y

salt2_y_wrms     = SALT2_RMSData[5][0] - Y_wRMS_GP[5][0]
err_salt2_y_wrms = np.sqrt(SALT2_RMSData[5][1]**2 + Y_wRMS_GP[5][1]**2)
ratio_y_wrms = salt2_y_wrms/err_salt2_y_wrms

salt2_y_rms     = SALT2_RMSData[5][2] - Y_wRMS_GP[5][2]
err_salt2_y_rms = np.sqrt(SALT2_RMSData[5][3]**2 + Y_wRMS_GP[5][3]**2)
ratio_y_rms = salt2_y_rms/err_salt2_y_rms

textfile_1.write('SALT2 - $Y$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    salt2_y, err_salt2_y, ratio_y,
    salt2_y_wrms, err_salt2_y_wrms, ratio_y_wrms,
    salt2_y_rms, err_salt2_y_rms, ratio_y_rms ))

textfile_2.write('SALTminusY           %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f \n'%(
    salt2_y, err_salt2_y, ratio_y,
    salt2_y_wrms, err_salt2_y_wrms, ratio_y_wrms,
    salt2_y_rms, err_salt2_y_rms, ratio_y_rms ))

#-------------------
snoopy_y     = Snoopy_RMSData[0][0] - Y_wRMS_GP[0][0]
err_snoopy_y = np.sqrt(Snoopy_RMSData[0][1]**2 + Y_wRMS_GP[0][1]**2)
ratio_y = snoopy_y/err_snoopy_y

snoopy_y_wrms     = Snoopy_RMSData[5][0] - Y_wRMS_GP[5][0]
err_snoopy_y_wrms = np.sqrt(Snoopy_RMSData[5][1]**2 + Y_wRMS_GP[5][1]**2)
ratio_y_wrms = snoopy_y_wrms/err_snoopy_y_wrms

snoopy_y_rms     = Snoopy_RMSData[5][2] - Y_wRMS_GP[5][2]
err_snoopy_y_rms = np.sqrt(Snoopy_RMSData[5][3]**2 + Y_wRMS_GP[5][3]**2)
ratio_y_rms = snoopy_y_rms/err_snoopy_y_rms

textfile_1.write('\snoopy{} - $Y$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    snoopy_y, err_snoopy_y, ratio_y,
    snoopy_y_wrms, err_snoopy_y_wrms, ratio_y_wrms,
    snoopy_y_rms, err_snoopy_y_rms, ratio_y_rms ))

textfile_2.write('SNOOPYminusY         %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f\n'%(
    snoopy_y, err_snoopy_y, ratio_y,
    snoopy_y_wrms, err_snoopy_y_wrms, ratio_y_wrms,
    snoopy_y_rms, err_snoopy_y_rms, ratio_y_rms ))

textfile_1.write(txt_hline)
#======================================

salt2_j     = SALT2_RMSData[0][0] - J_wRMS_GP[0][0]
err_salt2_j = np.sqrt(SALT2_RMSData[0][1]**2 + J_wRMS_GP[0][1]**2)
ratio_j = salt2_j/err_salt2_j

salt2_j_wrms     = SALT2_RMSData[5][0] - J_wRMS_GP[5][0]
err_salt2_j_wrms = np.sqrt(SALT2_RMSData[5][1]**2 + J_wRMS_GP[5][1]**2)
ratio_j_wrms = salt2_j_wrms/err_salt2_j_wrms

salt2_j_rms     = SALT2_RMSData[5][2] - J_wRMS_GP[5][2]
err_salt2_j_rms = np.sqrt(SALT2_RMSData[5][3]**2 + J_wRMS_GP[5][3]**2)
ratio_j_rms = salt2_j_rms/err_salt2_j_rms

textfile_1.write('SALT2 - $J$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    salt2_j, err_salt2_j, ratio_j,
    salt2_j_wrms, err_salt2_j_wrms, ratio_j_wrms,
    salt2_j_rms, err_salt2_j_rms, ratio_j_rms ))

textfile_2.write('SALTminusJ           %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f \n'%(
    salt2_j, err_salt2_j, ratio_j,
    salt2_j_wrms, err_salt2_j_wrms, ratio_j_wrms,
    salt2_j_rms, err_salt2_j_rms, ratio_j_rms ))

#-------------------
snoopy_j     = Snoopy_RMSData[0][0] - J_wRMS_GP[0][0]
err_snoopy_j = np.sqrt(Snoopy_RMSData[0][1]**2 + J_wRMS_GP[0][1]**2)
ratio_j = snoopy_j/err_snoopy_j

snoopy_j_wrms     = Snoopy_RMSData[5][0] - J_wRMS_GP[5][0]
err_snoopy_j_wrms = np.sqrt(Snoopy_RMSData[5][1]**2 + J_wRMS_GP[5][1]**2)
ratio_j_wrms = snoopy_j_wrms/err_snoopy_j_wrms

snoopy_j_rms     = Snoopy_RMSData[5][2] - J_wRMS_GP[5][2]
err_snoopy_j_rms = np.sqrt(Snoopy_RMSData[5][3]**2 + J_wRMS_GP[5][3]**2)
ratio_j_rms = snoopy_j_rms/err_snoopy_j_rms

textfile_1.write('\snoopy{} - $J$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    snoopy_j, err_snoopy_j, ratio_j,
    snoopy_j_wrms, err_snoopy_j_wrms, ratio_j_wrms,
    snoopy_j_rms, err_snoopy_j_rms, ratio_j_rms ))

textfile_2.write('SNOOPYminusJ         %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f\n'%(
    snoopy_j, err_snoopy_j, ratio_j,
    snoopy_j_wrms, err_snoopy_j_wrms, ratio_j_wrms,
    snoopy_j_rms, err_snoopy_j_rms, ratio_j_rms ))

textfile_1.write(txt_hline)
#======================================

salt2_h     = SALT2_RMSData[0][0] - H_wRMS_GP[0][0]
err_salt2_h = np.sqrt(SALT2_RMSData[0][1]**2 + H_wRMS_GP[0][1]**2)
ratio_h = salt2_h/err_salt2_h

salt2_h_wrms     = SALT2_RMSData[5][0] - H_wRMS_GP[5][0]
err_salt2_h_wrms = np.sqrt(SALT2_RMSData[5][1]**2 + H_wRMS_GP[5][1]**2)
ratio_h_wrms = salt2_h_wrms/err_salt2_h_wrms

salt2_h_rms     = SALT2_RMSData[5][2] - H_wRMS_GP[5][2]
err_salt2_h_rms = np.sqrt(SALT2_RMSData[5][3]**2 + H_wRMS_GP[5][3]**2)
ratio_h_rms = salt2_h_rms/err_salt2_h_rms

textfile_1.write('SALT2 - $H$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    salt2_h, err_salt2_h, ratio_h,
    salt2_h_wrms, err_salt2_h_wrms, ratio_h_wrms,
    salt2_h_rms, err_salt2_h_rms, ratio_h_rms ))

textfile_2.write('SALTminusH           %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f \n'%(
    salt2_h, err_salt2_h, ratio_h,
    salt2_h_wrms, err_salt2_h_wrms, ratio_h_wrms,
    salt2_h_rms, err_salt2_h_rms, ratio_h_rms ))

#-------------------
snoopy_h     = Snoopy_RMSData[0][0] - H_wRMS_GP[0][0]
err_snoopy_h = np.sqrt(Snoopy_RMSData[0][1]**2 + H_wRMS_GP[0][1]**2)
ratio_h = snoopy_h/err_snoopy_h

snoopy_h_wrms     = Snoopy_RMSData[5][0] - H_wRMS_GP[5][0]
err_snoopy_h_wrms = np.sqrt(Snoopy_RMSData[5][1]**2 + H_wRMS_GP[5][1]**2)
ratio_h_wrms = snoopy_h_wrms/err_snoopy_h_wrms

snoopy_h_rms     = Snoopy_RMSData[5][2] - H_wRMS_GP[5][2]
err_snoopy_h_rms = np.sqrt(Snoopy_RMSData[5][3]**2 + H_wRMS_GP[5][3]**2)
ratio_h_rms = snoopy_h_rms/err_snoopy_h_rms

textfile_1.write('\snoopy{} - $H$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    snoopy_h, err_snoopy_h, ratio_h,
    snoopy_h_wrms, err_snoopy_h_wrms, ratio_h_wrms,
    snoopy_h_rms, err_snoopy_h_rms, ratio_h_rms ))

textfile_2.write('SNOOPYminusH         %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f\n'%(
    snoopy_h, err_snoopy_h, ratio_h,
    snoopy_h_wrms, err_snoopy_h_wrms, ratio_h_wrms,
    snoopy_h_rms, err_snoopy_h_rms, ratio_h_rms ))

textfile_1.write(txt_hline)
#======================================

salt2_k     = SALT2_RMSData[0][0] - K_wRMS_GP[0][0]
err_salt2_k = np.sqrt(SALT2_RMSData[0][1]**2 + K_wRMS_GP[0][1]**2)
ratio_k = salt2_k/err_salt2_k

salt2_k_wrms     = SALT2_RMSData[5][0] - K_wRMS_GP[5][0]
err_salt2_k_wrms = np.sqrt(SALT2_RMSData[5][1]**2 + K_wRMS_GP[5][1]**2)
ratio_k_wrms = salt2_k_wrms/err_salt2_k_wrms

salt2_k_rms     = SALT2_RMSData[5][2] - K_wRMS_GP[5][2]
err_salt2_k_rms = np.sqrt(SALT2_RMSData[5][3]**2 + K_wRMS_GP[5][3]**2)
ratio_k_rms = salt2_k_rms/err_salt2_k_rms

textfile_1.write('SALT2 - $K_s$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    salt2_k, err_salt2_k, ratio_k,
    salt2_k_wrms, err_salt2_k_wrms, ratio_k_wrms,
    salt2_k_rms, err_salt2_k_rms, ratio_k_rms ))

textfile_2.write('SALTminusK           %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f \n'%(
    salt2_k, err_salt2_k, ratio_k,
    salt2_k_wrms, err_salt2_k_wrms, ratio_k_wrms,
    salt2_k_rms, err_salt2_k_rms, ratio_k_rms ))

#-------------------
snoopy_k     = Snoopy_RMSData[0][0] - K_wRMS_GP[0][0]
err_snoopy_k = np.sqrt(Snoopy_RMSData[0][1]**2 + K_wRMS_GP[0][1]**2)
ratio_k = snoopy_k/err_snoopy_k

snoopy_k_wrms     = Snoopy_RMSData[5][0] - K_wRMS_GP[5][0]
err_snoopy_k_wrms = np.sqrt(Snoopy_RMSData[5][1]**2 + K_wRMS_GP[5][1]**2)
ratio_k_wrms = snoopy_k_wrms/err_snoopy_k_wrms

snoopy_k_rms     = Snoopy_RMSData[5][2] - K_wRMS_GP[5][2]
err_snoopy_k_rms = np.sqrt(Snoopy_RMSData[5][3]**2 + K_wRMS_GP[5][3]**2)
ratio_k_rms = snoopy_k_rms/err_snoopy_k_rms

textfile_1.write('\snoopy{} - $K_s$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    snoopy_k, err_snoopy_k, ratio_k,
    snoopy_k_wrms, err_snoopy_k_wrms, ratio_k_wrms,
    snoopy_k_rms, err_snoopy_k_rms, ratio_k_rms ))

textfile_2.write('SNOOPYminusK         %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f\n'%(
    snoopy_k, err_snoopy_k, ratio_k,
    snoopy_k_wrms, err_snoopy_k_wrms, ratio_k_wrms,
    snoopy_k_rms, err_snoopy_k_rms, ratio_k_rms ))

textfile_1.write(txt_hline)
#======================================

salt2_anynir     = SALT2_RMSData[0][0] - AnyNIR_wRMS_GP[0][0]
err_salt2_anynir = np.sqrt(SALT2_RMSData[0][1]**2 + AnyNIR_wRMS_GP[0][1]**2)
ratio_anynir = salt2_anynir/err_salt2_anynir

salt2_anynir_wrms     = SALT2_RMSData[5][0] - AnyNIR_wRMS_GP[5][0]
err_salt2_anynir_wrms = np.sqrt(SALT2_RMSData[5][1]**2 + AnyNIR_wRMS_GP[5][1]**2)
ratio_anynir_wrms = salt2_anynir_wrms/err_salt2_anynir_wrms

salt2_anynir_rms     = SALT2_RMSData[5][2] - AnyNIR_wRMS_GP[5][2]
err_salt2_anynir_rms = np.sqrt(SALT2_RMSData[5][3]**2 + AnyNIR_wRMS_GP[5][3]**2)
ratio_anynir_rms = salt2_anynir_rms/err_salt2_anynir_rms

textfile_1.write('SALT2 - any $YJHK_s$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    salt2_anynir, err_salt2_anynir, ratio_anynir,
    salt2_anynir_wrms, err_salt2_anynir_wrms, ratio_anynir_wrms,
    salt2_anynir_rms, err_salt2_anynir_rms, ratio_anynir_rms ))

textfile_2.write('SALTminusAnyNIR      %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f \n'%(
    salt2_anynir, err_salt2_anynir, ratio_anynir,
    salt2_anynir_wrms, err_salt2_anynir_wrms, ratio_anynir_wrms,
    salt2_anynir_rms, err_salt2_anynir_rms, ratio_anynir_rms ))

#-------------------
snoopy_anynir     = Snoopy_RMSData[0][0] - AnyNIR_wRMS_GP[0][0]
err_snoopy_anynir = np.sqrt(Snoopy_RMSData[0][1]**2 + AnyNIR_wRMS_GP[0][1]**2)
ratio_anynir = snoopy_anynir/err_snoopy_anynir

snoopy_anynir_wrms     = Snoopy_RMSData[5][0] - AnyNIR_wRMS_GP[5][0]
err_snoopy_anynir_wrms = np.sqrt(Snoopy_RMSData[5][1]**2 + AnyNIR_wRMS_GP[5][1]**2)
ratio_anynir_wrms = snoopy_anynir_wrms/err_snoopy_anynir_wrms

snoopy_anynir_rms     = Snoopy_RMSData[5][2] - AnyNIR_wRMS_GP[5][2]
err_snoopy_anynir_rms = np.sqrt(Snoopy_RMSData[5][3]**2 + AnyNIR_wRMS_GP[5][3]**2)
ratio_anynir_rms = snoopy_anynir_rms/err_snoopy_anynir_rms

textfile_1.write('\snoopy{} - any $YJHK_s$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    snoopy_anynir, err_snoopy_anynir, ratio_anynir,
    snoopy_anynir_wrms, err_snoopy_anynir_wrms, ratio_anynir_wrms,
    snoopy_anynir_rms, err_snoopy_anynir_rms, ratio_anynir_rms ))

textfile_2.write('SNOOPYminusAnyNIR    %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f\n'%(
    snoopy_anynir, err_snoopy_anynir, ratio_anynir,
    snoopy_anynir_wrms, err_snoopy_anynir_wrms, ratio_anynir_wrms,
    snoopy_anynir_rms, err_snoopy_anynir_rms, ratio_anynir_rms ))

textfile_1.write(txt_hline)
#======================================

salt2_jh     = SALT2_RMSData[0][0] - JH_wRMS_GP[0][0]
err_salt2_jh = np.sqrt(SALT2_RMSData[0][1]**2 + JH_wRMS_GP[0][1]**2)
ratio_jh = salt2_jh/err_salt2_jh

salt2_jh_wrms     = SALT2_RMSData[5][0] - JH_wRMS_GP[5][0]
err_salt2_jh_wrms = np.sqrt(SALT2_RMSData[5][1]**2 + JH_wRMS_GP[5][1]**2)
ratio_jh_wrms = salt2_jh_wrms/err_salt2_jh_wrms

salt2_jh_rms     = SALT2_RMSData[5][2] - JH_wRMS_GP[5][2]
err_salt2_jh_rms = np.sqrt(SALT2_RMSData[5][3]**2 + JH_wRMS_GP[5][3]**2)
ratio_jh_rms = salt2_jh_rms/err_salt2_jh_rms

textfile_1.write('SALT2 - $JH$  & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    salt2_jh, err_salt2_jh, ratio_jh,
    salt2_jh_wrms, err_salt2_jh_wrms, ratio_jh_wrms,
    salt2_jh_rms, err_salt2_jh_rms, ratio_jh_rms ))

textfile_2.write('SALTminusJH          %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f \n'%(
    salt2_jh, err_salt2_jh, ratio_jh,
    salt2_jh_wrms, err_salt2_jh_wrms, ratio_jh_wrms,
    salt2_jh_rms, err_salt2_jh_rms, ratio_jh_rms ))

#-------------------
snoopy_jh     = Snoopy_RMSData[0][0] - JH_wRMS_GP[0][0]
err_snoopy_jh = np.sqrt(Snoopy_RMSData[0][1]**2 + JH_wRMS_GP[0][1]**2)
ratio_jh = snoopy_jh/err_snoopy_jh

snoopy_jh_wrms     = Snoopy_RMSData[5][0] - JH_wRMS_GP[5][0]
err_snoopy_jh_wrms = np.sqrt(Snoopy_RMSData[5][1]**2 + JH_wRMS_GP[5][1]**2)
ratio_jh_wrms = snoopy_jh_wrms/err_snoopy_jh_wrms

snoopy_jh_rms     = Snoopy_RMSData[5][2] - JH_wRMS_GP[5][2]
err_snoopy_jh_rms = np.sqrt(Snoopy_RMSData[5][3]**2 + JH_wRMS_GP[5][3]**2)
ratio_jh_rms = snoopy_jh_rms/err_snoopy_jh_rms

textfile_1.write('\snoopy{} - $JH$  & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    snoopy_jh, err_snoopy_jh, ratio_jh,
    snoopy_jh_wrms, err_snoopy_jh_wrms, ratio_jh_wrms,
    snoopy_jh_rms, err_snoopy_jh_rms, ratio_jh_rms ))

textfile_2.write('SNOOPYminusJH        %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f\n'%(
    snoopy_jh, err_snoopy_jh, ratio_jh,
    snoopy_jh_wrms, err_snoopy_jh_wrms, ratio_jh_wrms,
    snoopy_jh_rms, err_snoopy_jh_rms, ratio_jh_rms ))

textfile_1.write(txt_hline)
#======================================

salt2_yjh     = SALT2_RMSData[0][0] - YJH_wRMS_GP[0][0]
err_salt2_yjh = np.sqrt(SALT2_RMSData[0][1]**2 + YJH_wRMS_GP[0][1]**2)
ratio_yjh = salt2_yjh/err_salt2_yjh

salt2_yjh_wrms     = SALT2_RMSData[5][0] - YJH_wRMS_GP[5][0]
err_salt2_yjh_wrms = np.sqrt(SALT2_RMSData[5][1]**2 + YJH_wRMS_GP[5][1]**2)
ratio_yjh_wrms = salt2_yjh_wrms/err_salt2_yjh_wrms

salt2_yjh_rms     = SALT2_RMSData[5][2] - YJH_wRMS_GP[5][2]
err_salt2_yjh_rms = np.sqrt(SALT2_RMSData[5][3]**2 + YJH_wRMS_GP[5][3]**2)
ratio_yjh_rms = salt2_yjh_rms/err_salt2_yjh_rms

textfile_1.write('SALT2 - $YJH$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    salt2_yjh, err_salt2_yjh, ratio_yjh,
    salt2_yjh_wrms, err_salt2_yjh_wrms, ratio_yjh_wrms,
    salt2_yjh_rms, err_salt2_yjh_rms, ratio_yjh_rms ))

textfile_2.write('SALTminusYJH         %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f \n'%(
    salt2_yjh, err_salt2_yjh, ratio_yjh,
    salt2_yjh_wrms, err_salt2_yjh_wrms, ratio_yjh_wrms,
    salt2_yjh_rms, err_salt2_yjh_rms, ratio_yjh_rms ))

#-------------------
snoopy_yjh     = Snoopy_RMSData[0][0] - YJH_wRMS_GP[0][0]
err_snoopy_yjh = np.sqrt(Snoopy_RMSData[0][1]**2 + YJH_wRMS_GP[0][1]**2)
ratio_yjh = snoopy_yjh/err_snoopy_yjh

snoopy_yjh_wrms     = Snoopy_RMSData[5][0] - YJH_wRMS_GP[5][0]
err_snoopy_yjh_wrms = np.sqrt(Snoopy_RMSData[5][1]**2 + YJH_wRMS_GP[5][1]**2)
ratio_yjh_wrms = snoopy_yjh_wrms/err_snoopy_yjh_wrms

snoopy_yjh_rms     = Snoopy_RMSData[5][2] - YJH_wRMS_GP[5][2]
err_snoopy_yjh_rms = np.sqrt(Snoopy_RMSData[5][3]**2 + YJH_wRMS_GP[5][3]**2)
ratio_yjh_rms = snoopy_yjh_rms/err_snoopy_yjh_rms

textfile_1.write('\snoopy{} - $YJH$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    snoopy_yjh, err_snoopy_yjh, ratio_yjh,
    snoopy_yjh_wrms, err_snoopy_yjh_wrms, ratio_yjh_wrms,
    snoopy_yjh_rms, err_snoopy_yjh_rms, ratio_yjh_rms ))

textfile_2.write('SNOOPYminusYJH       %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f\n'%(
    snoopy_yjh, err_snoopy_yjh, ratio_yjh,
    snoopy_yjh_wrms, err_snoopy_yjh_wrms, ratio_yjh_wrms,
    snoopy_yjh_rms, err_snoopy_yjh_rms, ratio_yjh_rms ))

textfile_1.write(txt_hline)
#======================================

salt2_jhk     = SALT2_RMSData[0][0] - JHK_wRMS_GP[0][0]
err_salt2_jhk = np.sqrt(SALT2_RMSData[0][1]**2 + JHK_wRMS_GP[0][1]**2)
ratio_jhk = salt2_jhk/err_salt2_jhk

salt2_jhk_wrms     = SALT2_RMSData[5][0] - JHK_wRMS_GP[5][0]
err_salt2_jhk_wrms = np.sqrt(SALT2_RMSData[5][1]**2 + JHK_wRMS_GP[5][1]**2)
ratio_jhk_wrms = salt2_jhk_wrms/err_salt2_jhk_wrms

salt2_jhk_rms     = SALT2_RMSData[5][2] - JHK_wRMS_GP[5][2]
err_salt2_jhk_rms = np.sqrt(SALT2_RMSData[5][3]**2 + JHK_wRMS_GP[5][3]**2)
ratio_jhk_rms = salt2_jhk_rms/err_salt2_jhk_rms

textfile_1.write('%% SALT2 - $JHK_s$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    salt2_jhk, err_salt2_jhk, ratio_jhk,
    salt2_jhk_wrms, err_salt2_jhk_wrms, ratio_jhk_wrms,
    salt2_jhk_rms, err_salt2_jhk_rms, ratio_jhk_rms ))

textfile_2.write('# SALTminusJHK         %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f \n'%(
    salt2_jhk, err_salt2_jhk, ratio_jhk,
    salt2_jhk_wrms, err_salt2_jhk_wrms, ratio_jhk_wrms,
    salt2_jhk_rms, err_salt2_jhk_rms, ratio_jhk_rms ))

#-------------------
snoopy_jhk     = Snoopy_RMSData[0][0] - JHK_wRMS_GP[0][0]
err_snoopy_jhk = np.sqrt(Snoopy_RMSData[0][1]**2 + JHK_wRMS_GP[0][1]**2)
ratio_jhk = snoopy_jhk/err_snoopy_jhk

snoopy_jhk_wrms     = Snoopy_RMSData[5][0] - JHK_wRMS_GP[5][0]
err_snoopy_jhk_wrms = np.sqrt(Snoopy_RMSData[5][1]**2 + JHK_wRMS_GP[5][1]**2)
ratio_jhk_wrms = snoopy_jhk_wrms/err_snoopy_jhk_wrms

snoopy_jhk_rms     = Snoopy_RMSData[5][2] - JHK_wRMS_GP[5][2]
err_snoopy_jhk_rms = np.sqrt(Snoopy_RMSData[5][3]**2 + JHK_wRMS_GP[5][3]**2)
ratio_jhk_rms = snoopy_jhk_rms/err_snoopy_jhk_rms

textfile_1.write('%% \snoopy{} - $JHK_s$ & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f & $%.3f \\pm %.3f$ & %.1f \\\\ \n'%(
    snoopy_jhk, err_snoopy_jhk, ratio_jhk,
    snoopy_jhk_wrms, err_snoopy_jhk_wrms, ratio_jhk_wrms,
    snoopy_jhk_rms, err_snoopy_jhk_rms, ratio_jhk_rms ))

textfile_2.write('# SNOOPYminusJHK       %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f %7.3f %7.3f %5.1f\n'%(
    snoopy_jhk, err_snoopy_jhk, ratio_jhk,
    snoopy_jhk_wrms, err_snoopy_jhk_wrms, ratio_jhk_wrms,
    snoopy_jhk_rms, err_snoopy_jhk_rms, ratio_jhk_rms ))

textfile_1.write(txt_hline)
#======================================

# << --------- PART GENERATED WITH THE CODE BELOW ---
###################################################

textfile_1.close(); textfile_2.close();

textfile_1.close();textfile_1.close();textfile_1.close();
textfile_1.close();textfile_1.close();textfile_1.close();

textfile_2.close();textfile_2.close();textfile_2.close();
textfile_2.close();textfile_2.close();textfile_2.close();

#-----------------------------------------------------------------------------80

# ### Code to generate the NIR part of the python script written above.
#
# #### This is very useful!!

# Comment this cell if I'm NOT working on generating the python code to
# be used in the cell above.

# """
band_list = ['Y', 'J','H','K', 'AnyNIR', 'JH', 'YJH', 'JHK']
band_lowercase = ['y', 'j','h','k', 'anynir', 'jh', 'yjh', 'jhk']
band_text_list = ['$Y$', '$J$','$H$','$K_s$',
                  'any $YJHK_s$',
                  '$JH$ ', '$YJH$', '$JHK_s$']

######################################################
#    Create the text file to save the generated python script

textfile_1 = open(DirSavePythonOutput+'Code_Table_Latex_scatter_Opt_vs_NIR_.py','w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_Date_2   = '# On date: %s \n'%text_timenow
text_Author_2 = '# Data table created by: Arturo Avelino \n'
text_script_2 = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
text_line_2= '#'+'-'*60 + '\n'

textfile_1.write(text_line_2)
textfile_1.write(text_Author_2); textfile_1.write(text_Date_2);
textfile_1.write(text_script_2);
textfile_1.write(text_line_2)
textfile_1.write('# \n')

#--------------------------------------------------
for ii in range(len(band_list)):

    bandup = band_list[ii]
    bandlow = band_lowercase[ii]
    bandtxt = band_text_list[ii]

    ## Comment the results from 'JHK' because they produce no sense values.
    if bandup == 'JHK':
        comment1 = '%% '; comment2 = '# ';
    else:
        comment1 = ''; comment2 = '';

    textfile_1.write("salt2_%s     = SALT2_RMSData[0][0] - %s_wRMS_GP[0][0]\n"%(bandlow,bandup))
    textfile_1.write("err_salt2_%s = np.sqrt(SALT2_RMSData[0][1]**2 + %s_wRMS_GP[0][1]**2)\n"%(
        bandlow,bandup))
    textfile_1.write("ratio_%s = salt2_%s/err_salt2_%s\n"%(bandlow,bandlow,bandlow))
    textfile_1.write('\n')

    textfile_1.write("salt2_%s_wrms     = SALT2_RMSData[5][0] - %s_wRMS_GP[5][0]\n"%(bandlow,bandup))
    textfile_1.write("err_salt2_%s_wrms = np.sqrt(SALT2_RMSData[5][1]**2 + %s_wRMS_GP[5][1]**2)\n"%(
        bandlow,bandup))
    textfile_1.write("ratio_%s_wrms = salt2_%s_wrms/err_salt2_%s_wrms\n"%(bandlow,bandlow,bandlow))
    textfile_1.write('\n')

    textfile_1.write("salt2_%s_rms     = SALT2_RMSData[5][2] - %s_wRMS_GP[5][2]\n"%(bandlow,bandup))
    textfile_1.write("err_salt2_%s_rms = np.sqrt(SALT2_RMSData[5][3]**2 + %s_wRMS_GP[5][3]**2)\n"%(
        bandlow,bandup))
    textfile_1.write("ratio_%s_rms = salt2_%s_rms/err_salt2_%s_rms\n"%(bandlow,bandlow,bandlow))
    textfile_1.write('\n')

    textfile_1.write("textfile_1.write('%sSALT2 - %s & $%%.3f \\\\pm %%.3f$ & %%.1f & $%%.3f \\\\pm %%.3f$ & %%.1f \\\n"%(comment1,bandtxt))
    textfile_1.write("& $%.3f \\\\pm %.3f$ & %.1f \\\\\\\\ \\n'%(\n")
    textfile_1.write("    salt2_%s, err_salt2_%s, ratio_%s,\n"%(bandlow,bandlow,bandlow))
    textfile_1.write("    salt2_%s_wrms, err_salt2_%s_wrms, ratio_%s_wrms,\n"%(
        bandlow,bandlow,bandlow))
    textfile_1.write("    salt2_%s_rms, err_salt2_%s_rms, ratio_%s_rms ))\n"%(
        bandlow,bandlow,bandlow))
    textfile_1.write('\n')

    textfile_1.write("textfile_2.write('%sSALTminus%-10s %%7.3f %%7.3f %%5.1f %%7.3f %%7.3f %%5.1f %%7.3f %%7.3f %%5.1f \\n'%%(\n"%(comment2,bandup))
    textfile_1.write("    salt2_%s, err_salt2_%s, ratio_%s,\n"%(bandlow,bandlow,bandlow))
    textfile_1.write("    salt2_%s_wrms, err_salt2_%s_wrms, ratio_%s_wrms,\n"%(
        bandlow,bandlow,bandlow))
    textfile_1.write("    salt2_%s_rms, err_salt2_%s_rms, ratio_%s_rms ))\n"%(
        bandlow,bandlow,bandlow))
    textfile_1.write('\n')

    textfile_1.write("#--------------\n")
    textfile_1.write('\n')
    textfile_1.write("snoopy_%s     = Snoopy_RMSData[0][0] - %s_wRMS_GP[0][0]\n"%(bandlow,bandup))
    textfile_1.write("err_snoopy_%s = np.sqrt(Snoopy_RMSData[0][1]**2 + %s_wRMS_GP[0][1]**2)\n"%(
        bandlow,bandup))
    textfile_1.write("ratio_%s = snoopy_%s/err_snoopy_%s\n"%(bandlow,bandlow,bandlow))
    textfile_1.write('\n')

    textfile_1.write("snoopy_%s_wrms     = Snoopy_RMSData[5][0] - %s_wRMS_GP[5][0]\n"%(bandlow,bandup))
    textfile_1.write("err_snoopy_%s_wrms = np.sqrt(Snoopy_RMSData[5][1]**2 + %s_wRMS_GP[5][1]**2)\n"%(
        bandlow,bandup))
    textfile_1.write("ratio_%s_wrms = snoopy_%s_wrms/err_snoopy_%s_wrms\n"%(bandlow,bandlow,bandlow))
    textfile_1.write('\n')

    textfile_1.write("snoopy_%s_rms     = Snoopy_RMSData[5][2] - %s_wRMS_GP[5][2]\n"%(bandlow,bandup))
    textfile_1.write("err_snoopy_%s_rms = np.sqrt(Snoopy_RMSData[5][3]**2 + %s_wRMS_GP[5][3]**2)\n"%(
        bandlow,bandup))
    textfile_1.write("ratio_%s_rms = snoopy_%s_rms/err_snoopy_%s_rms\n"%(bandlow,bandlow,bandlow))
    textfile_1.write('\n')

    textfile_1.write("textfile_1.write('%s\snoopy{} - %s & $%%.3f \\\\pm %%.3f$ & %%.1f & $%%.3f \\\\pm %%.3f$ & %%.1f \\\n"%(comment1,bandtxt))
    textfile_1.write("& $%.3f \\\\pm %.3f$ & %.1f \\\\\\\\ \\n'%(\n")
    textfile_1.write("    snoopy_%s, err_snoopy_%s, ratio_%s,\n"%(bandlow,bandlow,bandlow))
    textfile_1.write("    snoopy_%s_wrms, err_snoopy_%s_wrms, ratio_%s_wrms,\n"%(
        bandlow,bandlow,bandlow))
    textfile_1.write("    snoopy_%s_rms, err_snoopy_%s_rms, ratio_%s_rms ))\n"%(
        bandlow,bandlow,bandlow))
    textfile_1.write('\n')

    textfile_1.write("textfile_2.write('%sSNOOPYminus%-8s %%7.3f %%7.3f %%5.1f %%7.3f %%7.3f %%5.1f %%7.3f %%7.3f %%5.1f\\n'%%(\n"%(comment2,bandup))
    textfile_1.write("    snoopy_%s, err_snoopy_%s, ratio_%s,\n"%(bandlow,bandlow,bandlow))
    textfile_1.write("    snoopy_%s_wrms, err_snoopy_%s_wrms, ratio_%s_wrms,\n"%(
        bandlow,bandlow,bandlow))
    textfile_1.write("    snoopy_%s_rms, err_snoopy_%s_rms, ratio_%s_rms ))\n"%(
        bandlow,bandlow,bandlow))
    textfile_1.write('\n')

    textfile_1.write('\n')
    textfile_1.write("textfile_1.write(txt_hline)\n")
    textfile_1.write("#======================================\n")
    textfile_1.write('\n')

#-----------------------------------------------------------------------------80
textfile_1.close();

if ScriptVersion == 'notebook':
    print "All done with no issues :)"

# """
0

textfile_1.close();textfile_1.close();textfile_1.close();
textfile_1.close();textfile_1.close();textfile_1.close();

##############################################################################80

# ## Latex Table: NIR-max vs B-max for GP intrinsic scatter (Andy's table)

band_list = ['Y', 'J','H','K', 'AnyNIR', 'JH', 'YJH', 'JHK']

#--------------------------------------------------------60
textfile_1 = open(DirSaveOutput+'Table_Latex_scatter_NIRmax_Bmax_.tex', 'w')
textfile_2 = open(DirSaveOutput+'Table_Latex_scatter_NIRmax_Bmax_.dat', 'w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_Date   = '%% On date: %s \n'%text_timenow

textfile_1.write(text_line)
textfile_1.write(text_Author); textfile_1.write(text_Date);
textfile_1.write(text_script);
textfile_1.write(text_line)
textfile_1.write(' \n')

textfile_2.write('#'+text_line)
textfile_2.write('#'+text_Author); textfile_2.write('#'+text_Date);
textfile_2.write('#'+text_script);
textfile_2.write('#'+text_line)

txt_hline = "\\hline \n"

# textfile_1.write("Optical $BVR$ Method - NIR band(s) & $\\Delta\\sigmaint$ & $n$-$\sigma$ \\\\ \n")
# textfile_1.write(txt_hline)

#--------------------------------------------------------60

# ---- PART GENERATED WITH THE CODE BELOW ---------->>

#-----------------------------------------------------------------------------80
# Data table created by: Arturo Avelino
# On date: 2018-12-18 (yyyy-mm-dd).
# Script used: 15_Latex_scatter.ipynb
#-----------------------------------------------------------------------------80
#
diff_y     = Y_wRMS_GP_tBmax[0][0] - Y_wRMS_GP[0][0]
err_diff_y = np.sqrt(Y_wRMS_GP_tBmax[0][1]**2 + Y_wRMS_GP[0][1]**2)
ratio_y = diff_y/err_diff_y
diff_y_wrms     = Y_wRMS_GP_tBmax[5][0] - Y_wRMS_GP[5][0]
err_diff_y_wrms = np.sqrt(Y_wRMS_GP_tBmax[5][1]**2 + Y_wRMS_GP[5][1]**2)
ratio_y_wrms = diff_y_wrms/err_diff_y_wrms
diff_y_rms     = Y_wRMS_GP_tBmax[5][2] - Y_wRMS_GP[5][2]
err_diff_y_rms = np.sqrt(Y_wRMS_GP_tBmax[5][3]**2 + Y_wRMS_GP[5][3]**2)
ratio_y_rms = diff_y_rms/err_diff_y_rms
textfile_1.write('$Y$ & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f \\\\ \n'%(
    diff_y, err_diff_y, ratio_y,
    diff_y_wrms, err_diff_y_wrms, ratio_y_wrms,
    diff_y_rms, err_diff_y_rms, ratio_y_rms ))
textfile_2.write('Y       %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f\n'%(
    diff_y, err_diff_y, ratio_y,
    diff_y_wrms, err_diff_y_wrms, ratio_y_wrms,
    diff_y_rms, err_diff_y_rms, ratio_y_rms ))

#-------------------
diff_j     = J_wRMS_GP_tBmax[0][0] - J_wRMS_GP[0][0]
err_diff_j = np.sqrt(J_wRMS_GP_tBmax[0][1]**2 + J_wRMS_GP[0][1]**2)
ratio_j = diff_j/err_diff_j
diff_j_wrms     = J_wRMS_GP_tBmax[5][0] - J_wRMS_GP[5][0]
err_diff_j_wrms = np.sqrt(J_wRMS_GP_tBmax[5][1]**2 + J_wRMS_GP[5][1]**2)
ratio_j_wrms = diff_j_wrms/err_diff_j_wrms
diff_j_rms     = J_wRMS_GP_tBmax[5][2] - J_wRMS_GP[5][2]
err_diff_j_rms = np.sqrt(J_wRMS_GP_tBmax[5][3]**2 + J_wRMS_GP[5][3]**2)
ratio_j_rms = diff_j_rms/err_diff_j_rms
textfile_1.write('$J$ & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f \\\\ \n'%(
    diff_j, err_diff_j, ratio_j,
    diff_j_wrms, err_diff_j_wrms, ratio_j_wrms,
    diff_j_rms, err_diff_j_rms, ratio_j_rms ))
textfile_2.write('J       %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f\n'%(
    diff_j, err_diff_j, ratio_j,
    diff_j_wrms, err_diff_j_wrms, ratio_j_wrms,
    diff_j_rms, err_diff_j_rms, ratio_j_rms ))

#-------------------
diff_h     = H_wRMS_GP_tBmax[0][0] - H_wRMS_GP[0][0]
err_diff_h = np.sqrt(H_wRMS_GP_tBmax[0][1]**2 + H_wRMS_GP[0][1]**2)
ratio_h = diff_h/err_diff_h
diff_h_wrms     = H_wRMS_GP_tBmax[5][0] - H_wRMS_GP[5][0]
err_diff_h_wrms = np.sqrt(H_wRMS_GP_tBmax[5][1]**2 + H_wRMS_GP[5][1]**2)
ratio_h_wrms = diff_h_wrms/err_diff_h_wrms
diff_h_rms     = H_wRMS_GP_tBmax[5][2] - H_wRMS_GP[5][2]
err_diff_h_rms = np.sqrt(H_wRMS_GP_tBmax[5][3]**2 + H_wRMS_GP[5][3]**2)
ratio_h_rms = diff_h_rms/err_diff_h_rms
textfile_1.write('$H$ & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f \\\\ \n'%(
    diff_h, err_diff_h, ratio_h,
    diff_h_wrms, err_diff_h_wrms, ratio_h_wrms,
    diff_h_rms, err_diff_h_rms, ratio_h_rms ))
textfile_2.write('H       %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f\n'%(
    diff_h, err_diff_h, ratio_h,
    diff_h_wrms, err_diff_h_wrms, ratio_h_wrms,
    diff_h_rms, err_diff_h_rms, ratio_h_rms ))

#-------------------
diff_k     = K_wRMS_GP_tBmax[0][0] - K_wRMS_GP[0][0]
err_diff_k = np.sqrt(K_wRMS_GP_tBmax[0][1]**2 + K_wRMS_GP[0][1]**2)
ratio_k = diff_k/err_diff_k
diff_k_wrms     = K_wRMS_GP_tBmax[5][0] - K_wRMS_GP[5][0]
err_diff_k_wrms = np.sqrt(K_wRMS_GP_tBmax[5][1]**2 + K_wRMS_GP[5][1]**2)
ratio_k_wrms = diff_k_wrms/err_diff_k_wrms
diff_k_rms     = K_wRMS_GP_tBmax[5][2] - K_wRMS_GP[5][2]
err_diff_k_rms = np.sqrt(K_wRMS_GP_tBmax[5][3]**2 + K_wRMS_GP[5][3]**2)
ratio_k_rms = diff_k_rms/err_diff_k_rms
textfile_1.write('$K_s$ & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f \\\\ \n'%(
    diff_k, err_diff_k, ratio_k,
    diff_k_wrms, err_diff_k_wrms, ratio_k_wrms,
    diff_k_rms, err_diff_k_rms, ratio_k_rms ))
textfile_2.write('K       %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f\n'%(
    diff_k, err_diff_k, ratio_k,
    diff_k_wrms, err_diff_k_wrms, ratio_k_wrms,
    diff_k_rms, err_diff_k_rms, ratio_k_rms ))

#-------------------
diff_anynir     = AnyNIR_wRMS_GP_tBmax[0][0] - AnyNIR_wRMS_GP[0][0]
err_diff_anynir = np.sqrt(AnyNIR_wRMS_GP_tBmax[0][1]**2 + AnyNIR_wRMS_GP[0][1]**2)
ratio_anynir = diff_anynir/err_diff_anynir
diff_anynir_wrms     = AnyNIR_wRMS_GP_tBmax[5][0] - AnyNIR_wRMS_GP[5][0]
err_diff_anynir_wrms = np.sqrt(AnyNIR_wRMS_GP_tBmax[5][1]**2 + AnyNIR_wRMS_GP[5][1]**2)
ratio_anynir_wrms = diff_anynir_wrms/err_diff_anynir_wrms
diff_anynir_rms     = AnyNIR_wRMS_GP_tBmax[5][2] - AnyNIR_wRMS_GP[5][2]
err_diff_anynir_rms = np.sqrt(AnyNIR_wRMS_GP_tBmax[5][3]**2 + AnyNIR_wRMS_GP[5][3]**2)
ratio_anynir_rms = diff_anynir_rms/err_diff_anynir_rms
textfile_1.write('any $YJHK_s$ & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f \\\\ \n'%(
    diff_anynir, err_diff_anynir, ratio_anynir,
    diff_anynir_wrms, err_diff_anynir_wrms, ratio_anynir_wrms,
    diff_anynir_rms, err_diff_anynir_rms, ratio_anynir_rms ))
textfile_2.write('AnyNIR  %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f\n'%(
    diff_anynir, err_diff_anynir, ratio_anynir,
    diff_anynir_wrms, err_diff_anynir_wrms, ratio_anynir_wrms,
    diff_anynir_rms, err_diff_anynir_rms, ratio_anynir_rms ))

#-------------------
diff_jh     = JH_wRMS_GP_tBmax[0][0] - JH_wRMS_GP[0][0]
err_diff_jh = np.sqrt(JH_wRMS_GP_tBmax[0][1]**2 + JH_wRMS_GP[0][1]**2)
ratio_jh = diff_jh/err_diff_jh
diff_jh_wrms     = JH_wRMS_GP_tBmax[5][0] - JH_wRMS_GP[5][0]
err_diff_jh_wrms = np.sqrt(JH_wRMS_GP_tBmax[5][1]**2 + JH_wRMS_GP[5][1]**2)
ratio_jh_wrms = diff_jh_wrms/err_diff_jh_wrms
diff_jh_rms     = JH_wRMS_GP_tBmax[5][2] - JH_wRMS_GP[5][2]
err_diff_jh_rms = np.sqrt(JH_wRMS_GP_tBmax[5][3]**2 + JH_wRMS_GP[5][3]**2)
ratio_jh_rms = diff_jh_rms/err_diff_jh_rms
textfile_1.write('$JH$  & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f \\\\ \n'%(
    diff_jh, err_diff_jh, ratio_jh,
    diff_jh_wrms, err_diff_jh_wrms, ratio_jh_wrms,
    diff_jh_rms, err_diff_jh_rms, ratio_jh_rms ))
textfile_2.write('JH      %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f\n'%(
    diff_jh, err_diff_jh, ratio_jh,
    diff_jh_wrms, err_diff_jh_wrms, ratio_jh_wrms,
    diff_jh_rms, err_diff_jh_rms, ratio_jh_rms ))

#-------------------
diff_yjh     = YJH_wRMS_GP_tBmax[0][0] - YJH_wRMS_GP[0][0]
err_diff_yjh = np.sqrt(YJH_wRMS_GP_tBmax[0][1]**2 + YJH_wRMS_GP[0][1]**2)
ratio_yjh = diff_yjh/err_diff_yjh
diff_yjh_wrms     = YJH_wRMS_GP_tBmax[5][0] - YJH_wRMS_GP[5][0]
err_diff_yjh_wrms = np.sqrt(YJH_wRMS_GP_tBmax[5][1]**2 + YJH_wRMS_GP[5][1]**2)
ratio_yjh_wrms = diff_yjh_wrms/err_diff_yjh_wrms
diff_yjh_rms     = YJH_wRMS_GP_tBmax[5][2] - YJH_wRMS_GP[5][2]
err_diff_yjh_rms = np.sqrt(YJH_wRMS_GP_tBmax[5][3]**2 + YJH_wRMS_GP[5][3]**2)
ratio_yjh_rms = diff_yjh_rms/err_diff_yjh_rms
textfile_1.write('$YJH$ & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f \\\\ \n'%(
    diff_yjh, err_diff_yjh, ratio_yjh,
    diff_yjh_wrms, err_diff_yjh_wrms, ratio_yjh_wrms,
    diff_yjh_rms, err_diff_yjh_rms, ratio_yjh_rms ))
textfile_2.write('YJH     %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f\n'%(
    diff_yjh, err_diff_yjh, ratio_yjh,
    diff_yjh_wrms, err_diff_yjh_wrms, ratio_yjh_wrms,
    diff_yjh_rms, err_diff_yjh_rms, ratio_yjh_rms ))

#-------------------
diff_jhk     = JHK_wRMS_GP_tBmax[0][0] - JHK_wRMS_GP[0][0]
err_diff_jhk = np.sqrt(JHK_wRMS_GP_tBmax[0][1]**2 + JHK_wRMS_GP[0][1]**2)
ratio_jhk = diff_jhk/err_diff_jhk
diff_jhk_wrms     = JHK_wRMS_GP_tBmax[5][0] - JHK_wRMS_GP[5][0]
err_diff_jhk_wrms = np.sqrt(JHK_wRMS_GP_tBmax[5][1]**2 + JHK_wRMS_GP[5][1]**2)
ratio_jhk_wrms = diff_jhk_wrms/err_diff_jhk_wrms
diff_jhk_rms     = JHK_wRMS_GP_tBmax[5][2] - JHK_wRMS_GP[5][2]
err_diff_jhk_rms = np.sqrt(JHK_wRMS_GP_tBmax[5][3]**2 + JHK_wRMS_GP[5][3]**2)
ratio_jhk_rms = diff_jhk_rms/err_diff_jhk_rms
textfile_1.write('%% $JHK_s$ & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f & $%.3f \\pm %.3f$ & %.2f \\\\ \n'%(
    diff_jhk, err_diff_jhk, ratio_jhk,
    diff_jhk_wrms, err_diff_jhk_wrms, ratio_jhk_wrms,
    diff_jhk_rms, err_diff_jhk_rms, ratio_jhk_rms ))
textfile_2.write('# JHK     %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f %7.3f %7.3f %7.2f\n'%(
    diff_jhk, err_diff_jhk, ratio_jhk,
    diff_jhk_wrms, err_diff_jhk_wrms, ratio_jhk_wrms,
    diff_jhk_rms, err_diff_jhk_rms, ratio_jhk_rms ))

#-------------------
# << --------- PART GENERATED WITH THE CODE BELOW ---
###################################################

textfile_1.close(); textfile_2.close();

textfile_1.close();textfile_1.close();textfile_1.close();
textfile_1.close();textfile_1.close();textfile_1.close();

textfile_2.close();textfile_2.close();textfile_2.close();
textfile_2.close();textfile_2.close();textfile_2.close();

#-----------------------------------------------------------------------------80

# ### Code to generate the NIR part of the python script written above.
#
# #### This is very useful!!

# Comment this cell if I'm NOT working on generating the python code to
# be used in the cell above.

# """
band_list = ['Y', 'J','H','K', 'AnyNIR', 'JH', 'YJH', 'JHK']
band_lowercase = ['y', 'j','h','k', 'anynir', 'jh', 'yjh', 'jhk']
band_text_list = ['$Y$', '$J$','$H$','$K_s$',
                  'any $YJHK_s$',
                  '$JH$ ', '$YJH$', '$JHK_s$']

######################################################
#    Create the text file to save the generated python script

textfile_1 = open(DirSavePythonOutput+'Code_Table_Latex_scatter_NIRmax_Bmax_.py','w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_Date_2   = '# On date: %s \n'%text_timenow
text_Author_2 = '# Data table created by: Arturo Avelino \n'
text_script_2 = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
text_line_2= '#'+'-'*60 + '\n'

textfile_1.write(text_line_2)
textfile_1.write(text_Author_2); textfile_1.write(text_Date_2);
textfile_1.write(text_script_2);
textfile_1.write(text_line_2)
textfile_1.write('# \n')

#-----------------------------------------------------------------------------80

for ii in range(len(band_list)):

    bandup = band_list[ii]
    bandlow = band_lowercase[ii]
    bandtxt = band_text_list[ii]

    ## Comment the results from 'JHK' because they produce no sense values.
    if bandup == 'JHK':
        comment1 = '%% '; comment2 = '# ';
    else:
        comment1 = ''; comment2 = '';

    textfile_1.write("diff_%s     = %s_wRMS_GP_tBmax[0][0] - %s_wRMS_GP[0][0]\n"%(
        bandlow,bandup,bandup))
    textfile_1.write("err_diff_%s = np.sqrt(%s_wRMS_GP_tBmax[0][1]**2 + %s_wRMS_GP[0][1]**2)\n"%(bandlow,bandup,bandup))
    textfile_1.write("ratio_%s = diff_%s/err_diff_%s\n"%(bandlow,bandlow,bandlow))

    textfile_1.write("diff_%s_wrms     = %s_wRMS_GP_tBmax[5][0] - %s_wRMS_GP[5][0]\n"%(bandlow,bandup,bandup))
    textfile_1.write("err_diff_%s_wrms = np.sqrt(%s_wRMS_GP_tBmax[5][1]**2 + %s_wRMS_GP[5][1]**2)\n"%(bandlow,bandup,bandup))
    textfile_1.write("ratio_%s_wrms = diff_%s_wrms/err_diff_%s_wrms\n"%(
        bandlow,bandlow,bandlow))

    textfile_1.write("diff_%s_rms     = %s_wRMS_GP_tBmax[5][2] - %s_wRMS_GP[5][2]\n"%(bandlow,bandup,bandup))
    textfile_1.write("err_diff_%s_rms = np.sqrt(%s_wRMS_GP_tBmax[5][3]**2 + %s_wRMS_GP[5][3]**2)\n"%(bandlow,bandup,bandup))
    textfile_1.write("ratio_%s_rms = diff_%s_rms/err_diff_%s_rms\n"%(
        bandlow,bandlow,bandlow))

    textfile_1.write("textfile_1.write('%s%s & $%%.3f \\\\pm %%.3f$ & %%.2f & $%%.3f \\\\pm %%.3f$ & %%.2f & \\\n"%(comment1,bandtxt))
    textfile_1.write("$%.3f \\\\pm %.3f$ & %.2f \\\\\\\\ \\n'%(\n")
    textfile_1.write("    diff_%s, err_diff_%s, ratio_%s,\n"%(
        bandlow,bandlow,bandlow))
    textfile_1.write("    diff_%s_wrms, err_diff_%s_wrms, ratio_%s_wrms, \n"%(
        bandlow,bandlow,bandlow))
    textfile_1.write("    diff_%s_rms, err_diff_%s_rms, ratio_%s_rms ))\n"%(
        bandlow,bandlow,bandlow))

    textfile_1.write("textfile_2.write('%s%-7s %%7.3f %%7.3f %%7.2f %%7.3f %%7.3f %%7.2f %%7.3f %%7.3f %%7.2f\\n'%%(\n"%(comment2,bandup))
    textfile_1.write("    diff_%s, err_diff_%s, ratio_%s,\n"%(
        bandlow,bandlow,bandlow))
    textfile_1.write("    diff_%s_wrms, err_diff_%s_wrms, ratio_%s_wrms, \n"%(
        bandlow,bandlow,bandlow))
    textfile_1.write("    diff_%s_rms, err_diff_%s_rms, ratio_%s_rms ))\n"%(
        bandlow,bandlow,bandlow))


    textfile_1.write('\n')
    # textfile_1.write("textfile_1.write(txt_hline)\n")
    textfile_1.write("#------------------------------------\n\n")

#-----------------------------------------------------------------------------80

textfile_1.close()

if ScriptVersion == 'notebook':
    print "All done with no issues :)"

# """
# 0

textfile_1.close();textfile_1.close();textfile_1.close();
textfile_1.close();textfile_1.close();textfile_1.close();

##############################################################################80

# # Latex macros
#
# Create the macros definitios used in the macros.tex file of the latex version of the paper

# Read scatter files

try:
    scatter_GPNIRmax_np = np.genfromtxt(DirSaveOutput+
                                'Table_Latex_scatter_GPNIRmax_Notes.dat',
                                dtype=['S10', 'S14', int,float,float,float,
                                      float,float,float,float,float])
except:
    scatter_GPNIRmax_np = np.genfromtxt(DirSaveOutput+
                                'Table_Latex_scatter_GPNIRmax.dat',
                                dtype=['S10', 'S14', int,float,float,float,
                                      float,float,float,float,float])

#-----------------------------------------------------------------------------80

try:
    scatter_GPBmax_np = np.genfromtxt(DirSaveOutput+
                                'Table_Latex_scatter_GPBmax_Notes.dat',
                                dtype=['S10', 'S14', int,float,float,float,
                                      float,float,float,float,float])
except:
    scatter_GPBmax_np = np.genfromtxt(DirSaveOutput+
                                'Table_Latex_scatter_GPBmax.dat',
                                dtype=['S10', 'S14', int,float,float,float,
                                      float,float,float,float,float])

#-----------------------------------------------------------------------------80

try:
    scatter_Temp_np = np.genfromtxt(DirSaveOutput+
                                'Table_Latex_scatter_Template_Notes.dat',
                                dtype=['S10', 'S14', int,float,float,float,
                                      float,float,float,float,float])
except:
    scatter_Temp_np = np.genfromtxt(DirSaveOutput+
                                'Table_Latex_scatter_Template.dat',
                                dtype=['S10', 'S14', int,float,float,float,
                                      float,float,float,float,float])
#-----------------------------------------------------------------------------80

try:
    scatter_OptvsNIR_np = np.genfromtxt(DirSaveOutput+
                                'Table_Latex_scatter_Opt_vs_NIR_Notes.dat',
                                dtype=['S20', float,float,float,float,
                                      float,float,float,float,float])
except:
    scatter_OptvsNIR_np = np.genfromtxt(DirSaveOutput+
                                'Table_Latex_scatter_Opt_vs_NIR_.dat',
                                dtype=['S20', float,float,float,float,
                                      float,float,float,float,float])

#-----------------------------------------------------------------------------80

try:
    scatter_NIRmaxBmax_np = np.genfromtxt(DirSaveOutput+
                                'Table_Latex_scatter_NIRmax_Bmax_Notes.dat',
                                dtype=['S12', float,float,float,float,
                                      float,float,float,float,float])
except:
    scatter_NIRmaxBmax_np = np.genfromtxt(DirSaveOutput+
                                'Table_Latex_scatter_NIRmax_Bmax_.dat',
                                dtype=['S12', float,float,float,float,
                                      float,float,float,float,float])

# Create the macros latex file.
#-----------------------------------------------------------------------------80

textfile_1 = open(DirSaveOutput+'macros_update.tex','w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd); %H:%M hrs.")
text_Date   = '%% On date: %s \n'%text_timenow
text_line_short = '%'+'-'*40 + '\n'

textfile_1.write(text_line_short)
textfile_1.write(text_Author); textfile_1.write(text_Date);
textfile_1.write(text_script);
textfile_1.write(text_line_short)
textfile_1.write(' \n')

textfile_1.write('%      Update the values defined in the macros file \n\n')

#--------------------------------------------------
textfile_1.write('%  Number of SNe in different situations \n')

textfile_1.write('\\newcommand{\\nsnIa}{%.0f} \n'%(AnyNIR_wRMS_Tp[3][1]))
textfile_1.write('\\newcommand{\\nY}{$%.0f$} \n'%(Y_wRMS_Tp[3][1]))
textfile_1.write('\\newcommand{\\nJ}{$%.0f$} \n'%(J_wRMS_Tp[3][1]))
textfile_1.write('\\newcommand{\\nH}{$%.0f$} \n'%(H_wRMS_Tp[3][1]))
textfile_1.write('\\newcommand{\\nK}{$%.0f$} \n'%(K_wRMS_Tp[3][1]))
textfile_1.write('\\newcommand{\\nKgp}{%.0f} \n'%(K_wRMS_GP[3][1]))

#-------------------
textfile_1.write('\n%  Bands definitions: \n')

textfile_1.write('\\newcommand{\\Y}{$Y$}\n')
textfile_1.write('\\newcommand{\\J}{$J$}\n')
textfile_1.write('\\newcommand{\\HH}{$H$}\n')
textfile_1.write('\\newcommand{\\K}{$K_s$}\n')
textfile_1.write('\\newcommand{\\AnyNIR}{any $YJHK_s$}\n')
textfile_1.write('\\newcommand{\\JH}{$JH$}\n')
textfile_1.write('\\newcommand{\\YJH}{$YJH$}\n')
textfile_1.write('\\newcommand{\\JHK}{$JHK_s$}\n')

textfile_1.write('\\newcommand{\\SALTminusY}{SALT2 $-Y$}\n')
textfile_1.write('\\newcommand{\\SALTminusJ}{SALT2 $-J$}\n')
textfile_1.write('\\newcommand{\\SALTminusH}{SALT2 $-H$}\n')
textfile_1.write('\\newcommand{\\SALTminusK}{SALT2 $-K_s$}\n')
textfile_1.write('\\newcommand{\\SALTminusAnyNIR}{SALT2 $-$ any $YJHK_s$}\n')
textfile_1.write('\\newcommand{\\SALTminusJH}{SALT2 $-JH$}\n')
textfile_1.write('\\newcommand{\\SALTminusYJH}{SALT2 $-YJH$}\n')
textfile_1.write('\\newcommand{\\SALTminusJHK}{SALT2 $-YJHK_s$}\n')

textfile_1.write('\\newcommand{\\SNOOPYminusY}{\snoopy$-Y$}\n')
textfile_1.write('\\newcommand{\\SNOOPYminusJ}{\snoopy$-J$}\n')
textfile_1.write('\\newcommand{\\SNOOPYminusH}{\snoopy$-H$}\n')
textfile_1.write('\\newcommand{\\SNOOPYminusK}{\snoopy$-K_s$}\n')
textfile_1.write('\\newcommand{\\SNOOPYminusAnyNIR}{\snoopy$-$any $YJHK_s$}\n')
textfile_1.write('\\newcommand{\\SNOOPYminusJH}{\snoopy$-JH$}\n')
textfile_1.write('\\newcommand{\\SNOOPYminusYJH}{\snoopy$-YJH$}\n')
textfile_1.write('\\newcommand{\\SNOOPYminusJHK}{\snoopy$-YJHK_s$}\n')

#-------------------
textfile_1.write('\n%  Intrinsic dispersion \n')

textfile_1.write('\\newcommand{\\nAnyYJHKgp}{%.0f} \n'%(AnyNIR_wRMS_GP[3][1]))

textfile_1.write('\\newcommand{\\sigmaIntAnyYJHKgp}{%.3f \\pm %.3f} \n'%(
    AnyNIR_wRMS_GP[0][0], AnyNIR_wRMS_GP[0][1]))

textfile_1.write('\\newcommand{\\sigmaIntAnyYJHKgpAtTBmaxGPsubsample}{%.3f \\pm %.3f} \n'%(
    AnyNIR_wRMS_GP_tBmax[0][0], AnyNIR_wRMS_GP_tBmax[0][1]))

textfile_1.write('\\newcommand{\\nAnyYJHKtemp}{%.0f} \n'%(AnyNIR_wRMS_Tp[3][1]))

textfile_1.write('\\newcommand{\\sigmaIntAnyYJHKTemp}{%.3f \\pm %.3f} \n'%(
    AnyNIR_wRMS_Tp[0][0], AnyNIR_wRMS_Tp[0][1]))

textfile_1.write('\\newcommand{\\sigmaIntAnyYJHKTempGPsubsample}{%.3f \\pm %.3f} \n'%(
    AnyNIR_wRMS_Tp_GPsubsample[0][0], AnyNIR_wRMS_Tp_GPsubsample[0][1]))

textfile_1.write('\\newcommand{\\sigmaIntSALT}{%.3f \\pm %.3f} \n'%(
    SALT2_RMSData[0][0], SALT2_RMSData[0][1]))

textfile_1.write('\\newcommand{\\sigmaIntSnoopy}{%.3f \\pm %.3f} \n'%(
    Snoopy_RMSData[0][0], Snoopy_RMSData[0][1]))

textfile_1.write('\\newcommand{\\sigmaIntHgp}{%.3f \\pm %.3f} \n'%(
    H_wRMS_GP[0][0], H_wRMS_GP[0][1]))

textfile_1.write('\\newcommand{\\sigmaIntYtemp}{%.3f \\pm %.3f} \n'%(
    Y_wRMS_Tp[0][0], Y_wRMS_Tp[0][1]))

#-------------------
textfile_1.write('\n%  wRMS \n')

textfile_1.write('\\newcommand{\\wrmsAnyYJHKgp}{%.3f \\pm %.3f} \n'%(
    AnyNIR_wRMS_GP[5][0], AnyNIR_wRMS_GP[5][1]))

textfile_1.write('\\newcommand{\\wrmsAnyYJHKgpAtTBmaxGPsubsample}{%.3f \\pm %.3f} \n'%(
    AnyNIR_wRMS_GP_tBmax[5][0], AnyNIR_wRMS_GP_tBmax[5][1]))

textfile_1.write('\\newcommand{\\wrmsAnyYJHKTempGPsubsample}{%.3f \\pm %.3f} \n'%(
    AnyNIR_wRMS_Tp_GPsubsample[5][0], AnyNIR_wRMS_Tp_GPsubsample[5][1]))

textfile_1.write('\\newcommand{\\wrmsSALT}{%.3f \\pm %.3f} \n'%(
    SALT2_RMSData[5][0], SALT2_RMSData[5][1]))

textfile_1.write('\\newcommand{\\wrmsSnoopy}{%.3f \\pm %.3f} \n'%(
    Snoopy_RMSData[5][0], Snoopy_RMSData[5][1]))

textfile_1.write('\\newcommand{\\wrmsHgp}{%.3f \\pm %.3f} \n'%(
    H_wRMS_GP[5][0], H_wRMS_GP[5][1]))

textfile_1.write('\\newcommand{\\wrmsYtemp}{%.3f \\pm %.3f} \n'%(
    Y_wRMS_Tp[5][0], Y_wRMS_Tp[5][1]))

##############################################################################80

textfile_1.write('\n%  RMS \n')

textfile_1.write('\\newcommand{\\rmsSALT}{%.3f \\pm %.3f} \n'%(
    SALT2_RMSData[5][2],SALT2_RMSData[5][3]))

textfile_1.write('\\newcommand{\\rmsSnoopy}{%.3f \\pm %.3f} \n'%(
    Snoopy_RMSData[5][2],Snoopy_RMSData[5][3]))

textfile_1.write('\\newcommand{\\rmsAnyYJHKgp}{%.3f \\pm %.3f} \n'%(
    AnyNIR_wRMS_GP[5][2], AnyNIR_wRMS_GP[5][3]))

textfile_1.write('\\newcommand{\\rmsAnyYJHKgpBmax}{%.3f \\pm %.3f} \n'%(
    AnyNIR_wRMS_GP_tBmax[5][2], AnyNIR_wRMS_GP_tBmax[5][3]))

textfile_1.write('\\newcommand{\\rmsAnyYJHKTempGPsample}{%.3f \\pm %.3f} \n'%(
    AnyNIR_wRMS_Tp_GPsubsample[5][2], AnyNIR_wRMS_Tp_GPsubsample[5][3]))

#-------------------
textfile_1.write("\n")
textfile_1.write('% Smallest value of the RMS in all NIR bands using the GP at NIR_max:\n')
textfile_1.write('\\newcommand{\\rmsGPNIRmaxSmallest}{%.3f \\pm %.3f} \n'%(
    min(scatter_GPNIRmax_np['f9']),
    scatter_GPNIRmax_np['f10'][np.argmin(scatter_GPNIRmax_np['f9'])] ))
textfile_1.write('% Band of the smallest value of the RMS in all NIR bands using the GP at NIR_max:\n')
textfile_1.write('\\newcommand{\\rmsGPNIRmaxSmallestBand}{\\%s} \n'%(
    scatter_GPNIRmax_np['f0'][np.argmin(scatter_GPNIRmax_np['f9'])] ))

textfile_1.write('% Largest value of the RMS in all NIR bands using the GP at NIR_max:\n')
textfile_1.write('\\newcommand{\\rmsGPNIRmaxLargest}{%.3f \\pm %.3f} \n'%(
    max(scatter_GPNIRmax_np['f9']),
    scatter_GPNIRmax_np['f10'][np.argmax(scatter_GPNIRmax_np['f9'])]  ))
textfile_1.write('% Band of the largest value of the RMS in all NIR bands using the GP at NIR_max:\n')
textfile_1.write('\\newcommand{\\rmsGPNIRmaxLargestBand}{\\%s} \n'%(
    scatter_GPNIRmax_np['f0'][np.argmax(scatter_GPNIRmax_np['f9'])] ))

#-------------------
textfile_1.write("\n")
textfile_1.write('% Smallest value of the RMS in individual bands using the GP at NIR_max:\n')
textfile_1.write('\\newcommand{\\rmsGPNIRmaxOneBandSmallest}{%.3f \\pm %.3f} \n'%(
    min(scatter_GPNIRmax_np['f9'][:4]),
    scatter_GPNIRmax_np['f10'][:4][np.argmin(scatter_GPNIRmax_np['f9'][:4])]  ))
textfile_1.write('% Band of the smallest value of the RMS in individual bands using the GP at NIR_max:\n')
textfile_1.write('\\newcommand{\\rmsGPNIRmaxOneBandSmallestBand}{\\%s} \n'%(
    scatter_GPNIRmax_np['f0'][:4][np.argmin(scatter_GPNIRmax_np['f9'][:4])] ))

textfile_1.write('% Largest value of the RMS in individual bands using the GP at NIR_max:\n')
textfile_1.write('\\newcommand{\\rmsGPNIRmaxOneBandLargest}{%.3f \\pm %.3f} \n'%(
    max(scatter_GPNIRmax_np['f9'][:4]),
    scatter_GPNIRmax_np['f10'][:4][np.argmax(scatter_GPNIRmax_np['f9'][:4])]  ))
textfile_1.write('% Band of the largest value of the RMS in individual bands using the GP at NIR_max:\n')
textfile_1.write('\\newcommand{\\rmsGPNIRmaxOneBandLargestBand}{\\%s} \n'%(
    scatter_GPNIRmax_np['f0'][:4][np.argmax(scatter_GPNIRmax_np['f9'][:4])] ))

#-------------------
textfile_1.write("\n")
textfile_1.write('% Smallest value of the RMS in all NIR bands using the template method:\n')
textfile_1.write('\\newcommand{\\rmsTempSmallest}{%.3f \\pm %.3f} \n'%(
    min(scatter_Temp_np['f9']),
    scatter_Temp_np['f10'][np.argmin(scatter_Temp_np['f9'])] ))
textfile_1.write('% Band of the smallest value of the RMS in all NIR bands using the template method:\n')
textfile_1.write('\\newcommand{\\rmsTempSmallestBand}{\\%s} \n'%(
    scatter_Temp_np['f0'][np.argmin(scatter_Temp_np['f9'])] ))

textfile_1.write('% Largest value of the RMS in all NIR bands using the template method:\n')
textfile_1.write('\\newcommand{\\rmsTempLargest}{%.3f \\pm %.3f} \n'%(
    max(scatter_Temp_np['f9']),
    scatter_Temp_np['f10'][np.argmax(scatter_Temp_np['f9'])] ))
textfile_1.write('% Band of the largest value of the RMS in all NIR bands using the template method:\n')
textfile_1.write('\\newcommand{\\rmsTempLargestBand}{\\%s} \n'%(
    scatter_Temp_np['f0'][np.argmax(scatter_Temp_np['f9'])] ))

#-------------------
textfile_1.write("\n")
textfile_1.write('% Smallest value of the RMS in individual bands using the template method:\n')
textfile_1.write('\\newcommand{\\rmsTempOneBandSmallest}{%.3f \\pm %.3f} \n'%(
    min(scatter_Temp_np['f9'][:4]),
    scatter_Temp_np['f10'][:4][np.argmin(scatter_Temp_np['f9'][:4])]  ))
textfile_1.write('% Band of the smallest value of the RMS in individual bands using the template method:\n')
textfile_1.write('\\newcommand{\\rmsTempOneBandSmallestBand}{\\%s} \n'%(
    scatter_Temp_np['f0'][:4][np.argmin(scatter_Temp_np['f9'][:4])] ))

textfile_1.write('% Largest value of the RMS in individual bands using the template method:\n')
textfile_1.write('\\newcommand{\\rmsTempOneBandLargest}{%.3f \\pm %.3f} \n'%(
    max(scatter_Temp_np['f9'][:4]),
    scatter_Temp_np['f10'][:4][np.argmax(scatter_Temp_np['f9'][:4])]  ))
textfile_1.write('% Band of the largest value of the RMS in individual bands using the template method:\n')
textfile_1.write('\\newcommand{\\rmsTempOneBandLargestBand}{\\%s} \n'%(
    scatter_Temp_np['f0'][:4][np.argmax(scatter_Temp_np['f9'][:4])] ))

#-------------------
textfile_1.write("\n")
textfile_1.write('% Smallest value of the RMS in all NIR bands using the GP at B_max:\n')
textfile_1.write('\\newcommand{\\rmsGPBmaxSmallest}{%.3f \\pm %.3f} \n'%(
    min(scatter_GPBmax_np['f9']),
    scatter_GPBmax_np['f10'][np.argmin(scatter_GPBmax_np['f9'])] ))
textfile_1.write('% Band of the smallest value of the RMS in all NIR bands using the GP at B_max:\n')
textfile_1.write('\\newcommand{\\rmsGPBmaxSmallestBand}{\\%s} \n'%(
    scatter_GPBmax_np['f0'][np.argmin(scatter_GPBmax_np['f9'])] ))

textfile_1.write('% Largest value of the RMS in all NIR bands using the GP at B_max:\n')
textfile_1.write('\\newcommand{\\rmsGPBmaxLargest}{%.3f \\pm %.3f} \n'%(
    max(scatter_GPBmax_np['f9']),
    scatter_GPBmax_np['f10'][np.argmax(scatter_GPBmax_np['f9'])] ))
textfile_1.write('% Band of the largest value of the RMS in all NIR bands using the GP at B_max:\n')
textfile_1.write('\\newcommand{\\rmsGPBmaxLargestBand}{\\%s} \n'%(
    scatter_GPBmax_np['f0'][np.argmax(scatter_GPBmax_np['f9'])] ))

#-------------------
textfile_1.write("\n")
textfile_1.write('% Smallest value of the RMS in individual bands using the GP at B_max:\n')
textfile_1.write('\\newcommand{\\rmsGPBmaxOneBandSmallest}{%.3f \\pm %.3f} \n'%(
    min(scatter_GPBmax_np['f9'][:4]),
    scatter_GPBmax_np['f10'][:4][np.argmin(scatter_GPBmax_np['f9'][:4])]  ))
textfile_1.write('% Band of the smallest value of the RMS in individual bands using the GP at B_max:\n')
textfile_1.write('\\newcommand{\\rmsGPBmaxOneBandSmallestBand}{\\%s} \n'%(
    scatter_GPBmax_np['f0'][:4][np.argmin(scatter_GPBmax_np['f9'][:4])] ))

textfile_1.write('% Largest value of the RMS in individual bands using the GP at B_max:\n')
textfile_1.write('\\newcommand{\\rmsGPBmaxOneBandLargest}{%.3f \\pm %.3f} \n'%(
    max(scatter_GPBmax_np['f9'][:4]),
    scatter_GPBmax_np['f10'][:4][np.argmax(scatter_GPBmax_np['f9'][:4])]  ))
textfile_1.write('% Band of the largest value of the RMS in individual bands using the GP at B_max:\n')
textfile_1.write('\\newcommand{\\rmsGPBmaxOneBandLargestBand}{\\%s} \n'%(
    scatter_GPBmax_np['f0'][:4][np.argmax(scatter_GPBmax_np['f9'][:4])] ))

#-------------------
# textfile_1.write("\n")
# textfile_1.write('% RMS of Y-band in the template method:\n')
# textfile_1.write('\\newcommand{\\rmsgpnirYonly}{%.3f \\pm %.3f} \n'%(
#     Y_wRMS_Tp[5][2], Y_wRMS_Tp[5][3]))

#-----------------------------------------------------------------------------80

textfile_1.write('\n%--- Comparing dispersion between Optical and NIR ---\n\n')

textfile_1.write("% Largest value of 'n-sigma' of SALT2/SNooPy - NIR, for intrinsic scatter:\n")
textfile_1.write('\\newcommand{\\optminusnirsigmabest}{%.1f}\n'%(
    max(scatter_OptvsNIR_np['f3']) ) )
textfile_1.write("% Band of the largest value of 'n-sigma' of SALT2/SNooPy - NIR, for intrinsic scatter:\n")
textfile_1.write('\\newcommand{\\optminusnirsigmabestBand}{\\%s}\n'%(
    scatter_OptvsNIR_np['f0'][np.argmax(scatter_OptvsNIR_np['f3'])] ) )

textfile_1.write("\n%     SALT2 - any YJHK\n\n")

textfile_1.write('\\newcommand{\\saltminusnirgp}{%.3f \\pm %.3f} \n'%(
    salt2_anynir, err_salt2_anynir))

textfile_1.write('\\newcommand{\\wrmsdiffSALTAnyYJHK}{%.3f \\pm %.3f} \n'%(
    salt2_anynir_wrms, err_salt2_anynir_wrms))

textfile_1.write('\\newcommand{\\saltminusnirgpAnyYJHKrms}{%.3f \\pm %.3f} \n'%(
    salt2_anynir_rms, err_salt2_anynir_rms))

textfile_1.write("\n")
textfile_1.write("% 'n-sigma' of SALT2 - anyYJHK, for RMS.\n")
textfile_1.write('\\newcommand{\\saltminusnirgpAnyYJHKrmsNS}{%.1f}\n'%(
    (salt2_anynir_rms/err_salt2_anynir_rms) ) )

#-------------------
textfile_1.write("\n%     SNooPy - any YJHK\n\n")

textfile_1.write('\\newcommand{\\snoopyminusnirgp}{%.3f \\pm %.3f} \n'%(
    snoopy_anynir, err_snoopy_anynir))

textfile_1.write('\\newcommand{\\wrmsdiffSnoopyAnyYJHK}{%.3f \\pm %.3f} \n'%(
    snoopy_anynir_wrms, err_snoopy_anynir_wrms))

textfile_1.write('\\newcommand{\\snoopyminusnirgpAnyYJHKrms}{%.3f \\pm %.3f} \n'%(
    snoopy_anynir_rms, err_snoopy_anynir_rms))

textfile_1.write("\n")
textfile_1.write("% 'n-sigma' of SNooPy - anyYJHK, for RMS.\n")
textfile_1.write('\\newcommand{\\snoopyminusnirgpAnyYJHKrmsNS}{%.1f}\n'%(
    (snoopy_anynir_rms/err_snoopy_anynir_rms) ) )

#-------------------
textfile_1.write('\n')
textfile_1.write("% Minimum and maximum values of 'n-sigma' of 'SALT2/SNooPy - NIR', for the intrinsic scatter,\n")
textfile_1.write("% where NIR = any YJHK, JH, YJH, only.\n")
textfile_1.write('\\newcommand{\\optminusnirsigma}{$\\sim %.1f$-$%.1f\\sigma$}\n'%(
    min(scatter_OptvsNIR_np['f3'][8:]), max(scatter_OptvsNIR_np['f3'][8:])) )

#-------------------
textfile_1.write('\n')
# Create a new array where I ignore K band:
array1_np = np.concatenate((scatter_OptvsNIR_np['f3'][:5],
                            scatter_OptvsNIR_np['f3'][8:]), axis=None)
textfile_1.write("% Minimum and maximum values of 'n-sigma' of 'SALT2/SNooPy - NIR', for the intrinsic scatter,\n")
textfile_1.write("% where NIR = Y, J, H, (K excluded), any YJHK, JH, YJH.\n")
textfile_1.write('\\newcommand{\\optminusnirsigmaGeneral}{$\\sim %.1f$-$%.1f\\sigma$}\n'%(
    min(array1_np), max(array1_np) ) )

#-------------------
# Create a new array where I ignore K band:
array2_np = np.concatenate((scatter_OptvsNIR_np['f9'][:5],
                            scatter_OptvsNIR_np['f9'][8:]), axis=None)
array2_bands_np = np.concatenate((scatter_OptvsNIR_np['f0'][:5],
                            scatter_OptvsNIR_np['f0'][8:]), axis=None)
array2_dRMS_np = np.concatenate((scatter_OptvsNIR_np['f7'][:5],
                            scatter_OptvsNIR_np['f7'][8:]), axis=None)
array2_edRMS_np = np.concatenate((scatter_OptvsNIR_np['f8'][:5],
                            scatter_OptvsNIR_np['f8'][8:]), axis=None)

textfile_1.write("\n")
textfile_1.write("% Smallest value of 'n-sigma' of 'SALT2/SNooPy - NIR' for RMS,\n")
textfile_1.write("% where NIR = Y, J, H, (K excluded), any YJHK, JH, YJH.\n")
textfile_1.write('\\newcommand{\\optminusnirRMSnsSmallest}{%.1f} \n'%(min(array2_np) ) )
textfile_1.write("% Band of the smallest value of  'n-sigma' of 'SALT2/SNooPy - NIR' for RMS:\n")
textfile_1.write('\\newcommand{\\optminusnirRMSnsSmallestBand}{\\%s} \n'%(
    array2_bands_np[np.argmin(array2_np)] ) )

textfile_1.write("\n")
textfile_1.write("% Largest value of 'n-sigma' of 'SALT2/SNooPy - NIR' for RMS,\n")
textfile_1.write("% where NIR = Y, J, H, (K excluded), any YJHK, JH, YJH.\n")
# Create a new array where I ignore K band:
textfile_1.write('\\newcommand{\\optminusnirRMSnsLargest}{%.1f} \n'%(max(array2_np) ) )
textfile_1.write("% Band of the largest value of  'n-sigma' of 'SALT2/SNooPy - NIR' for RMS:\n")
textfile_1.write('\\newcommand{\\optminusnirRMSnsLargestBand}{\\%s} \n'%(
    array2_bands_np[np.argmax(array2_np)] ) )
textfile_1.write("% Delta(RMS) of the largest value of 'n-sigma' of 'SALT2/SNooPy - NIR' for RMS,\n")
textfile_1.write('\\newcommand{\\optminusnirRMSdeltaLargest}{%.2f} \n'%(
    array2_dRMS_np[np.argmax(array2_np)] ) )

#-------------------
textfile_1.write("\n")
textfile_1.write("% Smallest value of 'n-sigma' of 'SALT2/SNooPy - NIR' for RMS,\n")
textfile_1.write("% where NIR =  any YJHK, JH, YJH only.\n")
# Create a new array where I ignore K band:
textfile_1.write('\\newcommand{\\optminusnirRMSnsSmallestCombined}{%.1f} \n'%(
    min(scatter_OptvsNIR_np['f9'][8:]) ) )
textfile_1.write("% Band of the Smallest value of  'n-sigma' of 'SALT2/SNooPy - NIR' for RMS:\n")
textfile_1.write('\\newcommand{\\optminusnirRMSnsSmallestCombinedBands}{\\%s} \n'%(
    scatter_OptvsNIR_np['f0'][8:][np.argmin(scatter_OptvsNIR_np['f9'][8:])] ) )

textfile_1.write("\n")
textfile_1.write("% Largest value of 'n-sigma' of 'SALT2/SNooPy - NIR' for RMS,\n")
textfile_1.write("% where NIR =  any YJHK, JH, YJH only.\n")
# Create a new array where I ignore K band:
textfile_1.write('\\newcommand{\\optminusnirRMSnsLargestCombined}{%.1f} \n'%(
    max(scatter_OptvsNIR_np['f9'][8:]) ) )
textfile_1.write("% Band of the largest value of  'n-sigma' of 'SALT2/SNooPy - NIR' for RMS:\n")
textfile_1.write('\\newcommand{\\optminusnirRMSnsLargestCombinedBands}{\\%s} \n'%(
    scatter_OptvsNIR_np['f0'][8:][np.argmax(scatter_OptvsNIR_np['f9'][8:])] ) )

#-------------------
textfile_1.write("\n")
textfile_1.write('\\newcommand{\\wrmsNIRmaxK}{%.3f} \n'%(K_wRMS_GP[5][0]))

#-------------------
textfile_1.write("\n")
textfile_1.write("% 'n-sigma' of SNooPy - J, for wRMS.\n")
textfile_1.write('\\newcommand{\\optminusnirsigmaJWorstwRMS}{%.1f} \n'%(
    (snoopy_j_wrms/err_snoopy_j_wrms) ) )

#-------------------
textfile_1.write("\n")
textfile_1.write("% 'n-sigma' of SNooPy - J, for RMS.\n")
textfile_1.write('\\newcommand{\\optminusnirsigmaJWorstRMS}{%.1f} \n'%(
    (snoopy_j_rms/err_snoopy_j_rms) ) )

#-------------------
textfile_1.write('\n')
textfile_1.write("% Minimum and maximum values of the intrinsic scatter in the GP NIR_max method only,\n")
textfile_1.write("% from all the NIR bands individually or combined.\n")
textfile_1.write('\\newcommand{\\nirgpbestscatter}{$\\sim %.2f$-$%.2f$} \n'%(
    min(scatter_GPNIRmax_np['f3']), max(scatter_GPNIRmax_np['f3'])  ) )

#-----------------------------------------------------------------------------80$$$$$$$$$$$

textfile_1.write("\n")
textfile_1.write("% Largest value of 'n-sigma' of 'Bmax - NIRmax' for intrinsic scatter,\n")
textfile_1.write("% where NIR = Y, J, H, K, any YJHK, JH, YJH.\n")
textfile_1.write('\\newcommand{\\BmaxminusNIRmaxLargestNS}{%.2f} \n'%(
    max(scatter_NIRmaxBmax_np['f3']) ) )
textfile_1.write("% Band of the largest value of 'Bmax - NIRmax' for intrinsic scatter:\n")
textfile_1.write('\\newcommand{\\BmaxminusNIRmaxLargestNSBand}{\\%s} \n'%(
    scatter_NIRmaxBmax_np['f0'][np.argmax(scatter_NIRmaxBmax_np['f3'])] ) )

#-------------------
textfile_1.write("\n")
textfile_1.write("% Largest value of 'n-sigma' of 'Bmax - NIRmax' for RMS,\n")
textfile_1.write("% where NIR = Y, J, H, K, any YJHK, JH, YJH.\n")
textfile_1.write('\\newcommand{\\BmaxminusNIRmaxLargestNSrms}{%.2f} \n'%(
    max(scatter_NIRmaxBmax_np['f9']) ) )
textfile_1.write("% Band of the largest value of 'Bmax - NIRmax' for RMS:\n")
textfile_1.write('\\newcommand{\\BmaxminusNIRmaxLargestNSrmsBand}{\\%s} \n'%(
    scatter_NIRmaxBmax_np['f0'][np.argmax(scatter_NIRmaxBmax_np['f9'])] ) )

textfile_1.write("% Delta(RMS) of the largest value of 'Bmax - NIRmax' for RMS,\n")
textfile_1.write('\\newcommand{\\BmaxminusNIRmaxLargestDeltarms}{%.2f} \n'%(
    scatter_NIRmaxBmax_np['f7'][np.argmax(scatter_NIRmaxBmax_np['f9'])] ) )

#-----------------------------------------------------------------------------80

#  textfile_1.write('\\newcommand{\\}{%.3f \\pm %.3f} \n'%())

textfile_1.close()

textfile_1.close();textfile_1.close();textfile_1.close();
textfile_1.close();textfile_1.close();textfile_1.close();

print "\n#    All done with no issues :)"

