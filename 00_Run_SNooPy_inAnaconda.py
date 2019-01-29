#!/usr/bin/env python
# coding: utf-8

# # Notebook to run SNooPy on python Anaconda

from snpy import *

# import pymysql as sql

# db_test = sqlbase.connect(host='localhost', user='root', passwd='Mimiraptus',db='SN')

sqlmod.setSQL('default')

s = sn('SN2004dt')

# # Fit with SNooPy a bunch of  LCs in batch
# # 3_Canopy_Fitting_v1_10.py

# The first path is for the python code to run, and the second for the location of the snoopy data.

#-----------------------------------------------------------------------------80

# ### CfA

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/3_Canopy_Fitting_v1_10.py"   "/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/AndyLCComp_2018_02/all/snoopy/2_Combine_Fit"')

# ### CSP

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/03_Canopy_Fitting_Main.py"   "/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/AndyLCComp_2018_02/all/snoopy/2_Combine_Fit/CSP/tmp_tests"')

# ### Others

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/Codes/Canopy/3_Canopy_Fitting_v1_10.py"   "/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/AndyLCComp_2018_02/all/snoopy/2_Combine_Fit"')

# ### RAISIN-1

# %run "/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/Codes/\
# Canopy/3_Canopy_Fitting_v1_10.py"   "/Users/arturo/Dropbox/Research/\
# Articulos/12_RAISINs/Data/RAISIN_1/Data/PS1/2018_01_30/DavidJones/snoopy_fit_v1/data"

# ### RAISIN-2

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/03_Canopy_Fitting_Main.py"   "/Users/arturo/Dropbox/Research/Articulos/14_RAISINs/Data/raisin12/data/2018_12_17/data_v2/DES/conda_test1"')

# # Canopy_Extract_DistanceMu_FromSnpy_v1_6_OK.py

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/Codes/Canopy/00_Canopy_Anaconda_Extract_DistanceMu_FromSnpy_OK.py"')

# # Canopy_Plot_AllLCData_WithoutFitting_v1_0.py

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/SoftwareResearch/Snoopy/Codes/Canopy/Canopy_Plot_AllLCData_WithoutFitting_v1_0.py"   "/Users/arturo/Dropbox/Research/Articulos/12_RAISINs/Data/RAISIN_2/Data/HST_DES/data"')

# # 6_Canopy_AbsAppMag_Phase_ForEachBand_v5_15.py
#
# The first path is for the python code to run.

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/06_Canopy_AbsAppMag_Phase_ForEachBand.py" J CfA')

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/06_Canopy_AbsAppMag_Phase_ForEachBand.py" J CSP')

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/06_Canopy_AbsAppMag_Phase_ForEachBand.py" J Others')

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/06_Canopy_AbsAppMag_Phase_ForEachBand.py" Y CfA')

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/06_Canopy_AbsAppMag_Phase_ForEachBand.py" Y CSP')

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/06_Canopy_AbsAppMag_Phase_ForEachBand.py" Y Others')

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/06_Canopy_AbsAppMag_Phase_ForEachBand.py" H CfA')

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/06_Canopy_AbsAppMag_Phase_ForEachBand.py" H CSP')

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/06_Canopy_AbsAppMag_Phase_ForEachBand.py" H Others')

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/06_Canopy_AbsAppMag_Phase_ForEachBand.py" K CfA')

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/06_Canopy_AbsAppMag_Phase_ForEachBand.py" K CSP')

get_ipython().run_line_magic('run', '"/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/github/LowzNIRFit/06_Canopy_AbsAppMag_Phase_ForEachBand.py" K Others')

