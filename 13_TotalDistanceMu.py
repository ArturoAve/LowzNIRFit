#!/usr/bin/env python
# coding: utf-8

# # Total distance modulus
#
# ### For Template and Gaussian-Processes Methods
#
# ### Computing the correlation among bands and using it to determine a single distance modulus using the 3 distance moduli in each band.
#
# Quantifying the correlation between bands by computing the population covariance matrix.
#
# The output are 2 datafiles:
#     - one similar in format to the 'DistanceMu_Good_AfterCutoffs_Main_.txt' files.
#     - one with latex format to easily put a table in the paper.
#
# #### Notes:
#
# - I have to plot the Hubble diagram of the output of this notebook by using the regular "11_DistanceMu_HubbleDiagram_vXX.ipynb" notebook.
#
#
# - If I need to remove an outlier in any of the "any YJHK", YJHK, YJH or JHK Gaussian-Process or template Hubble diagrams, then:
#
#
# 1. first I need to "comment" that SN in 'DistanceMu_Good_AfterCutoffs_Main_.txt' in all the bands,
# 2. then, rerun this "13_TotalDistanceMu_vXX.ipynb" iPython notebook to recompute the covariance matrices and their inverse matrix from the remaining SNe.
# 3. plot the Hubble diagram using "11_DistanceMu_HubbleDiagram_vXX.ipynb" notebook.
#

import numpy as np
from matplotlib import pyplot as plt
import os # To use command line like instructions
from scipy.stats.stats import pearsonr

#--------------------------------------------------------60
code_created_by = 'Arturo_Avelino'
# On date: 2017.01.10 (yyyy.mm.dd)
code_name = '13_TotalDistanceMu.ipynb'
version_code = '0.1.17'
last_update = '2019.01.28'
#--------------------------------------------------------60

##############################################################################80

# # USER

Method =  'Template_M'   # template method
# Method =  'GP_M'       # Gaussian-Process method at NIR max
# Method =  'GP_Bmax'       # Gaussian-Process method at B max

# Peculiar velocity uncertainty to use in the plots.
# This is just to -select- the appropiate folder containing
# the Hubble diagram that was previously computed using this value of vpecFix.
#  Options e.g., (150, 250, 300, 400). This has to be consistent wity
# the vpecFix value used for the individual NIR bands.

vpecFix = 150 # km/s.

MainDir = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/TheTemplates/'

if Method =='GP_M': SubFolderSave = 'GaussianProcess'
elif Method =='GP_Bmax': SubFolderSave = 'GP_Bmax'
elif Method =='Template_M': SubFolderSave = 'Template'

DirSaveOutput = MainDir+'AllBands/Plots/HubbleDiagram_vpec%s/'%(vpecFix)+SubFolderSave+'/'

# Minimum number of SNe to define and use the YJHK covariance matrix:
MinSNeForYJHK = 5  # 5 is the best option, given that for the
# GP HD I have 4 SNe only.

#--------------------------------------------------------60

#--- Fixed values ---

HoFix = 73.24 # Hubble constant (km/(s Mpc))

c = 299792.458  # Speed of light (km/s)
cc = 299792.458  # Speed of light (km/s)
OmMFix = 0.28 # Omega_Matter
OmLFix = 0.72 # Omega_Lambda
wFix = -1.0 # Dark energy EoS

# Minimum values for anything else:

DirSaveOutput

##############################################################################80

# # Automatic

#- Force the creation of the directory to save the outputs.
#- "If the subdirectory does not exist then create it"
import os # To use command line like instructions
if not os.path.exists(DirSaveOutput): os.makedirs(DirSaveOutput)

# NOTE: Below this cell is included automatically the
# case for 'DistanceMu_Good_AfterCutoffs_Main_Notes_.txt' in case it
# exist for each NIR band.

#-------------------
#    Y band

# Hubble diagram
if Method == 'Template_M':
    Y_HubbleDir = MainDir+'Y_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi3_EBVh0.4_Method7_MinData3_vpec%s_ok/plots_HD/'%vpecFix
    Y_HubbleDataFile = 'DistanceMu_Good_AfterCutoffs_Main_.txt'

elif Method == 'GP_M':
    Y_HubbleDir = MainDir+'Y_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_NIRmax/'%vpecFix
    Y_HubbleDataFile = 'DistanceMu_Good_AfterCutoffs_Main_.txt'

elif Method == 'GP_Bmax':
    Y_HubbleDir = MainDir+'Y_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_Bmax/'%vpecFix
    Y_HubbleDataFile = 'DistanceMu_Good_AfterCutoffs_Main_.txt'

#-------------------
#    J band

# Hubble diagram
if Method == 'Template_M':
    J_HubbleDir = MainDir+'J_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi_1e6_EBVh0.4_Method7_MinData3_vpec%s_ok/plots_HD/'%vpecFix
    J_HubbleDataFile = 'DistanceMu_Good_AfterCutoffs_Main_.txt'

elif Method == 'GP_M':
    J_HubbleDir = MainDir+'J_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_NIRmax/'%vpecFix
    J_HubbleDataFile = 'DistanceMu_Good_AfterCutoffs_Main_.txt'

elif Method == 'GP_Bmax':
    J_HubbleDir = MainDir+'J_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_Bmax/'%vpecFix
    J_HubbleDataFile = 'DistanceMu_Good_AfterCutoffs_Main_.txt'

#-------------------
#    H band

# Hubble diagram
if Method == 'Template_M':
    H_HubbleDir = MainDir+'H_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi_1e6_EBVh0.4_Method7_MinData3_vpec%s_ok/plots_HD/'%vpecFix
    H_HubbleDataFile = 'DistanceMu_Good_AfterCutoffs_Main_.txt'

elif Method == 'GP_M':
    H_HubbleDir = MainDir+'H_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_NIRmax/'%vpecFix
    H_HubbleDataFile = 'DistanceMu_Good_AfterCutoffs_Main_.txt'

elif Method == 'GP_Bmax':
    H_HubbleDir = MainDir+'H_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_Bmax/'%vpecFix
    H_HubbleDataFile = 'DistanceMu_Good_AfterCutoffs_Main_.txt'

#-------------------
#    K band

# Hubble diagram
if Method == 'Template_M':
    K_HubbleDir = MainDir+'K_band/Std_filters/4_HubbleDiagram_FlatPrior/AllSamples/Templ_AllSamples_z_gr_0/Phase-8_30_resid20_chi4_EBVh0.4_Method7_MinData3_vpec%s_ok/plots_HD/'%vpecFix
    K_HubbleDataFile = 'DistanceMu_Good_AfterCutoffs_Main_.txt'

elif Method == 'GP_M':
    K_HubbleDir = MainDir+'K_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_NIRmax/'%vpecFix
    K_HubbleDataFile = 'DistanceMu_Good_AfterCutoffs_Main_.txt'

elif Method == 'GP_Bmax':
    K_HubbleDir = MainDir+'K_band/Std_filters/4_HubbleDiagram_FlatPrior_GP/AllSamples/vpec%s_Bmax/'%vpecFix
    K_HubbleDataFile = 'DistanceMu_Good_AfterCutoffs_Main_.txt'

#--------------------------------------------------------60

# Reading the 'DistanceMu_Good_AfterCutoffs_Main_.txt' datatable files for each band

try:
    MuData_Y = np.genfromtxt(Y_HubbleDir+Y_HubbleDataFile[:-4]+'Notes_.txt',
                                 dtype=['S30',
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float])
except:
    MuData_Y = np.genfromtxt(Y_HubbleDir+Y_HubbleDataFile,
                                 dtype=['S30',
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float])

#--------------------------------------------------
try:
    MuData_J = np.genfromtxt(J_HubbleDir+J_HubbleDataFile[:-4]+'Notes_.txt',
                                 dtype=['S30',
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float])
except:
    MuData_J = np.genfromtxt(J_HubbleDir+J_HubbleDataFile,
                                 dtype=['S30',
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float])

#--------------------------------------------------
try:
    MuData_H = np.genfromtxt(H_HubbleDir+H_HubbleDataFile[:-4]+'Notes_.txt',
                                 dtype=['S30',
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float])

except:
    MuData_H = np.genfromtxt(H_HubbleDir+H_HubbleDataFile,
                                 dtype=['S30',
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float])

#--------------------------------------------------
try:
    MuData_K = np.genfromtxt(K_HubbleDir+K_HubbleDataFile[:-4]+'Notes_.txt',
                                 dtype=['S30',
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float])
except:
    MuData_K = np.genfromtxt(K_HubbleDir+K_HubbleDataFile,
                                 dtype=['S30',
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float,float,float,float,float,float,float,float,
                                       float,float,float])

# Create an array of sn NAMES but removing the last part of the name that
# refers to the specific band. This will be used to compare the SN names and
# find the SNe with observed in multiple bands.

MuData_Y_names = [];  sn=''
for sn in MuData_Y['f0']:
    MuData_Y_names += [sn[:-2]]

MuData_J_names = [];  sn=''
for sn in MuData_J['f0']:
    MuData_J_names += [sn[:-2]]

MuData_H_names = [];  sn=''
for sn in MuData_H['f0']:
    MuData_H_names += [sn[:-2]]

MuData_K_names = [];  sn=''
for sn in MuData_K['f0']:
    MuData_K_names += [sn[:-2]]

# ### Dictionaries

# Dictionaries (name: mu | error_mu | residual_mu | sample | z_CMB | error_zCMB | chi2dof)

Mu_Y_dict = {MuData_Y_names[i]: [ MuData_Y['f3'][i], MuData_Y['f4'][i], MuData_Y['f5'][i], MuData_Y['f7'][i],
                                  MuData_Y['f1'][i], MuData_Y['f2'][i], MuData_Y['f6'][i]  ]
             for i in range(len(MuData_Y_names))}

Mu_J_dict = {MuData_J_names[i]: [ MuData_J['f3'][i], MuData_J['f4'][i], MuData_J['f5'][i], MuData_J['f7'][i],
                                  MuData_J['f1'][i], MuData_J['f2'][i], MuData_J['f6'][i]  ]
             for i in range(len(MuData_J_names))}

Mu_H_dict = {MuData_H_names[i]: [ MuData_H['f3'][i], MuData_H['f4'][i], MuData_H['f5'][i], MuData_H['f7'][i],
                                  MuData_H['f1'][i], MuData_H['f2'][i], MuData_H['f6'][i]  ]
             for i in range(len(MuData_H_names))}

Mu_K_dict = {MuData_K_names[i]: [ MuData_K['f3'][i], MuData_K['f4'][i], MuData_K['f5'][i], MuData_K['f7'][i],
                                  MuData_K['f1'][i], MuData_K['f2'][i], MuData_K['f6'][i]  ]
             for i in range(len(MuData_K_names))}

#--------------------------------------------------------60

# ### $\mu_{\Lambda{\rm CDM}}$

from scipy.integrate import quad as intquad

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
ztest1 = 0.01

print 'Checking that the functions work well:', DistanceMu(ztest1, OmMFix, wFix, HoFix)
# Checking that the functions work well: 33.1141460988 # Ho=72
# Checking that the functions work well: 33.0773926577 # Ho=73.24

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

# #### Get the name of this ipython notebook
# To print it in the output text files as reference

# %%javascript
# var kernel = IPython.notebook.kernel;
# var thename = window.document.getElementById("notebook_name").innerHTML;
# var command = "NotebookName = " + "'"+thename+".ipynb"+"'";
# kernel.execute(command);

# print '#', (NotebookName)

# Get the current date and time
import datetime

# Read the time and date now
now = datetime.datetime.now()

##############################################################################80

# # Computing the population covariance matrix of
#
# # $\{\Delta \mu^Y_s, \Delta \mu^J_s, \Delta \mu^H_s \}$,  $\{\Delta \mu^J_s, \Delta \mu^H_s, \Delta \mu^K_s \}$ and $\{\Delta \mu^Y_s, \Delta \mu^J_s, \Delta \mu^H_s, \Delta \mu^K_s \}$, $\{\Delta \mu^J_s, \Delta \mu^H_s \}$,
#
# ### Arrays of residuals values

# Create an array of {residual_mu^Y_s, residual_mu^J_s, residual_mu^H_s}
# with SNe that are common to all bands

Residual_YJH = []
countInt = 0

for i in range(len(MuData_J_names)): # Loop over SNe in J band (it contains all the SNe)
    snName = MuData_J_names[i]

    # check that a given supernova was observed in YJH bands:
    if snName in MuData_H_names: # In H band
        if snName in MuData_Y_names: # In Y band
            Residual_YJH+=[ [Mu_Y_dict[snName][2], Mu_J_dict[snName][2], Mu_H_dict[snName][2]] ]
            countInt = countInt + 1

Residual_YJH_np =  np.array(Residual_YJH)
print '# Number of SNe common in the YJH bands:', countInt, ', for method: %s'%Method

# Number of SNe common in the YJH bands: 27 , for method: Template_M
# Number of SNe common in the YJH bands: 15 , for method: GP_M

# Create an array of {residual_mu^J_s, residual_mu^H_s}
# with SNe that are common to all bands

Residual_JH = []
countInt = 0

for i in range(len(MuData_J_names)): # Loop over SNe in J band (it contains all the SNe)
    snName = MuData_J_names[i]

    # check that a given supernova was observed in JH bands:
    if snName in MuData_H_names: # In H band
        Residual_JH+=[ [Mu_J_dict[snName][2], Mu_H_dict[snName][2]] ]
        countInt = countInt + 1

Residual_JH_np =  np.array(Residual_JH)
print '# Number of SNe common in the JH bands:', countInt, ', for method: %s'%Method

# Number of SNe common in the JH bands: 55 , for method: GP_Bmax

# Create an array of {residual_mu^J_s, residual_mu^H_s, residual_mu^K_s}
# with SNe that are common to all bands

Residual_JHK = []
countInt = 0

for i in range(len(MuData_J_names)): # Loop over SNe in J band (it contains all the SNe)
    snName = MuData_J_names[i]

    # check that a given supernova was observed in JHK bands:
    if snName in MuData_H_names: # In H band
        if snName in MuData_K_names: # In K band
            Residual_JHK+=[ [Mu_J_dict[snName][2], Mu_H_dict[snName][2], Mu_K_dict[snName][2] ] ]
            countInt = countInt + 1

Residual_JHK_np =  np.array(Residual_JHK)
print '# Number of SNe common in the JHK bands:', countInt, ', for method: %s'%Method

# Number of SNe common in the JHK bands: 40 , for method: Template_M
# Number of SNe common in the JHK bands: 19 , for method: GP_M

# Create an array of {residual_mu^Y_s, residual_mu^J_s, residual_mu^H_s, residual_mu^K_s}
# with SNe that are common to all bands

Residual_YJHK = []
count_YJHK = 0

for i in range(len(MuData_J_names)): # Loop over SNe in J band (it contains all the SNe)
    snName = MuData_J_names[i]

    # check that a given supernova was observed in YJHK bands:
    if snName in MuData_H_names: # In H band
        if snName in MuData_K_names: # In K band
            if snName in MuData_Y_names: # In Y band
                Residual_YJHK += [ [ Mu_Y_dict[snName][2], Mu_J_dict[snName][2],
                                    Mu_H_dict[snName][2], Mu_K_dict[snName][2] ] ]
                count_YJHK = count_YJHK + 1

# In Andy's compilation there are four SNe only with observations in YJHK band that
# passed my cutoffs for the Gaussian-Process Hubble diagram, then the covariance matrix
# for YJHK was constructed based on 4 data only;
# this makes the Cov_ResidualMu_YJHK very unreliable, so for these 4 SNe I'm going to ignore
# their K band (the most unstable band) and then use the Cov_ResidualMu_YJH for them instead,
# I mean, I treat those 4 SNe as if they had observations in the YJH bands only and then
# use the corresponding YJH covariance matrix for them. This hack is implemented
# in the main loop in the YJHK section as follows:
#    if count_YJHK >= 4: mu_K = Mu_K_dict[snName][0]; err_mu_K = Mu_K_dict[snName][1];
#    else: mu_K = -1.0 err_mu_K = -1.0;

Residual_YJHK_np =  np.array(Residual_YJHK)
print '# Number of SNe common in the YJHK bands:', count_YJHK, ', for method: %s'%Method
if count_YJHK < MinSNeForYJHK:
    print "# NOTE: YJHK covariance matrix from less or equal to 4 SNe, so I'll not define it."

# Number of SNe common in the YJHK bands: 8 , for method: Template_M

# Number of SNe common in the YJHK bands: 4 , for method: GP_M
# NOTE: YJHK covariance matrix from less than 4 SNe, so I'll not define it.

print '# Method: %s'%Method
Residual_YJHK_np

# Showing the residual arrays
print '-'*40
print '# Residual_YJH'
print '# Method:', Method
print Residual_YJH_np.T

# Showing the residual arrays
print '-'*40
print 'Residual_JHK'
print 'Method:', Method
print Residual_JHK_np.T

# Showing the residual arrays
print '-'*40
print 'Residual_YJHK'
print 'Method:', Method
print Residual_YJHK_np.T
if count_YJHK < MinSNeForYJHK:
    print "NOTE: YJHK covariance matrix from less or equal to 4 SNe, so I'll not define it."

#-----------------------------------------------------------------------------80

# ## Covariance matrix computation
#
# #### Definition of the Covariance matrices

# The covariance matrices

Cov_ResidualMu_YJH   = np.cov(Residual_YJH_np.T)
Cov_ResidualMu_JH   = np.cov(Residual_JH_np.T)
Cov_ResidualMu_JHK   = np.cov(Residual_JHK_np.T)
if count_YJHK >= MinSNeForYJHK: Cov_ResidualMu_YJHK  = np.cov(Residual_YJHK_np.T)
else: print "NOTE: YJHK covariance matrix from less or equal to 4 SNe, so I don't define it."

#--------------------------------------------------
# Latex file

txtfile_1 = open(DirSaveOutput+'Latex_CovM_JH.tex', 'w')

txtfile_1.write('% JH covariance matrix\n')
txtfile_1.write('%% Method: %s\n'%Method)

for ii in range(len(Cov_ResidualMu_JH)):
    txtfile_1.write('%.4f & %.4f \\\\\n'%(
        Cov_ResidualMu_JH[ii][0], Cov_ResidualMu_JH[ii][1]) )
txtfile_1.close();

print '-'*40
print 'The covariance matrix YJH:'
print 'Method:', Method
print Cov_ResidualMu_YJH

#--------------------------------------------------
# Latex file

txtfile_1 = open(DirSaveOutput+'Latex_CovM_YJH.tex', 'w')

txtfile_1.write('% YJH covariance matrix\n')
txtfile_1.write('%% Method: %s\n'%Method)

for ii in range(len(Cov_ResidualMu_YJH)):
    txtfile_1.write('%.4f & %.4f & %.4f \\\\\n'%(
        Cov_ResidualMu_YJH[ii][0], Cov_ResidualMu_YJH[ii][1],
        Cov_ResidualMu_YJH[ii][2]) )
txtfile_1.close();

print '-'*40
print 'The covariance matrix JHK:'
print 'Method:', Method
print Cov_ResidualMu_JHK

#--------------------------------------------------
# Latex file

txtfile_1 = open(DirSaveOutput+'Latex_CovM_JHK.tex', 'w')

txtfile_1.write('% JHK covariance matrix\n')
txtfile_1.write('%% Method: %s\n'%Method)

for ii in range(len(Cov_ResidualMu_JHK)):
    txtfile_1.write('%.4f & %.4f & %.4f \\\\\n'%(
        Cov_ResidualMu_JHK[ii][0], Cov_ResidualMu_JHK[ii][1],
        Cov_ResidualMu_JHK[ii][2]) )
txtfile_1.close();

print '-'*40
print 'The covariance matrix YJHK:'
print 'Method:', Method

#--------------------------------------------------
# Latex file

txtfile_1 = open(DirSaveOutput+'Latex_CovM_YJHK.tex', 'w')
txtfile_1.write('% YJHK covariance matrix\n')
txtfile_1.write('%% Method: %s\n'%Method)

if count_YJHK >= MinSNeForYJHK:
    print Cov_ResidualMu_YJHK
    print '-'*20
    print "% YJHK, Latex format:"
    for ii in range(len(Cov_ResidualMu_YJHK)):
        print "%.4f & %.4f & %.4f & %.4f \\\\"%(Cov_ResidualMu_YJHK[ii][0],
        Cov_ResidualMu_YJHK[ii][1], Cov_ResidualMu_YJHK[ii][2],
        Cov_ResidualMu_YJHK[ii][3])

    for ii in range(len(Cov_ResidualMu_YJHK)):
        txtfile_1.write('%.4f & %.4f & %.4f & %.4f \\\\\n'%(
        Cov_ResidualMu_YJHK[ii][0],
        Cov_ResidualMu_YJHK[ii][1], Cov_ResidualMu_YJHK[ii][2],
        Cov_ResidualMu_YJHK[ii][3]) )

else:
    txtfile_1.write("WARNING: There are less or equal to 4 SNe with YJHJ observations, so I'll not define the YJHK covariance matrix\n")

txtfile_1.close();

#--------------------------------------------------------60

# #### Definition of the inverse(Covariance matrix)

# Inverse of the total covariance matrix.
# In this case the inverse covariance matrix of 2 bands is NOT a submatrix of the
# inverse covariance matrix of 3 bands.

# There are not SNe in YHK bands that are not in J band, that's why it is
# not necessary to define the YH, HK, or YK covariance submatrices.

# These matrices are what I actually use in the main loop, i.e., I don't directly use
# the covariance matrices but their inverse.

# In Andy's compilation there are four SNe only with observations in YJHK band that
# passed my cutoffs for the Gaussian-Process Hubble diagram, then the covariance matrix
# for YJHK was constructed based on 4 data only;
# this makes the Cov_ResidualMu_YJHK very unreliable, so for these 4 SNe I'm going to ignore
# their K band (the most unstable band) and then use the Cov_ResidualMu_YJH for them instead,
# I mean, I treat those 4 SNe as if they had observations in the YJH bands only and then
# use the corresponding YJH covariance matrix for them. This hack is implemented
# in the main loop in the YJHK section as follows:
#    if count_YJHK >= 4: mu_K = Mu_K_dict[snName][0]; err_mu_K = Mu_K_dict[snName][1];
#    else: mu_K = -1.0 err_mu_K = -1.0;

InvCov_ResidualMu_YJH  = np.linalg.inv(Cov_ResidualMu_YJH) # YJH bands
InvCov_ResidualMu_JHK  = np.linalg.inv(Cov_ResidualMu_JHK) # JHK bands
InvCov_ResidualMu_YJ = np.linalg.inv(Cov_ResidualMu_YJH[:2,:2]) # YJ bands
InvCov_ResidualMu_JH = np.linalg.inv(Cov_ResidualMu_JH) # JH bands
InvCov_ResidualMu_JK = np.linalg.inv(Cov_ResidualMu_JHK[[0,2]][:,[0,2]]) # JK bands

if count_YJHK >= MinSNeForYJHK:
    InvCov_ResidualMu_YJHK = np.linalg.inv(Cov_ResidualMu_YJHK) # YJHK bands
else: print "NOTE: YJHK covariance matrix from less or equal to 4 SNe, so I don't define it."

print '-'*40
print 'InvCov_ResidualMu_YJH'
print 'Method:', Method
print InvCov_ResidualMu_YJH

print '-'*40
print 'InvCov_ResidualMu_JHK'
print 'Method:', Method
print InvCov_ResidualMu_JHK

print '-'*40
print 'InvCov_ResidualMu_YJHK'
print 'Method:', Method
if count_YJHK >= MinSNeForYJHK: print InvCov_ResidualMu_YJHK
else: print "YJHK covariance matrix from less or equal to 4 SNe, so I didn't define it."

##############################################################################80

# ### Just some cross checking of the matrices
#
# The quantities defined in this subsection are not used for any posterior computation

# Checking that the covariance matrix of YJ bands only is the submatrix
# of the 'np.cov(Residual_YJH)' covariance matrix

Cov_ResidualMu_YJ = np.cov(Residual_YJH_np[:,0:2].T)

print '-'*40
print 'Method:', Method
print ''
print 'The covariance matrix of YJ bands computed directly from the YJH residual matrix:'
print Cov_ResidualMu_YJ # Result: it IS a submatrix of 'np.cov(Residual_YJH)' as I expected.

print ''
print 'Comparing with the YJ submatrix of the YJH matrix: both (above) and this submatrices are the same:'
print Cov_ResidualMu_YJH[:2,:2]

# Checking that the inverse of InvCov_ResidualMu_YJHK is the equal to the
# covariance matrix Cov_ResidualMu_YJHK

if count_YJHK > MinSNeForYJHK:
    print '-'*40
    print 'Checking that inv(InvCov_ResidualMu_YJHK) should be equal to Cov_ResidualMu_YJHK.'
    print 'Method:', Method
    print ''
    print 'Inverse of InvCov_ResidualMu_YJHK:'
    print np.linalg.inv(InvCov_ResidualMu_YJHK) # Result:

    print ''
    print 'Covariance matrix:'
    print Cov_ResidualMu_YJHK
else:
    print "There are equal or less than 4 SNe, so the matrix for YJHK was not defined"

print '-'*40
print 'Method:', Method
print ''
print 'JH cov submatrix computed directly from the -YJH- residual matrix:'
print np.cov(Residual_YJH_np[:,1:3].T)

print ''
print 'JH cov submatrix computed directly from the -JHK- residual matrix:'
print np.cov(Residual_JHK_np[:,0:2].T)

# I notice that the inverse of the covariance matrix of the YJ submatrix
# is NOT a submatrix of the inverse covariance matrix of 3 bands.

print '-'*40
print 'Method:', Method
print InvCov_ResidualMu_YJ

"""
----------------------------------------
Method: Template_M
[[ 76.31565035 -44.4444541 ]
 [-44.4444541   97.51955589]]
"""
0

# I notice that the inverse of the covariance matrix of the JH submatrix
# is NOT a submatrix of the inverse covariance matrix of 3 bands.

print '-'*40
print 'Method:', Method
print InvCov_ResidualMu_JH

# Checking that the inverse of the inverse covariance matrix is
# the covariance matrix. OK.

print '-'*40
print 'Method:', Method
print np.linalg.inv(InvCov_ResidualMu_YJH)

# The correlation matrix

# This quantity is NOT needed in my computations. I just
# compute it here to visualize it and gain better understanding
# of my results.

Corr_ResidualMu_YJH = np.corrcoef(Residual_YJH_np.T)

print '-'*40
print 'Method:', Method
print 'The correlation matrix of YJH bands:'
print Corr_ResidualMu_YJH

##############################################################################80

# # Total distance modulus

# ## Any YJHK HD
# ### Main loop
#
# This cell produce also the final latex table that will be used in the paper. The "14_Latex_Tables_Figures_v2_5.ipynb" script doesn't do it.

# The MAIN loop: SNe with distance modulus in ANY band.
# (The loops in the cells below are copied from this loop cell)

# SNe with YJHK data, but if there is not data in the 4 bands, then use 3, or 2, or 1 band.

# Open the writting to a text file
textfile_1 = open(DirSaveOutput+'Table_TotalMu_AllBands_.txt', 'w')
textfile_2 = open(DirSaveOutput+'Table_TotalMu_AllBands_Latex_.txt', 'w')
textfile_3 = open(DirSaveOutput+'Table_TotalMu_AllBands_snNames_.txt', 'w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd).")
text_Date   = '# On date: %s \n'%text_timenow
text_Author = '# Data table created by: Arturo Avelino \n'
text_script = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
text_line = '#'+'-'*60 + '\n'

textfile_1.write(text_line)
textfile_1.write(text_Author); textfile_1.write(text_Date); textfile_1.write(text_script);
textfile_1.write(text_line)

textfile_3.write(text_line)
textfile_3.write(text_Author); textfile_3.write(text_Date); textfile_3.write(text_script);
textfile_3.write(text_line)

text_1 = '#   SN name                     z_CMB         error_zcmb       mu_total      error_mu       \
mu_residual    chi2_dof  Sample     mu_Y         error_mu_Y       mu_J          error_mu_J      mu_H           \
error_mu_H     mu_K            error_mu_K    Unused columns \n'
textfile_1.write(text_1);
textfile_3.write(text_1);
# textfile_2.write('% '+text_1)

# Some counters
count_YJHK = 0;
count_YJH = 0; count_YJ = 0; count_JH = 0;
count_JHK = 0; count_JK = 0;
count_J=0; count_Y=0;  count_H=0; count_K=0;
count_all = 0;

# Loop over SNe in J band (it contains almost all the SNe)
for i in range(len(MuData_J_names)):

    # Resetting variables
    snName = '';
    zcmbInt = 0; error_zcmbInt=0; chi2dofInt=0;
    mu_LCDMInt=0; mu_residualInt=0;
    SampleFlag = 0; SampleFlatText='';
    mu_Y = 0; mu_J = 0; mu_H = 0; mu_K = 0;
    err_mu_Y=0; err_mu_J=0; err_mu_H=0; err_mu_K=0;
    mu_np  = np.zeros(3); errorMu_np = np.zeros(3);
    var_np = np.zeros(3);
    weight_InvError2Mu_np=np.zeros(3);
    CovMatrix_Mu = np.zeros([3,3]); InvCovMatrix_Mu = np.zeros([3,3]);
    Numerator_Mu = 0; Denominator_Mu = 0;
    Numerator_Var = 0;
    muTotal = 0; error_muTotal = 0;
    JbandOnly = 0;

    #--------------------------------

    # Defining variables

    snName = MuData_J_names[i]

    # Begin textfile variables --->>
    # These variables are not used for computing, they are just to
    # be written in the text files
    zcmbInt = Mu_J_dict[snName][4]
    error_zcmbInt = Mu_J_dict[snName][5]
    chi2dofInt = Mu_J_dict[snName][6]
    mu_LCDMInt = DistanceMu(zcmbInt, OmMFix, wFix, HoFix)

    # Sample flag
    SampleFlag = Mu_J_dict[snName][3]
    if   SampleFlag ==1: SampleFlatText = 'CfA'
    elif SampleFlag ==2: SampleFlatText = 'CSP'
    elif SampleFlag ==3: SampleFlatText = 'Others'
    # <<---- End textfile variables


    if snName in MuData_H_names: # In H band
        if snName in MuData_K_names: # In K band
            if snName in MuData_Y_names: # In Y band

                # SNe with distance mu from YJHK bands:

                #--- Begin copy (this is the original) ---------->>
                count_YJHK = count_YJHK + 1
                count_all  = count_all  + 1
                print count_all, ', J-H-K-Y__%r'%count_YJHK, ',', snName

                # In Andy's compilation there are four SNe only with observations
                # in YJHK band that
                # passed my cutoffs for the Gaussian-Process Hubble diagram,
                # then the covariance matrix
                # for YJHK was constructed based on 4 data only;
                # this makes the Cov_ResidualMu_YJHK very unreliable, so for these
                # 4 SNe I'm going to ignore
                # their K band (the most unstable band) and then use the
                # Cov_ResidualMu_YJH for them instead,
                # I mean, I treat those 4 SNe as if they had observations in the
                # YJH bands only and then
                # use the corresponding YJH covariance matrix for them.

                # The distance mu and uncertainties from each NIR band.
                mu_Y = Mu_Y_dict[snName][0]; err_mu_Y = Mu_Y_dict[snName][1];
                mu_J = Mu_J_dict[snName][0]; err_mu_J = Mu_J_dict[snName][1];
                mu_H = Mu_H_dict[snName][0]; err_mu_H = Mu_H_dict[snName][1];
                if count_YJHK >= MinSNeForYJHK:
                    mu_K = Mu_K_dict[snName][0]; err_mu_K = Mu_K_dict[snName][1];
                else: mu_K = -1.000000000; err_mu_K = -1.000000000;

                # distance mu array
                mu_np = np.array([ mu_Y, mu_J, mu_H, mu_K  ])

                # variance of distance mu array
                var_np = np.array([err_mu_Y**2., err_mu_J**2., err_mu_H**2., err_mu_K**2.])

                #--- Computing the best estimated total distance modulus ---

                if count_YJHK >= MinSNeForYJHK:
                    for i1 in range(4):
                        Numerator_Mu = Numerator_Mu + mu_np[i1]*np.sum(InvCov_ResidualMu_YJHK[i1,:])
                        Denominator_Mu = Denominator_Mu + np.sum(InvCov_ResidualMu_YJHK[i1,:])
                else:
                    for i1 in range(3):
                        Numerator_Mu = Numerator_Mu + mu_np[i1]*np.sum(InvCov_ResidualMu_YJH[i1,:])
                        Denominator_Mu = Denominator_Mu + np.sum(InvCov_ResidualMu_YJH[i1,:])

                # --- Computing the VARIANCE of the total distance modulus ---
                if count_YJHK >= MinSNeForYJHK:
                    for i1 in range(4):
                        Numerator_Var = Numerator_Var + var_np[i1]*(np.sum(InvCov_ResidualMu_YJHK[i1,:])/Denominator_Mu)**2.
                else:
                    for i1 in range(3):
                        Numerator_Var = Numerator_Var + var_np[i1]*(np.sum(InvCov_ResidualMu_YJH[i1,:])/Denominator_Mu)**2.

                # <<-------- End copy (this is the original)---


            else: # SNe with distance mu from JHK bands only

                # Begin copy  ---------->>
                count_JHK = count_JHK + 1
                count_all  = count_all  + 1
                print count_all, ', J-H-K__%r'%count_JHK, ',', snName

                # The distance mu and uncertainties for the bands.
                mu_Y = -1.000000000;         err_mu_Y = -1.000000000;
                mu_J = Mu_J_dict[snName][0]; err_mu_J = Mu_J_dict[snName][1];
                mu_H = Mu_H_dict[snName][0]; err_mu_H = Mu_H_dict[snName][1];
                mu_K = Mu_K_dict[snName][0]; err_mu_K = Mu_K_dict[snName][1];

                # distance mu array
                mu_np = np.array([ mu_J, mu_H, mu_K ])
                # variance of distance mu array
                var_np = np.array([err_mu_J**2., err_mu_H**2., err_mu_K**2.])

                #--- Computing the best estimated total distance modulus ---
                for i2 in range(3):
                    Numerator_Mu = Numerator_Mu + mu_np[i2]*np.sum(InvCov_ResidualMu_JHK[i2,:])
                    Denominator_Mu = Denominator_Mu + np.sum(InvCov_ResidualMu_JHK[i2,:])

                # --- Computing the VARIANCE of the total distance modulus ---
                for i2 in range(3):
                    Numerator_Var = Numerator_Var + var_np[i2]*(np.sum(InvCov_ResidualMu_JHK[i2,:])/Denominator_Mu)**2.

                # <<-------- End copy


        elif snName in MuData_Y_names: # In Y band

            # SNe with distance mu from YJH bands.

            # Begin copy  ---------->>
            count_YJH = count_YJH + 1
            count_all  = count_all  + 1
            print count_all, ', J-H-Y__%r'%count_YJH, ',', snName

            # The distance mu and uncertainties for the bands.
            mu_Y = Mu_Y_dict[snName][0]; err_mu_Y = Mu_Y_dict[snName][1];
            mu_J = Mu_J_dict[snName][0]; err_mu_J = Mu_J_dict[snName][1];
            mu_H = Mu_H_dict[snName][0]; err_mu_H = Mu_H_dict[snName][1];
            mu_K = -1.000000000;         err_mu_K = -1.000000000;

            # distance mu array
            mu_np = np.array([ mu_Y, mu_J, mu_H  ])
            # variance of distance mu array
            var_np = np.array([err_mu_Y**2., err_mu_J**2., err_mu_H**2.])

            #--- Computing the best estimated total distance modulus ---
            for i3 in range(3):
                Numerator_Mu = Numerator_Mu + mu_np[i3]*np.sum(InvCov_ResidualMu_YJH[i3,:])
                Denominator_Mu = Denominator_Mu + np.sum(InvCov_ResidualMu_YJH[i3,:])

            # --- Computing the VARIANCE of the total distance modulus ---
            for i3 in range(3):
                Numerator_Var = Numerator_Var + var_np[i3]*( np.sum(InvCov_ResidualMu_YJH[i3,:])/Denominator_Mu )**2.

            # <<-------- End copy

        else: # SNe with distance mu from JH bands.

            # Begin copy  ---------->>
            count_JH = count_JH + 1
            count_all  = count_all  + 1
            print count_all, ', J-H__%r'%count_JH, ',', snName

            # The distance mu and uncertainties for the bands.
            mu_Y = -1.000000000;         err_mu_Y = -1.000000000;
            mu_J = Mu_J_dict[snName][0]; err_mu_J = Mu_J_dict[snName][1];
            mu_H = Mu_H_dict[snName][0]; err_mu_H = Mu_H_dict[snName][1];
            mu_K = -1.000000000;         err_mu_K = -1.000000000;

            # distance mu array
            mu_np = np.array([ mu_J, mu_H  ])
            # variance of distance mu array
            var_np = np.array([err_mu_J**2., err_mu_H**2.])

            #--- Computing the best estimated total distance modulus ---
            for i4 in range(2):
                Numerator_Mu = Numerator_Mu + mu_np[i4]*np.sum(InvCov_ResidualMu_JH[i4,:])
                Denominator_Mu = Denominator_Mu + np.sum(InvCov_ResidualMu_JH[i4,:])

            # --- Computing the VARIANCE of the total distance modulus ---
            for i4 in range(2):
                Numerator_Var = Numerator_Var + var_np[i4]*(np.sum(InvCov_ResidualMu_JH[i4,:])/Denominator_Mu )**2.

            # <<-------- End copy


    elif snName in MuData_Y_names: # In Y band

        # SNe with distance mu from YJ bands.

        #--- Begin copy  ---------->>
        count_YJ = count_YJ + 1
        count_all  = count_all  + 1
        print count_all, ', J-Y__%r'%count_YJ, ',', snName

        # The distance mu and uncertainties for the bands.
        mu_Y = Mu_Y_dict[snName][0]; err_mu_Y = Mu_Y_dict[snName][1];
        mu_J = Mu_J_dict[snName][0]; err_mu_J = Mu_J_dict[snName][1];
        mu_H = -1.000000000; err_mu_H = -1.000000000;
        mu_K = -1.000000000; err_mu_K = -1.000000000;

        # distance mu array
        mu_np = np.array([ mu_Y, mu_J  ])
        # variance of distance mu array
        var_np = np.array([err_mu_Y**2., err_mu_J**2.])

        #--- Computing the best estimated total distance modulus ---
        for i5 in range(2):
            Numerator_Mu = Numerator_Mu + mu_np[i5]*np.sum(InvCov_ResidualMu_YJ[i5,:])
            Denominator_Mu = Denominator_Mu + np.sum(InvCov_ResidualMu_YJ[i5,:])

        # --- Computing the VARIANCE of the total distance modulus ---
        for i5 in range(2):
            Numerator_Var = Numerator_Var + var_np[i5]*(np.sum(InvCov_ResidualMu_YJ[i5,:])/Denominator_Mu)**2.

        # <<-------- End copy ---


    elif snName in MuData_K_names: # In K band  # SNe with distance mu from JK bands.

        # Begin copy  ---------->>
        count_JK = count_JK + 1
        count_all  = count_all  + 1
        print count_all, ', J-K__%r'%count_JK, ',', snName

        # The distance mu and uncertainties for the bands.
        mu_Y = -1.000000000;         err_mu_Y = -1.000000000;
        mu_J = Mu_J_dict[snName][0]; err_mu_J = Mu_J_dict[snName][1];
        mu_H = -1.000000000;         err_mu_H = -1.000000000;
        mu_K = Mu_K_dict[snName][0]; err_mu_K = Mu_K_dict[snName][1];

        # distance mu array
        mu_np = np.array([ mu_J, mu_K  ])
        # variance of distance mu array
        var_np = np.array([err_mu_J**2., err_mu_K**2.])

        #--- Computing the best estimated total distance modulus ---
        for i6 in range(2):
            Numerator_Mu = Numerator_Mu + mu_np[i6]*np.sum(InvCov_ResidualMu_JK[i6,:])
            Denominator_Mu = Denominator_Mu + np.sum(InvCov_ResidualMu_JK[i6,:])

        # --- Computing the VARIANCE of the total distance modulus ---
        for i6 in range(2):
            Numerator_Var = Numerator_Var + var_np[i6]*(np.sum(InvCov_ResidualMu_JK[i6,:])/Denominator_Mu)**2.
        # <<-------- End copy


    else: # In J band only  # SNe with distance mu from J band only.

        # Begin copy  ---------->>
        count_J = count_J + 1
        count_all  = count_all  + 1
        print count_all, ', J__%r'%count_J, ',', snName

        # The distance mu and uncertainties for the bands.
        mu_Y = -1.000000000;         err_mu_Y = -1.000000000;
        mu_J = Mu_J_dict[snName][0]; err_mu_J = Mu_J_dict[snName][1];
        mu_H = -1.000000000;         err_mu_H = -1.000000000;
        mu_K = -1.000000000;         err_mu_K = -1.000000000;

        JbandOnly = 1 # Flag to indicate that this SNe only have data in J band.
        # <<-------- End copy

    #-----------------------------------------------
    # TOTAL distance mu and its uncertainty

    if JbandOnly == 1:
        muTotal = mu_J
        error_muTotal = err_mu_J  # Uncertainty in the total distance mu
    else:
        muTotal = Numerator_Mu/Denominator_Mu  # TOTAL distance mu
        #old. error_muTotal = np.sqrt(1/Denominator_Mu)  # Uncertainty in the total distance mu
        error_muTotal= np.sqrt(Numerator_Var) # Uncertainty in the total distance mu

    mu_residualInt = muTotal - mu_LCDMInt # Residual distance modulus

    #---------------------------------------------

    # For the latex text file:

    if mu_Y == -1.000000000: mu_Y_latex = '...'
    else: mu_Y_latex = '$%.2f \pm %.2f$'%(mu_Y, err_mu_Y)

    if mu_J == -1.000000000: mu_J_latex = '...'
    else: mu_J_latex = '$%.2f \pm %.2f$'%(mu_J, err_mu_J)

    if mu_H == -1.000000000: mu_H_latex = '...'
    else: mu_H_latex = '$%.2f \pm %.2f$'%(mu_H, err_mu_H)

    if mu_K == -1.000000000: mu_K_latex = '...'
    else: mu_K_latex = '$%.2f \pm %.2f$'%(mu_K, err_mu_K)

    #----------
    # Creating a string with the exact name of the SN

    if snName[7]=='_': name = snName[0:7] # To read correctly, e.g., "sn2011B"
    elif snName[7]!='_':
        if is_number(snName[7]): name = snName[0:15] # To read correctly, e.g., "snf20080514-002"
        else: name = snName[0:8]   # To read correctly, e.g., "sn1998bu"

    #------Write to the text file -------------------

    textfile_1.write('%-30s  %.10f  %.10f  %.10f  %.10f  %14.10f  1.00000000    %.0f    %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n'%
                     (snName, zcmbInt, error_zcmbInt, muTotal, error_muTotal, mu_residualInt, SampleFlag,
                      mu_Y, err_mu_Y, mu_J, err_mu_J, mu_H, err_mu_H, mu_K, err_mu_K ))
    # Latex file:
    textfile_2.write('SN%-10s & %s & %s & %s & %s & %s & $%.2f \pm %.3f$ \\\\ \n'%
                     (name[2:], SampleFlatText, mu_Y_latex, mu_J_latex, mu_H_latex, mu_K_latex,
                      muTotal, error_muTotal))

    textfile_3.write('%-30s  %.10f  %.10f  %.10f  %.10f  %14.10f  1.00000000    %.0f    %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n'%
                     (name, zcmbInt, error_zcmbInt, muTotal, error_muTotal, mu_residualInt, SampleFlag,
                      mu_Y, err_mu_Y, mu_J, err_mu_J, mu_H, err_mu_H, mu_K, err_mu_K ))

    #<<--------END main loop ------------
textfile_1.write(text_line)
textfile_2.write('%'+text_line)
textfile_3.write(text_line)

#------------ BEGIN LOOP 2: Y band ------------
# SNe that are not in J band HD list but should be included.

textfile_1.write("# Y band only. \n")
textfile_2.write('%'+"# Y band only. \n")
textfile_3.write("# Y band only. \n")
for j1 in range(len(MuData_Y_names)):

    snName = MuData_Y_names[j1]

    # If this SNe is not in J band, then include it in the text file:
    if snName not in MuData_J_names:
        # Begin textfile variables --->>
        # These variables are not used for computing, they are just to be written in the text files
        zcmbInt = Mu_Y_dict[snName][4]
        error_zcmbInt = Mu_Y_dict[snName][5]
        chi2dofInt = Mu_Y_dict[snName][6]
        mu_LCDMInt = DistanceMu(zcmbInt, OmMFix, wFix, HoFix)

        # Sample flag
        SampleFlag = Mu_Y_dict[snName][3]
        if   SampleFlag ==1: SampleFlatText = 'CfA'
        elif SampleFlag ==2: SampleFlatText = 'CSP'
        elif SampleFlag ==3: SampleFlatText = 'Others'
        # <<---- End textfile variables

        # The distance mu and uncertainties for the bands.
        mu_Y = Mu_Y_dict[snName][0]; err_mu_Y = Mu_Y_dict[snName][1];
        mu_J = -1.000000000;         err_mu_J = -1.000000000;
        mu_H = -1.000000000;         err_mu_H = -1.000000000;
        mu_K = -1.000000000;         err_mu_K = -1.000000000;

        # TOTAL distance mu and its uncertainty
        muTotal = mu_Y
        error_muTotal = err_mu_Y  # Uncertainty in the total distance mu
        mu_residualInt = muTotal - mu_LCDMInt # Residual distance modulus

        # For the latex text file:
        if mu_Y == -1.000000000: mu_Y_latex = '...'
        else: mu_Y_latex = '$%.2f \pm %.2f$'%(mu_Y, err_mu_Y)
        if mu_J == -1.000000000: mu_J_latex = '...'
        else: mu_J_latex = '$%.2f \pm %.2f$'%(mu_J, err_mu_J)
        if mu_H == -1.000000000: mu_H_latex = '...'
        else: mu_H_latex = '$%.2f \pm %.2f$'%(mu_H, err_mu_H)
        if mu_K == -1.000000000: mu_K_latex = '...'
        else: mu_K_latex = '$%.2f \pm %.2f$'%(mu_K, err_mu_K)

        #----------
        # Creating a string with the exact name of the SN
        if snName[7]=='_': name = snName[0:7] # To read correctly, e.g., "sn2011B"
        elif snName[7]!='_':
            if is_number(snName[7]): name = snName[0:15] # To read correctly, e.g., "snf20080514-002"
            else: name = snName[0:8]   # To read correctly, e.g., "sn1998bu"

        textfile_1.write('%-30s  %.10f  %.10f  %.10f  %.10f  %14.10f  1.00000000    %.0f    %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n'%
                     (snName, zcmbInt, error_zcmbInt, muTotal, error_muTotal, mu_residualInt, SampleFlag,
                      mu_Y, err_mu_Y, mu_J, err_mu_J, mu_H, err_mu_H, mu_K, err_mu_K ))
        # Latex file:
        textfile_2.write('SN%-10s & %s & %s & %s & %s & %s & $%.2f \pm %.3f$ \\\\ \n'%
                     (name[2:], SampleFlatText, mu_Y_latex, mu_J_latex, mu_H_latex, mu_K_latex,
                      muTotal, error_muTotal))

        textfile_3.write('%-30s  %.10f  %.10f  %.10f  %.10f  %14.10f  1.00000000    %.0f    %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n'%
                     (name, zcmbInt, error_zcmbInt, muTotal, error_muTotal, mu_residualInt, SampleFlag,
                      mu_Y, err_mu_Y, mu_J, err_mu_J, mu_H, err_mu_H, mu_K, err_mu_K ))

        count_all += 1; count_Y += 1
        print count_all, ', Y__%r'%count_Y, ',', snName

    #-------------END LOOP 2: Y band -----------------------------
textfile_1.write(text_line)
textfile_2.write('%'+text_line)
textfile_3.write(text_line)

#------------ BEGIN LOOP 3: K band ------------
# Copy/paste of "BEGIN LOOP 2: Y band" but just replacing:
# "_Y_" --> "_K_"
# "j1" --> "j2"

textfile_1.write("# K band only. \n")
textfile_2.write('%'+"# K band only. \n")
textfile_3.write("# K band only. \n")

# SNe that are not in J band HD list but should be included.
for j2 in range(len(MuData_K_names)):

    snName = MuData_K_names[j2]

    # If this SNe is not in J band, then include it in the text file:
    if snName not in MuData_J_names:
        # Begin textfile variables --->>
        # These variables are not used for computing, they are just to be written in the text files
        zcmbInt = Mu_K_dict[snName][4]
        error_zcmbInt = Mu_K_dict[snName][5]
        chi2dofInt = Mu_K_dict[snName][6]
        mu_LCDMInt = DistanceMu(zcmbInt, OmMFix, wFix, HoFix)

        # Sample flag
        SampleFlag = Mu_K_dict[snName][3]
        if   SampleFlag ==1: SampleFlatText = 'CfA'
        elif SampleFlag ==2: SampleFlatText = 'CSP'
        elif SampleFlag ==3: SampleFlatText = 'Others'
        # <<---- End textfile variables

        # The distance mu and uncertainties for the bands.
        mu_Y = -1.000000000;         err_mu_Y = -1.000000000;
        mu_J = -1.000000000;         err_mu_J = -1.000000000;
        mu_H = -1.000000000;         err_mu_H = -1.000000000;
        mu_K = Mu_K_dict[snName][0]; err_mu_K = Mu_K_dict[snName][1];

        # TOTAL distance mu and its uncertainty
        muTotal = mu_K
        error_muTotal = err_mu_K  # Uncertainty in the total distance mu
        mu_residualInt = muTotal - mu_LCDMInt # Residual distance modulus

        # For the latex text file:
        if mu_Y == -1.000000000: mu_Y_latex = '...'
        else: mu_Y_latex = '$%.2f \pm %.2f$'%(mu_Y, err_mu_Y)
        if mu_J == -1.000000000: mu_J_latex = '...'
        else: mu_J_latex = '$%.2f \pm %.2f$'%(mu_J, err_mu_J)
        if mu_H == -1.000000000: mu_H_latex = '...'
        else: mu_H_latex = '$%.2f \pm %.2f$'%(mu_H, err_mu_H)
        if mu_K == -1.000000000: mu_K_latex = '...'
        else: mu_K_latex = '$%.2f \pm %.2f$'%(mu_K, err_mu_K)

        #----------
        # Creating a string with the exact name of the SN
        if snName[7]=='_': name = snName[0:7] # To read correctly, e.g., "sn2011B"
        elif snName[7]!='_':
            if is_number(snName[7]): name = snName[0:15] # To read correctly, e.g., "snf20080514-002"
            else: name = snName[0:8]   # To read correctly, e.g., "sn1998bu"

        textfile_1.write('%-30s  %.10f  %.10f  %.10f  %.10f  %14.10f  1.00000000    %.0f    %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n'%
                     (snName, zcmbInt, error_zcmbInt, muTotal, error_muTotal, mu_residualInt, SampleFlag,
                      mu_Y, err_mu_Y, mu_J, err_mu_J, mu_H, err_mu_H, mu_K, err_mu_K ))
        # Latex file:
        textfile_2.write('SN%-10s & %s & %s & %s & %s & %s & $%.2f \pm %.3f$ \\\\ \n'%
                     (name[2:], SampleFlatText, mu_Y_latex, mu_J_latex, mu_H_latex, mu_K_latex,
                      muTotal, error_muTotal))

        textfile_3.write('%-30s  %.10f  %.10f  %.10f  %.10f  %14.10f  1.00000000    %.0f    %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n'%
                     (name, zcmbInt, error_zcmbInt, muTotal, error_muTotal, mu_residualInt, SampleFlag,
                      mu_Y, err_mu_Y, mu_J, err_mu_J, mu_H, err_mu_H, mu_K, err_mu_K ))

        count_all += 1; count_K += 1
        print count_all, ', K__%r'%count_K, ',', snName

    #-------------END LOOP 3: K band -----------------------------
textfile_1.write(text_line)
textfile_2.write('%'+text_line)
textfile_3.write(text_line)

#------------ BEGIN LOOP 4: H band ------------
# Copy/paste of "BEGIN LOOP 2: Y band" but just replacing:
# "_Y_" --> "_K_"
# "j1" --> "j3"

textfile_1.write("# H band only. \n")
textfile_2.write('%'+"# H band only. \n")
textfile_3.write("# H band only. \n")

for j3 in range(len(MuData_H_names)):

    snName = MuData_H_names[j3]

    # If this SNe is not in J band, then include it in the text file:
    if snName not in MuData_J_names:
        # Begin textfile variables --->>
        # These variables are not used for computing, they are just to be written in the text files
        zcmbInt = Mu_H_dict[snName][4]
        error_zcmbInt = Mu_H_dict[snName][5]
        chi2dofInt = Mu_H_dict[snName][6]
        mu_LCDMInt = DistanceMu(zcmbInt, OmMFix, wFix, HoFix)

        # Sample flag
        SampleFlag = Mu_H_dict[snName][3]
        if   SampleFlag ==1: SampleFlatText = 'CfA'
        elif SampleFlag ==2: SampleFlatText = 'CSP'
        elif SampleFlag ==3: SampleFlatText = 'Others'
        # <<---- End textfile variables

        # The distance mu and uncertainties for the bands.
        mu_Y = -1.000000000;         err_mu_Y = -1.000000000;
        mu_J = -1.000000000;         err_mu_J = -1.000000000;
        mu_H = Mu_H_dict[snName][0]; err_mu_H = Mu_H_dict[snName][1];
        mu_K = -1.000000000;         err_mu_K = -1.000000000;

        # TOTAL distance mu and its uncertainty
        muTotal = mu_H
        error_muTotal = err_mu_H  # Uncertainty in the total distance mu
        mu_residualInt = muTotal - mu_LCDMInt # Residual distance modulus

        # For the latex text file:
        if mu_Y == -1.000000000: mu_Y_latex = '...'
        else: mu_Y_latex = '$%.2f \pm %.2f$'%(mu_Y, err_mu_Y)
        if mu_J == -1.000000000: mu_J_latex = '...'
        else: mu_J_latex = '$%.2f \pm %.2f$'%(mu_J, err_mu_J)
        if mu_H == -1.000000000: mu_H_latex = '...'
        else: mu_H_latex = '$%.2f \pm %.2f$'%(mu_H, err_mu_H)
        if mu_K == -1.000000000: mu_K_latex = '...'
        else: mu_K_latex = '$%.2f \pm %.2f$'%(mu_K, err_mu_K)

        #----------
        # Creating a string with the exact name of the SN
        if snName[7]=='_': name = snName[0:7] # To read correctly, e.g., "sn2011B"
        elif snName[7]!='_':
            if is_number(snName[7]): name = snName[0:15] # To read correctly, e.g., "snf20080514-002"
            else: name = snName[0:8]   # To read correctly, e.g., "sn1998bu"

        textfile_1.write('%-30s  %.10f  %.10f  %.10f  %.10f  %14.10f  1.00000000    %.0f    %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n'%
                     (snName, zcmbInt, error_zcmbInt, muTotal, error_muTotal, mu_residualInt, SampleFlag,
                      mu_Y, err_mu_Y, mu_J, err_mu_J, mu_H, err_mu_H, mu_K, err_mu_K ))
        # Latex file:
        textfile_2.write('SN%-10s & %s & %s & %s & %s & %s & $%.2f \pm %.3f$ \\\\ \n'%
                     (name[2:], SampleFlatText, mu_Y_latex, mu_J_latex, mu_H_latex, mu_K_latex,
                      muTotal, error_muTotal))

        textfile_3.write('%-30s  %.10f  %.10f  %.10f  %.10f  %14.10f  1.00000000    %.0f    %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n'%
                     (name, zcmbInt, error_zcmbInt, muTotal, error_muTotal, mu_residualInt, SampleFlag,
                      mu_Y, err_mu_Y, mu_J, err_mu_J, mu_H, err_mu_H, mu_K, err_mu_K ))

        count_all += 1; count_H += 1
        print count_all, ', H__%r'%count_H, ',', snName

    #-------------END LOOP 4: H band -----------------------------

#-----------------------------------------------------------------------------80
# Close text files

textfile_1.write(text_line)
textfile_1.write("# %s SNe in this list. \n"%count_all)
textfile_1.write("# %s SNe with data in YJH  bands. \n"%count_YJH)
textfile_1.write("# %s SNe with data in JHK  bands. \n"%count_JHK)
textfile_1.write("# %s SNe with data in YJHK bands. \n"%count_YJHK)
textfile_1.write("# %s SNe with data in YJ   bands. \n"%count_YJ)
textfile_1.write("# %s SNe with data in JH   bands. \n"%count_JH)
textfile_1.write("# %s SNe with data in JK   bands. \n"%count_JK)
textfile_1.write("# %s SNe with data in J    bands only. \n"%count_J)
textfile_1.write("# %s SNe with data in Y    bands only. \n"%count_Y)
textfile_1.write("# %s SNe with data in K    bands only. \n"%count_K)
textfile_1.write("# %s SNe with data in H    bands only. \n"%count_H)

textfile_3.write(text_line)
textfile_3.write("# %s SNe in this list. \n"%count_all)
textfile_3.write("# %s SNe with data in YJH  bands. \n"%count_YJH)
textfile_3.write("# %s SNe with data in JHK  bands. \n"%count_JHK)
textfile_3.write("# %s SNe with data in YJHK bands. \n"%count_YJHK)
textfile_3.write("# %s SNe with data in YJ   bands. \n"%count_YJ)
textfile_3.write("# %s SNe with data in JH   bands. \n"%count_JH)
textfile_3.write("# %s SNe with data in JK   bands. \n"%count_JK)
textfile_3.write("# %s SNe with data in J    bands only. \n"%count_J)
textfile_3.write("# %s SNe with data in Y    bands only. \n"%count_Y)
textfile_3.write("# %s SNe with data in K    bands only. \n"%count_K)
textfile_3.write("# %s SNe with data in H    bands only. \n"%count_H)

textfile_1.close(); textfile_2.close(); textfile_3.close();

textfile_1.close(); textfile_2.close(); textfile_3.close();

##############################################################################80

# ### YJHK strict.
#
# This cell produce also the final latex table that will be used in the paper. The "14_Latex_Tables_Figures_v2_5.ipynb" script doesn't do it.

# Only if the YJHK cov matrix was defined:
if count_YJHK >= MinSNeForYJHK:

    # Loop 2: SNe with distance modulus in YJHK bands.
    # (Ths loop is a copy from the main loop, but with small modifications)

    # Open the writting to a text file
    textfile_1 = open(DirSaveOutput+'Table_TotalMu_YJHK_.txt', 'w')
    textfile_2 = open(DirSaveOutput+'Table_TotalMu_YJHK_Latex_.txt', 'w')

    now = datetime.datetime.now() # Read the time and date right now
    text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd)")
    text_Date   = '# On date: %s \n'%text_timenow
    text_Author = '# Data table created by: Arturo Avelino \n'
    text_script = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
    text_line = '#'+'-'*60 + '\n'

    textfile_1.write(text_line)
    textfile_1.write(text_Author); textfile_1.write(text_Date); textfile_1.write(text_script);
    textfile_1.write(text_line)

    text_1 = '#   SN name                     z_CMB         error_zcmb       mu_total      error_mu       \
    mu_residual    chi2_dof  Sample     mu_Y         error_mu_Y       mu_J          error_mu_J      mu_H           \
    error_mu_H     mu_K            error_mu_K    Unused columns \n'
    textfile_1.write(text_1);
    # textfile_2.write('% '+text_1)

    # Some counters
    count_YJHK = 0;
    count_YJH = 0; count_YJ = 0; count_JH = 0;
    count_JHK = 0; count_JK = 0;
    count_J=0; count_all = 0

    # Loop over SNe in J band (it contains almost all the SNe)
    for i in range(len(MuData_J_names)):

        # Resetting variables
        snName = '';
        zcmbInt = 0; error_zcmbInt=0; chi2dofInt=0;
        mu_LCDMInt=0; mu_residualInt=0;
        SampleFlag = 0; SampleFlatText='';
        mu_Y = 0; mu_J = 0; mu_H = 0; mu_K = 0;
        err_mu_Y=0; err_mu_J=0; err_mu_H=0; err_mu_K=0;
        mu_np  = np.zeros(3); errorMu_np = np.zeros(3);
        weight_InvError2Mu_np=np.zeros(3);
        var_np = np.zeros(3);
        CovMatrix_Mu = np.zeros([3,3]); InvCovMatrix_Mu = np.zeros([3,3]);
        Numerator_Mu = 0; Denominator_Mu = 0;
        Numerator_Var= 0;
        muTotal = 0; error_muTotal = 0;
        JbandOnly = 0;

        #--------------------------------

        # Defining variables

        snName = MuData_J_names[i]

        # Begin textfile variables --->>
        # These variables are not used for computing, they are just to be written in the text files
        zcmbInt = Mu_J_dict[snName][4]
        error_zcmbInt = Mu_J_dict[snName][5]
        chi2dofInt = Mu_J_dict[snName][6]
        mu_LCDMInt = DistanceMu(zcmbInt, OmMFix, wFix, HoFix)

        # Sample flag
        SampleFlag = Mu_J_dict[snName][3]
        if   SampleFlag ==1: SampleFlatText = 'CfA'
        elif SampleFlag ==2: SampleFlatText = 'CSP'
        elif SampleFlag ==3: SampleFlatText = 'Others'
        # <<---- End textfile variables

        if snName in MuData_H_names: # In H band
            if snName in MuData_K_names: # In K band
                if snName in MuData_Y_names: # In Y band  # SNe with distance mu from YJHK bands.

                    # Begin copy (this is the original) ---------->>
                    count_YJHK = count_YJHK + 1
                    count_all  = count_all  + 1
                    print count_all, ', J-H-K-Y__%r'%count_YJHK, ',', snName

                    # The distance mu and uncertainties for the bands.
                    mu_Y = Mu_Y_dict[snName][0]; err_mu_Y = Mu_Y_dict[snName][1];
                    mu_J = Mu_J_dict[snName][0]; err_mu_J = Mu_J_dict[snName][1];
                    mu_H = Mu_H_dict[snName][0]; err_mu_H = Mu_H_dict[snName][1];
                    mu_K = Mu_K_dict[snName][0]; err_mu_K = Mu_K_dict[snName][1];

                    # distance mu array
                    mu_np = np.array([ mu_Y, mu_J, mu_H, mu_K  ])
                    # variance of distance mu array
                    var_np = np.array([err_mu_Y**2., err_mu_J**2., err_mu_H**2., err_mu_K**2.])

                    #--- Computing the best estimated total distance modulus ---
                    for i in range(4):
                        Numerator_Mu = Numerator_Mu + mu_np[i]*np.sum(InvCov_ResidualMu_YJHK[i,:])
                        Denominator_Mu = Denominator_Mu + np.sum(InvCov_ResidualMu_YJHK[i,:])

                    # --- Computing the VARIANCE of the total distance modulus ---
                    for i in range(4):
                        Numerator_Var = Numerator_Var + var_np[i]*(np.sum(InvCov_ResidualMu_YJHK[i,:])/Denominator_Mu )**2.

                    # <<-------- End copy (this is the original)

                    #-------------------------------------------------------------

                    # TOTAL distance mu and its uncertainty

                    muTotal = Numerator_Mu/Denominator_Mu  # TOTAL distance mu

                    # Uncertainty in the total distance mu:
                    error_muTotal = np.sqrt(Numerator_Var)

                    mu_residualInt = muTotal - mu_LCDMInt # Residual distance modulus

                    #---------------------------------------------------------------------------------------

                    # For the latex text file:

                    if mu_Y == -1.000000000: mu_Y_latex = '...'
                    else: mu_Y_latex = '$%.2f \pm %.2f$'%(mu_Y, err_mu_Y)

                    if mu_J == -1.000000000: mu_J_latex = '...'
                    else: mu_J_latex = '$%.2f \pm %.2f$'%(mu_J, err_mu_J)

                    if mu_H == -1.000000000: mu_H_latex = '...'
                    else: mu_H_latex = '$%.2f \pm %.2f$'%(mu_H, err_mu_H)

                    if mu_K == -1.000000000: mu_K_latex = '...'
                    else: mu_K_latex = '$%.2f \pm %.2f$'%(mu_K, err_mu_K)

                    #----------
                    # Creating a string with the exact name of the SN

                    if snName[7]=='_': name = snName[0:7] # To read correctly "sn2011B"
                    elif snName[7]!='_':
                        if is_number(snName[7]): name = snName[0:15] # To read correctly "snf20080514-002"
                        else: name = snName[0:8]   # To read correctly "sn1998bu"

                    #---------------------------------------------------------------------------------------

                    textfile_1.write('%-30s  %.10f  %.10f  %.10f  %.10f  %14.10f  1.00000000    %.0f    %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n'%
                                     (snName, zcmbInt, error_zcmbInt, muTotal, error_muTotal, mu_residualInt, SampleFlag,
                                      mu_Y, err_mu_Y, mu_J, err_mu_J, mu_H, err_mu_H, mu_K, err_mu_K ))
                    # Latex file:
                    textfile_2.write('SN%-10s & %s & %s & %s & %s & %s & $%.2f \pm %.3f$ \\\\ \n'%
                                     (name[2:], SampleFlatText, mu_Y_latex, mu_J_latex, mu_H_latex, mu_K_latex,
                                      muTotal, error_muTotal))

    textfile_1.write(text_line)
    textfile_1.write("# %s SNe with data in YJHK bands \n"%count_YJHK)

    textfile_1.close()
    textfile_2.close()

else: print "# YJHK covariance matrix from less or equal to 4 SNe, \
so I didn't define it."

textfile_1.close(); textfile_2.close();

##############################################################################80

# ### YJH only
#
# This cell produce also the final latex table that will be used in the paper. The "14_Latex_Tables_Figures_v2_5.ipynb" script doesn't do it.

# Loop 3: SNe with distance modulus in YJH bands.
# (Ths loop is a copy from the main loop, but with small modifications)

# SNe with YJH data

# Open the writting to a text file
textfile_1 = open(DirSaveOutput+'Table_TotalMu_YJH_.txt', 'w')
textfile_2 = open(DirSaveOutput+'Table_TotalMu_YJH_Latex_.txt', 'w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd)")
text_Date   = '# On date: %s \n'%text_timenow
text_Author = '# Data table created by: Arturo Avelino \n'
text_script = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
text_line = '#'+'-'*60 + '\n'

textfile_1.write(text_line)
textfile_1.write(text_Author); textfile_1.write(text_Date); textfile_1.write(text_script);
textfile_1.write(text_line)

text_1 = '#   SN name                     z_CMB         error_zcmb       mu_total      error_mu       \
mu_residual    chi2_dof  Sample     mu_Y         error_mu_Y       mu_J          error_mu_J      mu_H           \
error_mu_H     mu_K            error_mu_K    Unused columns \n'
textfile_1.write(text_1);
# textfile_2.write('% '+text_1)

# Some counters
count_YJHK = 0;
count_YJH = 0; count_YJ = 0; count_JH = 0;
count_JHK = 0; count_JK = 0;
count_J=0; count_all = 0

# Loop over SNe in J band (it contains all the SNe)
for i in range(len(MuData_J_names)):

    # Resetting variables
    snName = '';
    zcmbInt = 0; error_zcmbInt=0; chi2dofInt=0;
    mu_LCDMInt=0; mu_residualInt=0;
    SampleFlag = 0; SampleFlatText='';
    mu_Y = 0; mu_J = 0; mu_H = 0; mu_K = 0;
    err_mu_Y=0; err_mu_J=0; err_mu_H=0; err_mu_K=0;
    mu_np  = np.zeros(3); errorMu_np = np.zeros(3);
    weight_InvError2Mu_np=np.zeros(3);
    CovMatrix_Mu = np.zeros([3,3]); InvCovMatrix_Mu = np.zeros([3,3]);
    var_np  = np.zeros(3);
    Numerator_Mu = 0; Denominator_Mu = 0
    Numerator_Var= 0;
    muTotal = 0; error_muTotal = 0;
    JbandOnly = 0;

    #--------------------------------

    # Defining variables

    snName = MuData_J_names[i]

    # Begin textfile variables --->>
    # These variables are not used for computing,
    # they are just to be written in the text files
    zcmbInt = Mu_J_dict[snName][4]
    error_zcmbInt = Mu_J_dict[snName][5]
    chi2dofInt = Mu_J_dict[snName][6]
    mu_LCDMInt = DistanceMu(zcmbInt, OmMFix, wFix, HoFix)

    # Sample flag
    SampleFlag = Mu_J_dict[snName][3]
    if   SampleFlag ==1: SampleFlatText = 'CfA'
    elif SampleFlag ==2: SampleFlatText = 'CSP'
    elif SampleFlag ==3: SampleFlatText = 'Others'
    # <<---- End textfile variables


    if snName in MuData_H_names: # In H band
        if snName in MuData_Y_names: # In Y band

            # SNe with distance mu from YJH bands.

            # Begin copy  ---------->>
            count_YJH = count_YJH + 1
            count_all  = count_all  + 1
            print count_all, ', J-H-Y__%r'%count_YJH, ',', snName

            # The distance mu and uncertainties for the bands.
            mu_Y = Mu_Y_dict[snName][0]; err_mu_Y = Mu_Y_dict[snName][1];
            mu_J = Mu_J_dict[snName][0]; err_mu_J = Mu_J_dict[snName][1];
            mu_H = Mu_H_dict[snName][0]; err_mu_H = Mu_H_dict[snName][1];
            mu_K = -1.000000000;         err_mu_K = -1.000000000;

            # distance mu array
            mu_np = np.array([ mu_Y, mu_J, mu_H  ])
            # variance of distance mu array
            var_np = np.array([err_mu_Y**2., err_mu_J**2., err_mu_H**2.])

            #--- Computing the best estimated total distance modulus ---
            for i in range(3):
                Numerator_Mu = Numerator_Mu + mu_np[i]*np.sum(InvCov_ResidualMu_YJH[i,:])
                Denominator_Mu = Denominator_Mu + np.sum(InvCov_ResidualMu_YJH[i,:])

            # --- Computing the VARIANCE of the total distance modulus ---
            for i in range(3):
                Numerator_Var = Numerator_Var + var_np[i]*(np.sum(InvCov_ResidualMu_YJH[i,:])/Denominator_Mu )**2.

            # <<-------- End copy

            #-------------------------------------------------------------
            # TOTAL distance mu and its uncertainty

            muTotal = Numerator_Mu/Denominator_Mu  # TOTAL distance mu

            # Uncertainty in the total distance mu
            error_muTotal = np.sqrt(Numerator_Var)

            mu_residualInt = muTotal - mu_LCDMInt # Residual distance modulus

            #---------------------------------------------------------------------------------------

            # For the latex text file:

            if mu_Y == -1.000000000: mu_Y_latex = '...'
            else: mu_Y_latex = '$%.2f \pm %.2f$'%(mu_Y, err_mu_Y)

            if mu_J == -1.000000000: mu_J_latex = '...'
            else: mu_J_latex = '$%.2f \pm %.2f$'%(mu_J, err_mu_J)

            if mu_H == -1.000000000: mu_H_latex = '...'
            else: mu_H_latex = '$%.2f \pm %.2f$'%(mu_H, err_mu_H)

            if mu_K == -1.000000000: mu_K_latex = '...'
            else: mu_K_latex = '$%.2f \pm %.2f$'%(mu_K, err_mu_K)

            #----------
            # Creating a string with the exact name of the SN

            if snName[7]=='_': name = snName[0:7] # To read correctly "sn2011B"
            elif snName[7]!='_':
                if is_number(snName[7]): name = snName[0:15] # To read correctly "snf20080514-002"
                else: name = snName[0:8]   # To read correctly "sn1998bu"

            #---------------------------------------------------------------------------------------

            textfile_1.write('%-30s  %.10f  %.10f  %.10f  %.10f  %14.10f  1.00000000    %.0f    %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n'%
                             (snName, zcmbInt, error_zcmbInt, muTotal, error_muTotal, mu_residualInt, SampleFlag,
                              mu_Y, err_mu_Y, mu_J, err_mu_J, mu_H, err_mu_H, mu_K, err_mu_K ))
            # Latex file:
            textfile_2.write('SN%-10s & %s & %s & %s & %s & %s & $%.2f \pm %.3f$ \\\\ \n'%
                             (name[2:], SampleFlatText, mu_Y_latex, mu_J_latex, mu_H_latex, mu_K_latex,
                              muTotal, error_muTotal))

textfile_1.write(text_line)
textfile_1.write("# %s SNe with data in YJH  bands \n"%count_YJH)

textfile_1.close()
textfile_2.close()

textfile_1.close(); textfile_2.close();

##############################################################################80

# ### JH only
#
# This cell produce also the final latex table that will be used in the paper.
# The "14_Latex_Tables_Figures_v2_5.ipynb" script doesn't do it.

# Loop 3: SNe with distance modulus in JH bands.
# (Ths loop is a copy from the main loop, but with small modifications)

# SNe with JH data

# Open the writting to a text file
textfile_1 = open(DirSaveOutput+'Table_TotalMu_JH_.txt', 'w')
textfile_2 = open(DirSaveOutput+'Table_TotalMu_JH_Latex_.txt', 'w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd)")
text_Date   = '# On date: %s \n'%text_timenow
text_Author = '# Data table created by: Arturo Avelino \n'
text_script = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
text_line = '#'+'-'*60 + '\n'

textfile_1.write(text_line)
textfile_1.write(text_Author); textfile_1.write(text_Date); textfile_1.write(text_script);
textfile_1.write(text_line)

text_1 = '#   SN name                     z_CMB         error_zcmb       mu_total      error_mu       \
mu_residual    chi2_dof  Sample     mu_Y         error_mu_Y       mu_J          error_mu_J      mu_H           \
error_mu_H     mu_K            error_mu_K    Unused columns \n'
textfile_1.write(text_1);
# textfile_2.write('% '+text_1)

# Some counters
count_YJHK = 0;
count_YJH = 0; count_YJ = 0; count_JH = 0;
count_JHK = 0; count_JK = 0;
count_J=0; count_all = 0

# Loop over SNe in J band (it contains all the SNe)
for i in range(len(MuData_J_names)):

    # Resetting variables
    snName = '';
    zcmbInt = 0; error_zcmbInt=0; chi2dofInt=0;
    mu_LCDMInt=0; mu_residualInt=0;
    SampleFlag = 0; SampleFlatText='';
    mu_Y = 0; mu_J = 0; mu_H = 0; mu_K = 0;
    err_mu_Y=0; err_mu_J=0; err_mu_H=0; err_mu_K=0;
    mu_np  = np.zeros(3); errorMu_np = np.zeros(3);
    weight_InvError2Mu_np=np.zeros(3);
    CovMatrix_Mu = np.zeros([3,3]); InvCovMatrix_Mu = np.zeros([3,3]);
    var_np  = np.zeros(3);
    Numerator_Mu = 0; Denominator_Mu = 0
    Numerator_Var= 0;
    muTotal = 0; error_muTotal = 0;
    JbandOnly = 0;

    #--------------------------------

    # Defining variables

    snName = MuData_J_names[i]

    # Begin textfile variables --->>
    # These variables are not used for computing,
    # they are just to be written in the text files
    zcmbInt = Mu_J_dict[snName][4]
    error_zcmbInt = Mu_J_dict[snName][5]
    chi2dofInt = Mu_J_dict[snName][6]
    mu_LCDMInt = DistanceMu(zcmbInt, OmMFix, wFix, HoFix)

    # Sample flag
    SampleFlag = Mu_J_dict[snName][3]
    if   SampleFlag ==1: SampleFlatText = 'CfA'
    elif SampleFlag ==2: SampleFlatText = 'CSP'
    elif SampleFlag ==3: SampleFlatText = 'Others'
    # <<---- End textfile variables


    if snName in MuData_H_names: # In H band

        # SNe with distance mu from YJH bands.

        # Begin copy  ---------->>
        count_JH = count_JH + 1
        count_all  = count_all  + 1
        print count_all, ', J-H__%r'%count_JH, ',', snName

        # The distance mu and uncertainties for the bands.
        mu_Y = -1.000000000;         err_mu_Y = -1.000000000;
        mu_J = Mu_J_dict[snName][0]; err_mu_J = Mu_J_dict[snName][1];
        mu_H = Mu_H_dict[snName][0]; err_mu_H = Mu_H_dict[snName][1];
        mu_K = -1.000000000;         err_mu_K = -1.000000000;

        # distance mu array
        mu_np = np.array([ mu_J, mu_H  ])
        # variance of distance mu array
        var_np = np.array([ err_mu_J**2., err_mu_H**2.])

        #--- Computing the best estimated total distance modulus ---
        for i in range(2):
            Numerator_Mu = Numerator_Mu + mu_np[i]*np.sum(InvCov_ResidualMu_JH[i,:])
            Denominator_Mu = Denominator_Mu + np.sum(InvCov_ResidualMu_JH[i,:])

        # --- Computing the VARIANCE of the total distance modulus ---
        for i in range(2):
            Numerator_Var = Numerator_Var + var_np[i]*(np.sum(InvCov_ResidualMu_JH[i,:])/Denominator_Mu )**2.

        # <<-------- End copy

        #-------------------------------------------------------------
        # TOTAL distance mu and its uncertainty

        muTotal = Numerator_Mu/Denominator_Mu  # TOTAL distance mu

        # Uncertainty in the total distance mu
        error_muTotal = np.sqrt(Numerator_Var)

        mu_residualInt = muTotal - mu_LCDMInt # Residual distance modulus

        #-------------------------------------------

        # For the latex text file:

        if mu_Y == -1.000000000: mu_Y_latex = '...'
        else: mu_Y_latex = '$%.2f \pm %.2f$'%(mu_Y, err_mu_Y)

        if mu_J == -1.000000000: mu_J_latex = '...'
        else: mu_J_latex = '$%.2f \pm %.2f$'%(mu_J, err_mu_J)

        if mu_H == -1.000000000: mu_H_latex = '...'
        else: mu_H_latex = '$%.2f \pm %.2f$'%(mu_H, err_mu_H)

        if mu_K == -1.000000000: mu_K_latex = '...'
        else: mu_K_latex = '$%.2f \pm %.2f$'%(mu_K, err_mu_K)

        #----------
        # Creating a string with the exact name of the SN

        if snName[7]=='_': name = snName[0:7] # To read correctly "sn2011B"
        elif snName[7]!='_':
            if is_number(snName[7]): name = snName[0:15] # To read correctly "snf20080514-002"
            else: name = snName[0:8]   # To read correctly "sn1998bu"

        #---------------------------------------------------------------------------------------

        textfile_1.write('%-30s  %.10f  %.10f  %.10f  %.10f  %14.10f  1.00000000    %.0f    %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n'%
                         (snName, zcmbInt, error_zcmbInt, muTotal, error_muTotal, mu_residualInt, SampleFlag,
                          mu_Y, err_mu_Y, mu_J, err_mu_J, mu_H, err_mu_H, mu_K, err_mu_K ))
        # Latex file:
        textfile_2.write('SN%-10s & %s & %s & %s & %s & %s & $%.2f \pm %.3f$ \\\\ \n'%
                         (name[2:], SampleFlatText, mu_Y_latex, mu_J_latex, mu_H_latex, mu_K_latex,
                          muTotal, error_muTotal))

textfile_1.write(text_line)
textfile_1.write("# %s SNe with data in JH  bands \n"%count_JH)

textfile_1.close()
textfile_2.close()

textfile_1.close(); textfile_2.close(); textfile_1.close(); textfile_2.close();

##############################################################################80

# ### JHK
#
# This cell produce also the final latex table that will be used in the paper. The "14_Latex_Tables_Figures_v2_5.ipynb" script doesn't do it.

# Loop 4: SNe with distance modulus in JHK bands.
# (This loop is a copy from the main loop, but with small modifications)

# Open the writting to a text file
textfile_1 = open(DirSaveOutput+'Table_TotalMu_JHK_.txt', 'w')
textfile_2 = open(DirSaveOutput+'Table_TotalMu_JHK_Latex_.txt', 'w')

now = datetime.datetime.now() # Read the time and date right now
text_timenow = now.strftime("%Y-%m-%d (yyyy-mm-dd)")
text_Date   = '# On date: %s \n'%text_timenow
text_Author = '# Data table created by: Arturo Avelino \n'
text_script = '# Script used: %s (version %s | last update: %s)\n'%(
        code_name, version_code, last_update)
text_line = '#'+'-'*60 + '\n'

textfile_1.write(text_line)
textfile_1.write(text_Author); textfile_1.write(text_Date); textfile_1.write(text_script);
textfile_1.write(text_line)

text_1 = '#   SN name                     z_CMB         error_zcmb       mu_total      error_mu       \
mu_residual    chi2_dof  Sample     mu_Y         error_mu_Y       mu_J          error_mu_J      mu_H           \
error_mu_H     mu_K            error_mu_K    Unused columns \n'
textfile_1.write(text_1);
# textfile_2.write('% '+text_1)

# Some counters
count_YJHK = 0;
count_YJH = 0; count_YJ = 0; count_JH = 0;
count_JHK = 0; count_JK = 0;
count_J=0; count_all = 0

for i in range(len(MuData_J_names)): # Loop over SNe in J band (it contains all the SNe)

    # Resetting variables
    snName = '';
    zcmbInt = 0; error_zcmbInt=0; chi2dofInt=0;
    mu_LCDMInt=0; mu_residualInt=0;
    SampleFlag = 0; SampleFlatText='';
    mu_Y = 0; mu_J = 0; mu_H = 0; mu_K = 0;
    err_mu_Y=0; err_mu_J=0; err_mu_H=0; err_mu_K=0;
    mu_np  = np.zeros(3); errorMu_np = np.zeros(3);
    var_np  = np.zeros(3);
    weight_InvError2Mu_np=np.zeros(3);
    CovMatrix_Mu = np.zeros([3,3]); InvCovMatrix_Mu = np.zeros([3,3]);
    Numerator_Mu = 0; Denominator_Mu = 0;
    Numerator_Var = 0;
    muTotal = 0; error_muTotal = 0;
    JbandOnly = 0;

    #--------------------------------

    # Defining variables

    snName = MuData_J_names[i]

    # Begin textfile variables --->>
    # These variables are not used for computing, they are just to be written in the text files
    zcmbInt = Mu_J_dict[snName][4]
    error_zcmbInt = Mu_J_dict[snName][5]
    chi2dofInt = Mu_J_dict[snName][6]
    mu_LCDMInt = DistanceMu(zcmbInt, OmMFix, wFix, HoFix)

    # Sample flag
    SampleFlag = Mu_J_dict[snName][3]
    if   SampleFlag ==1: SampleFlatText = 'CfA'
    elif SampleFlag ==2: SampleFlatText = 'CSP'
    elif SampleFlag ==3: SampleFlatText = 'Others'
    # <<---- End textfile variables


    if snName in MuData_H_names: # In H band
        if snName in MuData_K_names: # In K band # SNe with distance mu from JHK bands only

            # Begin copy  ---------->>
            count_JHK = count_JHK + 1
            count_all  = count_all  + 1
            print count_all, ', J-H-K__%r'%count_JHK, ',', snName

            # The distance mu and uncertainties for the bands.
            mu_Y = -1.000000000;         err_mu_Y = -1.000000000;
            mu_J = Mu_J_dict[snName][0]; err_mu_J = Mu_J_dict[snName][1];
            mu_H = Mu_H_dict[snName][0]; err_mu_H = Mu_H_dict[snName][1];
            mu_K = Mu_K_dict[snName][0]; err_mu_K = Mu_K_dict[snName][1];

            # distance mu array
            mu_np = np.array([ mu_J, mu_H, mu_K ])
            # variance of distance mu array
            var_np = np.array([err_mu_J**2., err_mu_H**2., err_mu_K**2.])

            #--- Computing the best estimated total distance modulus ---
            for i in range(3):
                Numerator_Mu = Numerator_Mu + mu_np[i]*np.sum(InvCov_ResidualMu_JHK[i,:])
                Denominator_Mu = Denominator_Mu + np.sum(InvCov_ResidualMu_JHK[i,:])

            # --- Computing the VARIANCE of the total distance modulus ---
            for i in range(3):
                Numerator_Var = Numerator_Var + var_np[i]*(np.sum(InvCov_ResidualMu_JHK[i,:])/Denominator_Mu )**2.

            # <<-------- End copy

            #-------------------------------------------------------------
            # TOTAL distance mu and its uncertainty

            muTotal = Numerator_Mu/Denominator_Mu  # TOTAL distance mu

            # Uncertainty in the total distance mu:
            error_muTotal = np.sqrt(Numerator_Var)

            mu_residualInt = muTotal - mu_LCDMInt # Residual distance modulus

            #---------------------------------------------------------------------------------------

            # For the latex text file:

            if mu_Y == -1.000000000: mu_Y_latex = '...'
            else: mu_Y_latex = '$%.2f \pm %.2f$'%(mu_Y, err_mu_Y)

            if mu_J == -1.000000000: mu_J_latex = '...'
            else: mu_J_latex = '$%.2f \pm %.2f$'%(mu_J, err_mu_J)

            if mu_H == -1.000000000: mu_H_latex = '...'
            else: mu_H_latex = '$%.2f \pm %.2f$'%(mu_H, err_mu_H)

            if mu_K == -1.000000000: mu_K_latex = '...'
            else: mu_K_latex = '$%.2f \pm %.2f$'%(mu_K, err_mu_K)

            #----------
            # Creating a string with the exact name of the SN

            if snName[7]=='_': name = snName[0:7] # To read correctly "sn2011B"
            elif snName[7]!='_':
                if is_number(snName[7]): name = snName[0:15] # To read correctly "snf20080514-002"
                else: name = snName[0:8]   # To read correctly "sn1998bu"

            #---------------------------------------------------------------------------------------

            textfile_1.write('%-30s  %.10f  %.10f  %.10f  %.10f  %14.10f  1.00000000    %.0f    %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  %.10f  %14.10f  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 \n'%
                             (snName, zcmbInt, error_zcmbInt, muTotal, error_muTotal, mu_residualInt, SampleFlag,
                              mu_Y, err_mu_Y, mu_J, err_mu_J, mu_H, err_mu_H, mu_K, err_mu_K ))
            # Latex file:
            textfile_2.write('SN%-10s & %s & %s & %s & %s & %s & $%.2f \pm %.3f$ \\\\ \n'%
                             (name[2:], SampleFlatText, mu_Y_latex, mu_J_latex, mu_H_latex, mu_K_latex,
                              muTotal, error_muTotal))

textfile_1.write(text_line)
textfile_1.write("# %s SNe with data in JHK  bands \n"%count_JHK)

textfile_1.close()
textfile_2.close()

textfile_1.close(); textfile_2.close();

print "# All done. "

