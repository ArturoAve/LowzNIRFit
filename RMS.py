# Functions to compute rms,  weighted RMS, and their uncertainty.

#--------------------------------------------------------60
code_created_by = 'Arturo_Avelino'
# On date: 2018.12.10 (yyyy.mm.dd)
code_name = ''
version_code = '0.1.2'
last_update = '2019.01.28'
#--------------------------------------------------------60

import numpy as np
#-------------------------------------------------------

def rms(x_np):
    """
    Function to compute the simple root mean square (rms).
    :param x_np : numpy array of data to determine their rms.
    :type x_np : ndarray.
    """
    N = len(x_np) # number of data
    rmsint = np.sqrt(np.sum(x_np**2.0)/N)
    return rmsint

#-------------------------------------------------------

def err_rms(x_np, sigma_np):
    """
    Function to compute the uncertainty on rms.
    :param x_np : numpy array of data to determine their rms.
    :type x_np : ndarray.
    :param sigma : numpy array of the uncertainties of the data.
    :type sigma : ndarray.
    """
    N = len(x_np) # number of data

    # This definition comes from setting w_s = 1 in:
    # Blondin, Mandel, Kirshner, 2011, "Do spectra improve distance
    #   measurements of Type Ia supernovae"
    # sigmaint = np.sqrt(np.sum(sigma_np**4.0)/(2.0 * N * np.sum(x**2.0)))

    # Using my own determination of err_rms
    sigmaint = np.sqrt(np.sum((x_np*sigma_np)**2.0)/(N * np.sum(x_np**2.0)))

    return sigmaint


#-------------------------------------------------------

def err_rms_boot(x_np, loopsize):
    """
    Function to compute the uncertainty on rms using bootstrap
    :param x_np : numpy array of data to determine their rms.
    :type x_np : ndarray.
    :param loopsize : loop size for bootstrap.
    :type sigma : integer.
    :random.seed = 12345
    """
    ndata = len(x_np) # size of data array
    rms_np = np.zeros(loopsize) # initialize
    # old. random_np = np.zeros(loopsize) # initialize
    np.random.seed(12345)

    #-----------------------------
    for i1 in range(loopsize):

        # Random array with values drawn from the original data
        random_np = np.random.choice(x_np, ndata)

        # Compute the RMS from the random array
        rmsint = np.sqrt(np.sum(random_np**2.0)/ndata)

        # Save the RMS to an array
        rms_np[i1] = rmsint

    #---------

    # Compute the standard deviation of the arrays of RMSs
    stdev_boot = np.std(rms_np)

    return stdev_boot

#-------------------------------------------------------

def wrms(x_np, w_np):
    """
    Function to compute the weighted rms (wrms).
    :param x : numpy array of data to determine their rms.
    :type x : ndarray.
    :param w : numpy array of the weights.
    :type w : ndarray.
    """
    numerator = np.sum(w_np * (x_np**2.0))
    denominator = np.sum(w_np)
    wrms_out = np.sqrt(numerator/denominator)
    return wrms_out

#-------------------------------------------------------

def err_wrms(x_np, w_np, sigma_np):
    """
    Function to compute the uncertainty on  weighted rms (wrms).
    :param x : numpy array of data to determine their rms.
    :type x: ndarray.
    :param w : numpy array of the weights.
    :type w : ndarray.
    :param sigma : numpy array of the uncertainties of the data.
    :type sigma : ndarray.
    """
    # Compute wRMS first:
    numerator = np.sum(w_np * (x_np**2.0))
    denominator = np.sum(w_np)
    wrms_out = np.sqrt(numerator/denominator)

    # Compute the uncertainty
    # This definition comes from Eq. (A.1) of Blondin, Mandel, Kirshner,
    # 2011, "Do spectra improve distance measurements of Type Ia supernovae"

    # varwRMS2 = ( 2.0/(np.sum(w))**2.0 ) * np.sum((w**2.0) * (sigma_np**4.0))
    # sigma_wrms = np.sqrt(varwRMS2)/(2.0*wrms_out)

    # Using my own determination of the expression to compute the uncertainty
    var = 1.0/np.sum(w_np)
    sigma_wrms = np.sqrt(var)

    return sigma_wrms

#-------------------------------------------------------

def err_wrms_boot(x_np, w_np, loopsize):
    """
    Function to compute the uncertainty on wrms using bootstrap
    :param x_np : numpy array of data to determine their wrms.
    :type x_np : ndarray.
    :param w : numpy array of the weights.
    :param loopsize : loop size for bootstrap.
    :type sigma : integer.
    :random.seed = 12345
    """
    ndata = len(x_np) # size of data array
    x_random_np = np.zeros(ndata) # initialize
    w_random_np = np.zeros(ndata) # initialize
    wrms_np = np.zeros(loopsize) # initialize
    np.random.seed(12345)

    #--------------------
    for i1 in range(loopsize):

        # loop over data
        for i2 in range(ndata):

            # Random integer
            randint = np.random.randint(0,ndata)

            x_random_np[i2] = x_np[randint]
            w_random_np[i2] = w_np[randint]

        # Compute the wRMS from the random arrays
        numerator = np.sum(w_random_np * (x_random_np**2.0))
        denominator = np.sum(w_random_np)
        wrms_out = np.sqrt(numerator/denominator)

        # Save the wRMS to an array
        wrms_np[i1] = wrms_out

    #---------

    # Compute the standard deviation of the arrays of RMSs
    stdev_boot = np.std(wrms_np)

    return stdev_boot

#-------------------------------------------------------

#       Test

# test1 = np.array([ 2.,  2.,  2.,  2.,  2.])
# test2 = np.ones(5)

# print rms(test1)

# print err_rms(test1, test2)
# 0.158113883008

