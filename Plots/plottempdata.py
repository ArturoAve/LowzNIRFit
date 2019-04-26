#!/usr/bin/env python
# coding: utf-8
# Anaconda python 2.7

# Function to plot the template + data for the low-z NIR paper.

# This is copy of my "plottempdata.py" file.

#--------------------------------------------------------60
code_created_by = 'Arturo_Avelino'
# On date: 2019.04.23 (yyyy.mm.dd)
code_name = 'plottempdata.py'
version_code = '0.1.3'
last_update = '2019.04.26'
#--------------------------------------------------------60

from matplotlib import pyplot as plt
import numpy as np

##############################################################################80

#        Plot templates with error OR standard deviation

def plot_templates(xt1, yt1, e_yt1, xt2, yt2, e_yt2, xt3, yt3, e_yt3,
    xt4, yt4, e_yt4,
    xlim=[], ylim=[],
    text_label_1='', text_label_2='', text_label_3='', text_label_4='',
    loc_label=[-5., 2.], color_label = 'k', fs=12,
    color_temp='g', alpha_temp = 0.5,
    color_temp_mean='k', alpha_temp_mean=1.0, ls_temp_mean='-',
    xlabel='x axis', ylabel='y axis', title = 'Templates',
    savefig = True,
    namesavefig = 'plot_templates_.png',
    resolution_dpi = 80,
    dir_save_output = ''   ):

    # Creating arrays with the required format to create the bands plots.
    xt1 = np.atleast_2d(xt1).T
    xt2 = np.atleast_2d(xt2).T
    xt3 = np.atleast_2d(xt3).T
    xt4 = np.atleast_2d(xt4).T

    #--------------------------------------------------------

    plt.clf()
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharex=True,
                                                figsize=(12, 5))
    #----  add a big axes, hide frame ---
    # Useful to put a global title and labels to the axes
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False,
                    right=False)

    # Close the space between the subplots
    plt.subplots_adjust(wspace=0.004, hspace = 0.005)

    #--------------------------------------------------------
    #     Template 1

    # Template: mean value
    ax1.plot(xt1, yt1, color=color_temp_mean, ls=ls_temp_mean, lw=2,
            alpha = alpha_temp_mean, label='plot 1')

    # Template: Standard deviation
    ax1.fill(np.concatenate([xt1, xt1[::-1]]),
            np.concatenate([yt1 - 1. * e_yt1,
                           (yt1 + 1. * e_yt1)[::-1]]),
            alpha=alpha_temp, fc=color_temp, ec='None', label='Standard deviation')

    ax1.text(loc_label[0], loc_label[1], text_label_1, fontsize=fs+2,
            color=color_label)

    ax1.grid(True, ls='--', alpha=0.3)

    if len(xlim) > 0: ax1.set_xlim(xlim[0], xlim[1])
    if len(ylim) > 0: ax1.set_ylim(ylim[0], ylim[1])

    ax1.set_ylabel(ylabel, fontsize=fs)
    ax1.tick_params(labelsize=(fs-2))

    #--------------------------------------------------------
    #     Template 2

    # Template: mean value
    ax2.plot(xt2, yt2, color=color_temp_mean, ls=ls_temp_mean, lw=2,
            alpha = alpha_temp_mean, label='plot 2')

    # Template: Standard deviation
    ax2.fill(np.concatenate([xt2, xt2[::-1]]),
            np.concatenate([yt2 - 1. * e_yt2,
                           (yt2 + 1. * e_yt2)[::-1]]),
            alpha=alpha_temp, fc=color_temp, ec='None', label='Standard deviation')

    ax2.text(loc_label[0], loc_label[1], text_label_2, fontsize=fs+2,
            color=color_label)

    ax2.grid(True, ls='--', alpha=0.3)

    if len(xlim) > 0: ax2.set_xlim(xlim[0], xlim[1])
    if len(ylim) > 0: ax2.set_ylim(ylim[0], ylim[1])

    # Hack to remove the y-axis numbers: make the axis numbers very small and
    # with white color.
    ax2.tick_params(axis='y', labelsize=1)
    for tl in ax2.get_yticklabels():
        tl.set_color('white')

    #-----------------------------------------------------
    #     Template 3

    # Template: mean value
    ax3.plot(xt3, yt3, color=color_temp_mean, ls=ls_temp_mean, lw=2,
            alpha = alpha_temp_mean, label='plot 3')

    # Template: Standard deviation
    ax3.fill(np.concatenate([xt3, xt3[::-1]]),
            np.concatenate([yt3 - 1. * e_yt3,
                           (yt3 + 1. * e_yt3)[::-1]]),
            alpha=alpha_temp, fc=color_temp, ec='None', label='Standard deviation')

    ax3.text(loc_label[0], loc_label[1], text_label_3, fontsize=fs+2,
            color=color_label)

    ax3.grid(True, ls='--', alpha=0.3)

    if len(xlim) > 0: ax3.set_xlim(xlim[0], xlim[1])
    if len(ylim) > 0: ax3.set_ylim(ylim[0], ylim[1])

    # Hack to remove the y-axis numbers: make the axis numbers very small and
    # with white color.
    ax3.tick_params(axis='y', labelsize=1)
    for tl in ax3.get_yticklabels():
        tl.set_color('white')

    #-----------------------------------------------------
    #     Template 4

    # Template: mean value
    ax4.plot(xt4, yt4, color=color_temp_mean, ls=ls_temp_mean, lw=2,
            alpha = alpha_temp_mean, label='plot 4')

    # Template: Standard deviation
    ax4.fill(np.concatenate([xt4, xt4[::-1]]),
            np.concatenate([yt4 - 1. * e_yt4,
                           (yt4 + 1. * e_yt4)[::-1]]),
            alpha=alpha_temp, fc=color_temp, ec='None', label='Standard deviation')

    ax4.text(loc_label[0], loc_label[1], text_label_4, fontsize=fs+2,
            color=color_label)

    ax4.grid(True, ls='--', alpha=0.3)

    if len(xlim) > 0: ax4.set_xlim(xlim[0], xlim[1])
    if len(ylim) > 0: ax4.set_ylim(ylim[0], ylim[1])

    # Hack to remove the y-axis numbers: make the axis numbers very small and
    # with white color.
    ax4.tick_params(axis='y', labelsize=1)
    for tl in ax4.get_yticklabels():
        tl.set_color('white')

    #-----------------------------------------------------

    # ax6.set_xlabel(r'Phase = (MJD - $T_{\rm Bmax}$)/(1+$z_{\rm hel}$)',
    #               fontsize=fs-2)

    plt.title(title, fontsize=fs)
    plt.xlabel(xlabel, fontsize=fs)

    # plt.tight_layout()

    if savefig:
        plt.savefig(dir_save_output+namesavefig, dpi=resolution_dpi)

    plt.close()

    return '# Plot template created.'

##############################################################################80

#                   Plot templates AND data

def templates_data(xt1, yt1, e_yt1, xt2, yt2, e_yt2, xt3, yt3, e_yt3,
    xt4, yt4, e_yt4,
    snname_1_np, appmag_1_np, snname_2_np, appmag_2_np,
    snname_3_np, appmag_3_np, snname_4_np, appmag_4_np,
    xlim=[], ylim=[],
    xlim_temp=[0, -1],

    dir_lc_1 = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/\
10Compute/TheTemplates/Y_band/Std_filters/2_Selection_FlatPrior_ok/\
AllSamples_appMag_vpec_0/Goods/',

    dir_lc_2 = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/\
10Compute/TheTemplates/J_band/Std_filters/2_Selection_FlatPrior_ok/\
AllSamples_appMag_vpec_0/Goods/',

    dir_lc_3 = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/\
10Compute/TheTemplates/H_band/Std_filters/2_Selection_FlatPrior_ok/\
AllSamples_appMag_vpec_0/Goods/',

    dir_lc_4 = '/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/\
10Compute/TheTemplates/K_band/Std_filters/2_Selection_FlatPrior_ok/\
AllSamples_appMag_vpec_0/Goods/',

    text_label_1='', text_label_2='', text_label_3='', text_label_4='',
    loc_label=[-5., 2.], color_label = 'k', fs=12,
    color_temp='g', alpha_temp = 0.5,
    color_temp_mean='k', alpha_temp_mean=0.5, ls_temp_mean='-',
    color_data='k', fmt='.', ms=6, elinewidth=1, capsize=2, alpha_data=0.6,
    xlabel='x axis', ylabel='y axis', title = 'Templates',
    alpha_grid = 0.5,
    savefig = True,
    namesavefig = 'plot_templates_data_.png',
    resolution_dpi = 80,
    dir_save_output = ''   ):

    """
    Function to create a plot with 4 horizontal panels, overlaying a band plot
    and error bar data on top of that in each panel.

    Arguments:
        (xti, yti, e_yti) : numpy arrays of to construct the band (i.e., template)
            plot. These are the (x, y, error_y) values of the band.
        (snname_i_np, appmag_i_np): numpy array of the filenames to plot its
            scatter data (i.e., the photometry) and their respective apparent
            magnitude values.
        (xlim, ylim): x and y axes limits of the plots. These values apply to
            all the panels.
        dir_lc_i: directory where the scatter data files (i.e., the photometry)
            are located.
        text_label_i: text to write on top of each panel to indicate what data is.
        loc_label: location where to print the text specified with "text_label_i".
        color_label; color of the text specified with "text_label_i".
        fs: font size of all the text: axes, title+2, "text_label_i".
        xlim_temp: Rows portion of the datafile to use to plot the band (i.e.,
            template) plots. Default is all the rows, i.e., [0,-1].
    """

    #--------------------------------------------------------
    # Creating arrays with the required format to create the bands plots.
    xt1 = np.atleast_2d(xt1[xlim_temp[0]:xlim_temp[1]] ).T
    xt2 = np.atleast_2d(xt2[xlim_temp[0]:xlim_temp[1]] ).T
    xt3 = np.atleast_2d(xt3[xlim_temp[0]:xlim_temp[1]] ).T
    xt4 = np.atleast_2d(xt4[xlim_temp[0]:xlim_temp[1]] ).T

    yt1 = yt1[xlim_temp[0]:xlim_temp[1]]
    yt2 = yt2[xlim_temp[0]:xlim_temp[1]]
    yt3 = yt3[xlim_temp[0]:xlim_temp[1]]
    yt4 = yt4[xlim_temp[0]:xlim_temp[1]]

    e_yt1 = e_yt1[xlim_temp[0]:xlim_temp[1]]
    e_yt2 = e_yt2[xlim_temp[0]:xlim_temp[1]]
    e_yt3 = e_yt3[xlim_temp[0]:xlim_temp[1]]
    e_yt4 = e_yt4[xlim_temp[0]:xlim_temp[1]]

    #--------------------------------------------------------
    plt.clf()
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharex=True,
                                                figsize=(12, 5))
    #----  add a big axes, hide frame ---
    # Useful to put a global title and labels to the axes
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False,
                    right=False)

    #---------------
    # Close the space between the subplots
    plt.subplots_adjust(wspace=0.004, hspace = 0.005)

    #########################################################
    #     TEMPLATE 1

    # Template: mean value
    ax1.plot(xt1, yt1, color=color_temp_mean,
            ls=ls_temp_mean, lw=2, alpha=alpha_temp_mean, label='plot 1')

    # Template: Standard deviation
    ax1.fill(np.concatenate([xt1, xt1[::-1]]),
            np.concatenate([yt1 - 1. * e_yt1,
                           (yt1 + 1. * e_yt1)[::-1]]),
            alpha=alpha_temp, fc=color_temp, ec='None', label='Standard deviation')

    #------------------------------
    #  Photometry 1

    # Print the LCs for all the SNe Ia.
    for i1 in range(len(snname_1_np)):
        # Define the file name of the SNe.
        file_lc = '%s.txt'%snname_1_np[i1]

        # Read the apparent magnitude at t_Bmax or any other time.
        appmag = appmag_1_np[i1]

        # Upload the data
        data_lc = np.genfromtxt(dir_lc_1+file_lc, skip_header=9,
                           usecols=[0, 3, 4])

        xx = data_lc[:,0]; yy= data_lc[:,1]-appmag; e_yy=data_lc[:,2];

        ax1.errorbar(xx, yy, yerr = e_yy,
                    fmt=fmt, color=color_data, ms=ms, alpha=alpha_data,
                    elinewidth=elinewidth, capsize=capsize)

    #------------------------------

    ax1.text(loc_label[0], loc_label[1], text_label_1, fontsize=fs+2,
            color=color_label)

    ax1.grid(True, ls='--', alpha=alpha_grid)

    if len(xlim) > 0: ax1.set_xlim(xlim[0], xlim[1])
    if len(ylim) > 0: ax1.set_ylim(ylim[0], ylim[1])

    # Size of the numbers in the axes.
    ax1.tick_params(labelsize=(fs-2))

    ax1.set_ylabel(ylabel, fontsize=fs)

    #########################################################
    #     TEMPLATE 2

    # Template: mean value
    ax2.plot(xt2, yt2, color=color_temp_mean, ls=ls_temp_mean, lw=2,
            alpha=alpha_temp_mean, label='plot 2')

    # Template: Standard deviation
    ax2.fill(np.concatenate([xt2, xt2[::-1]]),
            np.concatenate([yt2 - 1. * e_yt2,
                           (yt2 + 1. * e_yt2)[::-1]]),
            alpha=alpha_temp, fc=color_temp, ec='None', label='Standard deviation')

    #------------------------------
    #  Photometry 2

    # Print the LCs for all the SNe Ia.
    for i1 in range(len(snname_2_np)):
        # Define the file name of the SNe.
        file_lc = '%s.txt'%snname_2_np[i1]

        # Read the apparent magnitude at t_Bmax or any other time.
        appmag = appmag_2_np[i1]

        # Upload the data
        data_lc = np.genfromtxt(dir_lc_2+file_lc, skip_header=9,
                           usecols=[0, 3, 4])

        xx = data_lc[:,0]; yy= data_lc[:,1]-appmag; e_yy=data_lc[:,2];

        ax2.errorbar(xx, yy, yerr = e_yy,
                    fmt=fmt, color=color_data, ms=ms, alpha=alpha_data,
                    elinewidth=elinewidth, capsize=capsize)

    #------------------------------

    ax2.text(loc_label[0], loc_label[1], text_label_2, fontsize=fs+2,
            color=color_label)

    ax2.grid(True, ls='--', alpha=alpha_grid)

    if len(xlim) > 0: ax2.set_xlim(xlim[0], xlim[1])
    if len(ylim) > 0: ax2.set_ylim(ylim[0], ylim[1])

    # Size of the numbers in the axes.
    ax2.tick_params(labelsize=(fs-2))

    # Hack to remove the y-axis numbers: make the axis numbers very small and
    # with white color.
    ax2.tick_params(axis='y', labelsize=1)
    for tl in ax2.get_yticklabels():
        tl.set_color('white')

    #########################################################
    #     TEMPLATE 3

    # Template: mean value
    ax3.plot(xt3, yt3, color=color_temp_mean, ls=ls_temp_mean, lw=2,
            alpha=alpha_temp_mean, label='plot 3')

    # Template: Standard deviation
    ax3.fill(np.concatenate([xt3, xt3[::-1]]),
            np.concatenate([yt3 - 1. * e_yt3,
                           (yt3 + 1. * e_yt3)[::-1]]),
            alpha=alpha_temp, fc=color_temp, ec='None', label='Standard deviation')

    #------------------------------
    #  Photometry 3

    # Print the LCs for all the SNe Ia.
    for i1 in range(len(snname_3_np)):
        # Define the file name of the SNe.
        file_lc = '%s.txt'%snname_3_np[i1]

        # Read the apparent magnitude at t_Bmax or any other time.
        appmag = appmag_3_np[i1]

        # Upload the data
        data_lc = np.genfromtxt(dir_lc_3+file_lc, skip_header=9,
                           usecols=[0, 3, 4])

        xx = data_lc[:,0]; yy= data_lc[:,1]-appmag; e_yy=data_lc[:,2];

        ax3.errorbar(xx, yy, yerr = e_yy,
                    fmt=fmt, color=color_data, ms=ms, alpha=(alpha_data-0.1),
                    elinewidth=elinewidth, capsize=capsize)

    #------------------------------

    ax3.text(loc_label[0], loc_label[1], text_label_3, fontsize=fs+2,
            color=color_label)

    ax3.grid(True, ls='--', alpha=alpha_grid)

    if len(xlim) > 0: ax3.set_xlim(xlim[0], xlim[1])
    if len(ylim) > 0: ax3.set_ylim(ylim[0], ylim[1])

    # Size of the numbers in the axes.
    ax3.tick_params(labelsize=(fs-2))

    # Hack to remove the y-axis numbers: make the axis numbers very small and
    # with white color.
    ax3.tick_params(axis='y', labelsize=1)
    for tl in ax3.get_yticklabels():
        tl.set_color('white')

    #########################################################
    #     TEMPLATE 4

    # Template: mean value
    ax4.plot(xt4, yt4, color=color_temp_mean, ls=ls_temp_mean, lw=2,
            alpha=alpha_temp_mean, label='plot 4')

    # Template: Standard deviation
    ax4.fill(np.concatenate([xt4, xt4[::-1]]),
            np.concatenate([yt4 - 1. * e_yt4,
                           (yt4 + 1. * e_yt4)[::-1]]),
            alpha=alpha_temp, fc=color_temp, ec='None', label='Standard deviation')

    #------------------------------
    #  Photometry

    # Print the LCs for all the SNe Ia.
    for i1 in range(len(snname_4_np)):
        # Define the file name of the SNe.
        file_lc = '%s.txt'%snname_4_np[i1]

        # Read the apparent magnitude at t_Bmax or any other time.
        appmag = appmag_4_np[i1]

        # Upload the data
        data_lc = np.genfromtxt(dir_lc_4+file_lc, skip_header=9,
                           usecols=[0, 3, 4])

        xx = data_lc[:,0]; yy= data_lc[:,1]-appmag; e_yy=data_lc[:,2];

        ax4.errorbar(xx, yy, yerr = e_yy,
                    fmt=fmt, color=color_data, ms=ms, alpha=alpha_data,
                    elinewidth=elinewidth, capsize=capsize)

    #------------------------------

    ax4.text(loc_label[0], loc_label[1], text_label_4, fontsize=fs+2,
            color=color_label)

    ax4.grid(True, ls='--', alpha=alpha_grid)

    if len(xlim) > 0: ax4.set_xlim(xlim[0], xlim[1])
    if len(ylim) > 0: ax4.set_ylim(ylim[0], ylim[1])

    # Size of the numbers in the axes.
    ax4.tick_params(labelsize=(fs-2))

    # Hack to remove the y-axis numbers: make the axis numbers very small and
    # with white color.
    ax4.tick_params(axis='y', labelsize=1)
    for tl in ax4.get_yticklabels():
        tl.set_color('white')

    #########################################################

    # ax6.set_xlabel(r'Phase = (MJD - $T_{\rm Bmax}$)/(1+$z_{\rm hel}$)',
    #               fontsize=fs-2)

    plt.title(title, fontsize=fs)
    plt.xlabel(xlabel, fontsize=fs)

    # plt.tight_layout()

    if savefig:
        plt.savefig(dir_save_output+namesavefig, dpi=resolution_dpi)

    plt.close()

    return '# Plot template AND data created.'

##############################################################################80
