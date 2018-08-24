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
        # Create the variable "snName" containing the first 8 (or 7) letters of the SNe file name
        # I use "snName" to compare with the SNe names in 'SNeWithCepheidDistances.txt', so that
        # I will not compute a peculiar velocity uncertainty for those SNe.
        try:
            if   name[7] == '_': snName = name[:7] # To read correctly, e.g., "sn2011B_"
            elif name[7] != '_':
                if is_number(name[7]): snName = name[:15] # To read correctly, e.g., "snf20080514-002"
                else: snName = name[:8]  # To read correctly, e.g., "sn1998bu"
        except: snName = name[:6]  # To read correctly, e.g., "sn2011B"

        sigma_pecInt = 0 # Reset its value:
        # If this SNe has Cepheid distance, then use vpecFix=0 for it:
        if snName in ListSNeCepheid['f0']: sigma_pecInt = np.sqrt(sigma2_pec(zcmbInt, err_zcmb, 0))
        else: sigma_pecInt = np.sqrt(sigma2_pec(zcmbInt, err_zcmb, vpecFix))

        #---------------------------
        # If there are at least 3 photometric points in the app mag LC text file, then compute.
        if len(magData) >= InfoSN_NumLinesSkip + MinNumOfDataInLC:

            #---------------------------
            # Find out if there are at least one datum inside the interpolated phase range of the template
            NumDataInRange = 0
            for i in range(InfoSN_NumLinesSkip, len(magData)): # Loop over all the phases for a given SNe
                if magData[i,0] >= lowerPhase and magData[i,0] <= upperPhase:
                    NumDataInRange = NumDataInRange + 1

            #---------------------------
            # Fit the LCs with at least one data point within the phase range defined for the template.
            if NumDataInRange > 0:

                # Analytical minimization of chi2 to find the best estimate for distance modulus
                # and computing the uncertainty of distance modulus (analytic expression)

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

                for i in range(InfoSN_NumLinesSkip, len(magData)): # Loop over all the photometry for a given SNe
                    if magData[i,0] >= lowerPhase: # Ignore data outside the phase interpolation range:
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

                # Compute the residual and absolute magnitude value

                mu_LCDM = DistanceMu(magData[0,0], OmMFix, wFix, HoFix)
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
                # This data will be used to write down the "DistanceMu_All_BeforeCutoffs_.txt" text file

                textfile.write('%s%-30s   \
%.10f  %.10f  %.10f  %.10f  %13.10f  %12.10f        \
%1.f          %.10f  %.10f  \
%.10f  %.10f  \
%.10f  %.10f  %.6f  %.6f      0.00   \
%.10f  %.8f  %.10f  %.10f  \
%13.10f  %.10f  %10.7f  %.7f  %.10f  %.10f  %.10f  \
%.10f  %.10f    \
%.6f  %.6f  %.10f  %.10f # \n'%
                               (commentText, list_SNe[j][0][:-4],
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
                # This data will be used to write down the "DistanceMu_Good_AfterCutoffs_Main_.txt" text file

                if (zcmbInt < zcmbUpperLim and dm15Int > dm15LowerLim and dm15Int < dm15UpperLim and
                    EBVhost > EBVhostMin and EBVhost < EBVhostMax and EBV_MW < EBVMWLim and
                    chi2_dof_Int < chi2_dof_Max and mu_resid < residualMax):
                    textfileMain.write('%s%-30s   \
%.10f  %.10f  %.10f  %.10f  %13.10f  %12.10f        \
%1.f          %.10f  %.10f  \
%.10f  %.10f  \
%.10f  %.10f  %.6f  %.6f      0.00   \
%.10f  %.8f  %.10f  %.10f  \
%13.10f  %.10f  %10.7f  %.7f  %.10f  %.10f  %.10f  \
%.10f  %.10f    \
%.6f  %.6f  %.10f  %.10f # \n'%
                               (commentText, list_SNe[j][0][:-4],
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
                    # This data will be used to write down the "DistanceMu_Good_AfterCutoffs_z001_.txt"

                    if zcmbInt > 0.01:
                        textfil3.write('%s%-30s   \
%.10f  %.10f  %.10f  %.10f  %13.10f  %12.10f        \
%1.f          %.10f  %.10f  \
%.10f  %.10f  \
%.10f  %.10f  %.6f  %.6f      0.00   \
%.10f  %.8f  %.10f  %.10f  \
%13.10f  %.10f  %10.7f  %.7f  %.10f  %.10f  %.10f  \
%.10f  %.10f    \
%.6f  %.6f  %.10f  %.10f # \n'%
                               (commentText, list_SNe[j][0][:-4],
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
                    textfil5.write('%s%-30s   \
%.10f  %.10f  %.10f  %.10f  %13.10f  %12.10f        \
%1.f          %.10f  %.10f  \
%.10f  %.10f  \
%.10f  %.10f  %.6f  %.6f      0.00   \
%.10f  %.8f  %.10f  %.10f  \
%13.10f  %.10f  %10.7f  %.7f  %.10f  %.10f  %.10f  \
%.10f  %.10f    \
%.6f  %.6f  %.10f  %.10f # \n'%
                               (commentText, list_SNe[j][0][:-4],
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
                    plt.text(-7, y_lowerlim-0.7, 'Best fit: $m_{T_{Bmax}}$ = %r $\pm$ %r'%(round(appMag_TBmax,3),
                                                                              round(np.sqrt(error_appMag_TBmax),3)))

                plt.text(-7, y_lowerlim-0.5, '$\mu$ = %r $\pm$ %r'%(round(mu_photo_Analytic,3),
                                                                              round(np.sqrt(sigma2_mu),3)))
                plt.text(-7, y_lowerlim-0.3, '$z_{CMB}$ = %r'%round(magData[0,0],5))
                plt.text(-7, y_lowerlim-0.1, '$\Delta \mu$ = %r'%round(mu_resid,3))
                if Chi2dofPrint == True:
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
    textfil4.write('%6.2f       %5.2f      %.4f  %3.f     %.4f      %3.f  \n'%(lowerPhase, upperPhase,
                                                      rms, countGoodSNe, rms_z001, countGoodSNe_z001))

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

    print ' '
    print '--- All done smoothly, DistanceModuli_*.txt created ---'
    print 'Number of SNe before cutoffs = ', countSNe
    print 'Number of SNe after cutoffs = ', countGoodSNe
    print 'Number of SNe after cutoffs AND 0.01 < z < 0.06 = ', countGoodSNe_z001
    print text7, text_08_01, text_22