#  ADJUSTED SCRIPT FROM FORTRAN TO PYHTON FOR GISHYDRONXT
#  By: JAVER MARDONES
#  Last Modified: 02/13/2020
#
#  ---------------------------------------------------------------------------
#
#  Original description:
#
#  ---------------------------------------------------------------------------
#
#  Things remaining to do:
#  1) need to revise all regional "sig" equations (waiting on Will)
#  2) done - put in vmodel (rmse values from report - need to revise since model works with square of these numbers) (can be done now)
#  3) modify for new rat.dat file. (waiting on Will)
#
# Program to estimate flood frequency in Maryland
#

#
#  This program calculates the fixed region Flood Frequency based on the
#  regression equations developed by Wilbert Thomas as part of the study
#  conducted by Moglen, Thomas, and Cuneo (2006).
#
#  Program originally developed by Gary Tasker in 1999 to accompany the
#  regression equations developed by Dillow (1996).
#
#  Program modified by Glenn Moglen in 2010, and 2016 to accompany the fixed region
#  regression equations developed by Thomas and Moglen (2010, 2015).
#
#  Last Modified: 02/13/2020
#


# Import Libraries

import numpy as np
import os

####

#pythonaddins.MessageBox("TEST", "1")

def MLTPLY(X,Y,K1,K2):

# --------------------------------------------------------------
#  X IS K1 VECTOR
#  Y IS K1*K2 MATRIX
#  PROD = X*Y IS A K1*K2 MATRIX
# --------------------------------------------------------------
    PROD = []
    for k in range(K2):
        var_sum = 0
        for j in range(K1):
            try:
                var_sum = var_sum + X[j]*Y[k][j]
            except:
                var_sum = var_sum + X[j]*Y[j]
        PROD.append(var_sum)
    return PROD


def EYEARS(sig,zp,skew,vpi):

    #  compute equivalent years
    #
    #  Wilson-Hilferty approximation (Handbook of Hydrology, eqn 18.2.29)
    whkp=(2./skew)*(1.+(skew*zp/6.0)-skew**2/36.)**3-2.0/skew
    #  Handbook of hydrology eqn 18.4.11
    hard=1.0+skew*whkp+0.5*(1.+.75*skew**2)*whkp**2
    eqyrs = sig**2*(hard)/vpi
    return eqyrs


def GADJ(directory, area, staid, eqyrs, ip, region):
#
#  subroutine adjusts estimates for nearby gaging station
#  ref. Sauer, 1974, USGS WRI 52-73.
#

    rat_txt = open(os.path.join(directory, "tasker/thomasrat2020.txt"), "r")

    # Splits the element by "\n"
    rat_lines = rat_txt.readlines()
    rat_txt.close()
    rat_lines = [line[:-1] for line in rat_lines]

    gagepr = []
    sid = []
    agage = []
    r = []
    e = []
    for line in rat_lines:
        aux_line = line.split()
        gagepr.append(aux_line[0])
        sid.append(aux_line[1])
        agage.append(aux_line[2])
        r.append(aux_line[3::2])
        e.append(aux_line[4::2])

    radj=1.0
    yadj=0.0
    iflag_aux = 0
    for i in range(len(rat_lines)):
        if sid[i] == str(staid) and gagepr[i].lower() == region.lower():
            iflag_aux = 1
            delta = abs(area-float(agage[i]))
            iflag = 3
            if delta < .5*float(agage[i]):
                radj=float(r[i][ip])-delta*(float(r[i][ip])-1.)/(0.5*float(agage[i]))
                yadj=((float(e[i][ip])-eqyrs)/(float(r[i][ip])-1.0))*(radj-1.0)
            else:
                iflag = 1
        elif iflag_aux == 0:
            iflag=2

    return radj,yadj,iflag


def popvec(directory, region, npred):

    mrd_txt = open(os.path.join(directory, "tasker/mdregdat2020.txt"), "r")

    # Splits the element by "\n"
    mrd_lines = mrd_txt.readlines()
    mrd_txt.close()

    mrd_lines = [line[:-1] for line in mrd_lines]
    mrd_lines = [line.replace(',',' ') for line in mrd_lines]
    mrd_lines = [' '.join(line.split()) for line in mrd_lines]

    mrd1 = mrd_lines[0].split()
    mrd2 = mrd_lines[int(mrd1[1]) + 1].split()
    mrd3 = mrd_lines[int(mrd1[1]) + int(mrd2[1]) + 2].split()
    mrd4 = mrd_lines[int(mrd1[1]) + int(mrd2[1]) + int(mrd3[1]) + 3].split()

    n = 0    # Region A
    if region == 'B':
        n = int(mrd1[1]) + 1
    elif region == 'EC':
        n = int(mrd1[1]) + int(mrd2[1]) + 2
    elif region == 'P':
        n = int(mrd1[1]) + int(mrd2[1]) + int(mrd3[1]) + 3
    elif region == 'WC':
        n = int(mrd1[1]) + int(mrd2[1]) + int(mrd3[1]) + int(mrd4[1]) + 4

    stut = mrd_lines[n + 1].split()
    stut = [float(i) for i in stut[:4]]

    vmodel = mrd_lines[n + 2].split()
    vmodel = [float(i) for i in vmodel[:10]]

    vmodel = [x**2 for x in vmodel]

    bs = []
    xt = []
    for j in range(10):
        bs_read = mrd_lines[n + 3 + j].split()
        bs.append([float(i) for i in bs_read[:npred]])

        xt_aux = []
        for k in range(npred):
            xt_read = mrd_lines[n + 13 + k + npred*j].split()
            xt_aux.append([float(i) for i in xt_read[:npred]])
        xt.append(xt_aux)

    return stut,vmodel,bs,xt


def REGA(directory, optfolder, region, gageid, area, landslope):    # remaining: sigma, remove samax code

    skew = 0.39  #updated 6/16/04, checked on 8/20/10
    samax = [0.006611, 0.003836, 0.003164, 0.00317, 0.004217, 0.005448, 0.007263, 0.009328, 0.012089, 0.016071]
    ak = [-0.84163, -0.43224, 0.0, 0.84162, 1.28155, 1.75069, 2.05375, 2.32635, 2.57583, 2.87816]
    npred = 3
    region_name = 'Appalachian Plateaus and Allegheny Ridges'

    if gageid == 'False':
        gage_name = "No Adjustment"
    else:
        gage_name = gageid
    estim_par = [['Region', 'Area', 'Land Slope', 'Skew', 'Gage ID'], [region_name, round(area,2), round(landslope,3), round(skew,3), gage_name]]

    (stut,vmodel,bs,xt) = popvec(directory, region, npred)

    iflag = 0
    if not gageid == 'False':
        staid = gageid
        iflag = 3

    v = [1.0, np.log10(area), np.log10(landslope)]
    vt = v
    iwarn = 0
    ivpi = 0
    cu = []
    cl = []
    yhat_list = []
    sepc_list = []
    eqyrs_list = []
    sepred_list = []

    for ip in range(10):
        yhat = bs[ip][0] + bs[ip][1]*v[1] + bs[ip][2]*v[2]
        yhat = 10**yhat

        # Compute CI

        xtxi = []
        for i in range(npred):
            xtxi_aux = []
            for j in range(npred):
                xtxi_aux.append(xt[ip][j][i])
            xtxi.append(xtxi_aux)

        temp = MLTPLY(v,xtxi,npred,npred)
        temp2 = MLTPLY(temp,vt,npred,1)
        temp2 = float(temp2[0])

        if temp2 < 0:
            temp2 = 0
            ivpi = ivpi + 1

        vpi = vmodel[ip] + temp2
        sepred = np.sqrt(vpi)

        # compute equivalent years

        # Regional estimate of sigma
        sig = 0.2353
        eqyrs = EYEARS(sig,ak[ip],skew,vpi)

        if iflag == 3:
            (radj,yadj,iflag) = GADJ(directory, area,staid,eqyrs,ip,region)
            yhat = 10**(np.log10(yhat)) * radj

            vpi = vpi*eqyrs/(eqyrs+yadj)
            sepred = np.sqrt(vpi)
            eqyrs = eqyrs+yadj

        sepc = 100*np.sqrt(np.exp(vpi*5.302)-1)

        cu_aux = []
        cl_aux = []
        for iis in range(4):
            t = 10**(stut[iis]*vpi**.5)
            cu_aux.append(yhat*t)
            cl_aux.append(yhat/t)
        cu.append(cu_aux)
        cl.append(cl_aux)

        if temp2 > samax[ip]:
            iwarn = 1

        yhat_list.append(yhat)
        sepc_list.append(sepc)
        eqyrs_list.append(eqyrs)
        sepred_list.append(sepred)

    warning_message = ['']
    if iflag == 3:
        warning_message.append("Estimates adjusted for proximity to station %s" % (gageid))

    if iflag == 2 and not gageid == 'False':
        warning_message.append("Station %s unknown: NO ADJUSTMENT MADE" % (gageid))

    if iflag == 1:
        warning_message.append("Difference in drainage area for Station %s too great: NO ADUSTMENT MADE" % (gageid))

    if iwarn > 0:
        warning_message.append("WARNING -- Prediction beyond observed data")

    if ivpi > 0:
        warning_message.append("Warning: VPI is negative %s" % (ivpi))

    if area < 0.52 or area > 294.14:
        warning_message.append("WARNING - Drainage area out of range of observed data")

    if landslope < 0.066 or landslope > 0.227:
        warning_message.append("WARNING - Land slope out of range of observed data")

    response = {
        'region':region,
        'estim_par': estim_par,
        'warning_message': warning_message,
        'cl': [round(i,0) for i in cl],
        'cu': [round(i,0) for i in cu],
        'yhat_list': [round(i,0) for i in yhat_list],
        'sepc_list': [round(i,1) for i in sepc_list],
        'eqyrs_list': [round(i,2) for i in eqyrs_list],
        'sepred_list': [round(i,4) for i in sepred_list],
    }

    return(response)


def REGP(directory, optfolder, region, gageid, area, li, fc, impa):    # remaining: sigma, remove samax code

    skew = 0.48  #updated 02/13/20
    samax = [0.010764, 0.009292, 0.007926, 0.005814, 0.005044, 0.004724, 0.005105, 0.005785, 0.006569, 0.008877]
    ak = [-0.84163, -0.43224, 0.0, 0.84162, 1.28155, 1.75069, 2.05375, 2.32635, 2.57583, 2.87816]
    npred = 5
    region_name = 'Blue Ridge & Piedmont'

    if gageid == 'False':
        gage_name = "No Adjustment"
    else:
        gage_name = gageid
    estim_par = [['Region', 'Area', 'Lime', 'Forect Cover', 'Impervious Area', 'Skew', 'Gage ID'], [region_name, round(area,2), round(li,2), round(fc,2), round(impa,2), round(skew,3), gage_name]]

    (stut,vmodel,bs,xt) = popvec(directory, region, npred)

    iflag = 0
    if not gageid == 'False':
        staid = gageid
        iflag = 3

    v = [1.0, np.log10(area), np.log10(fc + 1.0), np.log10(impa + 1.0), np.log10(li + 1.0)]
    vt = v
    gv = [1.0, v[1], v[2], v[4]]
    gvt = gv
    iwarn = 0
    ivpi = 0
    cu = []
    cl = []
    yhat_list = []
    sepc_list = []
    eqyrs_list = []
    sepred_list = []

    for ip in range(8):
        yhat = bs[ip][0] + bs[ip][1]*v[1] + bs[ip][2]*v[2] + bs[ip][3]*v[3] + bs[ip][4]*v[4]
        yhat = 10**yhat

        # Compute CI

        xtxi = []
        for i in range(npred):
            xtxi_aux = []
            for j in range(npred):
                xtxi_aux.append(xt[ip][j][i])
            xtxi.append(xtxi_aux)

        temp = MLTPLY(v,xtxi,npred,npred)
        temp2 = MLTPLY(temp,vt,npred,1)
        temp2 = float(temp2[0])

        if temp2 < 0:
            temp2 = 0
            ivpi = ivpi + 1

        vpi = vmodel[ip] + temp2
        sepred = np.sqrt(vpi)

        # compute equivalent years

        # Regional estimate of sigma
        sig = 0.24862 - 0.05379*np.log10(area) + 0.09843*np.log10(fc+1) - 0.0297 *np.log10(impa+1)
        eqyrs = EYEARS(sig,ak[ip],skew,vpi)

        if iflag == 3:
            (radj,yadj,iflag) = GADJ(directory, area,staid,eqyrs,ip,region)
            yhat = 10**(np.log10(yhat)) * radj

            vpi = vpi*eqyrs/(eqyrs+yadj)
            sepred = np.sqrt(vpi)
            eqyrs = eqyrs+yadj

        sepc = 100*np.sqrt(np.exp(vpi*5.302)-1)

        cu_aux = []
        cl_aux = []
        for iis in range(4):
            t = 10**(stut[iis]*vpi**.5)
            cu_aux.append(yhat*t)
            cl_aux.append(yhat/t)
        cu.append(cu_aux)
        cl.append(cl_aux)

        if temp2 > samax[ip]:
            iwarn = 1

        yhat_list.append(yhat)
        sepc_list.append(sepc)
        eqyrs_list.append(eqyrs)
        sepred_list.append(sepred)

    #####################
    # beginning of IP = 9, 10 section
    #####################

    for ip in range(8,10):
        yhat = bs[ip][0] + bs[ip][1]*gv[1] + bs[ip][2]*gv[2] + bs[ip][3]*gv[3]
        yhat = 10**yhat

        # Compute CI

        gxtxi = []
        for i in range(npred-1):
            gxtxi_aux = []
            for j in range(npred-1):
                gxtxi_aux.append(xt[ip][j][i])
            gxtxi.append(gxtxi_aux)

        gtemp = MLTPLY(gv,gxtxi,npred-1,npred-1)
        temp2 = MLTPLY(gtemp,gvt,npred-1,1)
        temp2 = float(temp2[0])

        if temp2 < 0:
            print(ip)
            temp2 = 0
            ivpi = ivpi + 1

        vpi = vmodel[ip] + temp2
        sepred = np.sqrt(vpi)

        # compute equivalent years

        # Regional estimate of sigma
        sig = 0.3070
        eqyrs = EYEARS(sig,ak[ip],skew,vpi)

        if iflag == 3:
            (radj,yadj,iflag) = GADJ(directory, area,staid,eqyrs,ip,region)
            yhat = 10**(np.log10(yhat)) * radj

            vpi = vpi*eqyrs/(eqyrs+yadj)
            sepred = np.sqrt(vpi)
            eqyrs = eqyrs+yadj

        sepc = 100*np.sqrt(np.exp(vpi*5.302)-1)

        cu_aux = []
        cl_aux = []
        for iis in range(4):
            t = 10**(stut[iis]*vpi**.5)
            cu_aux.append(yhat*t)
            cl_aux.append(yhat/t)
        cu.append(cu_aux)
        cl.append(cl_aux)

        if temp2 > samax[ip]:
            iwarn = 1

        yhat_list.append(yhat)
        sepc_list.append(sepc)
        eqyrs_list.append(eqyrs)
        sepred_list.append(sepred)

    warning_message = []
    if iflag == 3:
        warning_message.append("Estimates adjusted for proximity to station %s" % (gageid))

    if iflag == 2 and not gageid == 'False':
        warning_message.append("Station %s unknown: NO ADJUSTMENT MADE" % (gageid))

    if iflag == 1:
        warning_message.append("Difference in drainage area for Station %s too great: NO ADUSTMENT MADE" % (gageid))

    if iwarn > 0:
        warning_message.append("WARNING -- Prediction beyond observed data")

    if ivpi > 0:
        warning_message.append("Warning: VPI is negative %s" % (ivpi))

    if area < 0.111 or area > 816.4:
        warning_message.append("WARNING - Drainage area out of range of observed data")

    if li < 0 or li > 81.7:
        warning_message.append("WARNING - Limestone percentage out of range of observed data")

    if fc < 0 or fc > 100:
        warning_message.append("WARNING - Forest cover out of range of observed data")

    if impa < 0 or impa > 53.3:
        warning_message.append("WARNING - Impervious Area percentage out of range of observed data")

    response = {
        'region':region,
        'estim_par': estim_par,
        'warning_message': warning_message,
        'cl': [[round(i, 0) for i in j] for j in cl],
        'cu': [[round(i, 0) for i in j] for j in cu],
        'yhat_list': [round(i,0) for i in yhat_list],
        'sepc_list': [round(i,1) for i in sepc_list],
        'eqyrs_list': [round(i,2) for i in eqyrs_list],
        'sepred_list': [round(i,4) for i in sepred_list],
    }

    return(response)



def REGWC(directory, optfolder, region, gageid, area, impa, hcd):    # remaining: sigma, remove samax code

    skew = 0.513  #updated 02/13/20
    samax = [0.096783, 0.085013, 0.071458, 0.093167, 0.11419, 0.14229, 0.174649, 0.21733, 0.269874, 0.355002]
    ak = [-0.84163, -0.43224, 0.0, 0.84162, 1.28155, 1.75069, 2.05375, 2.32635, 2.57583, 2.87816]
    npred = 4
    region_name = 'Western Coastal Plain'

    if gageid == 'False':
        gage_name = "No Adjustment"
    else:
        gage_name = gageid

    estim_par = [['Region', 'Area', 'Impervious Area', 'C & D Soils', 'Skew', 'Gage ID'], [region_name, round(area,2), round(impa,2), round(hcd,2), round(skew,3), gage_name]]

    (stut,vmodel,bs,xt) = popvec(directory, region, npred)

    iflag = 0
    if not gageid == 'False':
        staid = gageid
        iflag = 3

    v = [1.0, np.log10(area), np.log10(impa + 1.0), np.log10(hcd + 1.0)]
    vt = v
    iwarn = 0
    ivpi = 0
    cu = []
    cl = []
    yhat_list = []
    sepc_list = []
    eqyrs_list = []
    sepred_list = []

    for ip in range(10):
        yhat = bs[ip][0] + bs[ip][1]*v[1] + bs[ip][2]*v[2] + bs[ip][3]*v[3]
        yhat = 10**yhat

        # Compute CI

        xtxi = []
        for i in range(npred):
            xtxi_aux = []
            for j in range(npred):
                xtxi_aux.append(xt[ip][j][i])
            xtxi.append(xtxi_aux)

        temp = MLTPLY(v,xtxi,npred,npred)
        temp2 = MLTPLY(temp,vt,npred,1)
        temp2 = float(temp2[0])

        if temp2 < 0:
            temp2 = 0
            ivpi = ivpi + 1

        vpi = vmodel[ip] + temp2
        sepred = np.sqrt(vpi)

        # compute equivalent years

        # Regional estimate of sigma
        sig = 0.309
        eqyrs = EYEARS(sig,ak[ip],skew,vpi)

        if iflag == 3:
            (radj,yadj,iflag) = GADJ(directory, area,staid,eqyrs,ip,region)
            yhat = 10**(np.log10(yhat)) * radj

            vpi = vpi*eqyrs/(eqyrs+yadj)
            sepred = np.sqrt(vpi)
            eqyrs = eqyrs+yadj

        sepc = 100*np.sqrt(np.exp(vpi*5.302)-1)

        cu_aux = []
        cl_aux = []
        for iis in range(4):
            t = 10**(stut[iis]*vpi**.5)
            cu_aux.append(yhat*t)
            cl_aux.append(yhat/t)
        cu.append(cu_aux)
        cl.append(cl_aux)

        if temp2 > samax[ip]:
            iwarn = 1

        yhat_list.append(yhat)
        sepc_list.append(sepc)
        eqyrs_list.append(eqyrs)
        sepred_list.append(sepred)

    warning_message = []
    if iflag == 3:
        warning_message.append("Estimates adjusted for proximity to station %s" % (gageid))

    if iflag == 2 and not gageid == 'False':
        warning_message.append("Station %s unknown: NO ADJUSTMENT MADE" % (gageid))

    if iflag == 1:
        warning_message.append("Difference in drainage area for Station %s too great: NO ADUSTMENT MADE" % (gageid))

    if iwarn > 0:
        warning_message.append("WARNING -- Prediction beyond observed data")

    if ivpi > 0:
        warning_message.append("Warning: VPI is negative %s" "\n" % (ivpi))

    if area < 0.41 or area > 349.6:
        warning_message.append("WARNING - Drainage area out of range of observed data")

    if impa < 0 or impa > 36.8:
        warning_message.append("WARNING - Impervious Area percentage out of range of observed data")

    if hcd < 13 or hcd > 74.7:
        warning_message.append("WARNING - C&D-Soils out of range of observed data")

    response = {
        'region':region,
        'estim_par': estim_par,
        'warning_message': warning_message,
        'cl': [round(i,0) for i in cl],
        'cu': [round(i,0) for i in cu],
        'yhat_list': [round(i,0) for i in yhat_list],
        'sepc_list': [round(i,1) for i in sepc_list],
        'eqyrs_list': [round(i,2) for i in eqyrs_list],
        'sepred_list': [round(i,4) for i in sepred_list],
    }

    return(response)


def REGEC(directory, optfolder, region, gageid, area, landslope, ha):    # remaining: sigma, remove samax code

    skew = 0.484  #updated 02/13/20
    samax = [0.120639, 0.119933, 0.123117, 0.140115, 0.152867, 0.173699, 0.192408, 0.215096, 0.24263, 0.28509]
    ak = [-0.84163, -0.43224, 0.0, 0.84162, 1.28155, 1.75069, 2.05375, 2.32635, 2.57583, 2.87816]
    npred = 4
    region_name = 'Eastern Coastal Plain'

    if gageid == 'False':
        gage_name = "No Adjustment"
    else:
        gage_name = gageid

    estim_par = [['Region', 'Area', 'Land Soil', 'A Soil', 'Skew', 'Gage ID'], [region_name, round(area,2), round(landslope,3), round(ha,2), round(skew,3), gage_name]]

    (stut,vmodel,bs,xt) = popvec(directory, region, npred)

    iflag = 0
    if not gageid == 'False':
        staid = gageid
        iflag = 3

    v = [1.0, np.log10(area), np.log10(ha + 1.0), np.log10(landslope)]
    vt = v
    iwarn = 0
    ivpi = 0
    cu = []
    cl = []
    yhat_list = []
    sepc_list = []
    eqyrs_list = []
    sepred_list = []

    for ip in range(10):
        yhat = bs[ip][0] + bs[ip][1]*v[1] + bs[ip][2]*v[2] + bs[ip][3]*v[3]
        yhat = 10**yhat

        # Compute CI

        xtxi = []
        for i in range(npred):
            xtxi_aux = []
            for j in range(npred):
                xtxi_aux.append(xt[ip][j][i])
            xtxi.append(xtxi_aux)

        temp = MLTPLY(v,xtxi,npred,npred)
        temp2 = MLTPLY(temp,vt,npred,1)
        temp2 = float(temp2[0])

        if temp2 < 0:
            temp2 = 0
            ivpi = ivpi + 1

        vpi = vmodel[ip] + temp2
        sepred = np.sqrt(vpi)

        # compute equivalent years

        # Regional estimate of sigma
        sig = 0.295
        eqyrs = EYEARS(sig,ak[ip],skew,vpi)

        if iflag == 3:
            (radj,yadj,iflag) = GADJ(directory, area, staid, eqyrs, ip, region)
            yhat = 10**(np.log10(yhat)) * radj

            vpi = vpi*eqyrs/(eqyrs+yadj)
            sepred = np.sqrt(vpi)
            eqyrs = eqyrs+yadj

        sepc = 100*np.sqrt(np.exp(vpi*5.302)-1)

        cu_aux = []
        cl_aux = []
        for iis in range(4):
            t = 10**(stut[iis]*vpi**.5)
            cu_aux.append(yhat*t)
            cl_aux.append(yhat/t)
        cu.append(cu_aux)
        cl.append(cl_aux)

        if temp2 > samax[ip]:
            iwarn = 1

        yhat_list.append(yhat)
        sepc_list.append(sepc)
        eqyrs_list.append(eqyrs)
        sepred_list.append(sepred)

    warning_message = []
    if iflag == 3:
        warning_message.append("Estimates adjusted for proximity to station %s" % (gageid))

    if iflag == 2 and not gageid == 'False':
        warning_message.append("Station %s unknown: NO ADJUSTMENT MADE" % (staid))

    if iflag == 1:
        warning_message.append("Difference in drainage area for Station %s too great: NO ADUSTMENT MADE" % (gageid))

    if iwarn > 0:
        warning_message.append("WARNING -- Prediction beyond observed data")

    if ivpi > 0:
        warning_message.append("Warning: VPI is negative %s" % (ivpi))

    if area < 0.91 or area > 113.71:
        warning_message.append("WARNING - Drainage area out of range of observed data")

    if landslope < 0.002498 or landslope > 0.016:
        warning_message.append("WARNING - Land Slope out of range of observed data")

    if ha < 0 or ha > 78.8:
        warning_message.append("WARNING - A-Soils out of range of observed data")

    response = {
        'region':region,
        'estim_par': estim_par,
        'warning_message': warning_message,
        'cl': [round(i,0) for i in cl],
        'cu': [round(i,0) for i in cu],
        'yhat_list': [round(i,0) for i in yhat_list],
        'sepc_list': [round(i,1) for i in sepc_list],
        'eqyrs_list': [round(i,2) for i in eqyrs_list],
        'sepred_list': [round(i,4) for i in sepred_list],
    }

    return(response)

def RRE(directory, optfolder, intasker):

    region = intasker[0].upper()
    gageid = intasker[-1]

    if region == 'A':
        area = round(float(intasker[1]), 2)
        landslope = round(float(intasker[2]), 3)
        response = REGA(directory, optfolder, region, gageid, area, landslope)
    elif region == 'B' or region == 'P':
        area = round(float(intasker[1]), 2)
        fc = round(float(intasker[2]), 2)
        li = round(float(intasker[3]), 2)
        impa = round(float(intasker[4]), 2)
        response = REGP(directory, optfolder, region, gageid, area, li, fc, impa)
    elif region == 'WC':
        area = round(float(intasker[1]), 2)
        impa = round(float(intasker[2]), 2)
        hcd = round(float(intasker[3]), 2)
        response = REGWC(directory, optfolder, region, gageid, area, impa, hcd)
    elif region == 'EC':
        area = round(float(intasker[1]), 2)
        landslope = round(float(intasker[2]), 3)
        ha = round(float(intasker[3]), 2)
        response = REGEC(directory, optfolder, region, gageid, area, landslope, ha)

    return(response)
