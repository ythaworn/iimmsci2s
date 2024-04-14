#!/bin/bash
#
# Usage: iimmsci2s mode fout is_klt im_model n

wd=/Users/yuttapong/proj/geneflow/kl/iim

# CHOOSE IM model: 0 = IM, 1 = IIM, 2 = SC
# im_model=0
# im_model=1
im_model=2

# CHOOSE klt (infinite #sites) or klx (finite #sites, n)
# klt_opt=klt

klt_opt=klx
n=100

# set of parameters to estimate
# MtoI:
#  mode = 0: phi, tauR
#  mode = 1: phi, tauR, tauS
#  mode = 2: phi, thetaS, thetaR
#  mode = 3: phi, thetaS, thetaR, tauR, tauS
#  mode = 4: phi, thetaS = thetaR, tauR, tauS
#
# ItoM:
#  mode = 0: w, tauR
#  mode = 1: w, tauR, tauS (iim and sc only)
#  mode = 2: w, thetaT, thetaR
#  mode = 3: w, thetaT, thetaR, tauR, tauS (iim and sc only)
#  mode = 4: w, thetaT = thetaR, tauR, tauS (iim and sc only)
mode=3

fout=out.txt
flog=log.txt

# for klt
# ../iimmsci2s $mode $fout ${klt_opt} ${im_model} > $flog

# for klx
../iimmsci2s $mode $fout ${klt_opt} ${im_model} $n > $flog

# clean up
rm rub
