#!/bin/bash
#
# Thawornwattana et al (2024)

wd=/path/to/working/directory

numsites=(100 250 500 1000 2000 4000 8000 16000 64000)

# CHOOSE IM model: 0 = IM, 1 = IIM, 2 = SC
im_model=1

# CHOOSE klt (infinite #sites) or klx (finite #sites)
klt_opt=klt
nr=100

# klt_opt=klx
# nr=50

# specify which parameters to estimate
#  0 = w, tauR
#  1 = w, tauR, tauS
#  2 = w, thetaS, thetaR, tauR
#  3 = w, thetaS, thetaR, tauR, tauS
#  4 = w, thetaS = thetaR, tauR, tauS

for mode in {0..4}; do
  case $im_model in
    0)
      dout0=$wd/msci-im$doutflag-run2
      ;;
    1)
      dout0=$wd/msci-iim$doutflag-run2
      ;;
    2)
      dout0=$wd/msci-sc$doutflag-run2
      ;;
  esac

  case $mode in
    0)
      dout1=$dout0/M-tauR
      ;;
    1)
      dout1=$dout0/M-tauR-tauS
      ;;
    2)
      dout1=$dout0/M-thetaS-thetaR-tauR
      ;;
    3)
      dout1=$dout0/all
      ;;
    4)
      dout1=$dout0/sametheta
      ;;
    *)
      dout1=$dout0
      ;;
  esac

  dout=$dout1/${klt_opt}

  flog=$dout/log.txt

  echo "mode $mode, output $dout"
  mkdir -p $dout


  if [[ ${klt_opt} == "klt" ]]; then
    flag=klt
    echo $flag
    
    for r in $(seq 1 $nr); do
      fout=$dout/out-$flag-$r.txt
      ./iimmsci2s $mode $fout ${klt_opt} ${im_model} > $flog
    done

  else
    echo -n "klx: n ="
    for n in "${numsites[@]}"; do
      flag=klx-n$n
      echo -n " $n"

      drun=$dout/n$n
      mkdir -p $drun
      
      for r in $(seq 1 $nr); do
        ./iimmsci2s $mode $drun/out-$flag-$r.txt ${klt_opt} ${im_model} $n > $flog
      done
    done
    echo ""
  fi

  # remove log file
  rm $flog
done
