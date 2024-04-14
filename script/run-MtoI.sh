#!/bin/bash
#
# Huang et al (2022, Mol. Biol. Evol.)

wd=/path/to/working/directory

numsites=(100 250 500 1000 2000 4000 8000 16000 64000)

# CHOOSE klt (infinite #sites) or klx (finite #sites)
# klt_opt=klt
# nr=100

klt_opt=klx
nr=10

# specify which parameters to estimate
#  0 = phi, tauR
#  1 = phi, tauR, tauS
#  2 = phi, thetaS, thetaR, tauR
#  3 = phi, thetaS, thetaR, tauR, tauS
#  4 = phi, thetaS = thetaR, tauR, tauS

for mode in {0..4}; do
  dout=$wd/sc-large-M

  case $mode in
    0)
      dout=$dout/phi-tauR
      ;;
    1)
      dout=$dout/phi-tauR-tauS
      ;;
    2)
      dout=$dout/phi-thetaS-thetaR-tauR
      ;;
    3)
      dout=$dout/all
      ;;
    4)
      dout=$dout/sametheta
      ;;
    *)
      dout=$dout
      ;;
  esac

  dout=$dout/${klt_opt}

  flog=$dout/log.txt

  echo "mode $mode, output $dout"
  mkdir -p $dout

  if [[ ${klt_opt} == "klt" ]]; then
    flag=klt
    echo $flag
    for r in $(seq 1 $nr); do
      fout=$dout/out-$flag-$r.txt
      ./iimmsci2s $mode $fout ${klt_opt} > $flog
    done

  else
    echo -n "klx, n"
    for n in "${numsites[@]}"; do
      flag=klx-n$n
      echo -n " $n"

      drun=$dout/n$n
      mkdir -p $drun
      
      for r in $(seq 1 $nr); do
        ./iimmsci2s $mode $drun/out-$flag-$r.txt ${klt_opt} $n > $flog
      done
    done
    echo ""
  fi

  # remove log file
  rm $flog
done
