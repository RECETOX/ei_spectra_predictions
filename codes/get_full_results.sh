#!/bin/sh

#TODO: This script can still be optimized a bit
if [ -t 0 ]; then
  bin=$(grep 'BIN' $data $1 | awk '{ print $3 }')
fi

work_dir=$(pwd)

if [ ! -d "results" ]; then
  mkdir $work_dir/results
fi
rm -rf $work_dir/info_get_msp.log
rm -f $work_dir/results/simulated_spectra.msp

for dir in classes/*/*/*; do
  basedir=`basename $dir`
  inchikey=`dirname $dir | xargs basename`
  if [[ $basedir == "Spectra" ]]; then
    cd $work_dir/$dir
    tmp=${dir#*/}
    class_name=${tmp%%/*}
    if [ -f "tmpqcxms.res" ]; then
      if [ ! -f "result.jdx" ]; then
        $bin/plotms -f tmpqcxms.res 
      fi
      rm -f result.csv
    
      molname=`sed -n '2{p;q}'  $inchikey.xyz`
      kword=$(grep 'NPOINTS' result.jdx)
      num_peaks=$(sed 's/^[^=]*=//' <<< "$kword")
      echo `pwd`
      sed -n '/PEAK/,/END/{/PEAK/!{/END/!p}}' result.jdx > temp.dat
      awk '{print $1, $2}' temp.dat > tempa.dat
      sed "1s/^/NAME: $molname\nINCHIKEY: $inchikey\nNum Peaks: $num_peaks\n/" tempa.dat >> $work_dir/results/simulated_spectra.msp
      sed -i '$a\ ' $work_dir/results/simulated_spectra.msp
      rm temp.dat tempa.dat 
      echo "Get msp from class :" $class_name" --> molname :" $molname >> $work_dir/info_get_msp.log
    fi
  fi
    cd $work_dir
done

