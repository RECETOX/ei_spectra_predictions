#!/bin/sh

#TODO: This script can still be optimized a bit
if [ -t 0 ]
then
  bin=$(grep 'BIN' $data $1 | awk '{ print $3 }')
fi

work_dir=$(pwd)

if [ ! -d "results" ]; then
  mkdir $work_dir/results
fi

rm -f $work_dir/results/all_results.msp

for dir in classes/*/*/*; do
  basedir=`basename $dir`
  inchikey=`dirname $dir | xargs basename`
  if [[ $basedir == "Spectra" ]]; then
    cd $work_dir/$dir

    if [ -d "TMPQCXMS" ]; then

      rm -f $work_dir/$dir/tmpqcxms.res
      cd $work_dir/$dir/TMPQCXMS
      for i in TMP.*; do
        cd $work_dir/$dir/TMPQCXMS/$i
        if [ -f "qcxms.res" ]; then
          cat  qcxms.res >> $work_dir/$dir/tmpqcxms.res
        fi
      done

      cd $work_dir/$dir

      if [ ! -f "result.jdx" ]; then
        $bin/plotms -f tmpqcxms.res 
      fi
      rm -f result.csv
    
      molname=`sed -n '2{p;q}'  $inchikey.xyz`
      kword=$(grep 'NPOINTS' result.jdx)
      num_peaks=$(sed 's/^[^=]*=//' <<< "$kword")
      sed -n '/PEAK/,/END/{/PEAK/!{/END/!p}}' result.jdx > temp.dat
      awk '{print $1, $2}' temp.dat > tempa.dat
      sed "1s/^/NAME: $molname\nINCHIKEY: $inchikey\nNum Peaks: $num_peaks\n/" tempa.dat >> $work_dir/results/all_results.msp
      sed -i '$a\ ' $work_dir/results/all_results.msp
      rm temp.dat tempa.dat 

    fi
  fi
    cd $work_dir
done

