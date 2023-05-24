#!/bin/sh
if [ -t 0 ]
then
  bin=$(grep 'BIN' $data $1 | awk '{ print $3 }')
fi

work_dir=$(pwd)

if [ ! -d "results" ]; then
  mkdir $work_dir/results
fi

for dir in classes/*/*/*; do
  basedir=`basename $dir`
  molname=`dirname $dir | xargs basename`
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
        echo "$bin/plotms -f tmpqcxms.res"
        $bin/plotms -f tmpqcxms.res 
      fi
      rm -f $work_dir/$dir/result.csv
    
      if [ -f "result.jdx" ]; then
        kword=$(grep 'NPOINTS' result.jdx)
        num_peaks=$(sed 's/^[^=]*=//' <<< "$kword")
        sed -n '/PEAK/,/END/{/PEAK/!{/END/!p}}' result.jdx > temp.dat
        awk '{print $1, $2}' temp.dat > tempa.dat
        sed '1s/^/Name: '$molname'\nNum Peaks: '$num_peaks'\n/' tempa.dat > $work_dir/results/$(sed -e 's/.[^.]*$//' <<< "$molname.jdx")'.msp'
        rm temp.dat tempa.dat
      fi

    fi
  fi
    cd $work_dir
done

if [ ! -f "all_results.msp" ]; then
  cat $work_dir/results/*.msp >> $work_dir/results/all_results.msp
fi
