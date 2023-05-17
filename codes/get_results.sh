#!/bin/sh
date 
work_dir=$(pwd)
for dir in classes/*/*/*; do
  basedir=`basename $dir`
  molname=`dirname $dir | xargs basename`
  if [[ $basedir == "Spectra" ]]; then
    cd $work_dir/$dir
    if [ -d "TMPQCXMS" ]; then
      cd $work_dir/$dir/TMPQCXMS
      rm -f $work_dir/$dir/tmpqcxms.res
      for i in TMP.*; do
        cd $work_dir/$dir/TMPQCXMS/$i
        if [ -f "ready" ]; then
          if [ -f "qcxms.res" ]; then
            cat  qcxms.res >> $work_dir/$dir/tmpqcxms.res
          fi
        fi
      done
    fi
  fi
    cd $work_dir
done
date 