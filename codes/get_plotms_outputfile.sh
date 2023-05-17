#!/bin/sh
if [ -t 0 ]
then
  bin=$(grep 'BIN' $data $1 | awk '{ print $3 }')
fi

work_dir=$(pwd)
for dir in classes/*/*/*; do
  basedir=`basename $dir`
  molname=`dirname $dir | xargs basename`
  if [[ $basedir == "Spectra" ]]; then
    cd $work_dir/$dir
    rm -f $work_dir/$dir/result.csv
    rm -f $work_dir/$dir/result.jdx
    if [ -f "tmpqcxms.res" ]; then
      echo "$bin/plotms -f tmpqcxms.res"
      $bin/plotms -f tmpqcxms.res
    fi
  fi
    cd $work_dir
done


