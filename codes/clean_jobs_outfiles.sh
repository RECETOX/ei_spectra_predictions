#!/bin/sh


work_dir=$(pwd)
for dir in classes/*/*/*; do
   basedir=`basename $dir`
   molname=`dirname $dir | xargs basename`
   if [[ $basedir == "Spectra" ]]; then
      cd $work_dir/$dir 
      rm *.o*.*
    fi
    cd $work_dir
done

