#!/bin/sh

work_dir=$(pwd)
for dir in classes/*/*/*; do
   basedir=`basename $dir`
   molname=`dirname $dir | xargs basename`
   if [[ $basedir == "Optimization" ]]; then
      cd $dir
      if [ ! -f "$molname.log" ]; then
         echo "qsub  $dir/$molname.pbs"
         qsub  $molname.pbs
      elif [[ ! `grep "EQUILIBRIUM GEOMETRY LOCATED" $molname.log` ]]; then 
         echo "resubmit job: $dir/$molname.pbs"
         qsub  $molname.pbs
      fi
   fi
   cd $work_dir
done