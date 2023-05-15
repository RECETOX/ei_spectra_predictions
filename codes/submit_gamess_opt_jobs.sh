#!/bin/sh

work_dir=$(pwd)
for dir in classes/*/*/*; do
   basedir=`basename $dir`
   molname=`dirname $dir | xargs basename`
   if [[ $basedir == "Optimization" ]]; then
      cd $dir
      if [ ! -f "$molname.log" ]; then
         echo "qsub  $molname.pbs"
         qsub  $molname.pbs
      elif ! [[ `grep "EQUILIBRIUM GEOMETRY LOCATED" $molname.log` ]]; then
         if ! [[ `grep "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-" $molname.log` ]]; then 
            echo "resubmit job: $molname"
            qsub  $molname.pbs
         fi
      fi
   fi
   cd $work_dir
done