#!/bin/sh

if [ -t 0 ]; then
  keyword_ntraj=$(grep 'ntraj' $data $1 | awk '{ print $3 }')
fi

work_dir=$(pwd)
rm -rf $work_dir/info_explore_batch_jobs.log
var1=false

COUNTER_xyz=0
COUNTER_mol=0
COUNTER_TOTAL=0
for dir in classes/*/*/*; do
   basedir=`basename $dir`
   molname=`dirname $dir | xargs basename`
   COUNTER_TMP_UC=0
   COUNTER_TMP_C=0
   
   if [[ $basedir == "Spectra" ]]; then
      let COUNTER_TOTAL++
      cd $work_dir/$dir
      rm *.o*.*
      tmp=${dir#*/}
      class_name=${tmp%%/*}
      file_xyz=`ls $molname.xyz`

      if [ -f "$file_xyz" ]; then
        echo "found" $file_xyz "in class" $class_name": molname" $molname >> $work_dir/info_explore_batch_jobs.log
        let COUNTER_xyz++

        if [[ "$keyword_ntraj" ]]; then
          n_traj=$keyword_ntraj
        else 
          n_traj=$((`head -n 1 $file_xyz`*25))
        fi
        
        if [ -d "TMPQCXMS" ]; then
          echo "found TMPQCXMS in class" $class_name": molname" $molname >> $work_dir/info_explore_batch_jobs.log

          if [ -f "tmpqcxms.res" ]; then
            let COUNTER_mol++
            echo "successful SPECTRUM simulated in class" $class_name": molname" $molname >> $work_dir/info_explore_batch_jobs.log
            continue
          fi

          for i in `seq 1 $n_traj`; do
            cd $work_dir/$dir/TMPQCXMS/TMP.$i

            if [ -f "ready" ]; then
              let COUNTER_TMP_C++

            else
              # echo "NOT found "ready" file in TMP."$i >> $work_dir/info_explore_batch_jobs.log
              let COUNTER_TMP_UC++
            fi
          done
          echo "Number of completed TMP :" $COUNTER_TMP_C >> $work_dir/info_explore_batch_jobs.log
          echo "Number of uncompleted TMP :" $COUNTER_TMP_UC >> $work_dir/info_explore_batch_jobs.log
          
        else
          echo "Not found TMPQCXMS in class" $class_name": molname" $molname >> $work_dir/info_explore_batch_jobs.log
        fi

      else
        echo "NOT found xyz in class" $class_name": molname" $molname >> $work_dir/info_explore_batch_jobs.log
      fi
    fi
    cd $work_dir
    echo "==============="  >> $work_dir/info_explore_batch_jobs.log
done
echo "Number of xyz files found:" $COUNTER_xyz "missing:" $((COUNTER_TOTAL - COUNTER_xyz)) >> $work_dir/info_explore_batch_jobs.log
echo "Total number of mols:" $COUNTER_TOTAL >> $work_dir/info_explore_batch_jobs.log
echo "Number of mols ended successfully: " $COUNTER_mol >> $work_dir/info_explore_batch_jobs.log


