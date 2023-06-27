#!/bin/sh

if [ -t 0 ]; then
  walltime=$(grep 'WALLTIME' $data $1 | awk '{ print $3 }')
  user_email=$(grep 'USER_EMAIL' $data $1 | awk '{ print $3 }')
  bin=$(grep 'BIN' $data $1 | awk '{ print $3 }')
  keyword_ntraj=$(grep 'ntraj' $data $1 | awk '{ print $3 }')
fi

work_dir=$(pwd)
rm -rf $work_dir/info_submit_batch_jobs.log
for dir in classes/*/*/*; do
  not_ready=false
  basedir=`basename $dir`
  molname=`dirname $dir | xargs basename`
  if [[ $basedir == "Spectra" ]]; then
    cd $work_dir/$dir
    tmp=${dir#*/}
    class_name=${tmp%%/*}
    file_xyz=`ls $molname.xyz`
    if [ -f "$file_xyz" ]; then
      echo "found" $file_xyz "in class" $class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
      if [[ "$keyword_ntraj" ]]; then
        n_traj=$keyword_ntraj
      else 
        n_traj=$((`head -n 1 $file_xyz`*25))
      fi
        
      if [ -d "TMPQCXMS" ]; then
        echo "found TMPQCXMS in class" $class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
        for i in `seq 1 $n_traj`; do
          cd $work_dir/$dir/TMPQCXMS/TMP.$i

          if [ ! -f "ready" ]; then
            echo "NOT found ready file in TMP."$i "in class" $class_name": molname" $molname
            not_ready=true
            break        
          fi
        done

        if [ "$not_ready" = true ]; then
          echo "removed TMPQCXMS in class" $class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
          rm -rf $work_dir/$dir/TMPQCXMS
        fi
      fi

      if [ ! -d "TMPQCXMS" ]; then
        echo "Not found TMPQCXMS in class" $class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
        echo "submitted neutral MD job: "$class_name": molname" $molname
        echo "submitted neutral MD job: "$class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
        job1=$(qsub -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" $bin/neutral_run_md.pbs)

        echo "submitted prep production MD job: "$class_name": molname" $molname
        echo "submitted prep production MD job: "$class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
        job2=$(qsub -W depend=afterok:$job1 -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" $bin/prep_prod_run_ei_md.pbs)

        echo "submitted production run EI MD job: "$class_name": molname" $molname
        echo "submitted production run EI MD job: "$class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
        job3=$(qsub -W depend=afterok:$job2 -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" -J 1-$n_traj $bin/prod_run_ei_md.pbs)
      fi
    else
      echo "NOT found xyz in class" $class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
    fi
  fi
  cd $work_dir
  echo "===============" >> $work_dir/info_submit_batch_jobs.log
done
