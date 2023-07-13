#!/bin/sh

if [ -t 0 ]; then
  walltime=$(grep 'WALLTIME' $data $1 | awk '{ print $3 }')
  user_email=$(grep 'USER_EMAIL' $data $1 | awk '{ print $3 }')
  bin=$(grep 'BIN' $data $1 | awk '{ print $3 }')
  keyword_ntraj=$(grep 'ntraj' $data $1 | awk '{ print $3 }')
  threshold_traj=$(grep 'threshold_traj' $data $1 | awk '{ print $3 }')
  MEM=$(grep 'MEM' $data $1 | awk '{ print $3 }')
fi

work_dir=$(pwd)
ended_job=false
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
      if [[ "$keyword_ntraj" ]]; then
        n_traj=$keyword_ntraj
      else 
        n_traj=$((`head -n 1 $file_xyz`*25))
      fi
          
      if [ -f "missing_trajs.res" ]; then  
        readarray -t traj_array <  missing_trajs.res
        ended_job=true
        for i in ${traj_array[@]#*.}; do
          prod_job=$(qsub -N $molname -M $user_email -l walltime=$walltime -l select=1:ncpus=1:mem=$MEM:scratch_local=100gb -v "bin=$bin" -J $i-$((i+1)):2 $bin/prod_run_ei_md.pbs)
        done

      elif [ -f "tmpqcxms.res" ]; then
        echo "successful SPECTRUM simulated in class: " $class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
        continue
      else
        ended_job=true
        echo "submitted neutral MD job: "$class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
        neutral_job=$(qsub -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" $bin/neutral_run_md.pbs)

        echo "submitted prep production MD job: "$class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
        prep_job=$(qsub -W depend=afterok:$neutral_job -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" $bin/prep_prod_run_ei_md.pbs)
        
        echo "submitted production run EI MD job: "$class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
        prod_job=$(qsub -W depend=afterok:$prep_job -N $molname -M $user_email -l walltime=$walltime -l select=1:ncpus=1:mem=$MEM:scratch_local=100gb -v "bin=$bin" -J 1-$n_traj $bin/prod_run_ei_md.pbs)
      fi

      if [ "$ended_job"=true ];then
        echo "get tmpqcxms.res: "$class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
        getres_job=$(qsub -W depend=afterok:$prod_job -N $molname -M $user_email -l walltime=$walltime  -v "threshold_traj=$threshold_traj, n_traj=$n_traj" $bin/getres.pbs)
      fi

    else
      echo "NOT found xyz in class" $class_name": molname" $molname >> $work_dir/info_submit_batch_jobs.log
    fi
  fi
  cd $work_dir
  echo "==============="  >> $work_dir/info_submit_batch_jobs.log
done
