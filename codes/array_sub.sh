#!/bin/sh

if [ -t 0 ]
then
  walltime=$(grep 'WALLTIME' $data $1 | awk '{ print $3 }')
  user_email=$(grep 'USER_EMAIL' $data $1 | awk '{ print $3 }')
  bin=$(grep 'BIN' $data $1 | awk '{ print $3 }')
fi

work_dir=$(pwd)
for dir in classes/*/*/*; do
   basedir=`basename $dir`
   molname=`dirname $dir | xargs basename`
   if [[ $basedir == "Spectra" ]]; then
      cd $work_dir/$dir
      file_xyz=`ls $molname.xyz`
      file_in=`ls qcxms.in`
      if [ -f "$file_xyz" ]; then
        keyword_ntraj=`head -n 1 $file_in | awk '{ print $1 }'` 
        if [[ "$keyword_ntraj" ]]; then
          n_traj=`head -n 1 $file_in | awk '{ print $2 }'`
        else 
          n_traj=$((`head -n 1 $file_xyz`*25))
        fi
        echo 'submitted neutral MD job'
        job1=$(qsub -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" $bin/neutral_run_md.pbs)

        echo 'submitted prep production MD job'
        job2=$(qsub -W depend=afterok:$job1 -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" $bin/prep_prod_run_ei_md.pbs)

        echo 'submitted production run EI MD job'
        job3=$(qsub -W depend=afterok:$job2 -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" -J 1-$n_traj $bin/prod_run_ei_md.pbs)

      fi
    fi
    cd $work_dir
done

