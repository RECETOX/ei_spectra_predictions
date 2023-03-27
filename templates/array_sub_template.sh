#!/bin/sh
file_xyz=`basename $(ls *.xyz)`
file_in=`basename $(ls *.in)`

keyword_ntraj=`head -n 1 $file_in | awk '{ print $1 }'` 

if [[ "$keyword_ntraj" ]]; then
    n_traj=`head -n 1 $file_in | awk '{ print $2 }'`
else 
    n_traj=$((`head -n 1 $file_xyz`*25))
fi

molname=`basename $file_xyz .xyz`
walltime={{ WALLTIME }}
bin={{ BIN }}
user_email={{ USER_EMAIL }} 

echo 'submitted neutral MD job'
job1=$(qsub -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" $bin/neutral_run_md.pbs)

echo 'submitted prep production MD job'
job2=$(qsub -W depend=afterok:$job1 -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" $bin/prep_prod_run_ei_md.pbs)

echo 'submitted production run EI MD job'
job3=$(qsub -W depend=afterok:$job2 -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" -J 1-$n_traj $bin/prod_run_ei_md.pbs)