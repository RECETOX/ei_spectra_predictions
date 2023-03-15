#!/bin/sh
#molname=`basename $(pwd)`
#n_atoms=$(head -n 1 $file)

file=`basename $(ls *.xyz)`

molname=`basename $file .xyz`
walltime=02:00:00
bin=/storage/brno2/home/wrojasv/MassSpec/QC_Spectra_Predictions/local_run_test/bin_test
user_email=wrojasv@recetox.muni.cz

n_traj=`find TMPQCXMS/ -mindepth 1 -maxdepth 1 -type d | wc -l`

echo 'submitted neutral MD job'
job1=$(qsub -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" $bin/neutral_run_md.pbs)

echo 'submitted prep production MD job'
job2=$(qsub -W depend=afterok:$job1 -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" $bin/prep_prod_run_ei_md.pbs)

echo 'submitted production run EI MD job'
job3=$(qsub -W depend=afterok:$job2 -N $molname -M $user_email -l walltime=$walltime -v "bin=$bin" -J 1-$n_traj $bin/prod_run_ei_md.pbs)