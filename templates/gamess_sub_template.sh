#!/bin/sh
file_xyz=`basename $(ls *.inp)`

molname=`basename $file_xyz .inp`
walltime={{ WALLTIME }}
bin={{ BIN }}
user_email={{ USER_EMAIL }} 
ncpus={{ NCPUS }} 
mem={{ MEM }} 
scratch_local={{ SCRATCH_LOCAL }} 

echo 'submitted geometry optimiztion job'
job1=$(qsub -N $molname -M $user_email -l walltime=$walltime -l select=1:ncpus=$ncpus:mem=$mem:scratch_local=$scratch_local -v "bin=$bin" $bin/optimization_gamess.pbs)