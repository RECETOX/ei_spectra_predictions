molname=ch2
node=1
ncpus=1
mem=16gb
walltime=02:30:00  
user_email=wrojasv@recetox.muni.cz
bin=/storage/brno2/home/wrojasv/MassSpec/QC_Spectra_Predictions/bin

cat >${molname}.pbs<<EOF
#!/bin/sh

#PBS -l select=$node:ncpus=$ncpus:mem=$mem
#PBS -l walltime=$walltime 
#PBS -N $molname                           
#PBS -j oe
#PBS -m e
#PBS -M $user_email

export OMP_NUM_THREADS=1

cp -r $(echo '$PBS_O_WORKDIR')/$molname.xyz $(echo '$SCRATCHDIR') 
cd $(echo '$SCRATCHDIR') 

$bin/qcxms_v.5.2.1 -i $molname.xyz  >$molname.out 

cp -r $(echo '$SCRATCHDIR')/* $(echo '$PBS_O_WORKDIR')/

EOF

chmod +x ${molname}.pbs
qsub ${molname}.pbs