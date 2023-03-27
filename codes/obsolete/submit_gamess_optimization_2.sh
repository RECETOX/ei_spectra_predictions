#!/bin/sh
molname=ch2
node=1
ncpus=4
mem=16gb
walltime=02:30:00  
user_email=wrojasv@recetox.muni.cz

cat >${molname}.pbs<<EOF
#!/bin/sh

#PBS -l select=$node:ncpus=$ncpus:mem=$mem
#PBS -l walltime=$walltime 
#PBS -N $molname                           
#PBS -j oe
#PBS -m e
#PBS -M $user_email

export OMP_NUM_THREADS=1

module load gamess-sept2018

cp -r $(echo '$PBS_O_WORKDIR')/$molname.inp $(echo '$SCRATCHDIR') 
cd $(echo '$SCRATCHDIR')  

rungms $molname.inp 00 $(echo '$((PBS_NCPUS/2))') $(echo '$PBS_NCPUS') >$molname.log >&1 $molname.err

cp -r $(echo '$SCRATCHDIR')/$molname.log $(echo '$PBS_O_WORKDIR')/
EOF

chmod +x ${molname}.pbs
qsub ${molname}.pbs