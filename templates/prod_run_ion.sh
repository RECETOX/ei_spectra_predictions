#!/bin/sh

molname=ch2
node=1
ncpus=1
mem=16gb
walltime=02:30:00  
user_email=wrojasv@recetox.muni.cz
bin=/storage/brno2/home/wrojasv/MassSpec/QC_Spectra_Predictions/bin

DIR=TMPQCXMS

if [ ! -d "$DIR" ]; then 
    echo "run prep_ion_md first!"
    exit 
fi

cd $DIR        

for i in TMP.*
    do
        cd $i 
        rm -f ready 

        cat >${molname}_prod.pbs<<EOF
#!/bin/sh

#PBS -l select=$node:ncpus=$ncpus:mem=$mem
#PBS -l walltime=$walltime 
#PBS -N $molname                           
#PBS -j oe
#PBS -m e
#PBS -M $user_email

export OMP_NUM_THREADS=1

cp -rf $(echo '$PBS_O_WORKDIR')/qcxms.in  $(echo '$SCRATCHDIR') 
cp -rf $(echo '$PBS_O_WORKDIR')/start.xyz $(echo '$SCRATCHDIR') 
cp -rf $(echo '$PBS_O_WORKDIR')/qcxms.start $(echo '$SCRATCHDIR') 

cd $(echo '$SCRATCHDIR') 

$bin/qcxms_v.5.2.1 --prod >$molname.out 2>&1

cp -r $(echo '$SCRATCHDIR')/* $(echo '$PBS_O_WORKDIR')/
        
EOF

        chmod +x ${molname}_prod.pbs
        qsub ${molname}_prod.pbs

        cd ..
    done