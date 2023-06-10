#!/bin/sh

#PBS -l select=1:ncpus=1:mem=16gb:scratch_local=100gb
#PBS -l walltime=40:00:00
#PBS -M wudmir.rojas@recetox.muni.cz
#PBS -j oe
#PBS -m e

file=RECETOX_GC-EI-MS_20201028 
out_dirname=QC_SPECTRA_EI_RECETOX

bin=/storage/brno2/home/wrojasv/MassSpec/QC_Spectra_Predictions/source_code_git/ei_spectra_predictions
env_path=/auto/brno2/home/wrojasv/micromamba/envs/ei-spectra-predictions

module load micromamba/1.1.0 $env_path

micromamba activate 

echo 'starting job'
date 
cp -r $bin/templates  $SCRATCHDIR
cd $SCRATCHDIR

$env_path/bin/python $bin/codes/write_dirs_files_for_gamess_qcxms.py --sdf_filename $bin/data/$file.sdf   --params_filename $bin/codes/all_parameters.in  --project_dirname $out_dirname

cp -rf $SCRATCHDIR/$out_dirname $PBS_O_WORKDIR

cd $PBS_O_WORKDIR
echo 'finished job' 
date 
