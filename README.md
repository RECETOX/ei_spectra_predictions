# ei_spectra_predictions

## Spectra generation steps

### Handling input files for molecular optimization in the HPC

```
 ~/codes/write_dirs_files_for_gamess_qcxms.sh   
```

It submits the following command:
```
python codes/write_dirs_files_for_gamess_qcxms.py --sdf_filename data/RECETOX_GC-EI-MS_20201028.sdf   --params_filename codes/all_parameters.in --project_dirname /storage/du-cesnet/home/wrojasv/VO_metacentrum-tape_tape-archive/QC_SPECTRA_EI_RECETOX
```

### Submitting jobs for molecular optimization
```
 ~/codes/submit_gamess_opt_jobs.sh   ~/codes/all_parameters.in
```

### Extracting optimized molecules for spectra simulation
```
python codes/extract_coords_into_inp_xyz.py --sdf_filename data/RECETOX_GC-EI-MS_20201028.sdf   --params_filename codes/all_parameters.in --project_dirname /storage/du-cesnet/home/wrojasv/VO_metacentrum-tape_tape-archive/QC_SPECTRA_EI_RECETOX
```

### Handling batch jobs for spectra simulations
```
 ~/codes/submit_qcxms_batch_jobs.sh   ~/codes/all_parameters.in
```

### Collecting simulated spectra in MSP format
```
 ~/codes/get_full_results.sh   ~/codes/all_parameters.in
```

### Additional tools
Print resume of completed molecular spectra simulation and delete unnecessary files
```
 ~/codes/explore_qcxms_batch_jobs.sh   ~/codes/all_parameters.in
```