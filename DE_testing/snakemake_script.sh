# # snakemake job submission script that submits 100 jobs at a time
# if this is going really slow, can use '-P sage' to send jobs to sage instead
snakemake --keep-going -j 350 --cluster "qsub -N {params.job_name} -o {params.error_out_file}.out -e {params.error_out_file}.error -l h_rt={params.run_time} -pe serial {params.cores} -l h_vmem={params.memory}G"