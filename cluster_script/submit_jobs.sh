#!/bin/bash

config="../../misc/input.yml"

for i in {1..50} 
do
	job_name="MaSimCycling_${i}"
    output_name="cycling_out/cycling_${i}.txt"


	qsub -N $job_name -v input_files=$config,output_files=$output_name single_job.pbs

done
