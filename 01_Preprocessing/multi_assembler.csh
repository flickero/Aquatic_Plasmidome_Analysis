#!/bin/tcsh
#$ -o /gpfs0/tals/projects/Analysis/Oded_Project/multi_output_assembly_test
#$ -e /gpfs0/tals/projects/Analysis/Oded_Project/multi_error_assembly_test

# This is main submission script, main_job.sh
# Initialize counter
i=0

# Read each line from the args.csv file, skipping the first line (header)
tail -n +4 "/gpfs0/tals/projects/Analysis/Oded_Project/resources/args.csv" | head -n 8 | while IFS="," read -r study_accession sample_accession experiment_accession run_accession tax_id scientific_name experiment_title fastq_bytes fastq_ftp sra_ftp sample_title; do
    # Extract run_accession and fastq_ftp
    run_accession=$(echo "$run_accession" | tr -d '"')
    fastq_ftp=$(echo "$fastq_ftp" | tr -d '"')

    # Check if the variables are not empty before proceeding
    if [ -n "$run_accession" ] && [ -n "$fastq_ftp" ]; then
        echo "Run Accession: $run_accession"
        echo "Fastq FTP: $fastq_ftp"
        
        # Check the number of jobs in the queue
        r=$(qstat -u flickero | wc -l)
        while [ $r -ge 7 ]; do
            # If more than 11 jobs, wait and check again
            r=$(qstat -u flickero | wc -l)
            sleep 30
        done
    
        # Create a unique job name
        job_name="sample$i"
    
        # Submit the job using qsub
        qsub -N $job_name -cwd -q tals.q -V /gpfs0/tals/projects/Analysis/Oded_Project/scripts/assembly_runner.csh $run_accession $fastq_ftp 
        
        # Print the job name
        echo $job_name
    else
        echo "Error: Unable to extract run_accession or fastq_ftp from the line."
    fi

    # Increment counter
    let i=$i+1
done
