#!/bin/bash


for type in $( cat cancer_types.tsv )
do
    # Define the job script name
    job_script="job_beta_$type.sh"

    # Create the job script
    cat <<EOT > $job_script
#!/bin/bash

#SBATCH --job-name=$type   # Job name
#SBATCH --ntasks=1              # Number of tasks (equivalent to CPU cores)
#SBATCH --cpus-per-task=1       # Number of CPU cores required
#SBATCH --output=output_beta_$type.out  # Output file name
#SBATCH --error=error_beta_$type.err    # Error file name

# Load necessary modules (if required)
module load R/4.3.1-gfbf-2022b

# Run the R script with the parameters
R CMD BATCH --no-save --no-restore \
"--args cancer_type='${type}' chain='beta'" \
s2_cross_validation_single_CPU.R cross_validation_Rout/${type}_beta.Rout

EOT

    # Submit the job using bsub
    sbatch < $job_script

    # Optional: Pause between job submissions (if needed)
    sleep 1
done


