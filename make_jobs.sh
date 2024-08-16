#!/bin/bash

# Nested loops create AA and BB.  For each combination of AA and BB, a job is submitted.
# AA and BB are passed to the job script so that they can be used instead of
# $SLURM_ARRAY_TASK_ID, since we are avoiding use of arrays

for AA in `seq 1 5`; do
   for BB in `seq 1 100`; do
      sbatch -p shared --mem=4G -c 1 --time=3-00:00:00 $CHKPT_HOME/chkpt_job ./non_array_job.sh $AA $BB
   done
done
