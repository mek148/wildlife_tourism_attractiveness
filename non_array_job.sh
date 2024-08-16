#!/bin/bash

# This script takes two arguments, which are used to set $A and $B.

# Load software
module load r/4.1.2

A=$1
B=$2

# Run R script with a command line argument

R CMD BATCH "--args imputed_dataset$A $B" MCMCbird-models.R model$A_$B.Rout

