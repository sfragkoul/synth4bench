#!/bin/bash
set -e  # stop on error

# Define arrays for directories and variant callers
directories=(
  "/mnt/c/Users/sfragkoul/Desktop/synth_data/read_length/1000_50"
  "/mnt/c/Users/sfragkoul/Desktop/synth_data/read_length/1000_75"
  "/mnt/c/Users/sfragkoul/Desktop/synth_data/read_length/1000_100"
  "/mnt/c/Users/sfragkoul/Desktop/synth_data/read_length/1000_300"
  "/mnt/c/Users/sfragkoul/Desktop/synth_data/coverage_test/300_30_10"
  "/mnt/c/Users/sfragkoul/Desktop/synth_data/coverage_test/700_70_10"
  "/mnt/c/Users/sfragkoul/Desktop/synth_data/coverage_test/1000_100_10"
  "/mnt/c/Users/sfragkoul/Desktop/synth_data/coverage_test/3000_300_10"
  "/mnt/c/Users/sfragkoul/Desktop/synth_data/coverage_test/5000_500_10"
)

callers=("Mutect2" "Freebayes" "VarDict" "VarScan" "LoFreq")

# Loop over each directory and each caller
for dir in "${directories[@]}"; do
  for caller in "${callers[@]}"; do

    # print a separator and status line
    echo "========================================================"
    echo "Running caller: $caller"
    echo "  in directory: $dir"
    echo "--------------------------------------------------------"

    # Run the first R script command
    Rscript R/S4BR.R -c "$caller" -r 10 -w "$dir" -m Merged

    # Run the plot R script command
    Rscript R/S4BR_plot.R -c "$caller" -w "$dir" -m Merged

  done
done

echo "All jobs completed." 
