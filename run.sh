#!/bin/bash

# Name of the final log file; remove it if it already exists
FINAL_LOG="final_results.txt"
rm -f "$FINAL_LOG"

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
    # Log the current run details
    echo "========================================================" >> "$FINAL_LOG"
    echo "Folder: $dir | Caller: $caller" >> "$FINAL_LOG"
    echo "--------------------------------------------------------" >> "$FINAL_LOG"

    # Run the first R script command; redirect both stdout and stderr to the log file
    Rscript R/S4BR.R -c "$caller" -r 10 -w "$dir" -m Merged >> "$FINAL_LOG" 2>&1

    # Run the plot R script command; also redirect output
    Rscript R/S4BR_plot.R -c "$caller" -w "$dir" -m Merged >> "$FINAL_LOG" 2>&1

    # Optional: add a newline to separate runs
    echo "" >> "$FINAL_LOG"
  done
done
