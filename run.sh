#!/bin/bash
set -e  # stop on error

# Where to save timings
LOGFILE="run_times.txt"
# Header for the log file (overwrites existing file)
echo "Timestamp,Caller,Directory,Script,ElapsedSeconds" > "$LOGFILE"

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

    echo "========================================================"
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Starting caller: $caller"
    echo "               Directory: $dir"
    echo "--------------------------------------------------------"

    # --- Run the main analysis script and time it ---
    start_time=$(date +%s)
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ]   Running: R/S4BR.R"
    Rscript R/S4BR.R -c "$caller" -r 10 -w "$dir" -m Merged
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    now=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[ $now ]   Completed R/S4BR.R in ${elapsed}s"
    # append to log file
    echo "$now,$caller,$dir,S4BR.R,$elapsed" >> "$LOGFILE"
    echo "--------------------------------------------------------"

    # --- Run the plotting script and time it ---
    start_time=$(date +%s)
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ]   Running: R/S4BR_plot.R"
    Rscript R/S4BR_plot.R -c "$caller" -w "$dir" -m Merged
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    now=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[ $now ]   Completed S4BR_plot.R in ${elapsed}s"
    # append to log file
    echo "$now,$caller,$dir,S4BR_plot.R,$elapsed" >> "$LOGFILE"
    echo "--------------------------------------------------------"

    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Finished caller: $caller in $dir"
    echo "========================================================"
    echo

  done
done

echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] All jobs completed."
echo "Timings saved in $LOGFILE"
