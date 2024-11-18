Rscript R/S4BR.R -c Mutect2 -r 10 -v /mnt/c/Users/sfragkoul/Desktop/synth4bench/results -w /mnt/c/Users/sfragkoul/Desktop/synth4bench/results -m Merged
Rscript R/S4BR.R -c Freebayes -r 10 -v /mnt/c/Users/sfragkoul/Desktop/synth4bench/results -w /mnt/c/Users/sfragkoul/Desktop/synth4bench/results -m Merged
Rscript R/S4BR.R -c VarScan -r 10 -v /mnt/c/Users/sfragkoul/Desktop/synth4bench/results -w /mnt/c/Users/sfragkoul/Desktop/synth4bench/results -m Merged
Rscript R/S4BR.R -c VarDict -r 10 -v /mnt/c/Users/sfragkoul/Desktop/synth4bench/results -w /mnt/c/Users/sfragkoul/Desktop/synth4bench/results -m Merged
Rscript R/S4BR.R -c LoFreq -r 10 -v /mnt/c/Users/sfragkoul/Desktop/synth4bench/results -w /mnt/c/Users/sfragkoul/Desktop/synth4bench/results -m Merged


Rscript R/S4BR_plot.R -c Mutect2 -w /mnt/c/users/sfragkoul/Desktop/synth4bench/results -t /mnt/c/users/sfragkoul/Desktop/synth4bench/results -v /mnt/c/users/sfragkoul/Desktop/synth4bench/results -g /mnt/c/users/sfragkoul/Desktop/synth4bench/results -m Merged
Rscript R/S4BR_plot.R -c Freebayes -w /mnt/c/users/sfragkoul/Desktop/synth4bench/results -t /mnt/c/users/sfragkoul/Desktop/synth4bench/results -v /mnt/c/users/sfragkoul/Desktop/synth4bench/results -g /mnt/c/users/sfragkoul/Desktop/synth4bench/results -m Merged
Rscript R/S4BR_plot.R -c VarScan -w /mnt/c/users/sfragkoul/Desktop/synth4bench/results -t /mnt/c/users/sfragkoul/Desktop/synth4bench/results -v /mnt/c/users/sfragkoul/Desktop/synth4bench/results -g /mnt/c/users/sfragkoul/Desktop/synth4bench/results -m Merged
Rscript R/S4BR_plot.R -c VarDict -w /mnt/c/users/sfragkoul/Desktop/synth4bench/results -t /mnt/c/users/sfragkoul/Desktop/synth4bench/results -v /mnt/c/users/sfragkoul/Desktop/synth4bench/results -g /mnt/c/users/sfragkoul/Desktop/synth4bench/results -m Merged
Rscript R/S4BR_plot.R -c LoFreq -w /mnt/c/users/sfragkoul/Desktop/synth4bench/results -t /mnt/c/users/sfragkoul/Desktop/synth4bench/results -v /mnt/c/users/sfragkoul/Desktop/synth4bench/results -g /mnt/c/users/sfragkoul/Desktop/synth4bench/results -m Merged