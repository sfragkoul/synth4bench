#jvarkit sam2tsv /mnt/c/Users/sfragkoul/Desktop/synth4bench_Copy/results/Merged_auto.sorted.uniq.rg.bam -o  /mnt/c/Users/sfragkoul/Desktop/synth4bench_Copy/results/sam2tsv.tsv -R /mnt/d/sfragkoul/Synth_Data/Synthesizers/NEAT/testing/TP53/TP53.fasta

library(data.table)

data <- fread("sam2tsv.tsv")
unique(data$`REF-POS1`) |> as.numeric() |> summary()
unique(data$`CIGAR-OP`)

pos_interest <- data[which(`REF-POS1` == 31)]


