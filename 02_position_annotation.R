



rm(list = ls())
gc()

library(data.table)
library(stringr)

library(AnnotationHub)


df = "Ground_truth_vs_Mutect2.tsv" |> fread()

df2 = df[, c(
    "POS",
    
    "Ground Truth REF",
    "Ground Truth ALT",
    "Ground Truth DP",
    "Ground Truth AF",
    
    "Mutect2 REF",
    "Mutect2 ALT",
    "Mutect2 DP",
    "Mutect2 AF"
), with = FALSE]


df2 = df2[, by = c(
    "POS",
    "Ground Truth REF",
    "Ground Truth ALT",
    "Ground Truth DP",
    "Ground Truth AF"
), .(
    "Mutect2 REF" = `Mutect2 REF` |> tstrsplit(",") |> unlist(),
    "Mutect2 ALT" = `Mutect2 ALT` |> tstrsplit(",") |> unlist(),
    "Mutect2 DP"  = `Mutect2 DP` |> tstrsplit(",") |> unlist() |> as.integer(),
    "Mutect2 AF"  = `Mutect2 AF` |> tstrsplit(",") |> unlist() |> as.numeric()
)]







