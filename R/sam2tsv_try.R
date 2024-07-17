library(data.table)

data <- fread("sam2tsv.tsv")
unique(data$`REF-POS1`) |> as.numeric() |> summary()
unique(data$`CIGAR-OP`)

pos_interest <- data[which(`REF-POS1` == 31)]


