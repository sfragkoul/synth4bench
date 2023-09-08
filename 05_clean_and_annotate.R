#'
#'This R script takes the tsv file and adds annotation information based on 
#'the gene chromosomal positions.
#'
#'
#'Input: tsv file containing information regarding the ground truth variants
#'
#'Output: annotated tsv file containing information regarding the ground 
#'truth variants and annotation information
#'


rm(list = ls())
gc()

library(data.table)
library(stringr)
library(GenomicRanges)

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
    "Ground Truth DP",
    "Mutect2 REF",
    "Mutect2 ALT",
    "Mutect2 DP",
    "Mutect2 AF"
    
), .(
    "Ground Truth ALT" = `Ground Truth ALT` |> tstrsplit(",") |> unlist(),
    "Ground Truth AF"  = `Ground Truth AF` |> tstrsplit(",") |> unlist()
    # "Mutect2 REF" = `Mutect2 REF` |> tstrsplit(",") |> unlist(),
    # "Mutect2 ALT" = `Mutect2 ALT` |> tstrsplit(",") |> unlist(),
    # "Mutect2 DP"  = `Mutect2 DP` |> tstrsplit(",") |> unlist() |> as.integer(),
    # "Mutect2 AF"  = `Mutect2 AF` |> tstrsplit(",") |> unlist() |> as.numeric()
)]



mutect2_alt = str_split(df2$`Mutect2 ALT`, ",")
mutect2_af = str_split(df2$`Mutect2 AF`, ",")


cln = mapply(
    function(x, y, z) {
        
        index = which(y == x)
        
        return(
            c(y[index], z[index])
        )
        
    },
    
    df2$`Ground Truth ALT`, mutect2_alt, mutect2_af
)


df2$`Mutect2 ALT` = cln |> lapply(function(x) { return(x [1]) }) |> unlist()
df2$`Mutect2 AF` = cln |> lapply(function(x) { return(x [2]) }) |> unlist()

df2[which(is.na(`Mutect2 AF`))]$`Mutect2 DP` = NA
df2[which(is.na(`Mutect2 AF`))]$`Mutect2 REF` = NA

df2 = df2[, c(
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


# add position annotation ------------------------

mutation_gr <- GRanges(
    seqnames = Rle(rep("chr17", nrow(df2)), rep(1, nrow(df2))),
    ranges   = IRanges(start = df2$POS + 7668402, end = df2$POS + 7668402)
)

library(systemPipeR)
library(AnnotationHub)

ah <- AnnotationHub()

ah <- query(ah, c("hg38", "TxDb", "knownGene"))

db <- ah[["AH52260"]]

columns(db)

ann = genFeatures(db, reduce_ranges = FALSE)

overlaps = findOverlaps(mutation_gr, ann$exon)

df2$exonInfo = character()
df2[from(overlaps)]$exonInfo = "exon"

overlaps = findOverlaps(mutation_gr, ann$intron)

df2$intronInfo = character()
df2[from(overlaps)]$intronInfo = "intron"

rm(
    ah, ann, cln, df, mutation_gr, mutect2_af, mutect2_alt, overlaps, db
)

gc()

fwrite(
    df2, "Ground_truth_vs_Mutect2.clean.annotated.tsv",
    sep = "\t", row.names = FALSE, quote = FALSE
)

# 
# 
# common_v2 = melt(
#     df2[, c(
#         "POS", 
#         "Ground Truth DP", "Ground Truth AF", 
#         "Mutect2 DP", "Mutect2 AF"
#     ), with = FALSE], 
#     
#     id.vars = c("POS"),
#     value.factor = FALSE, variable.factor = FALSE,
#     value.name = "AF", variable.name = "type"
# )
# 
# 
# common_v2 = melt(
#     common[, c("POS", "Ground Truth AF", "Mutect2 AF"), with = FALSE], 
#     id.vars = c("POS"),
#     value.factor = FALSE, variable.factor = FALSE,
#     value.name = "AF", variable.name = "type"
# )
# 
# common_v2$variable = common_v2$type |> str_sub(-2, -1)
# 
# library(ggplot2)
# 
# ggplot(data = common_v2, aes(x = type, y = AF)) +
#     
#     geom_boxplot() +
#     
#     facet_wrap(vars(variable), scales = "free")
# 
# library(rstatix)
# 
# 
# wilcox_test(data = common_v2, AF ~ type, paired = TRUE)
# 
# #Fisher's F-test
# var.test(common$`Ground Truth AF`, common$`Mutect2 AF`)
# 
# 
# # library(patchwork)
# # 
# # gr1 = ggplot(data = df2, aes(x = POS, y = `Ground Truth ALT`)) +
# #     
# #     geom_point(aes(size = `Ground Truth AF`, color = `Ground Truth AF`)) +
# #     
# #     scale_size_continuous(limits = c(0, .7), labels = scales::percent) +
# #     scale_color_continuous(limits = c(0, .7)) +
# #     
# #     theme(
# #         legend.position = "bottom"
# #     )
# # 
# # gr2 = ggplot(data = df2[which(!is.na(`Mutect2 ALT`))], aes(x = POS, y = `Mutect2 ALT`)) +
# #     
# #     geom_point(aes(size = `Mutect2 AF`, color = `Mutect2 AF`)) +
# #     
# #     scale_size_continuous(limits = c(0, .7), labels = scales::percent) +
# #     scale_color_continuous(limits = c(0, .7)) +
# #     
# #     theme(
# #         legend.position = "bottom"
# #     )
# # 
# # 
# # gr1 / gr2
# 
# 
# 
# df2 = df[, c(
#     "POS",
#     
#     "Ground Truth REF",
#     "Ground Truth ALT",
#     "Ground Truth DP",
#     "Ground Truth AF",
#     
#     "Mutect2 REF",
#     "Mutect2 ALT",
#     "Mutect2 DP",
#     "Mutect2 AF"
# ), with = FALSE]
# 
# 
# df2 = df2[, by = c(
#     "POS",
#     "Ground Truth REF",
#     "Ground Truth ALT",
#     "Ground Truth DP",
#     "Ground Truth AF"
# ), .(
#     "Mutect2 REF" = `Mutect2 REF` |> tstrsplit(",") |> unlist(),
#     "Mutect2 ALT" = `Mutect2 ALT` |> tstrsplit(",") |> unlist(),
#     "Mutect2 DP"  = `Mutect2 DP` |> tstrsplit(",") |> unlist() |> as.integer(),
#     "Mutect2 AF"  = `Mutect2 AF` |> tstrsplit(",") |> unlist() |> as.numeric()
# )]
# 
# 
# 
# df2$`Ground Truth AF` = df2$`Ground Truth AF` / 100
# 
# df2 = df2[, c(
#     "POS",
#     "Ground Truth ALT",
#     "Ground Truth AF",
#     "Mutect2 ALT",
#     "Mutect2 AF"
# ), with = FALSE]
# 
# 
# df2 = df2[which(str_length(df2$`Mutect2 ALT`) == 1)]
# 
# df2 = list(
#     "Ground Truth" = df2[, 1:3],
#     "Mutect2"      = df2[, c(1, 4:5)]
# )
# 
# df2 = lapply(df2, function(x) { colnames(x) = c("POS", "ALT", "AF"); return(x)})
# 
# df2 = rbindlist(df2, idcol = "variable")
# 
# 
# library(ggforce)
# 
# ggplot() +
#     
#     geom_point(data = df2[which(variable == "Ground Truth")], 
#                aes(x = POS, y = ALT, size = AF, fill = variable), 
#                shape = 21, position = position_nudge(y = .15)) +
#     
#     geom_point(data = df2[which(variable == "Mutect2")], 
#                aes(x = POS, y = ALT, size = AF, fill = variable), 
#                shape = 21, position = position_nudge(y = -.15)) +
#     
#     scale_size_continuous(range = c(1, 10), breaks = c(.1, .3, .5, .6),
#                           labels = scales::percent)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
