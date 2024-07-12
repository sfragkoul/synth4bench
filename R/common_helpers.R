

gt_analysis <- function(runs, folder, merged_file) {
  
  nt_runs = list()
  
  for(r in runs) {
    
    a <- paste0(folder, "/", r, "/", r, "_report.tsv") |>
      readLines() |>
      str_split(pattern = "\t", simplify = TRUE) |>
      as.data.frame() |> 
      setDT()
    
    a$V1 = NULL
    # a$V3 = NULL
    a$V5 = NULL
    
    colnames(a) = c("POS", "REF", "DP", paste0("Nt_", 1:(ncol(a) - 3)))
    
    a = melt(
      a, id.vars = c("POS", "REF", "DP"),
      variable.factor = FALSE, value.factor = FALSE,
      variable.name = "Nt", value.name = "Count"
    )
    
    a = a[which(Count != "")]
    
    a$POS = as.numeric(a$POS)
    a$DP = as.numeric(a$DP)
    
    a$Nt = str_split_i(a$Count, "\\:", 1)
    
    a$Count = str_split_i(a$Count, "\\:", 2) |>
      as.numeric()
    
    a$Freq = round(100 * a$Count / a$DP, digits = 6)
    
    a = a[order(POS, -Count)]
    
    a = a[which(REF != a$Nt & Count != 0)]
    
    b = a[which(Nt %in% c("A", "C", "G", "T")), ]
    
    nt_runs[[ as.character(r) ]] = b
  }
  
  nt_runs = rbindlist(nt_runs, idcol = "Run")
  
  pos_of_interest = nt_runs[which(Freq == 100)]$POS |> unique()
  
  gt_runs = nt_runs[which(POS %in% pos_of_interest)]
  
  a <- paste0(folder, "/", merged_file , "_report.tsv") |> 
    readLines() |>
    str_split(pattern = "\t", simplify = TRUE) |> 
    as.data.frame() |> 
    setDT()
  
  a$V1 = NULL
  # a$V3 = NULL
  a$V5 = NULL
  
  colnames(a) = c("POS", "REF", "DP", paste0("Nt_", 1:(ncol(a) - 3)))
  
  a = melt(
    a, id.vars = c("POS", "REF", "DP"),
    variable.factor = FALSE, value.factor = FALSE,
    variable.name = "Nt", value.name = "Count"
  )
  
  a = a[which(Count != "")]
  
  a$POS = as.numeric(a$POS)
  a$DP = as.numeric(a$DP)
  
  a$Nt = str_split_i(a$Count, "\\:", 1)
  
  a$Count = str_split_i(a$Count, "\\:", 2) |>
    as.numeric()
  
  a$Freq = round(100 * a$Count / a$DP, digits = 6)
  
  a = a[order(POS, -Count)]
  
  a = a[which(REF != a$Nt & Count != 0)]
  
  b = a[which(Nt %in% c("A", "C", "G", "T")), ]
  
  
  merged_gt = b[which(POS %in% gt_runs$POS)]
  merged_gt = merged_gt[order(POS)]
  
  merged_gt$Freq = merged_gt$Freq / 100
  
  merged_gt = merged_gt[, by = .(POS, REF, DP), .(
    Nt = paste(Nt, collapse = ","),
    Count = paste(Count, collapse = ","),
    Freq = paste(round(Freq, digits = 3), collapse = ",")
  )]
  
  
  return(merged_gt)
  
}

read_vcf <- function(path, caller, gt, merged_file) {
    
    if(caller == "Freebayes") {
        
        vcf_df <- read_vcf_freebayes(path, gt, merged_file)
        
    } else if (caller == "Mutect2") {
        
        vcf_df <- read_vcf_mutect2(path, gt, merged_file)
        
    } else if (caller == "LoFreq") {
        
        vcf_df <- read_vcf_LoFreq(path, gt, merged_file)
        
    } else if (caller == "VarDict") {
        
        vcf_df <- read_vcf_VarDict(path, gt, merged_file)
        
    } else if (caller == "VarScan") {
        
        vcf_df <- read_vcf_VarScan(path, gt, merged_file)
        
    }
    
    return(vcf_df)
}

plot_synth4bench <- function(gt_comparison, vcf_path, gt_path, caller) {
    
    
    df = fread( gt_comparison )
    
    vcf_GT <- read.vcfR(gt_path, verbose = FALSE )
    
    vcf_caller <- read.vcfR(vcf_path, verbose = FALSE )
    
    if(caller == "freebayes") {
        
        plots <- plot_synth4bench_freebayes(df, vcf_GT, vcf_caller)
        
    } else if (caller == "mutect2") {
        
        plots <- plot_synth4bench_gatk(df, vcf_GT, vcf_caller)
        
    } else if (caller == "LoFreq") {
        
        plots <- plot_synth4bench_LoFreq(df, vcf_GT, vcf_caller)
        
    } else if (caller == "VarDict") {
        
        plots <- plot_synth4bench_VarDict(df, vcf_GT, vcf_caller)
        
    } else if (caller == "VarScan") {
        
        plots <- plot_synth4bench_VarScan(df, vcf_GT, vcf_caller)
        
    }
    
    
    return(plots)
    
}


report_varbp <- function(runs, folder, reference) {
    
    out = list()
    
    for (file in seq_len(runs)){
        
        out[[ file ]] <- list()
        
        message("started file: ", file)
        bam <- paste0(folder, file, "/", file, "_golden.bam")
        
        b <- scanBam(bam)
        b <- b[[1]]
        
        s = b$seq |> 
            as.character() |> 
            str_to_upper()
        
        
        r = read.fasta(reference)
        r = r[[1]] |> str_to_upper() |> paste(collapse = "")
        
        
        message("started main function")
        aln = mapply(function(x, y, z) {
            
            a = cigarRangesAlongReferenceSpace(x, with.ops = TRUE)[[1]] |> as.data.frame() |> setDT()
            b = cigarRangesAlongQuerySpace(x, with.ops = TRUE)[[1]] |> as.data.frame() |> setDT()
            
            b$POS = z - a$start + 1
            
            b$REF = r |> str_sub(
                start = z + a$start - 1,
                end   = z + a$end - 1
            )
            
            b$ALT = y |> str_sub(
                start = b$start,
                end   = b$end
            )
            
            return(b)
            
            
        }, b$cigar, s, b$pos, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        
        message("finished main function")
        
        aln = aln |> rbindlist(idcol = "SequenceID")
        
        r = aln$REF |> str_split("")
        a = aln$ALT |> str_split("")
        
        
        aln$ALT2 = mapply(function(x, y) {
            
            if( length(x) != length(y) ) return(NA)
            
            out = ifelse(x == y, "-", y) |> paste(collapse = "")
            
            return(out)
            
        }, r, a, SIMPLIFY = FALSE, USE.NAMES = FALSE) |> unlist()
        
        out[[ file ]][["report"]] <- aln
        
        # fwrite(
        #     aln, paste0(folder, file, "/", file, "_position_report.csv"),
        #     row.names = FALSE, quote = TRUE, sep = ","
        # )
        
        aln$ALT2 = ifelse(aln$names == "I", aln$ALT, aln$ALT2)
        aln$ALT = NULL
        
        
        aln = aln |> split(aln$names)
        
        aln2 = list()
        
        # alt ---------------------------------------------------------------------
        message("alt")
        
        x = aln$M
        
        ref = x$REF |> str_split("")
        alt = x$ALT2 |> str_split("")
        
        paste0(x$SequenceID, ":", x$start) |> duplicated() |> any()
        
        names(ref) = paste(x$SequenceID, x$POS, x$start, sep = ":")
        names(alt) = paste(x$SequenceID, x$POS, x$start, sep = ":")
        
        y = mapply(function(r, a) {
            
            out = data.table(
                "pos" = r |> length() |> seq_len(),
                "REF" = r, 
                "ALT" = a
            )
            
            out = out[which(ALT != "-")]
            
            return(out)
            
        }, ref, alt, SIMPLIFY = FALSE, USE.NAMES = TRUE)
        
        y = y |> rbindlist(idcol = "id")
        
        y$pos_chromosomal = as.numeric( str_split( y$id, "\\:", simplify = TRUE )[, 2] ) + y$pos - 1
        y$pos_read        = as.numeric( str_split( y$id, "\\:", simplify = TRUE )[, 3] ) + y$pos - 1
        y$id              = str_split( y$id, "\\:", simplify = TRUE )[, 1]
        
        y$pos = NULL
        
        y = y[, by = .(pos_chromosomal, pos_read, REF, ALT), .(
            `No. of sequences` = id |> unique() |> length()
        )]
        
        
        y = y[order(pos_chromosomal, pos_read, -`No. of sequences`)]
        
        aln2[[ "alt" ]] = y
        
        # deletions ---------------------------------------------------------------
        message("deletions")
        
        x = aln$D
        
        
        y = data.table(
            "pos_chromosomal" = x$POS,
            "pos_read"        = x$end,
            "REF"             = x$REF,
            "ALT"             = x$ALT2,
            "id"              = x$SequenceID
        )
        
        y = y[, by = .(pos_chromosomal, pos_read, REF, ALT), .(
            `No. of sequences` = id |> unique() |> length()
        )]
        
        y = y[order(pos_chromosomal, pos_read, -`No. of sequences`)]
        
        aln2[[ "deletion" ]] = y
        
        # insertions --------------------------------------------------------------
        message("insertions")
        
        x = aln$I
        
        
        y = data.table(
            "pos_chromosomal" = x$POS,
            "pos_read"        = x$start,
            "REF"             = x$REF,
            "ALT"             = x$ALT2,
            "id"              = x$SequenceID
        )
        
        y = y[, by = .(pos_chromosomal, pos_read, REF, ALT), .(
            `No. of sequences` = id |> unique() |> length()
        )]
        
        y = y[order(pos_chromosomal, pos_read, -`No. of sequences`)]
        
        aln2[[ "insertion" ]] = y
        
        
        # merge -------------------------------------------------------------------
        message("merge")
        
        aln2 = rbindlist(aln2, idcol = "names")
        
        out[[ file ]][["report2"]] <- aln2
        
        # fwrite(
        #     aln2, paste0(folder, file, "/", file, "_position_report2.csv"),
        #     row.names = FALSE, quote = TRUE, sep = ","
        # )
        
        paste0("finished file: ", file)
    }    
    
    return(out)
    
}


explore_mut_pos <- function(runs, folder, caller) {
    
    found = list()
    not_found = list()
    
    for (file in seq_len(runs)){
        #file= "1"
        message(paste0("working on file ", file))
        
        # input parameters --------------------------------------------------------
        qv = paste0(folder, "/Ground_truth_vs_", caller, ".clean_norm.tsv")
        gb = paste0(folder, file, "/", file, "_position_report2.csv")
        
        # golden alignments -------------------------------------------------------
        gb = gb |> fread()
        gb$key = paste(gb$pos_chromosomal, gb$REF, gb$ALT, sep = ":")
        
        # VarScan missing variants ----------------------
        
        qv = qv |> fread()
        qv$key = paste(qv$POS, qv$`Ground Truth REF`, qv$`Ground Truth ALT`, sep = ":")
        
        
        # filtering common mutations ----------------------------
        
        gb_filtered = gb[which(key %in% qv$key)]
        qv_filtered = qv[which(key %in% gb_filtered$key)]
        
        
        index = qv_filtered[[ paste0(caller, " AF") ]] |>
            is.na() |>
            which()
        
        q1 = qv_filtered[index]
        q2 = gb_filtered[which(key %in% q1$key)]
        
        q2$bin = floor(q2$pos_read / 10) + 1
        q2 = q2[order(q2$bin), ]
        
        q2$Found = "No"
        
        not_found[[file]] = q2
        
        # fwrite(
        #     q2, paste0(folder, file, "/", file, "_VarScan_read_pos_bins_not_found.tsv"),
        #     row.names = FALSE, quote = FALSE, sep = "\t"
        # )
        
        #boxplot-------------------------------------------------------------------
        
        #q2$bin = paste0("bin", q2$bin)
        
        #q2$bin = q2$bin |> factor(levels = unique(q2$bin))
        
        #gr1 = ggplot(data = q2) +
        
        # geom_point(aes(x = `bin`, y = `No. of sequences`, color = `bin`))
        
        #   geom_violin(aes(x = bin, y = `No. of sequences`, fill = bin)) + 
        
        #   ggtitle("Not Found")
        
        
        #ggsave(
        #    plot = gr1, filename = paste0(folder, file, "/", file, "_Box_plots_not_found.jpeg"),
        #    width = 8, height = 8, units = "in", dpi = 600
        #)    
        
        #FOUND---------------------------------------------------------------------
        
        index =  which(!is.na(qv_filtered[[ paste0(caller, " AF") ]]))
        
        q1 = qv_filtered[index]
        q2 = gb_filtered[which(key %in% q1$key)] 
        q2$bin = floor(q2$pos_read / 10) + 1
        
        q2 = q2[order(q2$bin), ]
        
        q2$Found = "Yes"
        
        found[[file]] = q2
        
        # fwrite(
        #     q2, paste0(folder, file, "/", file, "_VarScan_read_pos_bins_found.tsv"),
        #     row.names = FALSE, quote = FALSE, sep = "\t"
        # )
        
        
        
        #q2$bin = paste0("bin", q2$bin)
        
        #q2$bin = q2$bin |> factor(levels = unique(q2$bin))
        
        #gr1 = ggplot(data = q2) +
        
        #   geom_violin(aes(x = bin, y = `No. of sequences`, fill = bin)) + 
        
        #  ggtitle("Found")
        
        
        # ggsave(
        #    plot = gr1, filename = paste0(folder, file, "/", file, "_Box_plots_found.jpeg"),
        #     width = 8, height = 8, units = "in", dpi = 600
        # )  
    }
    
    found = found |> rbindlist(idcol = "Run")
    not_found = not_found |> rbindlist(idcol = "Run")
    
    out = rbind(found, not_found)
    
    return(out)
    
}
