
#TP SNVS-----------------------------------------------------------------------
gt_analysis <- function(runs, folder, merged_file) {
  
  nt_runs = list()
  
  for(r in runs) {
    #process reports.tsv files for individual files
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
  
  gt_runs = nt_runs[POS %in% pos_of_interest & Freq == "100"]
  
  #same process reports.tsv files for Merged file
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
  
  
  #merged_gt = b[which(POS %in% gt_runs$POS)]
  merged_gt <- merge(b, gt_runs, by = c("POS", "REF", "Nt"))
  colnames(merged_gt) = c("POS", "REF", "DP", "Nt", "Count", "Freq",
                          "Run", "DP Indiv", "Count Indiv", "Freq Indiv")
  merged_gt = merged_gt[order(POS)]
  
  merged_gt$Freq = merged_gt$Freq / 100
  
  # merged_gt = merged_gt[, by = .(POS, REF, DP), .(
  #   Nt = paste(Nt, collapse = ","),
  #   Count = paste(Count, collapse = ","),
  #   Freq = paste(round(Freq, digits = 3), collapse = ",")
  # )]
  
  
  return(merged_gt)
  
}

read_vcf_snvs_TP <- function(path, caller, gt, merged_file) {
    
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

plot_snvs_TP <- function(gt_snv_tp_comparison, vcf_path, gt_path, caller, merged_file) {
    
    df = fread(paste0(gt_snv_tp_comparison, "/", merged_file, "_", caller, "_snvs_TP.tsv"))
    
    vcf_GT <- read.vcfR(paste0(vcf_path, "/", merged_file, "_ground_truth_norm.vcf"), verbose = FALSE )
    
    vcf_caller <- read.vcfR(paste0(vcf_path, "/", merged_file, "_", caller, "_norm.vcf"), verbose = FALSE )
    
    if(caller == "Freebayes") {

        plots <- plot_snvs_TP_freebayes(df, vcf_GT, vcf_caller, merged_file)

    } else if (caller == "Mutect2") {

        plots <- plot_snvs_TP_gatk(df, vcf_GT, vcf_caller, merged_file)

    } else if (caller == "LoFreq") {

        plots <- plot_snvs_TP_LoFreq(df, vcf_GT, vcf_caller, merged_file)

    } else if (caller == "VarDict") {

        plots <- plot_snvs_TP_VarDict(df, vcf_GT, vcf_caller, merged_file)

    } else if (caller == "VarScan") {

        plots <- plot_snvs_TP_VarScan(df, vcf_GT, vcf_caller, merged_file)

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

#FP & FN SNVS------------------------------------------------------------------
#function "not in" def
`%ni%` <- Negate(`%in%`) 

load_gt_report <- function(path, merged_file) {
    #function to load Ground Truth bam-report 
    a <- paste0(path, "/", merged_file, "_report.tsv") |>
        readLines() |>
        str_split(pattern = "\t", simplify = TRUE) |>
        as.data.frame() |> 
        setDT()
    
    a$V1 = NULL
    a$V5 = NULL
    
    colnames(a) = c("POS", "REF", "DP", paste0("ALT_", 1:(ncol(a) - 3)))
    
    a = melt(
        a, id.vars = c("POS", "REF", "DP"),
        variable.factor = FALSE, value.factor = FALSE,
        variable.name = "ALT", value.name = "Count"
    )
    
    a = a[which(Count != "")]
    
    a$POS = as.numeric(a$POS)
    a$DP = as.numeric(a$DP)
    
    a$ALT = str_split_i(a$Count, "\\:", 1)
    
    a$Count = str_split_i(a$Count, "\\:", 2) |>
        as.numeric()
    
    a$Freq = round(a$Count / a$DP, digits = 6)
    
    a = a[order(POS, -Count)]
    
    a = a[which(REF != a$ALT & Count != 0)]
    
    # select SNVs
    a_snvs = a[which(ALT %in% c("A", "C", "G", "T")), ]
    #filter DEPTH>2
    a_snvs = a_snvs[which(a_snvs$Count >2), ]
    
    
    gt = list(
        all = a,
        snvs = a_snvs
        
    )
    return(gt)
}

load_gt_vcf <- function(path, merged_file){
    #function to load Ground Truth vcf
    ground_truth_vcf <- read.vcfR( paste0(path, "/",merged_file, 
                                          "_ground_truth_norm.vcf"),
                                   verbose = FALSE )
    
    ground_truth_vcf  = ground_truth_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    
    pick_gt = gt_snvs[which(gt_snvs$POS %in% ground_truth_vcf$POS)]
    pick_gt$mut = paste(pick_gt$POS, 
                        pick_gt$REF, 
                        pick_gt$ALT, sep = ":")
    return(pick_gt)
}

select_snvs <- function(df){
    # select SNVs from caller based on length of REF and ALT
    snvs = df[nchar(df$REF) == nchar(df$ALT)]
    snvs = snvs[which(nchar(snvs$REF) <2), ]
    snvs = snvs[which(nchar(snvs$ALT) <2), ]
    snvs$mut = paste(snvs$POS, snvs$REF, snvs$ALT, sep = ":")
    
    return(snvs)
}

define_fp <- function(caller, gt){
    #FP Variants
    fp_var = caller[which(caller$mut %ni% gt$mut)]
    
    return(fp_var)
}

define_fn <- function(caller, gt){
    #FN Variants
    fn_var = gt[which(gt$mut %ni% caller$mut)]
    fn_var$type = "FN"
    
    return(fn_var)
}

define_tp <- function(caller, gt){
    #FN Variants
    tp_var = caller[which(caller$mut %in% gt$mut)]
    tp_var$type = "TP"
    return(tp_var)
}

fn_dp_barplot <- function(q, caller){
    #FP DP plot
    df = q[, c(
        "POS", 
        "Ground Truth DP"
    ), with = FALSE] |>
        unique() |>
        melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
    
    #set color
    if(caller == "Freebayes") {
        
        color <- "#ae8d43"
        
    } else if (caller == "Mutect2") {
        
        color <- "#ae4364"
        
    } else if (caller == "LoFreq") {
        
        color <- "#c974ba"
        
    } else if (caller == "VarDict") {
        
        color <- "#8d43ae"
        
    } else if (caller == "VarScan") {
        
        color <- "#439aae"
    }
    
    
    o3=ggplot(data = df) +
        
        geom_point(aes(x = variable, y = value, fill = variable),
                   position = position_jitternormal(sd_x = .01, sd_y = 0),
                   shape = 21, stroke = .1, size = 2.5) +
        
        geom_boxplot(aes(x = variable, y = value, fill = variable),
                     width = .25, alpha = .5, outlier.shape = NA) +
        
        scale_fill_manual(
            values = c(
                "Ground Truth DP" = color
            )
        ) +
        
        scale_x_discrete(
            labels = c(paste0(caller, " FN Variants"))
        ) +
        
        scale_y_continuous(labels = scales::comma) +
        
        theme_minimal() +
        
        theme(
            legend.position = "none",
            
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", size = 13),
            
            axis.text.x = element_text(face = "bold", size = 13),
            axis.text.y = element_text(face = "bold", size = 13),
            
            axis.line = element_line(),
            axis.ticks = element_line(),
            
            panel.grid = element_blank(),
            
            plot.margin = margin(20, 20, 20, 20)
        ) +
        
        labs(
            y = "Coverage (No. of reads)"
        )
    return(o3)
    
}

fn_af_barplot <- function(q, caller){
    #FP AF plot
    df = q[, c(
        "POS",
        "Ground Truth AF"
    ), with = FALSE] |>
        unique() |>
        
        melt(id.vars = "POS", variable.factor = FALSE, value.factor = FALSE)
    
    #set color
    if(caller == "Freebayes") {
        
        color <- "#ae8d43"
        
    } else if (caller == "Mutect2") {
        
        color <- "#ae4364"
        
    } else if (caller == "LoFreq") {
        
        color <- "#c974ba"
        
    } else if (caller == "VarDict") {
        
        color <- "#8d43ae"
        
    } else if (caller == "VarScan") {
        
        color <- "#439aae"
    }
    
    o4 = ggplot(data = df[which(!is.na(value) & value != 0)]) +
        
        geom_point(aes(x = variable, y = value, fill = variable),
                   position = position_jitternormal(sd_x = .01, sd_y = 0),
                   shape = 21, stroke = .1, size = 2.5) +
        
        geom_boxplot(aes(x = variable, y = value, fill = variable),
                     width = .25, alpha = .5, outlier.shape = NA) +
        
        scale_fill_manual(
            values = c(
                "Ground Truth AF" = color
            )
        ) +
        
        scale_x_discrete(
            labels = c(paste0(caller, " FN Variants"))
        ) +
        
        scale_y_continuous(labels = scales::percent, trans = "log10") +
        
        theme_minimal() +
        
        theme(
            legend.position = "none",
            
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", size = 13),
            axis.text.x = element_text(face = "bold", size = 13),
            axis.text.y = element_text(face = "bold", size = 13),
            
            axis.line = element_line(),
            axis.ticks = element_line(),
            
            panel.grid = element_blank(),
            
            plot.margin = margin(20, 20, 20, 20)
        ) +
        
        labs(
            y = "Allele Frequency"
        )
    return(o4)
    
}

read_vcf_snvs_FP <- function(path, caller, merged_file, pick_gt, gt_all) {
    
    if(caller == "Freebayes") {
        
        fp_var <- final_fp_snvs_Freebayes(path, merged_file, pick_gt, gt_all)
        
    } else if (caller == "Mutect2") { 
        
        fp_var <- final_fp_snvs_gatk(path, merged_file, pick_gt, gt_all)
        
    } else if (caller == "LoFreq") {
        
        fp_var <- final_fp_snvs_LoFreq(path, merged_file, pick_gt, gt_all)
        
    } else if (caller == "VarDict") {
        
        fp_var <- final_fp_snvs_VarDict(path, merged_file, pick_gt, gt_all)
        
    } else if (caller == "VarScan") {
        
        fp_var <- final_fp_snvs_VarScan(path, merged_file, pick_gt, gt_all)
        
    }
    
    return(fp_var)
}

read_vcf_snvs_FN <- function(path, caller, merged_file, pick_gt) {
    
    if(caller == "Freebayes") {
        
        fn_var <- final_fn_snvs_Freebayes(path, merged_file, pick_gt)
        
    } else if (caller == "Mutect2") {
        
        fn_var <- final_fn_snvs_gatk(path, merged_file, pick_gt)
        
    } else if (caller == "LoFreq") {
        
        fn_var <- final_fn_snvs_LoFreq(path, merged_file, pick_gt)
        
    } else if (caller == "VarDict") {
        
        fn_var <- final_fn_snvs_VarDict(path, merged_file, pick_gt)
        
    } else if (caller == "VarScan") {
        
        fn_var <- final_fn_snvs_VarScan(path, merged_file, pick_gt)
        
    }
    
    return(fn_var)
}


plot_snvs_FP <- function(gt_comparison, caller, merged_file) {

    # Construct file path
    file_path <- paste0(gt_comparison, "/", merged_file, "_", caller, "_snvs_FP.tsv")

    # Check if file exists
    if (!file.exists(file_path)) {
        stop(paste("File does not exist:", file_path))
    }

    # Read the file
    df <- tryCatch(
        fread(file_path),
        error = function(e) {
            stop(paste("Error reading file:", file_path, "\n", e$message))
        }
    )

    # Check if the file is empty
    if (nrow(df) == 0) {
        warning(paste("File is empty:", file_path))
        # Return a placeholder plot or NULL
        return(ggplot() + labs(title = paste("No FP snvs data for", caller), x = NULL, y = NULL))
    }

    # Call specific plotting function based on the caller
    if (caller == "Freebayes") {
        fp_plot <- plot_snvs_FP_Freebayes(df, merged_file)
    } else if (caller == "Mutect2") {
        fp_plot <- plot_snvs_FP_gatk(df, merged_file)
    } else if (caller == "LoFreq") {
        fp_plot <- plot_snvs_FP_LoFreq(df, merged_file)
    } else if (caller == "VarDict") {
        fp_plot <- plot_snvs_FP_VarDict(df, merged_file)
    } else if (caller == "VarScan") {
        fp_plot <- plot_snvs_FP_VarScan(df, merged_file)
    } else {
        stop(paste("Unknown caller:", caller))
    }

    return(fp_plot)
}




plot_snvs_FN <- function(gt_comparison, caller, merged_file) {
    
    # Construct file path
    file_path <- paste0(gt_comparison, "/", merged_file, "_", caller, "_snvs_FN.tsv")
    
    # Check if file exists
    if (!file.exists(file_path)) {
        stop(paste("File does not exist:", file_path))
    }
    
    # Read the file
    df <- fread(file_path)
    
    # Check if the file is empty
    if (nrow(df) == 0) {
        warning(paste("File is empty:", file_path))
        # Return a placeholder plot or NULL
        return(ggplot() + labs(title = paste("No FN snvs data for", caller), x = NULL, y = NULL))
    }
    
    # Generate subplots if the file is not empty
    fn_plot1 <- fn_dp_barplot(df, caller)
    fn_plot2 <- fn_af_barplot(df, caller)
    
    # Combine the subplots
    fn_plot <- fn_plot1 + fn_plot2 +
        plot_layout(
            widths = c(1, 1)
        )
    
    return(fn_plot)
}


#INDELs------------------------------------------------------------------------

load_gt_report_indels <- function(path, merged_file) {
    #function to load Ground Truth bam-report 
    a <- paste0(path, "/", merged_file, "_report.tsv") |>
        readLines() |>
        str_split(pattern = "\t", simplify = TRUE) |>
        as.data.frame() |> 
        setDT()
    
    a$V1 = NULL
    a$V5 = NULL
    
    colnames(a) = c("POS", "REF", "DP", paste0("ALT_", 1:(ncol(a) - 3)))
    
    a = melt(
        a, id.vars = c("POS", "REF", "DP"),
        variable.factor = FALSE, value.factor = FALSE,
        variable.name = "ALT", value.name = "Count"
    )
    
    a = a[which(Count != "")]
    a$POS = as.numeric(a$POS)
    a$DP = as.numeric(a$DP)
    a$ALT = str_split_i(a$Count, "\\:", 1)
    a$Count = str_split_i(a$Count, "\\:", 2) |>
        as.numeric()
    a$Freq = round(a$Count / a$DP, digits = 6)
    a = a[order(POS, -Count)]
    a = a[which(REF != a$ALT & Count != 0)]
    
    # select indels
    a_indels = a[which(ALT %ni% c("A", "C", "G", "T")), ]
    #filter DEPTH>2
    a_indels = a_indels[which(a_indels$Count >2), ]
    
    gt = list(
        all = a,
        indels = a_indels
        
    )
    return(gt)
}

select_indels <- function(df){
    #function to select indels from caller based on length of REF and ALT
    
    #identify indels based on length
    indels = df[nchar(df$REF) != nchar(df$ALT)]
    indels$mut = paste(indels$POS, indels$REF, indels$ALT, sep = ":")
    
    return(indels)
}

load_gt_vcf_indels <- function(path, merged_file, gt_indels){
    #function to load Ground Truth vcf
    ground_truth_vcf <- read.vcfR( paste0(path, "/",merged_file, 
                                          "_ground_truth_norm.vcf"),
                                   verbose = FALSE )
    
    ground_truth_vcf  = ground_truth_vcf |> vcfR::getFIX() |> as.data.frame() |> setDT()
    
    pick_gt = gt_indels[which(gt_indels$POS %in% ground_truth_vcf$POS)]
    pick_gt$mut = paste(pick_gt$POS, 
                        pick_gt$REF, 
                        pick_gt$ALT, sep = ":")
    return(pick_gt)
}

standardize_indels <- function(dt) {
    #function to standardize indels
    setDT(dt)
    
    #deletions
    dt[grepl("^-", ALT), `:=` (
        ALT = substring(REF, 1, 1), 
        REF = paste0(REF, substring(ALT, 2)),
        POS = POS - 1  #Adjust POS for deletions
    )]
    
    #insertions
    dt[grepl("^\\+", ALT), ALT := paste0(REF, substring(ALT, 2))]
    
    dt$mut = paste(dt$POS, 
                   dt$REF, 
                   dt$ALT, sep = ":")
    return(dt)
}

gt_stdz_indels <- function(path, merged_file){
    gt_all = load_gt_report_indels(path, merged_file)$all
    gt_indels = load_gt_report_indels(path, merged_file)$indels
    pick_gt = load_gt_vcf_indels(path, merged_file, gt_indels)
    pick_gt_stdz = standardize_indels(pick_gt)
    return(pick_gt_stdz)
} 





call_tp_indels <- function(path, caller, merged_file, pick_gt_stdz) {
  
  if(caller == "Freebayes") {
    
    tp_indels <- final_tp_indels_Freebayes(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "Mutect2") {
    
    tp_indels <- final_tp_indels_gatk(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "LoFreq") {
    
    tp_indels <- final_tp_indels_LoFreq(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "VarDict") {
    
    tp_indels <- final_tp_indels_VarDict(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "VarScan") {
    
    tp_indels <- final_tp_indels_VarScan(path, merged_file, pick_gt_stdz)
    
  }
  
  return(tp_indels)
}



call_fn_indels <- function(path, caller, merged_file, pick_gt_stdz) {
  
  if(caller == "Freebayes") {
    
    fn_indels <- call_fn_indels_Freebayes(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "Mutect2") {
    
    fn_indels <- call_fn_indels_gatk(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "LoFreq") {
    
    fn_indels <- call_fn_indels_LoFreq(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "VarDict") {
    
    fn_indels <- call_fn_indels_VarDict(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "VarScan") {
    
    fn_indels <- call_fn_indels_VarScan(path, merged_file, pick_gt_stdz)
    
  }
  
  return(fn_indels)
}


call_fp_indels <- function(path, caller, merged_file, pick_gt_stdz) {
  
  if(caller == "Freebayes") {
    
    fp_indels <- call_fp_indels_Freebayes(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "Mutect2") {
    
    fp_indels <- call_fp_indels_gatk(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "LoFreq") {
    
    fp_indels <- call_fp_indels_LoFreq(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "VarDict") {
    
    fp_indels <- call_fp_indels_VarDict(path, merged_file, pick_gt_stdz)
    
  } else if (caller == "VarScan") {
    
    fp_indels <- call_fp_indels_VarScan(path, merged_file, pick_gt_stdz)
    
  }
  
  return(fp_indels)
}



plot_indels <- function(path, merged_file, caller) {
    
    if(caller == "Freebayes") {
        
        plots <- circular_plot_Freebayes(path, merged_file, caller)
        
    } else if (caller == "Mutect2") {
        
        plots <- circular_plot_gatk(path, merged_file, caller)
        
    } else if (caller == "LoFreq") {
        
        plots <- circular_plot_LoFreq(path, merged_file, caller)
        
    } else if (caller == "VarDict") {
        
        plots <- circular_plot_VarDict(path, merged_file, caller)
        
    } else if (caller == "VarScan") {
        
        plots <- circular_plot_VarScan(path, merged_file, caller)
        
    }
    
    return(plots)
    
}















