
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

#function "not in" def
`%ni%` <- Negate(`%in%`) 

define_fp <- function(caller, gt){
    #FP Variants
    fp_var = caller[which(caller$mut %ni% gt$mut)]
    fp_var$type = "FP"#################################
    
    return(fp_var)
}

define_fn <- function(caller, gt){
    #FN Variants
    fn_var = gt[which(gt$mut %ni% caller$mut)]
    fn_var = fn_var[,c("POS", "REF",  "ALT",  "DP", "AD", "Freq","mut" )]#####
    colnames(fn_var) = c("POS", "REF",  "ALT",  "DP", "AD", "AF","mut" )#####
    fn_var$type = "FN"
    
    return(fn_var)
}

define_tp <- function(caller, gt){
    #FN Variants
    tp_var = caller[which(caller$mut %in% gt$mut)]
    tp_var$type = "TP"
    return(tp_var)
}

