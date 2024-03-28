#'
#'A script to report all ground truth variants in each chromosomal position.
#'
#'Input: individual "golden" bam files produced by NEAT.
#'
#'Output: reports in tsv format with the variants in each chromosomal
#'position for each individual "golden" bam file
#'
#'Author: Nikos Pechlivanis(github:npechl) 
#'


rm(list = ls())
gc()

library(data.table)
library(stringr)

library(GenomicAlignments)
library(Rsamtools)
library(seqinr)


files = as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
folder = "1000_50/"

for (file in files){
    
    message("started file: ", file)
    bam <- paste0(folder, file, "/", file, "_golden.bam")
    ref <- "TP53.fasta"
    
    b <- scanBam(bam)
    b <- b[[1]]
    
    s = b$seq |> 
        as.character() |> 
        str_to_upper()
    
    
    r = read.fasta(ref)
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
    
    rm(r, s)
    gc()
    
    r = aln$REF |> str_split("")
    a = aln$ALT |> str_split("")
    
    
    aln$ALT2 = mapply(function(x, y) {
        
        if( length(x) != length(y) ) return(NA)
        
        out = ifelse(x == y, "-", y) |> paste(collapse = "")
        
        return(out)
        
    }, r, a, SIMPLIFY = FALSE, USE.NAMES = FALSE) |> unlist()
    
    
    # aln$REF = NULL
    
    rm(a, b, bam ,ref)
    gc()
    
    fwrite(
        aln, paste0(folder, file, "/", file, "_position_report.csv"),
        row.names = FALSE, quote = TRUE, sep = ","
    )
    
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
    
    fwrite(
        aln2, paste0(folder, file, "/", file, "_position_report2.csv"),
        row.names = FALSE, quote = TRUE, sep = ","
    )
    
    paste0("finished file: ", file)
}    






