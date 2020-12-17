#' Generate model input 
#'
#' This function prepares the input files for downstream analysis including annotating variants according to genes within which they fall. An RData object is created so that this function only needs to be run once.
#' @param ASE_file Specify file containing allele-specific expression information
#' @param legend_file Specify legend file containing SNP IDs and position info etc
#' @param haplotypes_file Provide haplotypes corresponding to SNPs in legend file, for given individuals found in samples file
#' @param samples_file containing cohort information
#' @param output_prefix Pre-fix to the name of the output file
#' @param species Select same species as cohort. Defaults to 'hsapiens' 
#' @param ensembl_version Specify version of Ensembl to download from. Defaults to NULL which is the most recent build 
#' @export
#' @examples
#' #' #Annotating the input files using the most recent human gene set 
#' Gen.input('path_to_ASE_file/ASEfile.txt.gz',
#'                     'path_to_legend_file/file.legendfile.gz',
#'                    'path_to_haplotypes_file/hapfile.hap.gz',
#'                     'path_to_samples_file/samples.txt')
#'                     
#' #Annotating the input files using a previous human gene set
#' Gen.input('path_to_ASE_file/ASEfile.txt.gz',
#'                     'path_to_legend_file/file.legendfile.gz',
#'                    'path_to_haplotypes_file/hapfile.hap.gz',
#'                     'path_to_samples_file/samples.txt', ensembl_version=78)
#'
#' #Annotating the input files using the most recent mouse gene set
#' Gen.input('path_to_ASE_file/ASEfile.txt.gz',
#'                     'path_to_legend_file/file.legendfile.gz',
#'                   'path_to_haplotypes_file/hapfile.hap.gz',
#'                     'path_to_samples_file/samples.txt',
#'                     'mmusculus')



# 1. What you put into the command line
Gen.input <- function(ASE_file, legend_file, haplotypes_file, samples_file, output_prefix, species = "hsapiens", ensembl_version = NULL) {
    # 2. For timing length of script
    Script.start.time <- Sys.time()
    # 3. Set your arguments. I will put a command that shows people what they have their arguments set to
    cat("\n\nARGUMENTS SET:\n")
    cat(paste(c("              ASE file:        ", ASE_file, "\n"), collapse = ""))
    cat(paste(c("              legend file:     ", legend_file, "\n"), collapse = ""))
    cat(paste(c("              haplotypes file: ", haplotypes_file, "\n"), collapse = ""))
    cat(paste(c("              Samples file:    ", samples_file, "\n"), collapse = ""))
    # cat(paste(c(' Number of tasks: ', NumTasks, '\n'), collapse=''))
    cat(paste(c("              species:         ", species, "\n"), collapse = ""))
    if (!is.null(ensembl_version)) {
        cat(paste(c("              ensembl version:  ", ensembl_version, "\n"), collapse = ""))
    } else {
        cat(paste(c("              ensembl version: ", "latest build", "\n"), collapse = ""))
    }
    
    dir.create("PackageTestWork", showWarnings = FALSE)
    dir.create("PackageTestWork/RDataFiles", showWarnings = FALSE)
    dir.create("PackageTestWork/ModelResults", showWarnings = FALSE)
    dir.create("PackageTestWork/ProgressFiles", showWarnings = FALSE)
    
    # 6. Load input files
    ASE <- readInputs(ASE_file, "allele counts")
    Samples <- readInputs(samples_file, "samples")
    LEG <- readInputs(legend_file, "legend")
    hap <- readInputs(haplotypes_file, "haplotypes")
    # ADD code to remove sites in ASE site not recorded as hets in hap file?
    
    if (nrow(LEG) != nrow(hap)) {
        stop("Haplotype and legend files do not contain the same number of rows")
    }
    if ((ncol(hap)/2) != nrow(Samples)) {
        stop("Must be two haplotypes per sample individual i.e. twice as many columns in the haplotype file as individuals in the samples file.")
    }
    dups <- duplicated(LEG$ID)
    if (sum(dups) > 0) {
        cat("Found", sum(dups), "duplicated IDs in the legend file. Removing these from the analysis.\n")
        LEG <- LEG[!dups, ]
        hap <- hap[!dups, ]
    }
    if (length(grep("chr", ASE$chr)) > 0) {
        cat("Removing leading chr from chromosome names\n")
        ASE$chr <- gsub("chr", "", ASE$chr)
    }
    colnames(hap) <- rep(Samples$Ind, each = 2)
    rownames(hap) <- LEG$ID
    
    colnames(LEG)[c(2:4)] <- c("end", "ref", "alt")
    ASE_vars <- merge(ASE, LEG, by = c("end", "ref", "alt"), all.x = FALSE, all.y = FALSE)
    ASE_vars <- merge(ASE_vars, Samples, by = "Ind", all.x = FALSE, all.y = FALSE)
    # ASE_vars<-merge(ASE_vars, counts, by.x='Ind', by.y='ID', all.x=FALSE, all.y=FALSE) some SNPs appear twice in ASE input when
    # overlapping two genes. may want to change this if doing gene level analysis
    ASE_vars <- ASE_vars[!duplicated(ASE_vars[, c("Ind", "ID")]), ]
    
    
    ensembl <- ensemblAttributes(species, ensembl_version)
    
    
    
    # To correct for overall level of expression at given haplotype for given individual
    refSums <- cbind(aggregate(refCount ~ Ind, data = ASE, FUN = sum), 0)
    altSums <- cbind(aggregate(altCount ~ Ind, data = ASE, FUN = sum), 1)
    colnames(refSums)[2:3] <- c("Reads", "All")
    colnames(altSums)[2:3] <- c("Reads", "All")
    counts <- rbind(refSums, altSums)
    
    # MAY NEED TO CHANGE COLUMN NAMES. COULD JUST DO NUMBERS AND SAY WHAT EACH NUMBER CORRESPONDS TO IN A HELP FILE
    ensembl$TSS <- ifelse(ensembl$strand == 1, ensembl$transcript_start, ensembl$transcript_end)
    ensemblmerge <- GRanges(seqnames = ensembl$chromosome_name, ranges = IRanges(start = ensembl$exon_chrom_start, end = ensembl$exon_chrom_end, 
        names = ensembl$ensembl_transcript_id), geneid = ensembl$ensembl_gene_id, TSS = ensembl$transcript_start, strand = ensembl$strand)
    
    ASE_varsmerge <- GRanges(seqnames = Rle(ASE_vars$chr), ranges = IRanges(start = ASE_vars$end, end = ASE_vars$end), names = ASE_vars$ID)
    
    merge <- mergeByOverlaps(ASE_varsmerge, ensemblmerge)
    # MAY NEED TO CHANGE COLUMN NAMES. COULD JUST DO NUMBERS AND SAY WHAT EACH NUMBER CORRESPONDS TO IN A HELP FILE
    
    
    # find out which snps are not found in my merge and so are without a known tss Pull them out of the metric file which also has their
    # positions
    snps_without_tss <- ASE_vars[which(!(ASE_vars$ID %in% merge$names)), ]
    
    # Add the TSS and geneid information to the ASE_vars df
    inGenes <- data.frame(ID = merge$names, TSS = merge$TSS, Gene = merge$geneid)
    outsideGenes <- data.frame(ID = snps_without_tss$ID, TSS = snps_without_tss$end, Gene = snps_without_tss$ID)
    inGenes <- inGenes[!duplicated(inGenes), ]
    outsideGenes <- outsideGenes[!duplicated(outsideGenes), ]
    annoDF <- rbind.data.frame(inGenes, outsideGenes)
    # This merge takes a while
    cat("Now merging\n")
    ASEGenes <- left_join(ASE_vars, annoDF, by = "ID")
    ################# removing redundancy here but possible this is going to cause problems if doing gene level analyses
    ASEGenes <- ASEGenes[!duplicated(ASEGenes[, c("Ind", "ID")]), ]
    ASEGenes$propRef <- ASEGenes$refCount/(ASEGenes$refCount + ASEGenes$altCount)
    ASEGenes$totalReads <- rowSums(ASEGenes[, c("refCount", "altCount")])
    ASEGenes$logRatio <- log2((ASEGenes$refCount + 1)/(ASEGenes$altCount + 1))
    #in theory user can input site with no reads so remove these
    ASEGenes<-ASEGenes[which(!is.na(ASEGenes$propRef)),]
    
    expectedRatio <- 0.5
    if (length(ASEGenes$propRef > 1000)) {
        expectedRatio <- median(ASEGenes$propRef)
    }
    cat("Setting expected proportion of reference reads to", expectedRatio, "in binomial test\n")
    ASEGenes$binomp <- mapply(applyBinom, ASEGenes$refCount, ASEGenes$refCount + ASEGenes$altCount, expectedRatio)
    cat("Merge complete\n")
    
    aseDat <- list(haps = hap, ASE = ASEGenes, leg = LEG, counts = counts, prefix = output_prefix, covar = Samples)
    output_file = paste(c(output_prefix, ".RData"), collapse = "")
    cat(paste(c("Saving ", output_file, "\n"), collapse = ""))
    save(aseDat, file = output_file)
    cat("Finished\n")
    return(aseDat)
}


applyBinom <- function(x, n, p) {
    binom.test(x, n, p)$p.value
}


readInputs <- function(thisFile, type) {
    if (file.exists(thisFile)) {
        if (type == "haplotypes") {
            thisData <- read.table(thisFile, header = F)
        } else {
            thisData <- read.table(thisFile, stringsAsFactors = F, header = T)
        }
        cols <- colnames(thisData)
        # lets change samples file so can take any number of different covariates do we need the ref and alt alleles in ASE file? Best to keep
        # required columns to a minimum also lets call allele counts rather than ASE
        if ((type == "samples") & (cols[1] != "Ind")) {
            print(sample_example)
            stop("First line of samples file does not start with Ind. Make sure a header line is specified as above")
            
        } else if ((type == "allele counts") & !all(colnames(ASE_example) %in% cols)) {
            print(ASE_example)
            stop("input file does not contain the correct column names. Allele counts file must be in the format shown above.")
            
        } else if ((type == "legend") & !all(colnames(LEG_example) %in% cols)) {
            
            print(LEG_example)
            stop("input file does not contain the correct column names. Legend file must be in the format shown above.")
            
        } else if (type == "haplotypes" & (length(which(thisData != 0 & thisData != 1)) > 0)) {
            stop("Error: Haplotype file contains values other than 0 or 1.")
        } else {
            cat(type, "file successfully loaded\n")
        }
    } else {
        stop(paste("Couldnt find specified ", type, " file at: ", getwd(), "/", thisFile, sep = ""))
    }
    return(thisData)
}

