ensemblConnect <- function(spec = "hsapiens", vers = NULL) {
    if (!is.null(vers)) {
        ensembl = useEnsembl(biomart = "ensembl", dataset = paste(spec, "_gene_ensembl", sep = ""), version = vers, host = "uswest.ensembl.org")
    } else {
        ensembl = useEnsembl(biomart = "ensembl", dataset = paste(spec, "_gene_ensembl", sep = ""), host = "uswest.ensembl.org")
    }
}

ensemblAttributes <- function(spec = "hsapiens", vers = NULL) {
    ens <- NULL
    ensAttr <- NULL
    attempt <- 0
    maxAttempts <- 20
    speciesID <- paste(spec, "_gene_ensembl", sep = "")
    while (is.null(ens) && attempt < maxAttempts) {
        attempt <- attempt + 1
        cat(paste("Trying to connect to Ensembl. Attempt", attempt, sep = " "), "\n")
        tryCatch({
            # Obtain dataframe with ensembl information
            ens <- ensemblConnect(spec, vers)
        }, error = function(err) {
            print("Connection failed.")
        }, finally = {
            if (!is.null(ens)) {
                cat("Connection successful. Getting gene info from Ensembl.\n")
                cat("This can take a few minutes.\n")
                ensAttr <- getBM(attributes = c("chromosome_name", "ensembl_transcript_id", "ensembl_gene_id", "start_position", "exon_chrom_start", 
                  "exon_chrom_end", "transcript_start", "transcript_end", "strand"), mart = ens)
                cat("Completed retrieving gene info.\n")
            } else if (attempt == maxAttempts) {
                stop("Connection unsuccessful. Giving up. Did you specify a valid species name? Options include hsapiens, mmusculus and others.\n")
            }
        })
    }
    return(ensAttr)
}
