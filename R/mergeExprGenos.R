mergeExprGenos <- function(thisId, inputObj, vars) {
    thisASE <- inputObj$ASE[which(inputObj$ASE$ID == thisId), ]
    all_0 <- cbind(thisASE[, c(colnames(inputObj$covar), "refCount")], 0)
    all_1 <- cbind(thisASE[, c(colnames(inputObj$covar), "altCount")], 1)
    colnames(all_0)[c(ncol(all_0) - 1, ncol(all_0))] <- c("Count", "All")
    colnames(all_1)[c(ncol(all_1) - 1, ncol(all_1))] <- c("Count", "All")
    thisVar <- rbind(all_0, all_1)
    # This shows you the allelic information for the variant you're currently on (i of hetCounts)
    thisVar <- merge(thisVar, inputObj$counts, all.x = FALSE, all.y = FALSE)
    idPairs <- colnames(inputObj$haps)[which(colnames(inputObj$haps) %in% all_0$Ind)]
    
    all <- t(inputObj$haps[which(rownames(inputObj$haps) %in% c(thisId, vars)), which(colnames(inputObj$haps) %in% all_0$Ind)])
    colnames(all)[which(colnames(all) == thisId)]<-"All"
    genos <- data.frame(Ind = idPairs, all, check.names = FALSE)
    #print('genos below') head(genos) colnames(genos)[which(colnames(genos)
    # == as.character(thisId))] <- 'All' colnames(genos)[1] <- 'Ind' Changed this to implement min-pvalue threshold
    exprVariants <- merge(thisVar, genos)
    if (dim(exprVariants)[2] > length(colnames(inputObj$covar)) + 4) {
        count_1 <- colSums(exprVariants[, (length(colnames(inputObj$covar)) + 4):dim(exprVariants)[2]])
        count_0 <- dim(exprVariants)[1] - count_1
        exprVar <- exprVariants[, (length(colnames(inputObj$covar)) + 4):dim(exprVariants)[2]][, which(count_0 > 0 & count_1 > 0)]
        exprVar <- cbind(exprVariants[, 1:(length(colnames(inputObj$covar)) + 3)], exprVar)
    } else if (dim(exprVariants)[2] == length(colnames(inputObj$covar)) + 4) {
        count_1 <- sum(exprVariants[, (length(colnames(inputObj$covar)) + 4):dim(exprVariants)[2]])
        count_0 <- dim(exprVariants)[1] - count_1
        if (count_0 > 0 & count_1 > 0) {
            exprVar <- exprVariants
        } else {
            stop("Specified variant is fixed in individuals")
        }
    } else {
        stop("No variant found. Was the ID specified correctly?")
    }
    if (dim(exprVar)[2] < (length(colnames(inputObj$covar)) + 4)) {
        stop("Variant not found or is fixed in individuals")
    }
    if(sum(count_0,count_1) < 171)
    {
      min_pval <- 1/(factorial(count_0 + count_1)/(factorial(count_0) * factorial(count_1)))
    } else {
      min_pval <- 0
    }
    exprVar <- unique(exprVar)
    return(list(exprVar = exprVar, min_pval = min_pval))
}
