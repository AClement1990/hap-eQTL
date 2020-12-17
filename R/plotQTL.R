#' Generate plots of expression values for given eQTL and expressed heterozygote variant
#'
#' This function takes the output of Run.model and produces various plots for the specified putative regulatory variant.
#' @param modelOut The output of the Run.model function
#' @param aseVar The id of the allele specific expression variant.
#' @param qtlVar The id of the putative regulatory variant.
#' @export
#' @examples
#' #' #Plot the output of Run.model.R stored in modelOut
#' plotQTL(aseDat,'rs116197803', 'rs175183')


plotQTL <- function(modelOut, aseVar, qtlVar, otherAll = FALSE) {
    qtlData <- mergeExprGenos(aseVar, modelOut, qtlVar)
    exprData <- qtlData$exprVar
    exprData<-exprData[with(exprData, order(Ind, All)),]
    exprData$RPM <- exprData$Count/(exprData$Reads/1e+06)
    exprData$alt_all<-exprData[,qtlVar][seq_len(length(exprData[,qtlVar])) + c(1,-1)]
    
    dodge <- position_dodge(width = 0.7)
    plots<-list()
    plots[[1]] <- ggplot(exprData, aes(x = as.factor(get(qtlVar)), y = RPM)) + geom_violin(trim = FALSE, position = dodge, aes(fill = as.factor(get(qtlVar)))) + 
        geom_boxplot(width = 0.2, position = dodge, fill = "white") + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) + 
        ylab(paste("RPM at ", aseVar, sep = "")) + xlab(paste("Allele at ", qtlVar, sep = "")) + theme_pubr() + stat_compare_means(label.x = 1.3) + 
        theme(legend.position = "none")
  
    if(isTRUE(otherAll))
    {
        #p1<-p1 
        plots[[2]] <- ggplot(exprData, aes(x = as.factor(alt_all), y = RPM)) + geom_violin(trim = FALSE, position = dodge, aes(fill = as.factor(get(qtlVar)))) + 
            geom_boxplot(width = 0.2, position = dodge, fill = "white") + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) + 
            ylab(paste("RPM at ", aseVar, sep = "")) + xlab(paste("Allele on other haplotype at ", qtlVar, sep = "")) + theme_pubr() + stat_compare_means(label.x = 1.3) + 
            theme(legend.position = "none") + facet_grid(. ~ as.factor(get(qtlVar)))
    }
    ggarrange(plotlist = plots, ncol=1)
    
    #scatDat<-spread(exprData[,c("Ind","All", "RPM")], key=All, value=RPM)
    #colnames(scatDat)<-c("Ind", "Allele0", "Allele1")
    #qtlGenos<-aggregate(rs140347 ~ Ind, data = exprData, paste, collapse = ",")
    #scatDat<-merge(scatDat, qtlGenos)
    #ggplot(scatDat, aes(x=Allele0, y=Allele1, colour=rs140347)) + geom_point()+ theme_pubr()+ 
    #  geom_smooth(method=lm, se=FALSE, fullrange=FALSE) +
    #  stat_ellipse(type = "norm")
    
    #ggscatter(scatDat, x = "Allele0", y = "Allele1",
    #          color = "rs175183", shape = "rs175183",
    #          ellipse = TRUE, mean.point = TRUE,
    #          star.plot = TRUE)
}
