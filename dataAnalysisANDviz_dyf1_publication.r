### this code replicates microarray and mass-spectrometry data analysis 
### and generates figures used in A.Segref manuscript: 
### "Thermosensation is linked /coupled to ubiquitin-dependent 
### protein turnover through insulin and calcineurin signaling"

# load library
library(ggrepel)
library(ggplot2)


### get a table with different IDs and annotation from the BioMart
mart_cel <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "celegans_gene_ensembl")
tx2gene <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "entrezgene_id" , "ensembl_gene_id","external_gene_name","transcript_biotype","gene_biotype"), mart = mart_cel)
# path to the directory containing results of DE analysis of the masspec data
wdir <- "Revision.data_CelSensing_20.10.2021/"
setwd(wdir)
# define function for making 1) sample distance heatmap and 2) PCA 
PCA_heatm_plot=function(count_table, samples=NA, groups , 
                       logtrans=F, title="PCA")
  {
  # sample - to use as labels of individual samples
  # groups - experimental groups, e.g. wild-type and mutant
  # logtrans - log transform the data
  
  ### Estimate and plot eucledean disances between samples
  # !!!! remember that not TPM but count data shall be normalized with rlog !!!!
  library("RColorBrewer")
  library("pheatmap")  
  require(ggplot2)
  
  # define sample names as column names if sample names are not provided
  if ( length(samples) == ncol(count_table)) colnames(count_table) <- samples else samples <- colnames(count_table)
  
  # select top 10 variable genes
  count_table <- as.matrix(count_table)
  count_table <- count_table[ which( rowVars(count_table) > quantile( rowVars(count_table) , 0.9)) ,  ]
  
  # log transform if necessary
  if( isTRUE(logtrans)) rld <- log( count_table + 1 ) else rld <- count_table
  rownames(rld)=rownames(count_table)
 
  # calculate distances
  euclDists <- dist( t( rld ) )
  euclDistsMatrix <- as.matrix( euclDists )
  rownames(euclDistsMatrix) <- names(count_table)
  colnames(euclDistsMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(euclDistsMatrix,
           clustering_distance_rows=euclDists,
           clustering_distance_cols=euclDists,
           col=colors , labels_row =  samples )
  
  
  ## PCA for count data rlog transform  using DEseq2 algorythm
  pr_comp_rnorm <- prcomp(t(rld),scale=F)
  #biplot(pr_comp_rnorm)
  XX = as.data.frame(pr_comp_rnorm$x)
  g <- qplot( x=PC1, y=PC2, data=XX,  colour=factor(groups), main=title) + 
    geom_point(size=4  ) 
  # labels and dots sizes
  g+theme_bw() + theme(axis.text=element_text(size=20),
                       axis.title=element_text(size=24,face="bold"), legend.text = element_text(size = 24),legend.title = element_text(size=24))+
    scale_size(range = c(5, 6))+theme(legend.key.size = unit(2, "cm"))+ labs(colour = "condition")+
    geom_text(aes(label=rownames(XX)),hjust=0.4, vjust=-0.8,size=5)
}


###================= Masspec data analysis of wild-type, dyf1, dyf1Unc13 and Dyf1Unc31 strains ================
library(DEP)

### DE analysis is performed by running DEP::run_app("LFQ") on proteinGroups.txt Masspec output, with default parameters
### Load raw data and results of DE analysis. 
# In total, 3922 proteins passed all filters and were analysed
# 68 proteins are significantly DE under qvalue < 0.1 cut-off
  {
    
    # load raw MasSpec intensities MA
  dyf1_intens <- read.table(sep = "\t",stringsAsFactors = F , header = T, 
                          file= paste( wdir, "analyis_replication/raw_data/masspec_LFQfiltered.tsv", sep = "" ) )
    # results of the DE analysis: pairwise comparison of dyf1 to all other genotypes
  dyf1_DEpair <- read.table(sep = "\t",stringsAsFactors = F , header = T, 
                            file=paste( wdir, "analyis_replication/allVSdyf1_masspecDEresults.tsv" , sep = "" ))
}

### Figure 3. a) PCA and volcano plots of DE results
  {
    ### plot samples distance heatmap and PCA 
    datTOplot <- as.matrix( dyf1_intens[,-1] )
    PCA_heatm_plot( datTOplot , groups = sub("_.*","", colnames(datTOplot)) , 
                            logtrans=F) 
    
  ### volcano plot
  datTOplot <- dyf1_DEpair
  # add a column of NAs
  datTOplot$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  datTOplot$diffexpressed[datTOplot$X563_ctrl_0_vs_X915_0_ratio < 0 & datTOplot$X563_ctrl_0_vs_X915_0_p.adj < 0.1] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  datTOplot$diffexpressed[datTOplot$X563_ctrl_0_vs_X915_0_ratio > 0 & datTOplot$X563_ctrl_0_vs_X915_0_p.adj < 0.1] <- "DOWN"
  
  
  
  # Now write down the name of genes beside the points...
  # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  datTOplot$delabel <- NA
  datTOplot$delabel[datTOplot$X563_ctrl_0_vs_X915_0_p.adj < 0.1 ] <- datTOplot$name[datTOplot$X563_ctrl_0_vs_X915_0_p.adj < 0.1 ]
  
  # reverse LFC to make wild-type the reference
  datTOplot$'-log2FC' <- -1*(datTOplot$X563_ctrl_0_vs_X915_0_ratio)
  datTOplot$'-log10(qvalue)' <- -log10(datTOplot$X563_ctrl_0_vs_X915_0_p.adj)
  
  
  # plot adding up all layers we have seen so far
  ggplot(data=datTOplot, aes( y=datTOplot$'-log10(qvalue)', 
                              x= datTOplot$'-log2FC' ,
                              col=diffexpressed, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_hline(yintercept=-log10(0.1), col="red") + 
    labs(y="-log10(qvalue)", x = "-log2FC", colour = "Diff. Expr") +
    scale_x_continuous( n.breaks = 8)+
    theme(axis.text=element_text(size=14 ), legend.title=element_text(size=16),
          axis.title=element_text(size=16), legend.text=element_text(size=16)) 
  
  
  }

### Figure 3. b) Heatmap of averaged MasSpec intensities for genes of interests values averaged across genotypes ###
  {
  library(RColorBrewer)
  sigl = 0.1
  XX <- dyf1_intens[,2:13]
  rownames(XX) <- dyf1_intens$name
  XX <- t (scale( t(XX ))) # scale and center the data to get normalised values
  # take the mean of 3 replicates of each genotype
  XX <- data.frame( wt = rowMeans(XX[,1:3]) , dyf1 = rowMeans(XX[,4:6]) , dyf1unc13 = rowMeans(XX[,7:9]) ,  dyf1unc31 = rowMeans(XX[,10:12]))
  
  # select genes significantly upregulated in Dyf1 compared to all other genotypes.
  set1 <-  subset( dyf1_DEpair$name, 
                   ( dyf1_DEpair$X563_ctrl_0_vs_X915_0_p.adj < sigl | 
                       dyf1_DEpair$X1058_0_vs_X915_0_p.adj < sigl | 
                       dyf1_DEpair$X1716_0_vs_X915_0_p.adj < sigl ) & 
                     (dyf1_DEpair$X563_ctrl_0_vs_X915_0_ratio < 0 & 
                        dyf1_DEpair$X1058_0_vs_X915_0_ratio < 0 &
                        dyf1_DEpair$X1716_0_vs_X915_0_ratio < 0 ) )
  # select 21 genes
  
  # subset the full expression table
  XX.m <- XX[which( rownames(XX) %in%  set1),]
  
  # make the plot  
  pheatmap::pheatmap( XX.m , col=rev(brewer.pal(n = 10, name = "RdYlBu") ) , 
                      treeheight_row = 0 , fontsize = 9,
                      scale="row" , clustering_method = "ward.D2" )
  
}
  
###=================  Microarray data analysis of wild-type and dyf-1 mutants ================
# ### processing and DE analysis of MA data 
#   {
#   library(gcrma)
#   library(limma)
#   library(oligo)
#   library("pd.elegene.1.0.st")
#   library("annotate")
#   
#   # first, download raw microarray data from GEO in CELL file format, save in dir cell_files
#   # list download files in the cell_files folder
#   targets <- list.files( path = paste( wdir,"cell_files", sep = ""),
#                                 full.names = T) 
#   rawData <- read.celfiles(targets)
#   ppData <- rma(rawData) # do preprocessing and normalization of the data 
#   library(affycoretools) # library to annotate affymetrix dataset
#   ppData <- annotateEset(ppData, "pd.elegene.1.0.st") # annotate Affymetrix geneset, matching probe IDs with genes etc
#   
#   
#   # make experiment design table 
#   wafer <- sub("-","",substr(sampleNames(rawData), 9, 11))
#   labels <- sub("N2","wild-type",sub("dyf","dyf-1(-)",sub("-","_",sub("_","",substr(sampleNames(rawData), 9, 13)))))
# 
#   
#   ## DE analysis
#   design <- cbind( WT=1 , MUvsWT=wafer=="dyf")
#   fit <- lmFit(ppData, design)
#   fit <- eBayes(fit)
#   dyf1_marrayTab <- topTable(fit, coef="MUvsWT", adjust="BH",p.value=1, number=Inf)
#   dyf1_marrayTab_sig <- topTable(fit, coef="MUvsWT", adjust="BH",p.value=0.01, number=Inf)
#   
#   }

### Load results of DE analysis and raw MA data 
# in total 26293 genes were analysed , 
# 2939 genes are significntly DE under qvalue < 0.01 cut-off
MAexprdata <- readRDS( paste( wdir, "analyis_replication/raw_data/dyf1_marray_rawRMA.rda" , sep = "") )
dyf1_marrayTab <- read.table( paste( wdir, "analyis_replication/dyf1_marray_DEall.csv" , sep = ""), 
                              header = T, row.names = 1 , sep = "\t")

### PCA and volcano plot of MA data 
  {
  
    ### PCA plot
    PCA_heatm_plot( MAexprdata , samples=labels , groups=design)
    
    ### Volcano plot
    {
      datTOplot <- dyf1_marrayTab
      # add a column of NAs
      datTOplot$diffexpressed <- "NO"
      # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
      datTOplot$diffexpressed[datTOplot$logFC > 0 & datTOplot$adj.P.Val < 0.1] <- "UP"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      datTOplot$diffexpressed[datTOplot$logFC < 0 & datTOplot$adj.P.Val < 0.1] <- "DOWN"
      
      
      
      # Now write down the name of genes beside the points...
      # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
      datTOplot$delabel <- NA
      datTOplot$delabel[datTOplot$adj.P.Val < 0.01 ] <- datTOplot$SYMBOL[datTOplot$adj.P.Val < 0.01 ]
      
      # reverse LFC to make wild-type the reference
      datTOplot$'-log10(qvalue)' <- -log10(datTOplot$adj.P.Val)
      
      
      # plot adding up all layers we have seen so far
      ggplot(data=datTOplot, aes( y=datTOplot$'-log10(qvalue)', 
                                  x= datTOplot$logFC ,
                                  col=diffexpressed, label=delabel)) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_hline(yintercept=-log10(0.1), col="red") + 
        labs(y="-log10(qvalue)", x = "-log2FC", colour = "Diff. Expr") +
        scale_x_continuous( n.breaks = 8)+
        theme(axis.text=element_text(size=14 ), legend.title=element_text(size=16),
              axis.title=element_text(size=16), legend.text=element_text(size=16)) 
      
    }
    
  }

### Figure 3. supplement 1A. plot mRNA expression lvls
# of ANMT-2 and C05D11.5 genes, using Marray data
  {
    # load results of MA DE analysis
  dyf1_marrayTab <- dyf1_marrayTab
  # select peptides of 2 genes of interest
  XX <- MAexprdata[which(rownames(MAexprdata) %in% dyf1_marrayTab$PROBEID[ which( dyf1_marrayTab$SYMBOL %in% c("F42A10.3", "C05D11.5"))]),]
  
  library(reshape2)
  XX.m <- melt(XX)
  XX.m$gtype <- c(rep("dyf1",8) , rep("wtype",8) )
  XX.m$gName <- rep(c("C05D11.5","anmt-2"),8)
  library(ggplot2)
  library(ggpubr)
  ggplot(XX.m, aes(x=gName,y=value,color=gtype) ) + geom_jitter(size=3, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8) ) +
    theme_classic() + stat_compare_means()
}

### Figure 4. a) microarray LFC heatmap of selected neuropeptide genes
  {
    # aggregate different probes for same genes
    dyf1_marrayTab_agg <- aggregate(.~SYMBOL , data=dyf1_marrayTab[,c(3,5,8,9)] , FUN=mean )
    # select genes
    neuropept_sel <- dyf1_marrayTab_agg[which(grepl("nlp|ins-|flp",dyf1_marrayTab_agg$SYMBOL,ignore.case = T) & dyf1_marrayTab_agg$adj.P.Val<0.01),]
    
    # ggplot
    library(ggplot2)
    library(scales)
  
  # plot
    tb <- neuropept_sel
  tb$SYMBOL <- as.character(tb$SYMBOL)
  # tb$SYMBOL[3] <- paste(as.character(tb$SYMBOL[3]),".",sep = "")
  tb$SYMBOL <- factor(  tb$SYMBOL , levels =  tb$SYMBOL[order(tb$logFC)], labels =   tb$SYMBOL[order(tb$logFC)] )
  ggplot(data = tb, aes(x = " " , y = SYMBOL )) +  geom_tile(aes(fill = logFC)) +
    theme_bw()+theme( panel.grid.major = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
                      panel.border = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                      axis.text.x = element_text(face="bold", angle=0,vjust=0.5,size=11),axis.text.y = element_text(size=11,face="italic",family = "arial"))+
    scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint=0,  guide = "colorbar", limits=c(min((tb$logFC), na.rm = T),max((tb$logFC), na.rm = T)), name="log2FC of\ngene expression")  # color for pcorr as corr
  
  
}


### SPIA pathway analysis
  {
  library("graphite")
  library("SPIA")
    options(connectionObserver = NULL)
    
  # download prepared KEGG and REACTOME pathways to run with graphite
    # KEGG: 117 entries, retrieved on 13-10-2020
    # Reactome : 1250 entries, retrieved on 13-10-2020
    { 
    path_REACT <- pathways("celegans", "reactome")
    path_REACT <- convertIdentifiers( path_REACT, "ENTREZID")
    prepareSPIA(path_REACT, "reactomeAll")
    
    path_KEGG <-pathways("celegans", "kegg")
    path_KEGG <- convertIdentifiers( path_KEGG, "ENTREZID")
    prepareSPIA(path_KEGG, "KEGGAll")
  }
  
  
  ### define background gene set
    ALL_genes <- unique( tx2gene$entrezgene_id )
    ALL_genes <- paste("ENTREZID", ALL_genes, sep = ":") # modify gene names 
  
  
  ### DE genes in entrez IDs
  x <- dyf1_marrayTab
  x <- merge( x[ which( x$adj.P.Val < 0.05), ], tx2gene , by.x ="SYMBOL", by.y = "external_gene_name",all = F)
  x <- aggregate(logFC~entrezgene_id, x ,sum)
  x <- setNames(x$logFC,paste("ENTREZID", x$entrezgene_id, sep = ":"))
  
  ### run SPIA
  SPIA_REACT <- runSPIA(de=x, all=ALL_genes, pathwaySetName="reactomeAll",verbose=F,nB=2000)
  SPIA_REACT_sig <- SPIA_REACT[ which(SPIA_REACT$pGFWER <= 0.05), ]
  SPIA_KEGG <- runSPIA( de=x, all=ALL_genes, pathwaySetName="KEGGAll",verbose=F,nB=2000)
  SPIA_KEGG_sig=SPIA_KEGG[ which(SPIA_KEGG$pGFWER <= 0.05), ]
  
  ### targeted SPIA analysis
  
  
}

### Fig.4 B) SPIA two-way evidence plot
  {
    
  # define a function
  spia2way_ggplot <- function( pathtype="",  SPIA_results = SPIA_REACT, 
                               combinemethod = "fisher", threshold = 0.05) 
    {
    # pathtype - plot title
    # combinemethod - method to combine p-values
    # threshold - threshold of SPIA pGFdr significance
    
    library(SPIA)
    library(ggplot2)
    library(ggrepel)
    # helper function from SPIA
    getP2<-function(pG, combine=combinemethod)
    {
      #given a pG returns two equal p-values such as   combfunc(p1,p2)=pG
      if(combine=="fisher"){
        ch=qchisq(pG,4,lower.tail = FALSE)
        return(sqrt(exp(-ch/2)))
      }
      
      if(combine=="norminv"){
        return(pnorm(qnorm(pG)*sqrt(2)/2))
        
      }
    }
    
    x <- SPIA_results
    x$pPERT[is.na(x$pPERT)] <- 1
    x$pNDE[is.na(x$pNDE)] <- 1
    
    
    tr1<-threshold/dim(na.omit(x))[1]
    trold=tr1
    tr2<-max(x[,"pG"][x[,"pGFdr"]<=threshold])
    if(tr2<=trold) tr2=trold*1.03
    
    # color points based on significance
    x$colourSig <- ifelse( x$pG < tr1, "red",  ifelse( x$pG < tr2, "blue","dark grey") )
    x$colourDir <- ifelse( (x$Status =="Activated" & x$pG < tr2 ), "red", ifelse( x$pG < tr2, "blue","dark grey") )
    
    x$NameSig <-  ifelse( x$pG < tr2, x$Name, "" )
    
    # make a plot
    PP <-  ggplot(x, aes(x=-log(pNDE),y=-log(pPERT)) )+ geom_point( aes(x =-log(pNDE), y=-log(pPERT),colour=colourDir , size = pSize ) , alpha=0.8 ) +   
      scale_size( range = c(1,20)) + labs(title = pathtype ) +
      scale_colour_manual(name = 'Pathway status', values = setNames(c("dark grey","blue","red") ,c("dark grey","blue","red")), labels=c('down-regulated','non-sig','up-regulated'))+
      geom_abline(aes(intercept = -log(getP2(tr1,combinemethod)^2), slope = -1), color = "green" ,lwd=1)+
      geom_abline(aes(intercept = -log(getP2(tr2,combinemethod)^2), slope = -1), color = "grey50" ,lwd=1) + 
      theme_minimal() + theme(plot.title = element_text(face="bold",size=20),axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"), 
                              legend.text = element_text(size = 20),legend.title = element_text(size=20),legend.key.size = unit(2, "cm"),strip.text.y = element_text(size = 20))+
      guides(colour = guide_legend(override.aes = list(size=20))) + # control size of of points in legend 
      geom_text_repel( size = 6 , point.padding =1 , aes(label = x$NameSig),box.padding   = 1, segment.color = 'grey50')  #ensure that labels are not jammed together
   
    
    print(PP)
  }
  
    # load the data 
    SPIA_REACTold <- readRDS("/home/tim_nevelsk/PROJECTS/celegans_sensingKavya/SPIA_REACT.dyf1MA.rda")
    SPIA_KEGGold <- readRDS("/home/tim_nevelsk/PROJECTS/celegans_sensingKavya/SPIA_KEGG.dyf1MA.rda")
    
    # render a plot 
    spia2way_ggplot( SPIA_results = SPIA_REACT , threshold=0.05 , 
                    pathtype="REACTOME" , combinemethod = "fisher" )
  
  }

### Figure 5. a) visualize Cel04020 KEGG pathway
  {
    # prepare data 
  x <- dyf1_marrayTab
  x <- merge( x , tx2gene , by.x ="SYMBOL", by.y = "external_gene_name",all = F)
  x <- aggregate(logFC~entrezgene_id, x ,sum)
  x <- setNames(x$logFC,x$entrezgene_id)
  
  # visualize
  library(pathview)
  pathview( gene.data=x, kegg.dir= wdir , pathway.id="04020", 
            species="cel", kegg.native=T, out.suffix = "dyf1_KEGG",
            bins=list(gene=20, cpd=10) )
  
}