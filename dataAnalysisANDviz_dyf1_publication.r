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
  
  # # select genes significantly upregulated in Dyf1 compared to all other genotypes.
  # set1 <-  subset( dyf1_DEpair$name, 
  #                  ( dyf1_DEpair$X563_ctrl_0_vs_X915_0_p.adj < sigl | 
  #                      dyf1_DEpair$X1058_0_vs_X915_0_p.adj < sigl | 
  #                      dyf1_DEpair$X1716_0_vs_X915_0_p.adj < sigl ) & 
  #                    (dyf1_DEpair$X563_ctrl_0_vs_X915_0_ratio < 0 & 
  #                       dyf1_DEpair$X1058_0_vs_X915_0_ratio < 0 &
  #                       dyf1_DEpair$X1716_0_vs_X915_0_ratio < 0 ) )
  
  # select proteins significantly upregulated in dyf1 relative to all other samples
  set1 <-  subset( dyf1_DEpair$name,
                   ( dyf1_DEpair$X563_ctrl_0_vs_X915_0_p.adj < sigl &  
                       dyf1_DEpair$X563_ctrl_0_vs_X915_0_ratio < 0) & 
                     ( (dyf1_DEpair$X1058_0_vs_X915_0_p.val < 0.05 & 
                          dyf1_DEpair$X1058_0_vs_X915_0_ratio < 0) | 
                         (dyf1_DEpair$X1716_0_vs_X915_0_p.val < 0.05 & 
                            dyf1_DEpair$X1716_0_vs_X915_0_ratio < 0) ) )
  
  

  # subset the full expression table
  XX.m <- XX[which( rownames(XX) %in%  set1),]
  
 
  # make the plot  
  pheatmap::pheatmap( XX.m , col=rev(brewer.pal(n = 10, name = "RdYlBu") ) , 
                      treeheight_row = 0 , fontsize = 14,
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
# of proteins of interest, using Marray data
  {
    set1_gnames <- set1
    # convert names to gene names available in microarray data 
    set1_gnames[set1_gnames %in% 
                  setdiff(set1_gnames, 
                          dyf1_marrayTab$SYMBOL) ] <- c("F42A10.3","MTCE.31",
      "F35E12.6","F33A8.4","C18E9.4","F20G2.1","F20G2.2", "Y50D4B.4")
    # select genes of interest from Microarray table
    # note that "MTCE.31","F35E12.6","F20G2.1","F20G2.2" are not available in microarray
    XX <- MAexprdata[ match( dyf1_marrayTab$PROBEID[ 
      which( dyf1_marrayTab$SYMBOL %in% set1_gnames)] , 
      rownames(MAexprdata) ),]

     
  library(reshape2)
  XX.m <- melt(XX)
  XX.m$gtype <- c(rep("dyf1", nrow(XX.m)/2) , rep("wtype", nrow(XX.m)/2) )
  XX.m$gName <- rep( dyf1_marrayTab$SYMBOL[ which( dyf1_marrayTab$SYMBOL %in% set1_gnames)] ,8)
  library(ggplot2)
  library(ggpubr)
  
  ggplot(XX.m, aes(x=gName,y=value,color=gtype) ) + geom_jitter(size=3, 
        position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8) ) +
    theme_classic() + 
    # remove test name results
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)))
  
  }
  


### Figure 4. b) microarray LFC heatmap of selected neuropeptide genes
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


### Figure 4. a) GSEA pathway analysis
  {
  ### extract neuronal system R-CEL-112316 and sensory perception R-CEL-9709957
    {
      library( data.tree)
      library( DiagrammeR )
      
      # load pathway hirerarchy
      react_hire <- read.table( paste( wdir, "ReactomePathwaysRelation.txt" , sep = ""), 
                                header = F,  sep = "\t")
      react_select <- react_hire[grep("R-CEL" , react_hire$V1 ),]
      
      # extract neuronal system tree
      name1 <- "R-CEL-112316"
      indexNS <- numeric()
      while( !isEmpty( name1 )){
        indexNS <- c( indexNS , which(  react_select$V1%in% name1 ))
        name1 <-  react_select$V2[ react_select$V1%in% name1]
      }
      react_NS <- react_select[ indexNS , ]
      # tree view
      react_NS_tree <- FromDataFrameNetwork( react_NS ) 
      
      # extract sensory perception tree
      name1 <- "R-CEL-9709957"
      indexSP <- numeric()
      while( !isEmpty( name1 )){
        indexSP <- c( indexSP , which(  react_select$V1%in% name1 ))
        name1 <-  react_select$V2[ react_select$V1%in% name1]
      }
      react_SP <- react_select[ indexSP , ]
      # tree view
      react_SP_tree <- FromDataFrameNetwork( react_SP ) 
      
      # convert paths IDs 
      React_pathID <- read.table( paste( wdir, "ReactomePathways.txt" , sep = ""), 
                                  header = F,  sep = "\t", quote = "")
      
      # get reactome IDs of SP and NS pathways
      react_SP.NS <- unique( c(react_SP$V1,react_SP$V2, react_NS$V1 , react_NS$V2 ))
     
    }  
 
  
  ### perform GSEA analysis with reactomePA
    {
    
    library( ReactomePA )
    library(enrichplot)
    options(connectionObserver = NULL)
    
    # get ranked list of genes
    x <- dyf1_marrayTab
    x <- merge( x, tx2gene , by.x ="SYMBOL", by.y = "external_gene_name",all = F)
    x <- aggregate( list( logFC= x$logFC , P.Value=x$P.Value), by= list(entrezgene_id=x$entrezgene_id) , mean)
    dyf1_marrayLFC <- setNames(x$logFC,x$entrezgene_id)
    dyf1_marrayLFC <- dyf1_marrayLFC[order(-dyf1_marrayLFC)]
    
    # perform GSEA
    path_REACT_GSEA <- gsePathway( geneList= dyf1_marrayLFC , organism = "celegans", pAdjustMethod="none",
                                   pvalueCutoff = 0.05 , minGSSize = 1 )
   
    # get a list of all fdr-corrected significant pathways for visualisation
    path_REACT_GSEAfdr0.01 <- gsePathway( geneList= dyf1_marrayLFC , organism = "celegans", pAdjustMethod="fdr",
                pvalueCutoff = 0.01 , minGSSize = 1 )
    
    # add LFC to the result table
    path_REACT_GSEAfdr0.01@result$Ave.log2FC <- sapply( 1:nrow(path_REACT_GSEAfdr0.01), 
                                                    function(ii){
                                                    gset <- unlist(stringr::str_split(path_REACT_GSEAfdr0.01@result$core_enrichment[ii],"/") )
                                                    mean( dyf1_marrayLFC[ names( dyf1_marrayLFC )%in% gset], na.rm =T)
                                                      } )
    # convert EntrezIDs to gene names
    path_REACT_GSEAfdr0.01@result$gNames <- Reduce( c, lapply( 1:nrow(path_REACT_GSEAfdr0.01), function(ii){
      paste( unique( tx2gene$external_gene_name[ match( unlist(stringr::str_split( 
        path_REACT_GSEAfdr0.01@result$core_enrichment[ii],"/") ) , 
        tx2gene$entrezgene_id)] ) , sep = "/", collapse = "/")
    }) )
    
    # save results
    write.table( path_REACT_GSEAfdr0.01, sep = "\t", row.names = F,
                quote = F, file = "path_REACT_GSEAfdr0.01.tsv")
   
    ### plot comprehensive overview of significantly enriched pathways
    # add reverse LFC for better vis with emapplot
    XX <- pairwise_termsim(path_REACT_GSEAfdr0.01)
    XX@result$aveRevLFC <- -1*sapply( 1:nrow(XX), function(ii){
                                gset <- unlist(stringr::str_split(XX@result$core_enrichment[ii],"/") )
                                mean( dyf1_marrayLFC [names( dyf1_marrayLFC )%in% gset], na.rm =T)
                                                          } )
    
    # shorten path names
    XX@result$Description <- sapply(XX@result$Description, wrap_text)
    # Figure 4. a) plot connection between pathways, color nodes based on ave. LFC
    emapplot(pairwise_termsim(XX), color="aveRevLFC",
             showCategory=nrow(XX),min_edge = 0.1, cex_label_category=0.5)
    

    ### select Sensory perception and Neuronal signaling pathways
    # parameter asis=T to return the original type of object 
    path_REACT_GSEA.SP.NS <- path_REACT_GSEAfdr0.01[(path_REACT_GSEAfdr0.01@result$Description %in%
                                               React_pathID$V2[ React_pathID$V1 %in%
                                                                  react_SP.NS ] ) , asis=T ]
    # add a column to the result table that shows if a pathway belongs to 
    # Sensory perception or Neuronal signaling categories
    path_REACT_GSEAfdr0.01@result$SP.NS <- ifelse( (path_REACT_GSEAfdr0.01@result$Description %in%
                                               React_pathID$V2[ React_pathID$V1 %in%
                                                                  react_SP.NS ] ) , "YES","NO")
    write.table( path_REACT_GSEAfdr0.01, sep = "\t", row.names = F,
                 quote = F, file = "path_REACT_GSEAfdr0.01.tsv")
    
    ### cnetplot: plot connection between termes and genes
    # wrap pathway names to improve visualization
    path_REACT_GSEA.SP.NS@result$Description <-  sapply( path_REACT_GSEA.SP.NS@result$Description , wrap_text )
    # make a plot
    PP <- cnetplot( path_REACT_GSEA.SP.NS, foldChange = dyf1_marrayLFC , cex_label_gene=0.6)
    # convert entrezID to gene names
    PP$layers[[4]]$data$name <- tx2gene$external_gene_name[ match( PP$layers[[4]]$data$name , tx2gene$entrezgene_id)]
    # plot
    print( PP )

    # # enrichemtn plots
    # lapply( 1:nrow(path_REACT_GSEA.SP.NS), function(ii){
    #   pdf(height = 6, width = 8, file = paste("/home/tim_nevelsk/PROJECTS/celegans_sensingKavya/GSEA/",
    #                                           rownames(path_REACT_GSEA.SP.NS@result)[ii],"_gsea.pdf",sep = ""))
    #   
    #   print( gseaplot( path_REACT_GSEA, title = (path_REACT_GSEA.SP.NS@result$Description[ii]),
    #                   geneSetID = rownames(path_REACT_GSEA.SP.NS@result)[ii]))
    #   dev.off()
    #   
    #   png(height = 400, width = 600, file = paste("/home/tim_nevelsk/PROJECTS/celegans_sensingKavya/GSEA/",
    #                                           rownames(path_REACT_GSEA.SP.NS@result)[ii],"_gsea.png",sep = ""))
    #   
    #   print( gseaplot( path_REACT_GSEA, title = (path_REACT_GSEA.SP.NS@result$Description[ii]),
    #                    geneSetID = rownames(path_REACT_GSEA.SP.NS@result)[ii]))
    #   dev.off()
    # })
  }
  
  
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