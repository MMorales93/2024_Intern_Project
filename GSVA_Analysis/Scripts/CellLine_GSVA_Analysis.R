#' ---
#'title: "RUVseq Data: GSVA analysis on cell lines gene expression (RNAseq)"
#'author: "Maritza Morales"
#'date: "`r Sys.Date()`"
#' output:
#'   pdf_document:
#'     highlight: zenburn
#'     fig_caption: yes
#'     toc: TRUE
#'     toc_depth: 3 
#'     includes:  
#'      in_header: Header.tex
#' ---
#+ fig.width=6, fig.height=6 

#+ echo=FALSE, warning=FALSE, message=FALSE
#+ echo=TRUE
#' # Background
#' * The input for this set GSVA analysis is the RUVseq normalized counts from RNASeq data. A total of nine databases are used for the following GSVA analysis. 
#+ echo=FALSE, warning=FALSE, message=FALSE

# Load libraries
library(GSVA)
library(pheatmap)
library(GSEABase)
options(stringsAsFactors=FALSE)
library(matrixStats)

#Sort out top 50 genes
fun.hsd <- function(gsva.mat, nn=50){ 
  rowsd.gsva <-  rowSds(gsva.mat) 
  rowsd.gsva <- rowsd.gsva[rev(order(rowsd.gsva))] 
  gsva.mat[rownames(gsva.mat) %in% names(rowsd.gsva)[1:c(nn)],] 
} 

#+ echo=FALSE, warning=FALSE, message=FALSE

#Databases

# Read in Hallmarks gene set collection
hallmarks <- getGmt("/Users/MMorales/Desktop/Databases/h.all.v2023.1.Hs.symbols.gmt", collectionType = BroadCollection(category = "c1"), geneIdType = SymbolIdentifier())
#Read in Canonical Pathways set collection
c2.cp <- getGmt("/Users/MMorales/Desktop/Databases/c2.cp.v2023.1.Hs.symbols.gmt", collectionType = BroadCollection(category = "c2"), geneIdType = SymbolIdentifier())
#Read in Chemical and Genetic Perturbations set collection
c2.cgp <- getGmt("/Users/MMorales/Desktop/Databases/c2.cgp.v2023.1.Hs.symbols.gmt", collectionType = BroadCollection(category = "c2"), geneIdType = SymbolIdentifier())
#Read in Cancer Gene Neighborhoods set collection
c4 <- getGmt("/Users/MMorales/Desktop/Databases/c4.cgn.v2023.1.Hs.symbols.gmt", collectionType = BroadCollection(category = "c4"), geneIdType = SymbolIdentifier())
#Read in Biological Processes set collection
c5.bp <- getGmt("/Users/MMorales/Desktop/Databases/c5.go.bp.v2023.1.Hs.symbols.gmt", collectionType = BroadCollection(category = "c5"), geneIdType = SymbolIdentifier())
#Read in Component Ontology set collection
c5.cc <- getGmt("/Users/MMorales/Desktop/Databases/c5.go.cc.v2023.1.Hs.symbols.gmt", collectionType = BroadCollection(category = "c5"), geneIdType = SymbolIdentifier())
#Read in Molecular Function set collection
c5.mf <- getGmt("/Users/MMorales/Desktop/Databases/c5.go.mf.v2023.1.Hs.symbols.gmt", collectionType = BroadCollection(category = "c5"), geneIdType = SymbolIdentifier())
#Read in Oncogenic Signature gene sets
c6 <- getGmt("/Users/MMorales/Desktop/Databases/c6.all.v2023.1.Hs.symbols.gmt", collectionType = BroadCollection(category = "c6"), geneIdType = SymbolIdentifier())
#Read in Cell Type Signature gene sets
c8 <- getGmt("/Users/MMorales/Desktop/Databases/c8.all.v2023.1.Hs.symbols.gmt", collectionType = BroadCollection(category = "c8"), geneIdType = SymbolIdentifier())

#+ echo=FALSE, warning=FALSE, message=FALSE

# Read in input matrix with expression values
# Row names will be gene names
expr <- read.csv("02.NormalizedCounts_58CellLines_STAR_GeneCounts_RNAseq_UQ_RUVg.csv", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
#Removing genes with no symbol
expr <- expr[!is.na(expr$Hugo_Symbol),]
#Removing gene that are duplicated
expr <- expr[!(duplicated(expr$Hugo_Symbol)),]
#Using the gene symbols as rownames and eliminating the gene symbol column
rownames(expr) <- expr$Hugo_Symbol
expr$Hugo_Symbol <- NULL

#+ echo=FALSE, warning=FALSE, message=FALSE

# Convert dataframe to matrix
expr.mat <- as.matrix(expr)

cacheDir <- system.file("extdata", package="GSVA")
cachePrefix <- "cache4vignette_"

#Hallmarks
#' # Input: RUVseq normalized counts
#' ## Hallmark gene sets 
#+ echo=FALSE, warning=FALSE, message=FALSE, results='hide'
hm.gsva <- gsva(expr.mat, hallmarks, method = "gsva",min.sz=10, max.sz=200, verbose=TRUE,  kcdf="Poisson", abs.ranking=FALSE)

#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="GSVA score heatmap for Hallmark gene sets (n=50). The samples are ordered by hierarchical clustering of GSVA scores.", fig.height=10, fig.width=9, fig.align = "center"

#Heatmap top 50 genes
hm.es <- fun.hsd(gsva.mat = hm.gsva, nn = 50)

pheatmap(hm.gsva, fontsize_row = 7, fontsize_col = 7)

# Saving the file 
write.csv(hm.gsva, "GSVAscores_58_CellLines_Input_RUVseqNormalizedCounts_Hallmark.csv")


#Canonical Pathways
#' # Input: RUVseq normalized counts
#' ## Canonical Pathways
#+ echo=FALSE, warning=FALSE, message=FALSE, results='hide'
c2cp.gsva <- gsva(expr.mat, c2.cp, method = "gsva",min.sz=10, max.sz=200, verbose=TRUE,  kcdf="Poisson", abs.ranking=FALSE)

#Saving the file
write.csv(c2cp.gsva, "GSVAscores_58_CellLines_Input_RUVseqNormalizedCounts_C2_Canonical Pathways.csv")

#Heatmap top 50 genes
c2cp.es <- fun.hsd(gsva.mat = c2cp.gsva, nn = 50)

#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="GSVA score heatmap for Canonical Pathways gene sets (n=50). The samples are ordered by hierarchical clustering of GSVA scores.", fig.height=10, fig.width=9, fig.align = "center"
pheatmap(c2cp.es, fontsize_row = 6, fontsize_col = 6)


#Chemical and Genetic Perturbations
#' # Input: RUVseq normalized counts
#' ## Chemical and Genetic Perturbations 
#+ echo=FALSE, warning=FALSE, message=FALSE, results='hide'
c2cgp.gsva <- gsva(expr.mat, c2.cgp, method = "gsva",min.sz=10, max.sz=200, verbose=TRUE,  kcdf="Poisson", abs.ranking=FALSE)

#Saving the file
write.csv(c2cgp.gsva, "GSVAscores_58_CellLines_Input_RUVseqNormalizedCounts_C2_ChemicalAndGenticPerturbations.csv")

#Heatmap top 50 genes
c2cgp.es <- fun.hsd(gsva.mat = c2cgp.gsva, nn = 50)

#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="GSVA score heatmap for Chemical and Genetic Perturbations gene sets (n=50). The samples are ordered by hierarchical clustering of GSVA scores.", fig.height=10, fig.width=9, fig.align = "center"
pheatmap(c2cgp.es, fontsize_row = 6, fontsize_col = 6)


#Cancer Gene Neighborhoods
#' # Input: RUVseq normalized counts
#' ## Cancer Gene Neighborhood gene sets 
#+ echo=FALSE, warning=FALSE, message=FALSE, results='hide'
c4.gsva <- gsva(expr.mat, c4, method = "gsva",min.sz=10, max.sz=200, verbose=TRUE,  kcdf="Poisson", abs.ranking=FALSE)

# Saving the file 
write.csv(c4.gsva, "GSVAscores_58_CellLines_Input_RUVseqNormalizedCounts_C4_Cancer Gene Neighborhoods.csv")

#Heatmap top 50 genes
c4.es <- fun.hsd(gsva.mat = c4.gsva, nn = 50)

#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="GSVA score heatmap for Cancer Gene Neighborhoods gene sets (n=50). The samples are ordered by hierarchical clustering of GSVA scores.", fig.height=10, fig.width=9, fig.align = "center"
pheatmap(c4.es, fontsize_row = 7, fontsize_col = 6)


#Biological Processes
#' # Input: RUVseq normalized counts
#' ## Biological Processes gene sets 
#+ echo=FALSE, warning=FALSE, message=FALSE, results='hide'
c5bp.gsva <- gsva(expr.mat, c5.bp, method = "gsva",min.sz=10, max.sz=200, verbose=TRUE,  kcdf="Poisson", abs.ranking=FALSE)

# Saving the file 
write.csv(c5bp.gsva, "GSVAscores_58_CellLines_Input_RUVseqNormalizedCounts_C5_Biological_Processes.csv")

#Heatmap top 50 genes
c5bp.es <- fun.hsd(gsva.mat = c5bp.gsva, nn = 50)

#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="GSVA score heatmap for Biological Processes gene sets (n=50). The samples are ordered by hierarchical clustering of GSVA scores.", fig.height=10, fig.width=9, fig.align = "center"
pheatmap(c5cc.es, fontsize = 5)


#Component Ontology
#' # Input: RUVseq normalized counts
#' ## Component Ontology gene sets 
#+ echo=FALSE, warning=FALSE, message=FALSE, results='hide'
c5cc.gsva <- gsva(expr.mat, c5.cc, method = "gsva",min.sz=10, max.sz=200, verbose=TRUE,  kcdf="Poisson", abs.ranking=FALSE)

# Saving the file 
write.csv(c5cc.gsva, "GSVAscores_58_CellLines_Input_RUVseqNormalizedCounts_C5_Component_Ontology.csv")

#Heatmap top 50 genes
c5cc.es <- fun.hsd(gsva.mat = c5cc.gsva, nn = 50)

#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="GSVA score heatmap for Component Ontology gene sets (n=50). The samples are ordered by hierarchical clustering of GSVA scores.", fig.height=10, fig.width=9, fig.align = "center"
pheatmap(c5cc.es, fontsize = 5)


#Molecular Function
#' # Input: RUVseq normalized counts
#' ## Molecular Function gene sets 
#+ echo=FALSE, warning=FALSE, message=FALSE, results='hide'
c5mf.gsva <- gsva(expr.mat, c5.mf, method = "gsva",min.sz=10, max.sz=200, verbose=TRUE,  kcdf="Poisson", abs.ranking=FALSE)

# Saving the file 
write.csv(c5mf.gsva, "GSVAscores_58_CellLines_Input_RUVseqNormalizedCounts_C5_Molecular_Function.csv")

#Heatmap top 50 genes
c5mf.es <- fun.hsd(gsva.mat = c5mf.gsva, nn = 50)

#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="GSVA score heatmap for Molecular Function gene sets (n=50). The samples are ordered by hierarchical clustering of GSVA scores.", fig.height=10, fig.width=9, fig.align = "center"
pheatmap(c5mf.es, fontsize_row = 5, fontsize_col = 5)


#Oncogenic Signature
#' # Input: RUVseq normalized counts
#' ## Oncogenic Signature gene sets 
#+ echo=FALSE, warning=FALSE, message=FALSE, results='hide'
c6.gsva <- gsva(expr.mat, c6, method = "gsva",min.sz=10, max.sz=200, verbose=TRUE,  kcdf="Poisson", abs.ranking=FALSE)

# Saving the file 
write.csv(c6.gsva, "GSVAscores_58_CellLines_Input_RUVseqNormalizedCounts_C6_Oncogenic Signature.csv")

#Heatmap top 50 genes
c6.es <- fun.hsd(gsva.mat = c6.gsva, nn = 50)

#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="GSVA score heatmap for Oncogenic Signature gene sets (n=50). The samples are ordered by hierarchical clustering of GSVA scores.", fig.height=10, fig.width=9, fig.align = "center"
pheatmap(c6.es, fontsize = 6)


#Cell Type
#' # Input: RUVseq normalized counts
#' ## Cell Type Signature gene sets 
#+ echo=FALSE, warning=FALSE, message=FALSE, results='hide'
c8.gsva <- gsva(expr.mat, c8, method = "gsva",min.sz=10, max.sz=200, verbose=TRUE,  kcdf="Poisson", abs.ranking=FALSE)

# Saving the file 
write.csv(c8.gsva, "GSVAscores_58_CellLines_Input_RUVseqNormalizedCounts_C8_Cell Type.csv")

#Heatmap top 50 genes
c8.es <- fun.hsd(gsva.mat = c8.gsva, nn = 50)

#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="GSVA score heatmap for Cell Type Signature gene sets (n=50). The samples are ordered by hierarchical clustering of GSVA scores.", fig.height=10, fig.width=9, fig.align = "center"
pheatmap(c8.es, fontsize = 5)