#' ---
#'title: "Comparison between DESeq2 Normalization and RUVg Normalization"
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
#' # Background Information
#' * Upon starting in the Duval Lab, I was given RNA seq data that had already need normalized with using two methods: 
#'    * DESeq2
#'    * RUVg
#' * I was tasked with comparing the 2 normalization methods via the following methods: 
#'    * Basic Statistics 
#'      * Mean Gene Expression
#'      * Median Gene Expression
#'      * Max Gene Expression
#'    * Correlation Plots
#'      * Correlation plots were generated using GSVA score data for both normalization methods
#' * This task served two purposes: 
#'    1. The lab had wanted to compare both methods, but due to the current project load when I started this project was put on the back burner. My arrival meant I could take some of those projects that had been deprioritized and work on them, helping to decrease the pending project work load.  
#'    2. This project allowed me to become familiar and comfortable with utilizing R and RStudio. Prior to the start of the program I had extremely limited coding experience. I was very up front with Dr. Duval about this when she interviewed me. This project was her way of allowing me the time I needed to became more familar with the skills I needed to be successful in the lab. 
#+ echo=FALSE, warning=FALSE, message=FALSE
#Load Packages
library(tidyverse) 
library(corrplot)
library(EDASeq)
library(conflicted) 
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::filter) 
conflicts_prefer(dplyr::count) 

#Load Data
meta <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Input/CancerType_Canine_CellLines.csv")
# Listing cancer types as factors 
x <- as.factor(meta$Cancer_Type)
b <- as.factor(meta$Batch)
#DESeq2 Data----
des_norm <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Input/04.NormalizedCounts_58CellLines_STAR_GeneCounts_RNAseq_Deseq2.csv", row.names = 1, check.names = FALSE)
##GSVA ----
hm_des_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/DEgsva/GSVAscores_58_CellLines_Input_DESeqNormalizedCounts_Hallmark.csv", row.names = 1)
cp_des_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/DEgsva/GSVAscores_58_CellLines_Input_DESeqNormalizedCounts_C2_Canonical Pathways.csv", row.names = 1)
cp_des_gsva <- cp_des_gsva[!grepl("SARS", rownames(cp_des_gsva)), ]
cgp_des_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/DEgsva/GSVAscores_58_CellLines_Input_DESeqNormalizedCounts_C2_ChemicalAndGenticPerturbations.csv", row.names = 1)
cgp_des_gsva <- cgp_des_gsva[!grepl("SARS", rownames(cgp_des_gsva)), ]
cgn_des_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/DEgsva/GSVAscores_58_CellLines_Input_DESeqNormalizedCounts_C4_Cancer Gene Neighborhoods.csv", row.names = 1)
bp_des_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/DEgsva/GSVAscores_58_CellLines_Input_DESeqNormalizedCounts_C5_Biological_Processes.csv", row.names = 1)
cc_des_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/DEgsva/GSVAscores_58_CellLines_Input_DESeqNormalizedCounts_C5_Component_Ontology.csv", row.names = 1)
mf_des_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/DEgsva/GSVAscores_58_CellLines_Input_DESeqNormalizedCounts_C5_Molecular_Function.csv", row.names = 1)
c6_des_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/DEgsva/GSVAscores_58_CellLines_Input_DESeqNormalizedCounts_C6_Oncogenic Signature.csv", row.names = 1)
c8_des_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/DEgsva/GSVAscores_58_CellLines_Input_DESeqNormalizedCounts_C8_Cell Type.csv", row.names = 1)

#RUVg Data ----
ruv_norm <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Input/02.NormalizedCounts_58CellLines_STAR_GeneCounts_RNAseq_UQ_RUVg.csv", row.names = 1, check.names = FALSE)
##GSVA ----
hm_ruv_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/RUVgsva/GSVAscores_58_CellLines_Input_RUVSeqNormalizedCounts_Hallmark.csv", row.names = 1)
cp_ruv_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/RUVgsva/GSVAscores_58_CellLines_Input_RUVSeqNormalizedCounts_C2_Canonical Pathways.csv", row.names = 1)
cp_ruv_gsva <- cp_ruv_gsva[!grepl("SARS", rownames(cp_ruv_gsva)), ]
cgp_ruv_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/RUVgsva/GSVAscores_58_CellLines_Input_RUVSeqNormalizedCounts_C2_ChemicalAndGenticPerturbations.csv", row.names = 1)
cgp_ruv_gsva <- cgp_ruv_gsva[!grepl("SARS", rownames(cgp_ruv_gsva)), ]
cgn_ruv_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/RUVgsva/GSVAscores_58_CellLines_Input_RUVSeqNormalizedCounts_C4_Cancer Gene Neighborhoods.csv", row.names = 1)
bp_ruv_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/RUVgsva/GSVAscores_58_CellLines_Input_RUVSeqNormalizedCounts_C5_Biological_Processes.csv", row.names = 1)
cc_ruv_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/RUVgsva/GSVAscores_58_CellLines_Input_RUVSeqNormalizedCounts_C5_Component_Ontology.csv", row.names = 1)
mf_ruv_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/RUVgsva/GSVAscores_58_CellLines_Input_RUVSeqNormalizedCounts_C5_Molecular_Function.csv", row.names = 1)
c6_ruv_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/RUVgsva/GSVAscores_58_CellLines_Input_RUVSeqNormalizedCounts_C6_Oncogenic Signature.csv", row.names = 1)
c8_ruv_gsva <- read.csv("/Users/MMorales/Desktop/Duval_Lab/GSVACellLines8.22.23/2023.08.21_GSVA_Analysis_58_Canine_CellLines/58_Cell_Line/Output/RUVgsva/GSVAscores_58_CellLines_Input_RUVSeqNormalizedCounts_C8_Cell Type.csv", row.names = 1)

#Cell Line names to ensure that all cell lines are spelled the same way across all data frames
cell_lines <- c("DEN-HSA", "Cindy", "SB", "FACC-19-CHP25", "DH82", "Nike", "FACC-21-CHS72", "1771", "CLBL1", "CLL1390", "OSW", "CMT12", "BRMCT", "C2", "CML-10C2", "CML-6M", "Parks", "17CM98", "Jones", "FACC-19-CMN18", "D-17", "MacKinley", "Abrams", "FACC-18-COS03", "FACC-19-COS12", "FACC-19-COS17", "FACC-19-COS32", "Gracie", "HMPOS", "Moresco", "OS2-4", "OSA8", "Yamane", "FACC-21-COS124", "FACC-22-COS150", "FACC-21-COS73", "FACC-21-COS74", "FACC-21-COS77", "FACC-21-COS79", "FACC-21-COS81", "FACC-21-COS83", "FACC-22-CPNST161", "FACC-21-CPU97", "FACC-21-CPU98", "STSA-1", "FACC-19-CSTS36", "FACC-22-STS132", "FACC-22-STS138", "FACC-19-CPU26", "FACC-19-CSC27", "FACC-22-CSC165", "CTAC", "FACC-21-CTHY80", "Angus", "Bliley", "Kinsey", "Tyler1", "Tyler2")

colnames(hm_des_gsva) <- cell_lines
colnames(hm_ruv_gsva) <- cell_lines
colnames(cp_des_gsva) <- cell_lines
colnames(cp_ruv_gsva) <- cell_lines
colnames(cgp_des_gsva) <- cell_lines
colnames(cgp_ruv_gsva) <- cell_lines
colnames(cgn_des_gsva) <- cell_lines
colnames(cgn_ruv_gsva) <- cell_lines
colnames(bp_des_gsva) <- cell_lines
colnames(bp_ruv_gsva) <- cell_lines
colnames(cc_des_gsva) <- cell_lines
colnames(cc_ruv_gsva) <- cell_lines
colnames(mf_des_gsva) <- cell_lines
colnames(mf_ruv_gsva) <- cell_lines
colnames(c6_des_gsva) <- cell_lines
colnames(c6_ruv_gsva) <- cell_lines
colnames(c8_des_gsva) <- cell_lines
colnames(c8_ruv_gsva) <- cell_lines

#+ echo=FALSE, warning=FALSE, message=FALSE
## Bar Graph ----
#Generate count chart of how many cell lines per cancer type
Chart <- meta %>%
  count(Cancer_Type) %>%
  group_by(Cancer_Type)

#Define colors to cancer type
col <- c("Hemangiosarcoma"= "#A6CEE3","Hepatocellular Carcinoma"= "#1F78B4","Histiocytic Sarcoma"= "#B2DF8A" ,"Leukemia / Lymphoma"= "#33A02C","Mammary Carcinoma"=  "#FB9A99","Mast Cell Tumor"= "#E31A1C","Melanoma"= "#FDBF6F","Meningioma"= "#FF7F00" ,"Osteosarcoma"= "#CAB2D6","Soft Tissue Sarcoma"=  "#6A3D9A","Squamous cell carcinoma"= "cyan2","Thyroid Carcinoma"= "#B15928" ,"Transitional Cell Carcinoma" ="gray60",  "Pulmonary adenocarcinoma" = "#009E73" , "Pulmonary carcinoma"="#1E8E99", "Peripheral Nerve Sheath Tumor" ="#CC5800")

#' ## Cell Lines Per Cancer
#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="Bar plot illustrating the distribution of 58 cancer cell lines across 16 different cancer types. Osteosarcoma has the highest number of cell lines (21), highlighting its prevalence in cancer research. This plot can be used to identify cancer types with readily available cell lines for further study.", fig.align = "center"
ggplot(data = Chart, aes(x = reorder(Cancer_Type, +n), y = n, fill = Cancer_Type)) +
  theme_classic() +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5, size = 2.75) +
  labs(x = "Cancer Type",
       y = "Number of Cell Lines") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.title = element_text(hjust = 0.40, size = 9),
        legend.key.size = unit(0.8, "line"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA)) +
  scale_color_manual(values = col,
                     limits = names(col),
                     aesthetics = c("colour", "fill")) +
  guides(fill=guide_legend(ncol = 1,
                           title = "Cancer Type"))

## Correlation ----
#Combine files
gsvaDES <- bind_rows(hm_des_gsva, cp_des_gsva, cgp_des_gsva, cgn_des_gsva, bp_des_gsva, cc_des_gsva, mf_des_gsva, c6_des_gsva, c8_des_gsva)

gsvaRUV <- bind_rows(hm_ruv_gsva, cp_ruv_gsva, cgp_ruv_gsva, cgn_ruv_gsva, bp_ruv_gsva, cc_ruv_gsva, mf_ruv_gsva, c6_ruv_gsva, c8_ruv_gsva)

#Combined correlation on both normalization methods
corrmat_combined <- cor(gsvaDES, gsvaRUV, method = "pearson")
#' ## Correlation
#' ### Correlation Plot 
#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="Heatmap visualization of the correlation coefficients between RUVg and DESeq2 normalization methods for gene set variation analysis (GSVA) data. Each cell in the matrix represents the correlation coefficient for a pair of cancer cell lines, with RUVg normalization on the y-axis and DESeq2 normalization on the x-axis. The color intensity reflects the strength and direction of the correlation (typically red for positive, blue for negative). This heatmap allows for a detailed examination of the relationships between the two normalization methods across all 58 cancer cell lines, potentially revealing systematic biases or trends depending on the cancer type.", fig.height=8, fig.align = "center"
#Create correlation plot using corrplot on the correlation matrix with both normalization methods
corrplot(corrmat_combined, method = "color", tl.cex = 0.5, mar = c(1, 1, 1, 1)) 
mtext("DESeq2", side = 3, line = -2, cex = 1)
mtext("RUVseq", side = 2, line = 3, cex = 1)

#' ### Histogram 
#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="Histogram visualization comparing RUVg and DESeq2 normalization methods for gene set variation analysis (GSVA) data. The x-axis represents the correlation coefficient values, indicating the strength and direction of the association between the two normalization methods. The y-axis shows the frequency of each correlation value across all 58 cancer cell lines. This histogram allows for a comprehensive analysis of the distribution of correlation coefficients, highlighting potential biases or trends favoring one normalization method over the other for GSVA analysis in this cancer cell line dataset.", fig.height=8, fig.align = "center"
hist(corrmat_combined)

#Pull values of matching cell lines from correlation matrix 
corr_value <- diag(corrmat_combined)
corr.df <- data.frame(corr_value)
corr.df <- data.frame("corr_value" = round(corr_value, digits = 2))

#Merge cancer types into working data frame
corr.df <- merge(meta, corr.df, by.x="Cell_Line_ID", by.y = 0) 

# Colors for cancer type 
hcol = list(Histology = c("Hemangiosarcoma"= "#A6CEE3","Hepatocellular Carcinoma"= "#1F78B4","Histiocytic Sarcoma"= "#B2DF8A" ,"Leukemia / Lymphoma"= "#33A02C","Mammary Carcinoma"=  "#FB9A99","Mast Cell Tumor"= "#E31A1C","Melanoma"= "#FDBF6F","Meningioma"= "#FF7F00" ,"Osteosarcoma"= "#CAB2D6","Soft Tissue Sarcoma"=  "#6A3D9A","Squamous cell carcinoma"= "cyan2","Thyroid Carcinoma"= "#B15928" ,"Transitional Cell Carcinoma" ="gray60",  "Pulmonary adenocarcinoma" = "#009E73" , "Pulmonary carcinoma"="#1E8E99", "Peripheral Nerve Sheath Tumor" ="#CC5800"))

# Define PCA colors 
colors <- hcol$Histology

#' ###  Correlation Bar Plot
#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="Bar plot visualizing the correlation coefficients between RUVseq and DESeq normalization methods for gene set variation analysis (GSVA) on 58 cancer cell lines. Each bar is colored according to the cancer type it represents, allowing for quick identification of trends between normalization methods and specific cancer types. This visualization facilitates the comparison of normalization techniques, highlighting which method yields more consistent results (higher correlation) for GSVA analysis within different cancer types.", fig.align="center"

ggplot(data = corr.df, aes(x = reorder(Cell_Line_ID, corr_value), y = corr_value, fill = Cancer_Type)) +
  theme_classic() +
  geom_bar(stat = "identity") +
  geom_text(aes(label = corr_value), angle = 90, vjust = 0.5, hjust = 1.2, size = 3, color = "black") +
  labs(x = "Cell Line",
       y = "Correlation Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2, size = 6),
        axis.text.y = element_text(size = 7)) +
  theme(legend.position = "bottom",  
        legend.direction = "horizontal",  
        legend.box = "horizontal",
        legend.text = element_text(size = 6),
        legend.box.margin=margin(-2,-2,-2,-2),
        legend.key.size = unit(0.5, "line"),
        plot.margin = margin(1,1,1,1)) + #Change to zero and add the following code: 
  #panel.background = element_rect(fill = "transparent", colour = NA),
  #plot.background = element_rect(fill = "transparent", colour = NA),
  #legend.background = element_rect(fill = "transparent", colour = NA),
  #legend.box.background = element_rect(fill = "transparent", colour = NA)
  scale_fill_manual(values = colors) +
  guides(fill = guide_legend(ncol = 4,
                             title.position = "top"))

#+ echo=FALSE, warning=FALSE, message=FALSE
#Remove Gene names
des_norm2 <- des_norm
des_norm2$Hugo_Symbol <- NULL
## Stats ---- 
### DESeq2----
#' ## DESeq2
#' ### PCA Plot 
#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="PCA plot showing the expression values of the DESeq normalized data organized by the cancer type.", fig.height=10, fig.width=9, fig.align = "center"
#Plot PCA 
plotPCA(as.matrix(des_norm2), col=colors[x], cex=1.2, labels=FALSE, pch=c(19, 17, 15)[b])
legend("topleft", legend= unique(x), fill = unique(colors[x]), xpd=NA, bty="n" , cex = 1)
legend("topright", legend= c("Batch 1", "Batch 2", "Batch 3"), pch=c(19, 17, 15), xpd=NA, bty="n" , cex = 1)

#' ### Mean Box Plot
#+ echo=FALSE, warning=FALSE, message=FALSE
#' \mbox{}
#+ echo=FALSE, warning=FALSE, message=FALSE

#### Mean ----
#Calculate mean 
df.mean <- des_norm2 %>%   
  summarise_all(mean)
#Tranpose data
df.mean <- data.frame("Mean_expression" = t(df.mean)) 

#Merge values with cancer type 
df.mean <- merge(meta, df.mean, by.x="Cell_Line_ID", by.y = 0)  
df.mean <-  df.mean[order(df.mean$Cancer_Type),]  

#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="Box plot of the mean expression values calculated using the DESeq normalized data with a jitter. The jitter also illustrates the batch numbers that mean values were calculated from. Entire plot is organized by the cancer type.", fig.height=7.5, fig.align = "center"
ggplot(df.mean, aes(x = Mean_expression, y = Cancer_Type)) +  
  theme_classic() +
  geom_boxplot() + 
  labs(x= "Mean Expression", y= "Cancer Type") + 
  geom_jitter(shape=19, position=position_jitter(0.2),
              aes(color = factor(Batch))) +
  scale_color_manual(name = "Batch Number",
                     values = c("#d55e00", "#56b4e9", "#009e73"))


#' \newpage
#+ echo=FALSE, warning=FALSE, message=FALSE
#### Median ----
#Calculate Median 
df.median <- des_norm2 %>%   
  summarise_all(median) 
#Tranpose data
df.median <- data.frame("Median_expression" = t(df.median)) 

#Merge values with Cancer Type 
df.median <- merge(meta, df.median, by.x="Cell_Line_ID", by.y = 0)  
df.median <-  df.median[order(df.median$Cancer_Type),]

#' \newpage
#' ### Median Box Plot  
#' \mbox{}
#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="Box plot of the median expression values calculated using the DESeq normalized data with a jitter. The jitter also illustrates the batch numbers that median values were calculated from. Entire plot is organized by the cancer type.", fig.align = "center", fig.height=7.5
ggplot(df.median, aes(x=Median_expression, y=Cancer_Type)) +  
  theme_classic() +
  geom_boxplot() +   
  labs(x= "Median Expression", y= "Cancer Type") +  
  geom_jitter(shape=19, position = position_jitter(0.2),
              aes(color = factor(Batch))) +
  scale_color_manual(name = "Batch Number",
                     values = c("#d55e00", "#56b4e9", "#009e73"))

#' \newpage
#+ echo=FALSE, warning=FALSE, message=FALSE
#### Max ----
#Calculate max for cell lines
df.max1 <- des_norm %>% 
  summarise_all(max)
#Remove Hugo_Symbol column
df.max1$Hugo_Symbol <- NULL
#Transpose data
df.max1 <- data.frame("Max_expression" = t(df.max1))
#Retrieve row names of max values 
df.max <- des_norm2 %>% 
  summarise_all(~ row.names(des_norm2)[which.max(.)]) %>% 
  pivot_longer(everything())
#Add Gene names to working dataframe  
df.max1$Gene_Name <- df.max$value
#Create data frame of row names and Hugo_Symbols from master data set 
id_sym <-  data.frame(row.names(des_norm))
id_sym$sym <- des_norm$Hugo_Symbol
id_sym$value <- id_sym$row.names.des_norm.
id_sym$row.names.des_norm. <- NULL
#Create reference variable to sort id_sym
ref <- df.max$value
#Filter out only necessary variables into new reference dataframe 
df.max3 <- filter(id_sym, value %in% ref)
#Add Gene IDs to max value data frame
df.max1$Gene_ID <- df.max3$sym[match(df.max1$Gene_Name, df.max3$value)]
#Merge and order by Cancer types
df.max1 <- merge(meta, df.max1, by.x="Cell_Line_ID", by.y = 0) 
df.max1 <-  df.max1[order(df.max1$Cancer_Type),]
#Remove unnecessary columns
df.max1$Batch <- NULL
df.max1$Gene_Name <- NULL
#Reorder columns so that Gene_ID is before Max value
df.max1 <- df.max1[,c(1, 2, 4, 3)]
#Generate bar graph showing number of occurrences of genes with max expression
#Generate chart of max count grouped by cancer type
max_count <- df.max1 %>% 
  count(Gene_ID, Cancer_Type)
group_total <- max_count %>% 
  mutate(sum(n))

max_count <- max_count %>%
  mutate(Gene_ID = factor(Gene_ID, levels = unique(Gene_ID[order(n, Gene_ID)])))

max_count <- max_count %>%
  group_by(Cancer_Type, Gene_ID) %>%
  summarise(n = sum(n)) %>%
  ungroup()

#' ### Max Gene Expression
#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="This bar graph visualizes the distribution of genes with the highest expression level across the 58 cell lines. While nearly 20,000 genes were analyzed, only 15 emerged as the most highly expressed gene in at least one cell line. Each bar represents a gene, and its color indicates the specific cancer(s) where that gene showed the highest expression. Refer to Appendix Table 19 for detailed descriptions of these 15 genes."
ggplot(data = max_count, aes(x = reorder(Gene_ID, n, FUN = sum), y = n, fill = Cancer_Type, group = interaction(Gene_ID, Cancer_Type))) +
  theme_classic() +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Gene ID",
       y = "Number of Occurrences",
       title = "Number Occurrences of Genes in Max Value Grouped by Cancer Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.text = element_text(size = 7),
        legend.title = element_text(hjust = 0.40, size = 9),
        legend.key.size = unit(0.8, "line")) +
  scale_color_manual(values = col,
                     limits = names(col),
                     aesthetics = c("colour", "fill")) +
  geom_text(aes(y = n, label = n), stat = "identity", position = position_stack(vjust = 0.5), size = 2.5, color = "black") +
  guides(fill=guide_legend(ncol = 1,
                           title = "Cancer Type"))

#+ echo=FALSE, warning=FALSE, message=FALSE
#Remove Gene names
ruv_norm2 <- ruv_norm
ruv_norm2$Hugo_Symbol <- NULL
### RUVg----
#' ## RUVg
#' ### PCA Plot 
#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="PCA plot showing the expression values of the RUVg normalized data organized by the cancer type.", fig.height=10, fig.width=9, fig.align = "center"
#Plot PCA 
plotPCA(as.matrix(ruv_norm2), col=colors[x], cex=1.2, labels=FALSE, pch=c(19, 17, 15)[b])
legend("topright", legend= unique(x), fill = unique(colors[x]), xpd=NA, bty="n" , cex = 1)
legend("topright", legend= c("Batch 1", "Batch 2", "Batch 3"), pch=c(19, 17, 15), xpd=NA, bty="n" , cex = 1)
#+ echo=FALSE, warning=FALSE, message=FALSE
#### Mean ----
#Calculate mean 
df.mean <- ruv_norm2 %>%   
  summarise_all(mean)
#Tranpose data
df.mean <- data.frame("Mean_expression" = t(df.mean)) 

#Merge values with cancer type 
df.mean <- merge(meta, df.mean, by.x="Cell_Line_ID", by.y = 0)  
df.mean <-  df.mean[order(df.mean$Cancer_Type),]  

#' ### Mean Box Plot 
#' \mbox{} 
#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="Box plot of the mean expression values calculated using the RUVg normalized data with a jitter. The jitter also illustrates the batch numbers that mean values were calculated from. Entire plot is organized by the cancer type.", fig.height=7.5, fig.align = "center" 
ggplot(df.mean, aes(x = Mean_expression, y = Cancer_Type)) +  
  theme_classic() +
  geom_boxplot() + 
  labs(x= "Mean Expression", y= "Cancer Type") + 
  geom_jitter(shape=19, position=position_jitter(0.2),
              aes(color = factor(Batch))) +
  scale_color_manual(name = "Batch Number",
                     values = c("#d55e00", "#56b4e9", "#009e73"))

#### Median ----
#Calculate Median 
df.median <- des_norm2 %>%   
  summarise_all(median) 
#Tranpose data
df.median <- data.frame("Median_expression" = t(df.median)) 


#Merge values with Cancer Type 
df.median <- merge(meta, df.median, by.x="Cell_Line_ID", by.y = 0)  
df.median <-  df.median[order(df.median$Cancer_Type),]

#' ### Median Box Plot  
#' \mbox{}
#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="Box plot of the median expression values calculated using the DESeq normalized data with a jitter. The jitter also illustrates the batch numbers that median values were calculated from. Entire plot is organized by the cancer type.", fig.height=7.5, fig.align = "center" 
ggplot(df.median, aes(x=Median_expression, y=Cancer_Type)) +  
  theme_classic() +
  geom_boxplot() +   
  labs(x= "Median Expression", y= "Cancer Type") +  
  geom_jitter(shape=19, position = position_jitter(0.2),
              aes(color = factor(Batch))) +
  scale_color_manual(name = "Batch Number",
                     values = c("#d55e00", "#56b4e9", "#009e73"))

#### Max ----
#Calculate max for cell lines
df.max1 <- ruv_norm %>% 
  summarise_all(max)
#Remove Hugo_Symbol column
df.max1$Hugo_Symbol <- NULL
#Transpose data
df.max1 <- data.frame("Max_expression" = t(df.max1))
#Retrieve row names of max values 
df.max <- ruv_norm2 %>% 
  summarise_all(~ row.names(ruv_norm2)[which.max(.)]) %>% 
  pivot_longer(everything())
#Add Gene names to working dataframe  
df.max1$Gene_Name <- df.max$value
#Create data frame of row names and Hugo_Symbols from master data set 
id_sym <-  data.frame(row.names(des_norm))
id_sym$sym <- des_norm$Hugo_Symbol
id_sym$value <- id_sym$row.names.des_norm.
id_sym$row.names.des_norm. <- NULL
#Create reference variable to sort id_sym
ref <- df.max$value
#Filter out only necessary variables into new reference dataframe 
df.max3 <- filter(id_sym, value %in% ref)
#Add Gene IDs to max value data frame
df.max1$Gene_ID <- df.max3$sym[match(df.max1$Gene_Name, df.max3$value)]
#Merge and order by Cancer types
df.max1 <- merge(meta, df.max1, by.x="Cell_Line_ID", by.y = 0) 
df.max1 <-  df.max1[order(df.max1$Cancer_Type),]
#Remove unnecessary columns
df.max1$Batch <- NULL
df.max1$Gene_Name <- NULL

#Generate chart of max count grouped by cancer type
max_count <- df.max1 %>% 
  count(Gene_ID, Cancer_Type)
group_total <- max_count %>% 
  mutate(sum(n))

max_count <- max_count %>%
  mutate(Gene_ID = factor(Gene_ID, levels = unique(Gene_ID[order(n, Gene_ID)])))

max_count <- max_count %>%
  group_by(Cancer_Type, Gene_ID) %>%
  summarise(n = sum(n)) %>%
  ungroup()

#' ### Max Gene Expression
#+ echo=FALSE, warning=FALSE, message=FALSE, fig.cap="This bar graph visualizes the distribution of genes with the highest expression level across the 58 cell lines. While nearly 20,000 genes were analyzed, only 18 emerged as the most highly expressed gene in at least one cell line. Each bar represents a gene, and its color indicates the specific cancer(s) where that gene showed the highest expression. Refer to Appendix Table 20 for detailed descriptions of these 18 genes."
ggplot(data = max_count, aes(x = reorder(Gene_ID, n, FUN = sum), y = n, fill = Cancer_Type, group = interaction(Gene_ID, Cancer_Type))) +
  theme_classic() +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Gene ID",
       y = "Number of Occurrences",
       title = "Number Occurrences of Genes in Max Value Grouped by Cancer Type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.text = element_text(size = 7),
        legend.title = element_text(hjust = 0.40, size = 9),
        legend.key.size = unit(0.8, "line")) +
  scale_color_manual(values = col,
                     limits = names(col),
                     aesthetics = c("colour", "fill")) +
  geom_text(aes(y = n, label = n), stat = "identity", position = position_stack(vjust = 0.5), size = 2.5, color = "black") +
  guides(fill=guide_legend(ncol = 1,
                           title = "Cancer Type"))