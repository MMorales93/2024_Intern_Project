# Gene Descriptions

# This script uses gene symbols as the input and pulls gene descriptions from MyGene.info and places them in a table that can be rendered into a pdf with other scripts. I use this code most when I need to create a quick reference for genes of interest

#Libraries needed for script
library(gprofiler2)
library(mygene)
library(flextable)
library(tidyr)

#Generate Descriptions

deg_genes_up <- deg_genes[deg_genes$sig=='up',] #Subset line by condition in this case significance or you can create a list of gene names 'genes <- c("TMEM131L", "APC", "DNAJB1", "LARP4B", "CSF1")'

#Convert gene symbols to Entrez Gene ID
deg_up_info <- gconvert(deg_genes_up$Hugo_Symbol, organism = "hsapiens", target = "ENTREZGENE_ACC")
#Get gene descriptions 
deg_up_sum <- getGenes(deg_up_info$target, fields = "summary")
#Pull gene symbol and gene summary into seperate dataframe to be placed into a table 
deg_up_df <- data.frame(Hugo_Symbol = deg_up_info$input,
                        Gene_Summary = deg_up_sum@listData[["summary"]])

#Display descriptions in table
flextable(deg_up_df) %>%
  fontsize(part = "body", size = 6) %>%
  fontsize(part = "header", size = 8) %>%
  set_table_properties(layout ="fixed", align = "center") %>%
  flextable::width(width = 5, j = 2, unit = "in") %>%
  align(part = "header", align = "center") %>%
  hline( part = "body")