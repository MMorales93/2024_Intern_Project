# Script to run several simple correlations in a lop and pull out significant correlations

#Load Data ----
#Load expression or gsva data
hall <- read.csv("input/GSVAscores_Osteosarcoma_Input_RUVseqNormalizedCounts_Hallmark.csv", row.names = 1)
#Load drug data
drug <- read.csv("input/Alisertib.csv", row.names = 1)
#Or of you have a data frame with all data combined just load that
all_data <- read.csv("Olaparib_All_Data.csv", row.names = 1)

#Create a data frame that p-values will populate 
rows <- colnames(all_data)[5:783] #Use genes or GeneSet names as row names
cols <- colnames(drug) #Drug sensitvity data will be column names 
#Could also use
cols <- colnames(all_data)[1:4]

check_list <- data.frame(matrix(nrow = 779, ncol = 4))
row.names(check_list) <- rows
colnames(check_list) <- cols

#For Loop to run simple correlations
for(i in 5:ncol(all_data)){#Make sure to start correlations with genes and/or GeneSets and not drug sensitvity parameters
  
  column <- names(all_data[i])
  
  #GR50.CRE
  avz <- cor.test(all_data$GR50.CRE, all_data[, i], method = "p")
  check_list[column, "GR50.CRE"] <- avz$p.value
  
  #Res.vs.Sens
  avz <- wilcox.test(all_data[,i] ~ Res.vs.Sens, data = all_data, exact = F) 
  check_list[column, "Res.vs.Sens"] <- avz$p.value[1]
  
  #GR.RESPONSE
  avz <- cor.test(all_data$GR.Response, all_data[, i], method = "p")
  check_list[column, "GR.Response"] <- avz$p.value
  
  #Fraction.Dead.Cells
  avz <- cor.test(all_data$Fraction.Dead.Cells, all_data[, i], method = "p")
  check_list[column, "Fraction.Dead.Cells"] <- avz$p.value
}

#Filter out significant genes and/or GeneSets
#Filters out rows so long as 1 of the drug sensitvity parameters is significant p<0.05
filtered_df <- check_list[apply(check_list, 1, function(x) any(x < 0.05)), ]
#Filters out rows so long as all 4 of the drug sensitvity parameters is significant p<0.05
filtered_df_1 <- check_list[apply(check_list, 1, function(x) all(x < 0.05)), ]