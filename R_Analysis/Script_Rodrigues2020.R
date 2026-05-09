


# Imports
library(DESeq2)
library('dplyr')

# Constants
RAW_COUNTS_FOLDER <- 'C:\\Users\\Ana\\Documents\\Estágio\\Rodrigues2020\\'

# Functions
import_featureCounts <- function(counts_file, clinical_data, te_class = FALSE) {
  
  counts_df <- read.table(counts_file, skip = 1, header = T)
  total_cols = length(colnames(counts_df))
  colnames(counts_df)[7:total_cols] <- sapply(strsplit(colnames(counts_df)[-c(1:6)], "\\."), '[', 9)
  colnames(counts_df)[7:total_cols] <- sub("_Aligned", "", colnames(counts_df)[7:total_cols])
  counts_df <- counts_df[, -c(2:6)] #Fica só os gene ids e as counts
  
  srr_names <- colnames(counts_df)[-1] # Todos os nomes das colunas exceto os gene ids
  colnames(counts_df)[-1] <- clinical_data[["Sample_Name"]][match(srr_names, clinical_data$Sample_Name)] #Colocar na tabela final, todos os srr presentes na clinical data
  print(colnames(counts_df))
  
  if (te_class == TRUE) {
    counts_df <- counts_df[startsWith(counts_df$Geneid, "TE_"), ] #Se for uma TEClass, vai apenas buscar os que começam por TE
  }
  
  rownames(counts_df) <- counts_df$Geneid
  counts_df$Geneid <- NULL
  
  return(counts_df)
  
}

# Ver o que faz o save_expression_name
dea_deseq <- function(counts_df, clinical_info, save_expression = TRUE, save_expression_name = "Expression.csv") {
  print(clinical_info)
  counts_df <- counts_df[, clinical_info[["Sample_Name"]]] #Ter a certeza que os sample names no sample info são iguais aos do counts df

  dds <- DESeqDataSetFromMatrix(countData = counts_df, colData = clinical_info, design = ~ Condition )
  keep <- rowSums(counts(dds)) > 10 # Certifica-se que temos um bom número de counts para determinado gene!
  dds <- dds[keep,] # O keep é um vetor com valores TRUE ou FALSE para saber se todos os genes têm acima de 10 reads! Só entram no dds se forem true!
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("Condition", "RTT", "WT"), alpha = 0.05)
  print(summary(res))
  
  res <- res[order(res$padj),]
  
  if (save_expression) {
    print(colnames(counts_df))
    write.csv(counts(dds, normalize = F),
              paste0(RAW_COUNTS_FOLDER, save_expression_name), #Cria o nome do ficheiro final
              quote = FALSE)
  }
  
  return(as.data.frame(res)) #Retorna o resultado numa data frame
  
  
}

collapse_duplicate_lanes <- function(clinical_data) { # Forma uma tibble é basically uma data frame
  
  clinical_data <- clinical_data %>%
    group_by(Sample_Name) %>%
    
    summarise (
      Run = paste(unique(Run), collapse = '/'),
      Condition = first(Condition),
      Cells = first(Cells),
      Cells_Abv = first(Cells_Abv),
      Sample_Description = paste(unique(Sample_Description), collapse = '/')
    ) %>%
    
    ungroup()
  
  
  
  return(as.data.frame(clinical_data))
    
}

############ MAIN PROCESS #################

clinical_data <- read.csv(paste0(RAW_COUNTS_FOLDER, 'Clinical_Data_Rodrigues20.csv'), header = TRUE, sep = ";")

clinical_data <- clinical_data[clinical_data$Cells_Abv == 'Neu',]
clinical_data <- collapse_duplicate_lanes(clinical_data)
clinical_data <- clinical_data[c('Run', 'Sample_Name', 'Condition')]

print(clinical_data)

# Process individual TE counts
te_counts <- import_featureCounts(paste0(RAW_COUNTS_FOLDER, 'Counts_TEIndividual.txt'), clinical_data)
res_dea <- dea_deseq(te_counts, clinical_data, save_expression = TRUE, 'Expression_TE_DESeq_Neu.csv')
#Saves the result, porque é que guardo aqui novamente?
write.csv(res_dea, file = paste0(RAW_COUNTS_FOLDER, 'DEA_TE_DESeq_Neu.csv'), quote = FALSE)


# Process gene counts
gene_counts <- import_featureCounts(paste0(RAW_COUNTS_FOLDER, 'Counts_Genes_Unique.txt'), clinical_data)
res_dea <- dea_deseq(gene_counts, clinical_data, save_expression = TRUE, "Expression_Genes_DESeq_Neu.csv")
write.csv(res_dea, file = paste0(RAW_COUNTS_FOLDER, 'DEA_Genes_DESeq_Neu.csv'), quote = FALSE)


# Process subfamily TE counts
te_class_counts <- import_featureCounts(paste0(RAW_COUNTS_FOLDER, 'Counts_TEClass.txt'), clinical_data, te_class = TRUE)
res_dea <- dea_deseq(te_class_counts, clinical_data, save_expression = TRUE, 'Expression_TEClass_DESeq_Neu.csv')
write.csv(res_dea, file = paste0(RAW_COUNTS_FOLDER, 'DEA_TEClass_DESeq_Neu.csv'), quote = FALSE)









