


# Imports
library(DESeq2)

# Constants
RAW_COUNTS_FOLDER <- 'C:\\Users\\Ana\\Documents\\Estágio\\Aldinger2020\\'

# Functions
import_featureCounts <- function(counts_file, clinical_data, te_class = FALSE) {
  
  counts_df <- read.table(counts_file, skip = 1, header = T)
  total_cols = length(colnames(counts_df))
  colnames(counts_df)[7:total_cols] <- sapply(strsplit(colnames(counts_df)[-c(1:6)], "\\."), '[', 9)
  colnames(counts_df)[7:total_cols] <- sub("/mnt/comics-data/comics02/RTT_files/Aldinger2020/03_Uniquelly_Mapped/", "", colnames(counts_df)[7:total_cols])
  colnames(counts_df)[7:total_cols] <- sub("_Aligned", "", colnames(counts_df)[7:total_cols])
  counts_df <- counts_df[, -c(2:6)] #Fica só os gene ids e as counts
  srr_names <- colnames(counts_df)[-1] # Todos os nomes das colunas exceto os gene ids
  
  matches <- match(srr_names, clinical_data$Run)
  
  counts_df <- counts_df[, c(1, which(!is.na(matches)) + 1)] # Incluir apenas os que não são na
  
  srr_names_fixed <- colnames(counts_df)[-1]
  colnames(counts_df)[-1] <- clinical_data[["Sample.Name"]][match(srr_names_fixed, clinical_data$Run)] #Colocar na tabela final, todos os srr presentes na clinical data
  
  if (te_class == TRUE) {
    counts_df <- counts_df[startsWith(counts_df$Geneid, "TE_"), ] #Se for uma TEClass, vai apenas buscar os que começam por TE
  }
  
  rownames(counts_df) <- counts_df$Geneid
  counts_df$Geneid <- NULL
  return(counts_df)
  
}

# Ver o que faz o save_expression_name
dea_deseq <- function(counts_df, clinical_info, save_expression = TRUE, save_expression_name = "Expression.csv") {
  
  counts_df <- counts_df[, clinical_info[["Sample.Name"]]] #Ter a certeza que os sample names no sample info são iguais aos do counts df
  dds <- DESeqDataSetFromMatrix(countData = counts_df, colData = clinical_info, design = ~ Condition )
  
  keep <- rowSums(counts(dds)) > 10 # Certifica-se que temos um bom número de counts para determinado gene!
  dds <- dds[keep,] # O keep é um vetor com valores TRUE ou FALSE para saber se todos os genes têm acima de 10 reads! Só entram no dds se forem true!
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("Condition", "RTT", "Ctrl"), alpha = 0.05)
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

############ MAIN PROCESS #################

clinical_data <- read.csv(paste0(RAW_COUNTS_FOLDER, 'SraRunTable_Aldinger20.csv'))
clinical_data <- clinical_data[clinical_data$source_cortex_name == 'Cingluate cortex',]
clinical_data <- clinical_data[c('Run', 'Sample.Name', 'Condition')]


# Process individual TE counts
te_counts <- import_featureCounts(paste0(RAW_COUNTS_FOLDER, 'Counts_TEIndividual.txt'), clinical_data)
res_dea <- dea_deseq(te_counts, clinical_data, save_expression = TRUE, 'Expression_TE_DESeq_Cingluate_cortex.csv')
#Saves the result, porque é que guardo aqui novamente?
write.csv(res_dea, file = paste0(RAW_COUNTS_FOLDER, 'DEA_TE_DESeq_Cingluate_cortex.csv'), quote = FALSE)


# Process gene counts
gene_counts <- import_featureCounts(paste0(RAW_COUNTS_FOLDER, 'Counts_Genes_Unique.txt'), clinical_data)
res_dea <- dea_deseq(gene_counts, clinical_data, save_expression = TRUE, "Expression_Genes_DESeq_Cingluate_cortex.csv")
write.csv(res_dea, file = paste0(RAW_COUNTS_FOLDER, 'DEA_Genes_DESeq_Cingluate_cortex.csv'), quote = FALSE)


# Process subfamily TE counts
te_class_counts <- import_featureCounts(paste0(RAW_COUNTS_FOLDER, 'Counts_TEClass.txt'), clinical_data, te_class = TRUE)
res_dea <- dea_deseq(te_class_counts, clinical_data, save_expression = TRUE, 'Expression_TEClass_DESeq_Cingluate_cortex.csv')
write.csv(res_dea, file = paste0(RAW_COUNTS_FOLDER, 'DEA_TEClass_DESeq_Cingluate_cortex.csv'), quote = FALSE)









