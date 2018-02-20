#Keep tracking the running time of this program
starttime <- proc.time()
library(data.table)
library(googlesheets)
library(stringr)
library(tidyverse)
library(ggpubr)
library(plyr)
library(glmnet)
#source code contains function to read googlesheet named as master and convert Sample ID
source('~/shenz/gray_brca/r/basics.R')

#Read the maf file containing all mutations profile of each sample (BRCA mutations)
maf = fread('~/shenz/gray_brca/joe/final-filtered-facets.vep.maf')

#Read the maf file containing all mutations profile of each (BRCA WT control)
wt_maf <- fread('~/shenz/gray_brca/joe/Proj_93017_data_mutations_extended_010318.vep.maf')

#Read the decomposed file processed from maf to obtain all mutational signatures
decomp = fread('~/shenz/gray_brca/joe/sigs.out')

#Read the file contains the total number of mutations
sample_n =fread('~/shenz/gray_brca/joe/spectrum.txt')

#Read the maf file that contains annotated germline and somatic variants information
gml <- fread('~/shenz/gray_brca/joe/gml-som-variants-facets-annotated.txt')

gml$sample <- sapply(gml$sample, sample_conv)

master$CMO_ID <- sapply(master$CMO_ID, sample_conv)

master <- master %>% filter(is.na(Sample_QC))
  
subset_master_gml <- master %>% filter(master$CMO_ID %in% gml$sample ) 

subset_master_non_gml <- master %>% filter(!(master$CMO_ID %in% gml$sample)) 

germline_wt_cohort <- subset_master_gml %>% filter(`12_245_PartC_Consented` == 'YES' & (BRCA_Type == 'Somatic' | BRCA_Type == 'Somatic,Somatic'))

germline_unknown_cohort <- subset_master_gml %>% filter(`12_245_PartC_Consented` == 'NO'& (BRCA_Type == 'Somatic' | BRCA_Type == 'Somatic,Somatic'))

germline_mutation_cohort <-gml %>% filter(BRCA_Type == 'Germline' | BRCA_Type == 'Germline,Somatic')

germline_wt_dataframe <- germline_wt_cohort  %>% select('sample' = CMO_ID) %>% mutate(Germline_mutation = 'WT') 

germline_unknown_dataframe <- germline_unknown_cohort %>% select('sample' = CMO_ID)%>% mutate(Germline_mutation = 'Unknown') 

germline_mutation_dataframe <- germline_mutation_cohort %>% 
  select(sample) %>% 
  mutate(Germline_mutation = paste(germline_mutation_cohort$HGVSp_Short, germline_mutation_cohort$Hugo_Symbol, sep = ' ')) 

tot_germline_dataframe <- rbind(germline_wt_dataframe, germline_unknown_dataframe, germline_mutation_dataframe) %>%
  group_by(sample)  %>% unique() %>% transmute(Germline_mutation = paste(c(Germline_mutation), collapse = ',')) %>% unique() %>% arrange(sample)

somatic_mutation_search <- function(sample_ID) {
  maf[(grepl('BRCA1', maf$Hugo_Symbol) | grepl('BRCA2', maf$Hugo_Symbol)) & grepl(sample_ID, maf$Tumor_Sample_Barcode), ]
}

somatic_mutation_cohort_maf <- lapply(gml$sample, somatic_mutation_search) %>% 
  do.call(rbind, .)  

somatic_mutation_dataframe_maf <- somatic_mutation_cohort_maf %>% select('sample' = Tumor_Sample_Barcode) %>% 
  mutate(Somatic_mutation = paste(somatic_mutation_cohort_maf$HGVSp_Short, somatic_mutation_cohort_maf$Hugo_Symbol, sep = ' '))

somatic_mutation_cohort_gml <- gml %>% filter(BRCA_Type == 'Somatic' | BRCA_Type == 'Germline,Somatic')

somatic_mutation_dataframe_gml <- somatic_mutation_cohort_gml %>% select(sample) %>%
  mutate(Somatic_mutation = paste(somatic_mutation_cohort_gml$HGVSp_Short, somatic_mutation_cohort_gml$Hugo_Symbol, sep = ' '))

combined_somatic_mutation_dataframe <- rbind(somatic_mutation_dataframe_gml,somatic_mutation_dataframe_maf) %>%
  group_by(sample) %>% unique() %>% transmute(Somatic_mutation = paste(c(Somatic_mutation), collapse = ',')) %>% unique() %>% arrange(sample)

somatic_zygosity_subset <- gml %>% select(sample,Tumor_zygosity = zygosity_flag) 

merged_dataset_gml <- merge.data.frame(tot_germline_dataframe, combined_somatic_mutation_dataframe,by.x = 'sample', all= TRUE) %>%
  merge.data.frame(somatic_zygosity_subset, by.x = 'sample', all = TRUE)

somatic_hom_deletion_subset <- subset_master_non_gml %>% filter(BRCA_Type == 'Somatic' & Variant_Type == 'deletion') 

somatic_hom_deletion_df <- somatic_hom_deletion_subset %>% select('sample' = CMO_ID) %>%  
  mutate(Germline_mutation = NA) %>% 
  mutate(Somatic_mutation = paste(somatic_hom_deletion_subset$HGVSp_Short, somatic_hom_deletion_subset$Hugo_Symbol, sep = ' ')) %>%
  mutate(Tumor_zygosity = 'Homozygous Deletion')
         
WT_tumor_df <- subset_master_non_gml %>% filter(BRCA_Type == 'WT') %>% select('sample' = CMO_ID) %>% 
  mutate(Germline_mutation = 'WT', Somatic_mutation = 'WT', Tumor_zygosity = 'None')

germline_deletion_subset <- subset_master_non_gml %>% filter(BRCA_Type == 'Germline' & Variant_Type == 'deletion')

germline_deletion_somatic_mut <- germline_deletion_subset$CMO_ID %>% lapply(somatic_mutation_search) 

for (index in 1 : length(germline_deletion_somatic_mut)) {
  if (!is_empty(germline_deletion_somatic_mut[[index]]$Tumor_Sample_Barcode)) {
  names(germline_deletion_somatic_mut)[index] <- germline_deletion_somatic_mut[[index]]$Tumor_Sample_Barcode 
}
}

match_index <- match(germline_deletion_subset$CMO_ID, names(germline_deletion_somatic_mut))

germline_deletion_df <- germline_deletion_subset %>% select('sample' = CMO_ID) %>% 
  mutate(Germline_mutation = paste(germline_deletion_subset$HGVSp_Short, germline_deletion_subset$Hugo_Symbol, sep = ' ' )) %>%
  mutate(Somatic_mutation = ifelse (!is.na(match_index), paste(germline_deletion_somatic_mut[[match_index[!is.na(match_index)]]]$HGVSc, germline_deletion_somatic_mut[[match_index[!is.na(match_index)]]]$Hugo_Symbol, sep = ' '), NA) )%>%
  mutate(Tumor_zygosity = 'Germline Deletion Rare Case')

merged_dataset_non_gml <- rbind(somatic_hom_deletion_df, WT_tumor_df, germline_deletion_df)

Final_merged_table <- unique(rbind(merged_dataset_gml,merged_dataset_non_gml))

print('Finish Organizing Table')
#----------------------------------------------------------------------------------------------------------------------
#Part1.5 clean the BRCA WT tumor samples and add to the Final_merged_table

#clear BRCA mutant samples from BRCA wt maf file
avoid_vect <- c('P-0020647-T01-WES','P-0021156-T01-WES','P-0021154-T01-WES', 'P-0009423-T02-WES',
                'P-0021254-T01-WES','P-0021100-T01-WES','P-0020719-T01-WES','P-0020399-T01-WES',
                'P-0020831-T01-WES','P-0020791-T01-WES','P-0020841-T01-WES','P-0014885-T02-WES',
                'P-0020633-T01-WES','P-0020639-T01-WES','P-0020673-T01-WES','P-0020668-T01-WES')

wt_maf_cleared <- wt_maf %>% filter(!(Tumor_Sample_Barcode %in% avoid_vect))

#create a BRCA WT samples dataframe
BRCA_wt_df <- data.frame('sample' = unique(wt_maf_cleared$Tumor_Sample_Barcode)) %>% mutate(Germline_mutation = 'WT') %>%
  mutate(Somatic_mutation = 'WT') %>% mutate(Tumor_zygosity = 'None')

#add BRCA WT dataframe to the Final_merged_table
Final_merged_table_with_WT <- rbind(Final_merged_table, BRCA_wt_df)

print('Add WT Sample to table')
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Part2 starts , add a new column to further annotate the tumor based on the somatic and germline mutation pattern according to the two-hits hypothesis
pathogenic_mutation_maf_match <- function(sample_ID, symbol) {
  patho_vect <- c("Frame_Shift_Ins", "Nonsense_Mutation", "Frame_Shift_Del", "Translation_Start_Site", "Nonstop_Mutation", "Splice_Site")
  
  selected_row <- maf %>% filter(Tumor_Sample_Barcode == sample_ID & grepl(symbol,Hugo_Symbol))
  
  variant_vect <- selected_row$Variant_Classification
  number <- 0
  if (length(variant_vect) > 0) {
    for (index in 1: length(variant_vect)){
      if (variant_vect[index] %in% patho_vect == TRUE){
        number = number + 1
      }
    } 
  }
  return (number)
}

#Germline Mutation + LOH (Biallelic Inactivation)
Annotated_merged_table <- Final_merged_table_with_WT %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 0 & Tumor_zygosity == 'AI_LOH_ALT'& grepl('BRCA', Germline_mutation), 'Biallelic Inactivation', 'Others'))

#Germline Deletion
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'Germline Deletion Rare Case', 'Germline Deletion',Tumor_type_based_mutation_pattern))

#Mono allelic Germline mutation 
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'None' & grepl('BRCA', Germline_mutation) & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 0, 'Monoallelic Germline Inactivation', Tumor_type_based_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'AI_LOH_REF' & grepl('BRCA', Germline_mutation) & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 0, 'Monoallelic Germline Inactivation', Tumor_type_based_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'AI_LOH_REF' & grepl('BRCA', Germline_mutation) & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 1, 'Biallelic Inactivation', Tumor_type_based_mutation_pattern))

# Wildtype Tumor
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'None' & Germline_mutation == 'WT' & Somatic_mutation == 'WT', 'BRCA Wildtype Tumor',Tumor_type_based_mutation_pattern))


#Monoallelic Somatic Mutation (according to Maf and annotated table)
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'None' & Germline_mutation == 'WT' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 1 , 'Monoallelic Somatic Inactivation', Tumor_type_based_mutation_pattern))

Annotated_merged_table <-  Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Germline_mutation == 'Unknown' & Tumor_zygosity == 'None'& sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 1, 'Monoallelic Somatic Inactivation', Tumor_type_based_mutation_pattern))


#Monoallelic Somatic Mutation (according to Maf and annotated table) + LOH (Biallelic Inactivation)
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'AI_LOH_ALT'  & Germline_mutation == 'WT' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 1, 'Biallelic Inactivation', Tumor_type_based_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'AI_LOH_ALT'  & Germline_mutation == 'Unknown' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 1, 'Biallelic Inactivation', Tumor_type_based_mutation_pattern))


#Annotate Subclonal Tumor
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'subclonal' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') >= 1, 'Subclonal BRCA Inactivation', Tumor_type_based_mutation_pattern))


#Annotate AI_ALT Tumor
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'AI_ALT' &  Germline_mutation == 'Unknown' &sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 1, 'Monoallelic Somatic Inactivation', Tumor_type_based_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'AI_ALT' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 0 & grepl('BRCA', Germline_mutation), 'Monoallelic Germline Inactivation', Tumor_type_based_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'AI_ALT' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 0 & !grepl('BRCA', Germline_mutation), 'BRCA Wildtype Tumor', Tumor_type_based_mutation_pattern))

#Annotate AI_REF Tumor
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'AI_REF' &  Germline_mutation == 'Unknown' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 1, 'Monoallelic Somatic Inactivation', Tumor_type_based_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'AI_REF' & is.na(Somatic_mutation), 'Monoallelic Germline Inactivation', Tumor_type_based_mutation_pattern))

#annotate somatic homozygous deletion cases
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'Homozygous Deletion' & grepl('BRCA', Somatic_mutation), 'Biallelic Inactivation',Tumor_type_based_mutation_pattern))

#Unusual ALT_LOH cases (with more than 1 pathogenic mutation, G and S both have mutations)
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'AI_LOH_ALT' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') > 1 & Germline_mutation == 'WT', 'Biallelic Inactivation', Tumor_type_based_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'AI_LOH_ALT' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') > 1 & Germline_mutation == 'Unknown', 'Biallelic Inactivation', Tumor_type_based_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'AI_LOH_ALT' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 1 & grepl('BRCA', Germline_mutation), 'Biallelic Inactivation', Tumor_type_based_mutation_pattern))

#Unusual None cases (Germline and Somatic mutations leading to biallelic inactivation)
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_zygosity == 'None' & grepl('BRCA', Germline_mutation)  & sapply(sample, pathogenic_mutation_maf_match, 'BRCA') == 1 , 'Biallelic Inactivation', Tumor_type_based_mutation_pattern))

print('Finish broad classification of mutation type')
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Part 2.5 Based on the classification of of Part 2, further classify the mutation events according
#to BRCA1 or BRCA2

#detail classification of Biallelic Inactivation

#Homodels
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Biallelic Inactivation' & Somatic_mutation =='Deletion BRCA1', 'Homodels BRCA1', 'Others'))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Biallelic Inactivation' & Somatic_mutation =='Deletion BRCA2', 'Homodels BRCA2', Detailed_mutation_pattern))

#Germline + LOH
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Biallelic Inactivation' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA1') == 0  & grepl('BRCA1', Germline_mutation) & Tumor_zygosity == 'AI_LOH_ALT', 'Germline BRCA1 + LOH', Detailed_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Biallelic Inactivation' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA2') == 0 & grepl('BRCA2', Germline_mutation) & Tumor_zygosity == 'AI_LOH_ALT', 'Germline BRCA2 + LOH', Detailed_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Biallelic Inactivation' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA2') == 1 & grepl('BRCA2', Germline_mutation) & Tumor_zygosity == 'AI_LOH_ALT', 'Germline BRCA2 + LOH', Detailed_mutation_pattern))

#Somatic + LOH
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Biallelic Inactivation' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA1') >= 1 & !grepl('BRCA', Germline_mutation) & Tumor_zygosity == 'AI_LOH_ALT', 'Somatic BRCA1 + LOH', Detailed_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Biallelic Inactivation' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA2') >= 1 & !grepl('BRCA', Germline_mutation) & Tumor_zygosity == 'AI_LOH_ALT', 'Somatic BRCA2 + LOH', Detailed_mutation_pattern))

#Germline mutation + second mutation
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Biallelic Inactivation' & sapply(sample, pathogenic_mutation_maf_match, 'BRCA2') == 1 & grepl('BRCA2', Germline_mutation) & Tumor_zygosity != 'AI_LOH_ALT', 'Germline BRCA2 + second mutation', Detailed_mutation_pattern))


#Detail classification of Mono allelic Germline mutation
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Monoallelic Germline Inactivation' & grepl('BRCA1', Germline_mutation) & sapply(sample, pathogenic_mutation_maf_match, 'BRCA1') == 0, 'Monoallelic BRCA1 Germline Inactivation', Detailed_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Monoallelic Germline Inactivation' & grepl('BRCA2', Germline_mutation) & sapply(sample, pathogenic_mutation_maf_match, 'BRCA2') == 0, 'Monoallelic BRCA2 Germline Inactivation', Detailed_mutation_pattern))

#Detail classification of Mono allelic Somatic mutation
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Monoallelic Somatic Inactivation'  & sapply(sample, pathogenic_mutation_maf_match, 'BRCA1') >=1, 'Monoallelic BRCA1 Somatic Inactivation', Detailed_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Monoallelic Somatic Inactivation'  & sapply(sample, pathogenic_mutation_maf_match, 'BRCA2') >=1, 'Monoallelic BRCA2 Somatic Inactivation', Detailed_mutation_pattern))

#Detail classification of subclonal tumor sample
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Subclonal BRCA Inactivation'  & sapply(sample, pathogenic_mutation_maf_match, 'BRCA1') >= 1, 'Subclonal BRCA1 Inactivation', Detailed_mutation_pattern))

Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Subclonal BRCA Inactivation'  & sapply(sample, pathogenic_mutation_maf_match, 'BRCA2') >= 1, 'Subclonal BRCA2 Inactivation', Detailed_mutation_pattern))

#Detail classification of BRCA WT tumor sample
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'BRCA Wildtype Tumor', 'BRCA Wildtype Tumor', Detailed_mutation_pattern))

#Detail classification of Germline deletion tumor sample
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Detailed_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Germline Deletion' & grepl('BRCA1', Germline_mutation), 'BRCA1 Germline Deletion', Detailed_mutation_pattern))

#change the Tumor_type_based_mutation_pattern to factor and make each name into multiple lines
Annotated_merged_table <- Annotated_merged_table %>% 
  mutate(Tumor_type_based_mutation_pattern = factor(Tumor_type_based_mutation_pattern, levels = c('BRCA Wildtype Tumor', 'Biallelic Inactivation', 'Monoallelic Germline Inactivation', 'Monoallelic Somatic Inactivation',  'Germline Deletion', 'Subclonal BRCA Inactivation'))) %>%
  arrange(Tumor_type_based_mutation_pattern)         

levels(Annotated_merged_table$Tumor_type_based_mutation_pattern) <- gsub(' ', '\n', levels(Annotated_merged_table$Tumor_type_based_mutation_pattern))

print('Finish detail classification of mutation type')

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Part 3 starts here, draw plots and look at the association between different tumor types (based on mutations) and mutational signature 3

# read a file contains all mutational signatures for BRCA WT Tumors
WT_decomp <-fread('~/ifs/res/taylorlab/jonssonp/Proj_93017/mutsig/Proj_93017_data_mutations_extended_010318.spectrum.mutsig.out')

#read a file contains total number of mutations for BRCA WT Tumors
WT_sample_n <-fread('~/ifs/res/taylorlab/jonssonp/Proj_93017/mutsig/Proj_93017_data_mutations_extended_010318.spectrum.out')

# Create a new dataframe to add all 64 possible mutations to get total mutations
sample_n_total <- sample_n %>% select('sample' = Tumor_Sample_Barcode) %>% mutate(total = rowSums(sample_n[,-1]))

#Create a new dataframe to add all 64 possible mutations to get total mutations for BRCA WT tumor
WT_sample_n_total <- WT_sample_n %>% select('sample' = Tumor_Sample_Barcode) %>% mutate(total = rowSums(WT_sample_n[,-1]))

#combine total number mutation df together
total_number_mutation <- rbind(sample_n_total, WT_sample_n_total)

# Create a new dataframe that extract signature 3 info from decomp dataframe
decomp_sig_3 <- decomp %>% filter (Signature == 'Signature.3') %>% select ('sample' = Tumor_Sample_Barcode,mean)

# Create a new dataframe that extract signature 3 info from WT decomp dataframe
WT_decomp_sig_3 <- WT_decomp %>% filter (Signature == 'Signature.3') %>% select ('sample' = Tumor_Sample_Barcode,mean)

#combine total signature 3 df
total_decomp_sig_3 <- rbind(decomp_sig_3, WT_decomp_sig_3)

#Combine previous two dataframe according to sample ID (Tumor Sample Barcode), modify the mean and total to one column of signature 3 mutation activity (number).
combine_sample_decomp_sig_3 <- join(total_decomp_sig_3,total_number_mutation,by = 'sample', type = 'inner') 

combine_sample_decomp_sig_3 <- combine_sample_decomp_sig_3 %>% select(sample) %>% 
  mutate(signature_3_activity = combine_sample_decomp_sig_3$mean* combine_sample_decomp_sig_3$total) 

Annotated_merged_table <- Annotated_merged_table %>% filter(!(sample =='s_C_001261_T001_d' & Tumor_zygosity == 'subclonal'))

# combine annotated dataframe and signature 3 activity dataframe
combine_annotate_sig_3 <- join(Annotated_merged_table, combine_sample_decomp_sig_3, by = 'sample', type = 'inner')

#Make a function to delete the 'BRCA' in Tumor_type_based_mutation_pattern column
delete_brca <- function(df)  {
  name_list <- levels(df$Tumor_type_based_mutation_pattern)
  new_name_list <- vector()
  for (name in name_list) {
    if (grepl('BRCA\n', name)){
      name <- gsub('BRCA\n', '',name)
    }
    new_name_list <- c(new_name_list, name)
  }
  mapvalues(df$Tumor_type_based_mutation_pattern, from = name_list, to = new_name_list)
}

combine_annotate_sig_3$Tumor_type_based_mutation_pattern = delete_brca(combine_annotate_sig_3)
#create violin plot to look at the association of mutational events and signature 3 mutation

signature_3_tumor_asso_plot <- ggplot(combine_annotate_sig_3, aes(x = Tumor_type_based_mutation_pattern, y = signature_3_activity)) + 
  geom_violin() + 
  theme_bw()+
  theme (axis.text.x = element_text(face = 'bold', size = 12), axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none') +
  ylab('Signature 3 Activity') +
  theme(axis.title = element_text(face = 'bold'))+
  ylim(0,510)
  

signature_3_tumor_asso_plot <- signature_3_tumor_asso_plot + geom_jitter(aes(x = Tumor_type_based_mutation_pattern, y = signature_3_activity, color = Tumor_type_based_mutation_pattern), width = 0.15)+
  geom_boxplot(width = 0.15, outlier.shape = NA)
  
#use ggpubr package to calculate p value via Wilcoxon rank-sum test and add into graph
signature_3_tumor_asso_plot + stat_compare_means(label = 'p.signif', label.y = 480,  method = 'wilcox.test', method.args = list(alternative = 'greater'), ref.group = 'Wildtype\nTumor',fontface = 'bold', size = 5,
                                                 symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))) 

# Save the plot into R plot file
ggsave ("~/shenz/script/R_plot/BRCA_Signature_3_Plot/Whole/Mutation events signature 3 activity correlation plot.tiff", width = 8, height = 5)

print('Finish the plot of mutation events and signature 3 activity')
#-----------------------------------------------------------------------------------------------------------------
#Part 3.5 This part is going to look at the association between mutation events and the percentage of signature 3 , to 
# eliminate the effects of differential mutational burdens

combine_annotate_decomp_sig_3_percent <- join(Annotated_merged_table, total_decomp_sig_3, by = 'sample', type = 'inner')

#delete BRCA from column
combine_annotate_decomp_sig_3_percent$Tumor_type_based_mutation_pattern <- delete_brca(combine_annotate_decomp_sig_3_percent)

#create Violin plot and plot P value
signature_3_percent_tumor_asso_plot <- ggplot(combine_annotate_decomp_sig_3_percent, aes(x = Tumor_type_based_mutation_pattern, y = mean)) + 
  geom_violin() + 
  theme_bw()+
  theme (axis.text.x = element_text(face = 'bold', size = 12), axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none') +
  ylab('Signature 3 Percentage')+
  theme(axis.title = element_text(face = 'bold'))+
  ylim(0,1.2)

signature_3_percent_tumor_asso_plot <- signature_3_percent_tumor_asso_plot + geom_jitter(  aes(x = Tumor_type_based_mutation_pattern, y = mean, color = Tumor_type_based_mutation_pattern), width = 0.15)+
  geom_boxplot(width = 0.15, outlier.shape = NA) 

#use ggpubr package to calculate p value via Wilcoxon rank-sum test and add into graph
signature_3_percent_tumor_asso_plot + stat_compare_means(label = 'p.signif', label.y = 1.15, ref.group = 'Wildtype\nTumor', method = 'wilcox.test', method.args = list(alternative = 'greater'), fontface = 'bold', size = 5,
                                                         symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")))

# Save the plot into R plot file
ggsave ("~/shenz/script/R_plot/BRCA_Signature_3_Plot/Whole/Mutation events signature 3 percentage correlation plot.tiff", width = 8, height = 5)

print('Finish the plot of mutation events and signature 3 percentage')


#------------------------------------------------------------------------------------------------------------------
#Part4 sort the table according to different cancer types and plot the correlations with signature 3 activity
#extract cancer type related info
master_edit <- master %>% select('sample' = CMO_ID, Cancer_Type) 

#read the file contaning cancer types of BRCA WT Cancer
WT_master <- fread('~/ifs/res/taylorlab/jonssonp/Proj_93017/Proj_93017_data_clinical_sample_010318.txt')

#extract cancer type related info for BRCA WT tumor
WT_master <- WT_master %>% select('sample' = SAMPLE_ID,'Cancer_Type' = CANCER_TYPE)

#combine to cancer type table
cancer_type_master <- rbind(WT_master, master_edit)

#merge with signature 3 dataframe
combine_annotate_sig_3_type <- join(combine_annotate_sig_3, cancer_type_master, by = 'sample', type = 'inner')

# group together the lower incidence cancer (less than 15) as "rare cases" (not rare in tumor frequency, rare in this cohort)

combine_annotate_sig_3_type_regroup <- combine_annotate_sig_3_type %>% group_by(Cancer_Type) %>% dplyr::mutate (Cases = n()) %>% ungroup()

combine_annotate_sig_3_type_regroup <- combine_annotate_sig_3_type_regroup %>% mutate(Cancer_Type = ifelse(Cases < 15, 'Rare Case Cancer', Cancer_Type))

combine_annotate_sig_3_type_regroup[,"signature_3_activity"] = log(combine_annotate_sig_3_type_regroup[,"signature_3_activity"],2)

#reorder the column and make the names multiline
combine_annotate_sig_3_type_regroup <- combine_annotate_sig_3_type_regroup %>% 
  mutate(Cancer_Type = factor(Cancer_Type, levels = c('Ovarian Cancer', 
                                                      'Breast Cancer', 
                                                      'Bladder Cancer', 
                                                      'Colorectal Cancer',  
                                                      'Glioma', 
                                                      'Non-Small Cell Lung Cancer', 
                                                      'Prostate Cancer', 
                                                      'Endometrial Cancer',
                                                      'Melanoma', 
                                                      'Hepatobiliary Cancer', 
                                                      'Rare Case Cancer', 
                                                      'Cancer of Unknown Primary',
                                                      'Pancreatic Cancer',
                                                      'Soft Tissue Sarcoma'))) %>%
  arrange(Cancer_Type)         

levels(combine_annotate_sig_3_type_regroup$Cancer_Type) <- gsub(' ', '\n', levels(combine_annotate_sig_3_type_regroup$Cancer_Type))

#plot the assoaciation between cancer types and signature 3 activity
signature_3_activity_cancer_type_asso_plot <- ggplot(combine_annotate_sig_3_type_regroup, aes(x = Cancer_Type, y = signature_3_activity)) + 
  geom_violin()+
  theme_bw()+
  theme (axis.text.x = element_text(face = 'bold', size = 10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none') +
  xlab('Cancer Type') +
  ylab('Signature 3 Activity (log2)')+
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter( aes(x = Cancer_Type, y = signature_3_activity, color = Cancer_Type), width = 0.15)+
  geom_boxplot(width = 0.2, outlier.shape = NA)+
  ylim(0,12)

#use ggpubr package to calculate p value via Wilcoxon rank-sum test and add into graph
signature_3_activity_cancer_type_asso_plot + stat_compare_means(ref.group = 'Ovarian\nCancer', label = 'p.signif', label.y = 12, method = 'wilcox.test', method.args = list(alternative = 'less'), size = 5, fontface = 'bold')

# Save the plot into R plot file
ggsave ("~/shenz/script/R_plot/BRCA_Signature_3_Plot/Cancer_Type_Signature_3/Cancer type signature 3 activity correlation plot.tiff", width = 12, height = 7)

print('Finish the plot of cancer type and signature 3 activity correlation plot')

#--------------------------------------------------------------------------------------------------------------------------------
#Part4.5 plot the correlation between Cancer Types and Signature 3 Percentage.

#merge cancer type and mutation signature 3 percent dataframe
combine_type_sig_3_percent <- join(combine_annotate_decomp_sig_3_percent, cancer_type_master, by = 'sample', type = 'inner')

# group together the lower incidence cancer (less than 15) as "rare cases" (not rare in tumor frequency, rare in this cohort)

combine_type_sig_3_percent <- combine_type_sig_3_percent %>% group_by(Cancer_Type) %>% dplyr::mutate (Cases = n()) %>% ungroup()

combine_type_sig_3_percent <- combine_type_sig_3_percent %>% dplyr::mutate(Cancer_Type = ifelse(Cases < 15, 'Rare Case Cancer', Cancer_Type))

#reorder the column and make the names multiline
combine_type_sig_3_percent <- combine_type_sig_3_percent %>% 
  mutate(Cancer_Type = factor(Cancer_Type, levels = c('Ovarian Cancer', 
                                                      'Breast Cancer', 
                                                      'Bladder Cancer', 
                                                      'Colorectal Cancer',  
                                                      'Glioma', 
                                                      'Non-Small Cell Lung Cancer', 
                                                      'Prostate Cancer', 
                                                      'Endometrial Cancer',
                                                      'Melanoma', 
                                                      'Hepatobiliary Cancer', 
                                                      'Rare Case Cancer', 
                                                      'Cancer of Unknown Primary',
                                                      'Pancreatic Cancer',
                                                      'Soft Tissue Sarcoma'))) %>%
  arrange(Cancer_Type)         

levels(combine_type_sig_3_percent$Cancer_Type) <- gsub(' ', '\n', levels(combine_type_sig_3_percent$Cancer_Type))

#plot the assoaciation between cancer types and signature 3 activity
signature_3_percent_cancer_type_asso_plot <- ggplot(combine_type_sig_3_percent, aes(x = Cancer_Type, y = mean)) + 
  geom_violin()+
  theme_bw()+
  theme (axis.text.x = element_text(face = 'bold', size = 10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none') +
  xlab('Cancer Type') +
  ylab('Signature 3 Percent')+
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter( aes(x = Cancer_Type, y = mean, color = Cancer_Type), width = 0.15)+
  geom_boxplot(width = 0.2, outlier.shape = NA)+
  ylim(0,1)

#use ggpubr package to calculate p value via Wilcoxon rank-sum test and add into graph
signature_3_percent_cancer_type_asso_plot + stat_compare_means(ref.group = 'Ovarian\nCancer', label = 'p.signif', label.y = 1, method = 'wilcox.test', method.args = list(alternative = 'less'), size = 5, fontface = 'bold')

# Save the plot into R plot file
ggsave ("~/shenz/script/R_plot/BRCA_Signature_3_Plot/Cancer_Type_Signature_3/Cancer type signature 3 Percentage correlation plot.tiff", width = 12, height = 7)

print('Finish the plot of cancer type and signature 3 Percentage correlation plot')

#--------------------------------------------------------------------------------------------------------------------------------
#Part5 plot the correlation between Mutational Events and mutational burdens.(also color the hypermuated point)

#Add a new mutational burdens column to the original annotated table

mutation_burden_combine <- join(Annotated_merged_table, total_number_mutation, by = 'sample', type = 'inner')

colnames(mutation_burden_combine)[7] <- 'Total_Mutational_Burdens'

mutation_burden_combine[,'Total_Mutational_Burdens'] <- log(mutation_burden_combine[,'Total_Mutational_Burdens'],2) 

#Hypermutated part
Mutation_medium <- median(mutation_burden_combine$Total_Mutational_Burdens)
Mutation_IQR <- IQR(mutation_burden_combine$Total_Mutational_Burdens)

Mutation_cutoff <- Mutation_medium + 1.5 * Mutation_IQR

mutation_burden_combine_hypermutation <- mutation_burden_combine %>% filter(Total_Mutational_Burdens >= Mutation_cutoff)

#read the msi score text
msi <- fread('~/ifs/res/taylorlab/impact_msisensor/data_msi.txt')

#read the Sample ID mapping file
mapping_ID <- fread('~/ifs/res/taylorlab/jonssonp/Proj_93017/Proj_93017_sample_map_010318.txt')

mapping_ID$SAMPLE_ID <- sapply(mapping_ID$SAMPLE_ID , sample_conv)

#make a file contains all matching IDs of BRCA mutated and WT tumors
total_matching_ID <- master %>% select ('SAMPLE_ID' = CMO_ID, DMP_ID)

mapping_ID_subset <- mapping_ID %>% select (SAMPLE_ID, DMP_ID)

total_matching_ID <- rbind(total_matching_ID, mapping_ID_subset)

#Sample ID converter
WES_converter <- function(sample_ID_list) {
  sample_DMP <- vector()
  for (sample_ID in sample_ID_list) {
  if (grepl('WES',sample_ID)){
    sample_ID<- gsub('WES', 'IM6', sample_ID)
    sample_DMP <- c(sample_DMP, sample_ID)
  }
  
  if(grepl('s_C_', sample_ID)){
    index <- match(sample_ID, total_matching_ID$SAMPLE_ID)
    sample_ID <- total_matching_ID$DMP_ID[index]
    sample_DMP <- c(sample_DMP,sample_ID)
  }
  }
  sample_DMP
}

#msi score finder
msi_finder <- function(sample_ID_vector) {
  msi_score_vector <- vector()
  for (sample_ID in sample_ID_vector) {
    index <- match(sample_ID, msi$SAMPLE_ID)
    msi_score <- msi$MSI_SCORE[index]
    msi_score_vector <- c(msi_score_vector, msi_score)
  }
  msi_score_vector
}

wt_maf_revised <- wt_maf %>% select(Tumor_Sample_Barcode,Hugo_Symbol, HGVSp_Short, Variant_Classification)
maf_revised <- maf %>% select(Tumor_Sample_Barcode,Hugo_Symbol, HGVSp_Short, Variant_Classification)

maf_revised <- rbind(maf_revised, wt_maf_revised)

#POLE finder
POLE_finder <- function(sample_ID_vector) {
  POLE_vector <- vector()
  subset_df <- maf_revised %>% filter (grepl('POLE', Hugo_Symbol)  & (grepl('P286', HGVSp_Short) | grepl('V411', HGVSp_Short)))
  for (sample_ID in sample_ID_vector) {
    index <- match(sample_ID, subset_df$Tumor_Sample_Barcode)
    if (!is.na(index)){
    Somatic_mutation<- subset_df$HGVSp_Short[index]
    POLE_vector <- c(POLE_vector, Somatic_mutation)
    }
    else {
      POLE_vector <- c(POLE_vector, 'POLE Somatic WT')
    }
  }
  POLE_vector
}




#extract MMR signature info
MMR_decomp <- rbind(decomp, WT_decomp) %>% filter (Signature == 'Signature.6' | Signature == 'Signature.15'|
                                                       Signature == 'Signature.20'|Signature == 'Signature.26') %>%
  select(Tumor_Sample_Barcode, mean) %>% group_by(Tumor_Sample_Barcode) %>% dplyr::summarise(MMR_signature_percentage = sum(mean) )


#make a function to extract MMR signature
MMR_signature_extracter <- function(sample_ID_list) {
  MMR_SCORE <- vector()
  for (sample_ID in sample_ID_list){
    subset <- MMR_decomp %>% filter(Tumor_Sample_Barcode == sample_ID)
    MMR_SCORE <- c(MMR_SCORE, subset$MMR_signature_percentage)
  }
  MMR_SCORE
}

#extract APOBEC signature info
APOBEC_decomp <- rbind(decomp, WT_decomp) %>% filter (Signature == 'Signature.2' | Signature == 'Signature.13') %>%
  select(Tumor_Sample_Barcode, mean) %>% group_by(Tumor_Sample_Barcode) %>% dplyr::summarise(APOBEC_signature_percentage = sum(mean) )

#make a function to extract APOBEC signature
APOBEC_signature_extracter <- function(sample_ID_list) {
  APOBEC_SCORE <- vector()
  for (sample_ID in sample_ID_list){
    subset <- APOBEC_decomp %>% filter(Tumor_Sample_Barcode == sample_ID)
    APOBEC_SCORE <- c(APOBEC_SCORE, subset$APOBEC_signature_percentage)
  }
  APOBEC_SCORE
}

#extract APOBEC signature info
POLE_decomp <- rbind(decomp, WT_decomp) %>% filter (Signature == 'Signature.10') %>%
  select(Tumor_Sample_Barcode, 'POLE_signature_percentage' = mean)

#make a function to extract POLE signature
POLE_signature_extracter <- function(sample_ID_list) {
  POLE_SCORE <- vector()
  for (sample_ID in sample_ID_list){
    subset <- POLE_decomp %>% filter(Tumor_Sample_Barcode == sample_ID)
    POLE_SCORE <- c(POLE_SCORE, subset$POLE_signature_percentage)
  }
  POLE_SCORE
}

# Add a column of DMP_ID conversion, msi_SCORE, POLE and APOBEC germline mutation, finally give proposed mechanism
mutation_burden_combine_hypermutation_convert <- mutation_burden_combine_hypermutation %>% 
  select(sample) %>% mutate(DMP_ID = WES_converter(sample)) %>% mutate(MSI_SCORE = msi_finder(DMP_ID)) %>%
  mutate(MSI_SCORE = ifelse(is.na(MSI_SCORE), 0.04, MSI_SCORE)) %>%
  mutate(POLE_Somatic_Mutation = POLE_finder(sample)) %>%
  mutate(APOBEC_Score = APOBEC_signature_extracter(sample))%>%
  mutate(MMR_Signature_Percentage = MMR_signature_extracter(sample))%>%
  mutate(POLE_Score = POLE_signature_extracter(sample))%>%
  mutate(Proposed_mechanism = ifelse(MSI_SCORE >= Mutation_cutoff & MMR_Signature_Percentage >= 0.5, 'MSI High', 'Unknown Significance')) %>%
  mutate(Proposed_mechanism = ifelse(APOBEC_Score >= 0.5, 'APOBEC Mutation', Proposed_mechanism)) %>%
  mutate(Proposed_mechanism = ifelse(POLE_Score >= 0.5, 'POLE Mutation', Proposed_mechanism))
  

mutation_burden_combine_hypermutation_convert_subset <- mutation_burden_combine_hypermutation_convert %>% 
  select(sample,Proposed_mechanism)

#look at the proposed mechanism of hypermutated tumor (add to the previous plot)
mutation_burden_combine_hyper <- join(mutation_burden_combine, mutation_burden_combine_hypermutation_convert_subset, by = 'sample') %>%
  mutate(Proposed_mechanism = ifelse(is.na(Proposed_mechanism), 'nonhypermutated tumor',Proposed_mechanism )) 

#plot the correlation between mutational events and mutational burdens
mutation_events_burden_asso_plot <- ggplot(mutation_burden_combine_hyper, aes(x = Tumor_type_based_mutation_pattern, y = Total_Mutational_Burdens))+
  geom_violin()+
  theme_bw()+
  theme(axis.text.x = element_text(face = 'bold', size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none')+
  labs(x = 'Mutational Events', y = 'Total Mutational Burdens (log2)' )+
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter(aes(x = Tumor_type_based_mutation_pattern, y = Total_Mutational_Burdens, color = Proposed_mechanism), width = 0.15, size = 4)+
  scale_color_hue(l=45, c=100)+
  geom_boxplot(width = 0.15, outlier.shape = NA)+
  ylim(3.5,15)

#calculate significance(P value) among Biallelic Inactivation, Monoallelic Germline and Monoallelic Somatic via Wilcox-Rank-Sum test. Add on
#to the plot
mutation_events_burden_asso_plot + stat_compare_means(label = 'p.signif', label.y = 15, ref.group = 'BRCA\nWildtype\nTumor', method = 'wilcox.test', method.args = list(alternative = 'greater'),fontface = 'bold', size = 5)

#Save the plot to Rplot file
ggsave('~/shenz/script/R_plot/BRCA_Mutational_Burdens/Mutational Burdens versus Mutational Events Correlation Graph.tiff', width = 8, height = 5)

print('Finish the plot of mutational events and mutation burden correlation plot (hypermutated highlight)')

#----------------------------------------------------------------------------------------------------------------------------------------------
#Part6 read the file containing all microhomology information for each sample, generate a plot of the correlation between
#microhomology and mutational events

#read the micrgohomology maf into R
microhomo <- fread('~/shenz/gray_brca/joe/final-filtered-facets-microhomology.maf')

#read the microhomology maf of BRCA WT tumor into R
WT_microhomo <- fread('~/shenz/gray_brca/joe/Proj_93017_data_mutations_extended_010318_microhomology_output.vep.maf')

#set a cutoff for microhomology length
microhomo_filtered <- microhomo %>% group_by(Sample) %>% filter(Homology_Length > 10) %>% dplyr::summarise(number_of_homology = n())

wt_microhomo_filtered <- WT_microhomo %>% group_by(Sample) %>% filter(Homology_Length > 10) %>% dplyr::summarise(number_of_homology = n())

#calculate total deletions
maf_total_INDEL <- maf %>% group_by(Tumor_Sample_Barcode) %>% filter(Variant_Type == 'DEL' | Variant_Type == 'INS') %>% dplyr::summarise(Total_Indel = n())

wt_maf_total_INDEL <- wt_maf %>% group_by(Tumor_Sample_Barcode) %>% filter(Variant_Type == 'DEL' | Variant_Type == 'INS') %>% dplyr::summarise(Total_Indel = n())

#set up a df of microhomology percent in total deletions
microhomo_percent<- merge(microhomo_filtered, maf_total_INDEL, by.x = 'Sample', by.y = 'Tumor_Sample_Barcode', all = FALSE)

microhomo_percent <- microhomo_percent %>% select ('sample' = Sample) %>% mutate(Microhomology_Percent = microhomo_percent$number_of_homology / microhomo_percent$Total_Indel)

#set up a df of microhomology percent in total deletions for BRCA WT tumor sample
wt_microhomo_percent<- merge(wt_microhomo_filtered, wt_maf_total_INDEL, by.x = 'Sample', by.y = 'Tumor_Sample_Barcode', all = FALSE) 

wt_microhomo_percent <- wt_microhomo_percent %>% select ('sample' = Sample) %>% mutate(Microhomology_Percent = wt_microhomo_percent$number_of_homology / wt_microhomo_percent$Total_Indel )

#first bind two microhomo-percent table together, then merge with annotated table
total_microhomo_percent <- rbind(microhomo_percent, wt_microhomo_percent)

microhomo_anno_comb <- join(Annotated_merged_table, total_microhomo_percent, by = 'sample', type = 'inner')

#plot the correlation between mutational events and number of microhomology
mutation_events_microhomo_asso_plot <- ggplot(microhomo_anno_comb, aes(x = Tumor_type_based_mutation_pattern, y = Microhomology_Percent))+
  geom_violin()+
  theme_bw()+
  theme(axis.text.x = element_text(face = 'bold', size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none')+
  labs(x = 'Mutational Events', y = 'Microhomology / Total Indels (bp > 10)' )+
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter(aes(x = Tumor_type_based_mutation_pattern, y = Microhomology_Percent, color = Tumor_type_based_mutation_pattern ), width = 0.15)+
  geom_boxplot(width = 0.15, outlier.shape = NA)+
  ylim(0,0.6)

#calculate significance(P value) among Biallelic Inactivation, Monoallelic Germline and Monoallelic Somatic via Wilcox-Rank-Sum test. Add on
#to the plot
mutation_events_microhomo_asso_plot + stat_compare_means(label = 'p.signif', label.y = 0.6, method = 'wilcox.test', method.args = list(alternative = 'greater'), ref.group = 'BRCA\nWildtype\nTumor', size = 5, fontface = 'bold')

#Save the plot to Rplot file
ggsave('~/shenz/script/R_plot/BRCA_Microhomology/Mutational Events versus Microhomology Percent Correlation Graph.tiff', width = 8, height = 5)

print('Finish the plot of Mutation events and Microhomology percentage correlation plot')


#--------------------------------------------------------------------------------------------------------------------
#Part7 Look at and plot the association between NtAI score and mutational events

#read the file that contains the NtAI score
NtAI_score <- fread('~/ifs/res/taylorlab/jonssonp/gray_brca/hrd_scores_wes/ntai-scores-em.txt')

#read the file that contains the NtAI score for BRCA WT patient
wt_NtAI_score <- fread('~/ifs/res/taylorlab/jonssonp/Proj_93017/hrd_scores/ntai-scores-em.txt')


#make a simple function to convert the ID in HRD scores to ID used in wt_maf df 
sample_ID_match <- function(sample_ID) {
  if (TRUE %in% grepl(sample_ID, mapping_ID$SAMPLE_ID)){
    index <- grep(sample_ID, mapping_ID$SAMPLE_ID)
    matching_ID <- mapping_ID$DMP_ID[index]
    matching_ID <- gsub('IM6', 'WES', matching_ID)
  }
  else {
        matching_ID <- NA
  }
  matching_ID
}

wt_NtAI_score$Sample <- sapply(wt_NtAI_score$Sample, sample_ID_match)

wt_NtAI_score <- wt_NtAI_score %>% filter(!is.na(Sample))

NtAI<- NtAI_score %>% select('sample' = Sample,NtAI)

wt_NtAI <- wt_NtAI_score %>% select('sample' = Sample, NtAI)

#combine two NtAI df in to one total NtAI df 
total_NtAI <- rbind(NtAI, wt_NtAI)

#merge two dataframe
anno_combine_NtAI <- join(Annotated_merged_table, total_NtAI, by = 'sample', type = 'inner')

#delete BRCA from the variable names
anno_combine_NtAI$Tumor_type_based_mutation_pattern = delete_brca(anno_combine_NtAI)

#plot the correlation between NtAI score and mutational events
mutation_events_NtAI_asso_plot <- ggplot(anno_combine_NtAI, aes(x = Tumor_type_based_mutation_pattern, y = NtAI))+
  geom_violin()+ 
  theme_bw()+
  theme(axis.text.x = element_text(face = 'bold', size = 12), axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = 'none')+
  labs(y = 'NTAI Score')+
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter(aes(x = Tumor_type_based_mutation_pattern, y = NtAI, color = Tumor_type_based_mutation_pattern), width = 0.15)+
  geom_boxplot(width = 0.15, outlier.shape = NA)+
  ylim(0, 38)

#calculate significance(P value) among Biallelic Inactivation, Monoallelic Germline and Monoallelic Somatic via Wilcox-Rank-Sum test. Add on
#to the plot
mutation_events_NtAI_asso_plot + stat_compare_means(label = 'p.signif', label.y = 38, method = 'wilcox.test', method.args = list(alternative = 'greater'), ref.group = 'Wildtype\nTumor', fontface = 'bold', size = 5,
                                                    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")))

#Save the plot to Rplot file
ggsave('~/shenz/script/R_plot/BRCA_Genomic_Scar/Whole/Mutational Events versus NtAI Score Correlation Graph.tiff', width = 8, height = 5)


print('Finish the plot of Mutation events and NTAI Score correlation plot')

#--------------------------------------------------------------------------------------------------------------------
#Part7 Look at and plot the association between LST score and mutational events

#read file containing information of LST Scores
LST_Scores <- fread('~/ifs/res/taylorlab/jonssonp/gray_brca/hrd_scores_wes/lst-scores-em.txt')

#read the file contains LST scores for BRCA WT tumor
WT_LST_Scores  <- fread('~/ifs/res/taylorlab/jonssonp/Proj_93017/hrd_scores/lst-scores-em.txt')

WT_LST_Scores$Sample <- sapply(WT_LST_Scores$Sample, sample_ID_match)

WT_LST_Scores <- WT_LST_Scores %>% filter(!is.na(Sample))

#combine to LST Score df into one df
total_LST <- rbind(LST_Scores, WT_LST_Scores)

colnames(total_LST)[1] <- 'sample'

#merge two dataframe
anno_combine_LST <- join(Annotated_merged_table, total_LST, by = 'sample', type = 'inner')

#plot the correlation between LST score and mutational events
mutation_events_LST_asso_plot <- ggplot(anno_combine_LST, aes(x = Tumor_type_based_mutation_pattern, y = LST))+
  geom_violin()+ 
  theme_bw()+
  theme(axis.text.x = element_text(face = 'bold', size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = 'none')+
  labs(x = 'Mutational Events', y = 'LST Score')+
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter(aes(x = Tumor_type_based_mutation_pattern, y = LST, color = Tumor_type_based_mutation_pattern), width = 0.15)+
  geom_boxplot(width = 0.15, outlier.shape = NA)+
  ylim(0, 52)

#calculate significance(P value) among Biallelic Inactivation, Monoallelic Germline and Monoallelic Somatic via Wilcox-Rank-Sum test. Add on
#to the plot

mutation_events_LST_asso_plot + stat_compare_means(label = 'p.signif', label.y = 52, method = 'wilcox.test', method.args = list(alternative = 'greater'), ref.group = 'BRCA\nWildtype\nTumor', fontface = 'bold', size = 5)

#Save the plot to Rplot file
ggsave('~/shenz/script/R_plot/BRCA_Genomic_Scar/Whole/Mutational Events versus LST Score Correlation Graph.tiff', width = 8, height = 5)


print('Finish the plot of Mutation events and LST Score correlation plot')


#--------------------------------------------------------------------------------------------------------------------
#Part8 Look at and plot the association between HRD-LOH score and mutational events

#read file containing information of HRD_LOH Scores
HRD_LOH_Scores <- fread('~/ifs/res/taylorlab/jonssonp/gray_brca/hrd_scores_wes/hrdloh-scores-em.txt')

#read the file containing HRD_LOH Scores for BRCA WT tumor
WT_HRD_LOH <- fread('~/ifs/res/taylorlab/jonssonp/Proj_93017/hrd_scores/hrdloh-scores-em.txt')

WT_HRD_LOH$Sample <- sapply(WT_HRD_LOH$Sample, sample_ID_match)

WT_HRD_LOH <- WT_HRD_LOH %>% filter(!is.na(Sample))

#combine two HRD_LOH df into one dataframe
total_HRD_LOH <- rbind(HRD_LOH_Scores, WT_HRD_LOH)
colnames(total_HRD_LOH)[1] <- 'sample'

#merge two dataframe
anno_combine_HRD <- join(Annotated_merged_table, total_HRD_LOH, by = 'sample', type = 'inner')

#plot the correlation between HRD-LOH score and mutational events
mutation_events_HRD_LOH_asso_plot <- ggplot(anno_combine_HRD, aes(x = Tumor_type_based_mutation_pattern, y = HRD_LOH))+
  geom_violin()+ 
  theme_bw()+
  theme(axis.text.x = element_text(face = 'bold', size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = 'none')+
  labs(x = 'Mutational Events', y = 'HRD-LOH Score')+
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter(aes(x = Tumor_type_based_mutation_pattern, y = HRD_LOH, color = Tumor_type_based_mutation_pattern), width = 0.15)+
  geom_boxplot(width = 0.15, outlier.shape = NA)+
  ylim(0, 38)

#calculate significance(P value) among Biallelic Inactivation, Monoallelic Germline and Monoallelic Somatic via Wilcox-Rank-Sum test. Add on
#to the plot
mutation_events_HRD_LOH_asso_plot + stat_compare_means(label = 'p.signif', label.y = 38, method = 'wilcox.test', method.args = list(alternative = 'greater'), ref.group = 'BRCA\nWildtype\nTumor', fontface = 'bold', size = 5)

#Save the plot to Rplot file
ggsave('~/shenz/script/R_plot/BRCA_Genomic_Scar/Whole/Mutational Events versus HRD-LOH Score Correlation Graph.tiff', width = 8, height = 5)


print('Finish the plot of Mutation events and HRD_LOH Score correlation plot')

#--------------------------------------------------------------------------------------------------------------------
#Part9 Look at and plot the association between number of Indels and mutational events

#combine previous established dataframes to create a total Indel df
total_indel <- rbind(maf_total_INDEL, wt_maf_total_INDEL)
colnames(total_indel)[1] <- 'sample'

#merge two dataframes together
anno_tot_INDELS_combine <- join(Annotated_merged_table, total_indel, by = 'sample', type = 'inner')

anno_tot_INDELS_combine[,"Total_Indel"] <- log(anno_tot_INDELS_combine[,"Total_Indel"],2)

colnames(anno_tot_INDELS_combine)[7] <- 'Total_Indel_log2'

#plot the correlation between number of Indels and mutational events
mutation_events_indel_asso_plot <- ggplot(anno_tot_INDELS_combine, aes(x = Tumor_type_based_mutation_pattern, y = Total_Indel_log2))+
  geom_violin()+ 
  theme_bw()+
  theme(axis.text.x = element_text(face = 'bold', size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = 'none')+
  labs(x = 'Mutational Events', y = 'Total Indels(log2)')+
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter(aes(x = Tumor_type_based_mutation_pattern, y = Total_Indel_log2, color = Tumor_type_based_mutation_pattern), width = 0.15)+
  geom_boxplot(width = 0.15, outlier.shape = NA)+
  ylim(0, 10)
#calculate significance(P value) among Biallelic Inactivation, Monoallelic Germline and Monoallelic Somatic via Wilcox-Rank-Sum test. Add on
#to the plot

mutation_events_indel_asso_plot + stat_compare_means(label = 'p.signif', label.y = 10, method = 'wilcox.test', method.args = list(alternative = 'greater'), ref.group = 'BRCA\nWildtype\nTumor', size = 5, fontface = 'bold')

#Save the plot to Rplot file
ggsave('~/shenz/script/R_plot/BRCA_total_indel/Mutational Events versus Indel numbers Correlation Graph.tiff', width = 8, height = 5)

print('Finish the plot of Mutation events and Indel correlation plot')
#--------------------------------------------------------------------------------------------------------------------
#Part10 Add background graph one to show the overall distributions of Mutational events in the Gray_BRCA cohort
Annotated_merged_table_bg1 <- Annotated_merged_table %>% mutate(Tumor_type_based_mutation_pattern = factor(Tumor_type_based_mutation_pattern, levels = rev(levels(Tumor_type_based_mutation_pattern))))

#Delete BRCA from column
Annotated_merged_table_bg1$Tumor_type_based_mutation_pattern <- delete_brca(Annotated_merged_table_bg1)

Mutational_events_distribution_plot <- ggplot(Annotated_merged_table_bg1, aes(x= Tumor_type_based_mutation_pattern)) +
  geom_bar(stat = 'count', aes(fill = Tumor_type_based_mutation_pattern))+
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_text(face = 'bold', size = 12), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = 'none')+
  labs(y = 'Number of Samples')+
  theme(axis.title = element_text(face = 'bold', size = 12))+
  ylim(0, 430)+
  geom_text(stat = 'count', aes(label =..count..), hjust = -0.4, size = 6, fontface = 'bold')

#save the file to Rplot file
ggsave("~/shenz/script/R_plot/Presentation_Background_Graph/Mutation_Events_Distribution/Mutational events number count (Overall) background graph one.tiff", width = 8, height = 6)

Annotated_merged_table_bg1_bi_subset <- Annotated_merged_table_bg1 %>% filter(Tumor_type_based_mutation_pattern == 'Biallelic\nInactivation') 


Mutational_events_distribution_plot <- ggplot(Annotated_merged_table_bg1_bi_subset, aes(x= Detailed_mutation_pattern)) +
  geom_bar(stat = 'count', aes(fill = Detailed_mutation_pattern))+
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_text(face = 'bold', size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = 'none')+
  labs(x = '\nMutational Events', y = 'Number of sample')+
  theme(axis.title = element_text(face = 'bold', size = 12))+
  ylim(0, 60)+
  geom_text(stat = 'count', aes(label =..count..), hjust = -0.4, size = 6, fontface = 'bold')

#save the file to Rplot file
ggsave("~/shenz/script/R_plot/Presentation_Background_Graph/Mutation_Events_Distribution/Mutational events number count background graph one (Biallelic Inactivation).tiff", width = 8, height = 6)

Annotated_merged_table_bg1_rest <- Annotated_merged_table_bg1 %>% filter(Tumor_type_based_mutation_pattern != 'Biallelic\nInactivation' & Tumor_type_based_mutation_pattern != 'BRCA\nWildtype\nTumor') %>%
  mutate(Detailed_mutation_pattern = factor(Detailed_mutation_pattern, levels = rev(unique(Detailed_mutation_pattern))))

Mutational_events_distribution_plot <- ggplot(Annotated_merged_table_bg1_rest, aes(x= Detailed_mutation_pattern)) +
  geom_bar(stat = 'count', aes(fill = Detailed_mutation_pattern))+
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_text(face = 'bold', size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = 'none')+
  labs(x = 'Mutational Events', y = 'Number of sample')+
  theme(axis.title = element_text(face = 'bold', size = 12))+
  ylim(0, 35)+
  geom_text(stat = 'count', aes(label =..count..), hjust = -0.4, size = 6, fontface = 'bold')

#save the file to Rplot file
ggsave("~/shenz/script/R_plot/Presentation_Background_Graph/Mutation_Events_Distribution/Mutational events number count background graph one (Non-Biallelic).tiff", width = 8, height = 6)

print('Finish the plot of Mutation Events Overall distrinution (bg1), one overall and two specialized')

#--------------------------------------------------------------------------------------------------------------------
#Part11 Add background graph two to show the overall distributions of Cancer types in the Gray_BRCA cohort

#add cancer type information to annotated table
Annotated_merged_table_bg2 <- join(Annotated_merged_table, cancer_type_master , by = 'sample', type = 'inner') %>%
  group_by(Cancer_Type)%>%
  dplyr::mutate(cases = n()) %>%
  ungroup()%>%
  mutate(Cancer_Type = ifelse(cases < 15, 'Others', Cancer_Type))   %>%
  arrange(cases,Cancer_Type) %>%
  mutate(Cancer_Type = factor(Cancer_Type,levels = unique(Cancer_Type)))


Cancer_type_distribution_plot <- ggplot(Annotated_merged_table_bg2, aes(x = Cancer_Type)) +
  geom_bar(stat = 'count', aes(fill = Cancer_Type))+
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_text(face = 'bold', size = 12), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor =  element_blank(), legend.position = 'none')+
  labs( y = 'Number of Samples')+
  theme(axis.title = element_text(face = 'bold', size = 12))+
  geom_text(stat = 'count', aes(label =..count..), hjust = -0.4, size = 6, fontface = 'bold')+
  ylim(0, 130)

#save the file to Rplot file
ggsave("~/shenz/script/R_plot/Presentation_Background_Graph/Cancer_Type_Distribution/Cancer Types number count background graph two.tiff", width = 8, height = 6)

print('Finish the plot of Cancer Types Overall distrinution (bg2)')

#---------------------------------------------------------------------------------------------------------------------
# Part 12 plot the corelation between mutational events and signature 3 activity in prostate, ovarian and breast cancer
combine_annotate_sig_3_OBP <- join(combine_annotate_sig_3, cancer_type_master, by = 'sample', type = 'inner') %>%
  filter(Cancer_Type == 'Breast Cancer' | Cancer_Type == 'Prostate Cancer' | Cancer_Type == 'Ovarian Cancer')

#create violin plot to look at the association of mutational events and signature 3 mutation

signature_3_tumor_asso_plot_OBP <- ggplot(combine_annotate_sig_3_OBP, aes(x = Tumor_type_based_mutation_pattern, y = signature_3_activity)) + 
  geom_violin() + 
  theme_bw()+
  theme (axis.text.x = element_text(face = 'bold', size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none') +
  xlab('Mutation Events') +
  ylab('Signature 3 Acvtivity') +
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter(aes(x = Tumor_type_based_mutation_pattern, y = signature_3_activity, color = Tumor_type_based_mutation_pattern), width = 0.15)+
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ylim(0,470)

#statistical comparison using wilcox
signature_3_tumor_asso_plot_OBP + stat_compare_means(label = 'p.signif', label.y = 460, method = 'wilcox.test', method.args = list(alternative = 'greater'), ref.group = 'BRCA\nWildtype\nTumor',fontface = 'bold', size = 5)

ggsave('~/shenz/script/R_plot/BRCA_Signature_3_Plot/Ovarian_Breast_Prostate/Mutational Events and Signature 3 Activity Association Plot in Ovarian, Prostate and Breast Cancer.tiff', width = 8, height = 5)

print('Finish Mutational Events and Signature 3 Activity Association Plot in OBP')

#---------------------------------------------------------------------------------------------------------------------
# Part 13 plot the corelation between mutational events and signature 3 percentage in prostate, ovarian and breast cancer

combine_annotate_decomp_sig_3_percent_OBP <- join(combine_annotate_decomp_sig_3_percent, cancer_type_master, by = 'sample', type = 'inner') %>%
  filter(Cancer_Type == 'Breast Cancer' | Cancer_Type == 'Prostate Cancer' | Cancer_Type == 'Ovarian Cancer')

#create violin plot to look at the association of mutational events and signature 3 percentage in OBP

signature_3_tumor_asso_plot_OBP <- ggplot(combine_annotate_decomp_sig_3_percent_OBP, aes(x = Tumor_type_based_mutation_pattern, y = mean)) + 
  geom_violin() + 
  theme_bw()+
  theme (axis.text.x = element_text(face = 'bold', size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none') +
  xlab('Mutation Events') +
  ylab('Signature 3 Percentages') +
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter(aes(x = Tumor_type_based_mutation_pattern, y = mean, color = Tumor_type_based_mutation_pattern), width = 0.15)+
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ylim(0,1)

#statistical comparison using wilcox
signature_3_tumor_asso_plot_OBP + stat_compare_means(label = 'p.signif', label.y = 1, method = 'wilcox.test', method.args = list(alternative = 'greater'), ref.group = 'BRCA\nWildtype\nTumor', fontface = 'bold', size = 5)

ggsave('~/shenz/script/R_plot/BRCA_Signature_3_Plot/Ovarian_Breast_Prostate/Mutational Events and Signature 3 Percentage Association Plot in Ovarian, Prostate and Breast Cancer.tiff', width = 8, height = 5)

print('Finish Mutational Events and Signature 3 Percentage Association Plot in OBP')

#---------------------------------------------------------------------------------------------------------------------
# Part 14 plot the corelation between mutational events and NtAI score in prostate, ovarian and breast cancer

anno_combine_NtAI_OBP <- join(anno_combine_NtAI, cancer_type_master, by = 'sample', type = 'inner') %>%
  filter(Cancer_Type == 'Breast Cancer' | Cancer_Type == 'Prostate Cancer' | Cancer_Type == 'Ovarian Cancer')

#plot the correlation between NtAI score and mutational events in OBP
mutation_events_NtAI_asso_plot_OBP <- ggplot(anno_combine_NtAI_OBP, aes(x = Tumor_type_based_mutation_pattern, y = NtAI))+
  geom_violin()+ 
  theme_bw()+
  theme(axis.text.x = element_text(face = 'bold', size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = 'none')+
  labs(x = 'Mutational Events', y = 'NTAI Score')+
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter(aes(x = Tumor_type_based_mutation_pattern, y = NtAI, color = Tumor_type_based_mutation_pattern), width = 0.15)+
  geom_boxplot(width = 0.15, outlier.shape = NA)+
  ylim(0, 40)

#calculate significance(P value) using wilcox compare to BRCA Wildtype Tumor
mutation_events_NtAI_asso_plot_OBP + stat_compare_means(label = 'p.signif', label.y = 39, method = 'wilcox.test', method.args = list(alternative = 'greater'), ref.group = 'BRCA\nWildtype\nTumor', fontface = 'bold', size = 5)

ggsave('~/shenz/script/R_plot/BRCA_Genomic_Scar/Ovarian_Breast_Prostate/Mutational Events and NtAI Score Association Plot in Ovarian, Prostate and Breast Cancer.tiff', width = 8, height = 5)

print('Finish Mutational Events and NtAI Score Association Plot in OBP')

#---------------------------------------------------------------------------------------------------------------------
# Part 15 plot the corelation between mutational events and LST score in prostate, ovarian and breast cancer

anno_combine_LST_OPB <- join(anno_combine_LST, cancer_type_master, by = 'sample', type = 'inner') %>%
  filter(Cancer_Type == 'Breast Cancer' | Cancer_Type == 'Prostate Cancer' | Cancer_Type == 'Ovarian Cancer')

#plot the correlation between LST score and mutational events in OBP
mutation_events_LST_asso_plot_OBP <- ggplot(anno_combine_LST_OPB, aes(x = Tumor_type_based_mutation_pattern, y = LST))+
  geom_violin()+ 
  theme_bw()+
  theme(axis.text.x = element_text(face = 'bold', size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = 'none')+
  labs(x = 'Mutational Events', y = 'LST Score')+
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter(aes(x = Tumor_type_based_mutation_pattern, y = LST, color = Tumor_type_based_mutation_pattern), width = 0.15)+
  geom_boxplot(width = 0.15, outlier.shape = NA)+
  ylim(0, 50)

#calculate significance(P value) using wilcox compare to BRCA Wildtype Tumor
mutation_events_LST_asso_plot_OBP + stat_compare_means(label = 'p.signif', label.y = 50, method = 'wilcox.test', method.args = list(alternative = 'greater'), ref.group = 'BRCA\nWildtype\nTumor', fontface = 'bold', size = 5)

ggsave('~/shenz/script/R_plot/BRCA_Genomic_Scar/Ovarian_Breast_Prostate/Mutational Events and LST Score Association Plot in Ovarian, Prostate and Breast Cancer.tiff', width = 8, height = 5)

print('Finish Mutational Events and LST Score Association Plot in OBP')

#---------------------------------------------------------------------------------------------------------------------
# Part 16 plot the corelation between mutational events and HRD-LOH score in prostate, ovarian and breast cancer

anno_combine_HRD_OBP <- join(anno_combine_HRD, cancer_type_master, by = 'sample', type = 'inner') %>%
  filter(Cancer_Type == 'Breast Cancer' | Cancer_Type == 'Prostate Cancer' | Cancer_Type == 'Ovarian Cancer')

#plot the correlation between HRD-LOH score and mutational events in OBP
mutation_events_HRD_asso_plot_OBP <- ggplot(anno_combine_HRD_OBP, aes(x = Tumor_type_based_mutation_pattern, y = HRD_LOH))+
  geom_violin()+ 
  theme_bw()+
  theme(axis.text.x = element_text(face = 'bold', size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = 'none')+
  labs(x = 'Mutational Events', y = 'HRD-LOH Score')+
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter(aes(x = Tumor_type_based_mutation_pattern, y = HRD_LOH, color = Tumor_type_based_mutation_pattern), width = 0.15)+
  geom_boxplot(width = 0.15, outlier.shape = NA)+
  ylim(0, 38)

#calculate significance(P value) using wilcox compare to BRCA Wildtype Tumor
mutation_events_HRD_asso_plot_OBP + stat_compare_means(label = 'p.signif', label.y = 38, method = 'wilcox.test', method.args = list(alternative = 'greater'), ref.group = 'BRCA\nWildtype\nTumor', fontface = 'bold', size = 5)

ggsave('~/shenz/script/R_plot/BRCA_Genomic_Scar/Ovarian_Breast_Prostate/Mutational Events and HRD-LOH Score Association Plot in Ovarian, Prostate and Breast Cancer.tiff', width = 8, height = 5)

print('Finish Mutational Events and HRD-LOH Score Association Plot in OBP')

#------------------------------------------------------------------------------------------------------------
# Part17 HRDetection Score (both the whole patient cohort and OBP selected patient)
#combine all previous parameters to the annotated table ready to do HRD detection

HRD_detect_tbl <- combine_annotate_decomp_sig_3_percent  

names(HRD_detect_tbl)[7] <- 'signature_3_percentage'

#add Cancer Type info
HRD_detect_tbl <- join(HRD_detect_tbl, cancer_type_master, by = 'sample', type = 'inner')

# Create a new dataframe that extract signature 8 info from decomp dataframe
decomp_sig_8 <- decomp %>% filter (Signature == 'Signature.8') %>% select ('sample' = Tumor_Sample_Barcode,'signature_8_percentage' = mean)

# Create a new dataframe that extract signature 8 info from WT decomp dataframe
WT_decomp_sig_8 <- WT_decomp %>% filter (Signature == 'Signature.8') %>% select ('sample' = Tumor_Sample_Barcode, 'signature_8_percentage' = mean)

#combine total signature 3 df
total_decomp_sig_8 <- rbind(decomp_sig_8, WT_decomp_sig_8)

#add signature 8 to the table
HRD_detect_tbl <- join(HRD_detect_tbl, total_decomp_sig_8, by = 'sample', type = 'inner')

#add microhomology info
HRD_detect_tbl <- join(HRD_detect_tbl, select(microhomo_anno_comb, sample, Microhomology_Percent), by = 'sample', type = 'inner')

#add NtAI info
HRD_detect_tbl <- join(HRD_detect_tbl, total_NtAI, by = 'sample', type = 'inner')

#add LST info
HRD_detect_tbl <- join(HRD_detect_tbl, total_LST, by = 'sample', type = 'inner')

#add HRD-LOH info
HRD_detect_tbl <- join(HRD_detect_tbl, total_HRD_LOH, by = 'sample', type = 'inner')

#add total Indel info
HRD_detect_tbl <- join(HRD_detect_tbl, select(anno_tot_INDELS_combine, sample, Total_Indel_log2), by = 'sample', type = 'inner')

#add Mutation burden info
HRD_detect_tbl <- join(HRD_detect_tbl, select(mutation_burden_combine, sample, Total_Mutational_Burdens), by = 'sample', type = 'inner')

HRD_detect_tbl_OBP <- HRD_detect_tbl %>% filter(Cancer_Type == 'Prostate Cancer'|Cancer_Type =='Breast Cancer'|Cancer_Type =='Ovarian Cancer')

print('Finish organizing the HRD_detect table')
#----------------------------------------------------------------------------------------------------------------------------
#Utilize HRDetect algorithm to give a better representation of HRD in this patient cohort

set.seed(12345)
#Sample BRCA positive Tumor samples 
brca_pos = filter(HRD_detect_tbl, Tumor_type_based_mutation_pattern == 'Biallelic\nInactivation') %>% 
  sample_n(22) %>% 
  select(sample) %>% 
  mutate(brca = 1)

#Sample BRCA negative Tumor samples
brca_neg = filter(HRD_detect_tbl, Tumor_type_based_mutation_pattern == 'BRCA\nWildtype\nTumor') %>% 
  sample_n(100) %>% 
  select(sample) %>% 
  mutate(brca = 0)

#combine two training set and scale them
train_set = bind_rows(brca_pos,
                      brca_neg) %>%
  left_join(., select(HRD_detect_tbl, sample, Total_Mutational_Burdens, Total_Indel_log2, signature_3_percentage, signature_8_percentage, LST, HRD_LOH, NtAI, Microhomology_Percent)) %>% 
  mutate_at(vars(matches('Total_Mutational_Burdens|Total_Indel_log2|signature_3_percentage|signature_8_percentage|LST|HRD_LOH|NtAI|Microhomology_Percent')), funs(scale(.)))

#Start Training
y = as.factor(train_set$brca)
x = as.matrix(select(train_set, matches('Total_Mutational_Burdens|Total_Indel_log2|signature_3_percentage|signature_8_percentage|LST|HRD_LOH|NtAI|Microhomology_Percent')))
glmmod = glmnet(x, y, alpha = 1, family = 'binomial', lower.limits = 0) # alpha --> lasso only, no ridge penalty
plot(glmmod, xvar = 'lambda', label = T)
cv_glmmod = cv.glmnet(x, y , alpha = 1, family = 'binomial', nfolds = 10)
plot(cv_glmmod, label = T)
lmin = cv_glmmod$lambda.min

#Start testing
test_set = filter(HRD_detect_tbl, !(sample %in% train_set$sample)) %>%
  select(sample, Tumor_type_based_mutation_pattern, Detailed_mutation_pattern, Total_Mutational_Burdens, Total_Indel_log2, signature_3_percentage, signature_8_percentage, LST, HRD_LOH, NtAI, Microhomology_Percent) %>% 
  mutate_at(vars(matches('Total_Mutational_Burdens|Total_Indel_log2|signature_3_percentage|signature_8_percentage|LST|HRD_LOH|NtAI|Microhomology_Percent')), funs(scale(.)))

newx = as.matrix(select(test_set, matches('Total_Mutational_Burdens|Total_Indel_log2|signature_3_percentage|signature_8_percentage|LST|HRD_LOH|NtAI|Microhomology_Percent')))
test_pred = predict(cv_glmmod, newx, s = lmin, type = 'response')
colnames(test_pred) = 'HRD_detection_score'

Final_set <- cbind(test_set, test_pred) %>% arrange(HRD_detection_score) 

Final_set$sample <- factor(Final_set$sample, levels = c(Final_set$sample))

# Plot the density plot to show the HRDetction score
HRD_detect_mutational_events_asso <- ggplot(Final_set, aes(x = sample, y = HRD_detection_score))+
  geom_bar(stat = 'identity', aes(color = Tumor_type_based_mutation_pattern, fill = Tumor_type_based_mutation_pattern))+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = 'none')+
  ylim(0,0.9)+
  labs(x = '411 Cancer Samples', y = 'HRDetect Scores')+
  theme(axis.title = element_text(size = 12, face = 'bold'))+
  theme(axis.text.y = element_text(size = 12, face = 'bold'))+
  scale_color_manual(values = c('gray', 'yellow', 'blue', 'brown', 'purple', 'black'))+
  scale_fill_manual(values = c('gray', 'yellow', 'blue', 'brown', 'purple', 'black'))


ggsave('~/shenz/script/R_plot/BRCA_HRDetect/Whole/BRCA Mutational Events HRDetction Association Plot.tiff', width = 10, height = 7)

print('Finish Mutational Events and HRDetection Association Plot')

#----------------------------------------------------------------------------------------------------------------------------
#Utilize HRDetect algorithm to give a better representation of HRD in this patient cohort of Ovarian, Breast and Prostate cancer

#Start testing
test_set_OBP = filter(HRD_detect_tbl_OBP, !(sample %in% train_set$sample)) %>%
  select(sample, Tumor_type_based_mutation_pattern, Detailed_mutation_pattern, Total_Mutational_Burdens, Total_Indel_log2, signature_3_percentage, signature_8_percentage, LST, HRD_LOH, NtAI, Microhomology_Percent) %>% 
  mutate_at(vars(matches('Total_Mutational_Burdens|Total_Indel_log2|signature_3_percentage|signature_8_percentage|LST|HRD_LOH|NtAI|Microhomology_Percent')), funs(scale(.)))

newx = as.matrix(select(test_set_OBP, matches('Total_Mutational_Burdens|Total_Indel_log2|signature_3_percentage|signature_8_percentage|LST|HRD_LOH|NtAI|Microhomology_Percent')))
test_pred_OBP = predict(cv_glmmod, newx, s = lmin, type = 'response')
colnames(test_pred_OBP) = 'HRD_detection_score'

Final_set_OBP <- cbind(test_set_OBP, test_pred_OBP) %>% arrange(HRD_detection_score) 

Final_set_OBP$sample <- factor(Final_set_OBP$sample, levels = c(Final_set_OBP$sample))

# Plot the density plot to show the HRDetction score
HRD_detect_mutational_events_asso <- ggplot(Final_set_OBP, aes(x = sample, y = HRD_detection_score))+
  geom_bar(stat = 'identity', aes(fill = Tumor_type_based_mutation_pattern, color = Tumor_type_based_mutation_pattern))+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = 'none')+
  ylim(0,0.4)+
  labs(x = '146 Breast, Ovarian and Prostate Cancer Samples ', y = 'HRDetect Scores')+
  theme(axis.title = element_text(size = 12, face = 'bold'))+
  theme(axis.text.y = element_text(size = 12, face = 'bold'))+
  scale_fill_manual(values = c('gray', 'yellow', 'blue', 'brown', 'purple', 'black'))+
  scale_color_manual(values = c('gray', 'yellow', 'blue', 'brown', 'purple', 'black'))

ggsave('~/shenz/script/R_plot/BRCA_HRDetect/Ovarian_Breast_Prostate/BRCA Mutational Events HRDetction Association Plot in Breast, Prostate and Ovarian Cancer.tiff', width = 10, height = 7)

print('Finish Mutational Events and HRDetection Association Plot OBP')

#--------------------------------------------------------------------------------------------------------------------------
# Part 17 ONLY subset the Mono somatic and germline inactivation samples from the HRD-detect_tbl for philip

mono_subset_tbl <- HRD_detect_tbl %>% filter (grepl('Monoallelic', Tumor_type_based_mutation_pattern))

#save as a text file
write.table(mono_subset_tbl, file = '~/shenz/script/R_file/Monoallelic Germline and Somatic Subset.txt', sep = '\t')

#calculate the ruuning time of this program
print(proc.time()-starttime)
#---------------------------------------------------------------------------------------------------------------------------
# Part 18 start to plot the association to signature 3 (differentiate monoallelic Germline/Somatic and biallelic somatic/germline)
Further_combine_annotate_sig_3 <- combine_annotate_sig_3 %>% mutate_if(is.factor, as.character) %>%
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Biallelic\nInactivation' & grepl('Germline', Detailed_mutation_pattern), 'Germline\nBiallelic\nInactivation', Tumor_type_based_mutation_pattern)) %>%
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Biallelic\nInactivation' & grepl('Somatic', Detailed_mutation_pattern), 'Somatic\nBiallelic\nInactivation', Tumor_type_based_mutation_pattern)) %>%
  mutate(Tumor_type_based_mutation_pattern = ifelse(Tumor_type_based_mutation_pattern == 'Biallelic\nInactivation' & grepl('Homodels', Detailed_mutation_pattern), 'Homozygous\nDeletion', Tumor_type_based_mutation_pattern)) %>%
  mutate(Tumor_type_based_mutation_pattern = factor(Tumor_type_based_mutation_pattern, levels = c("BRCA\nWildtype\nTumor", 'Germline\nBiallelic\nInactivation', 'Somatic\nBiallelic\nInactivation', 'Homozygous\nDeletion', "Monoallelic\nGermline\nInactivation", "Monoallelic\nSomatic\nInactivation", "Germline\nDeletion", "Subclonal\nBRCA\nInactivation")))%>%
  arrange(Tumor_type_based_mutation_pattern)


#Make a violin plot to show the difference
mutation_events_sig3_asso_plot_detailed <- ggplot(Further_combine_annotate_sig_3, aes(x = Tumor_type_based_mutation_pattern, y = signature_3_activity))+
  geom_violin()+
  theme_bw()+
  theme(axis.text.x = element_text(face = 'bold', size = 10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none')+
  labs(x = 'Mutational Events', y = 'Signature 3 Activity' )+
  theme(axis.title = element_text(face = 'bold'))+
  geom_jitter(aes(x = Tumor_type_based_mutation_pattern, y = signature_3_activity, color = Tumor_type_based_mutation_pattern), width = 0.15)+
  geom_boxplot(width = 0.15, outlier.shape = NA)+
  ylim(0,480)

#calculate significance(P value) among Biallelic Inactivation, Monoallelic Germline and Monoallelic Somatic via Wilcox-Rank-Sum test. Add on
#to the plot
mutation_events_sig3_asso_plot_detailed + stat_compare_means(label = 'p.signif', label.y = 480, ref.group = 'BRCA\nWildtype\nTumor', method = 'wilcox.test', method.args = list(alternative = 'greater'),fontface = 'bold', size = 5)

ggsave('~/shenz/script/R_plot/BRCA_Signature_3_Plot/Other_Plots/BRCA Mutational Events (Detailed) Signature 3 Association Plot.tiff', width = 8, height = 5)

print('Finish Mutational Events and Signature 3 Acitivity Association Plot (Detailed Version)')

#-----------------------------------------------------------------------------------------------------------------------------
#generate spearman correlation plot of NTAI against LST score

#Add NtAI, LST and HRD-LOH score into Annotated table to combine together
HRD_Score_Combine_df <-join(Annotated_merged_table, total_NtAI, by = 'sample', type = 'inner') 

HRD_Score_Combine_df <- join(HRD_Score_Combine_df, total_LST, by = 'sample', type = 'inner') 

HRD_Score_Combine_df <- join(HRD_Score_Combine_df, total_HRD_LOH, by = 'sample', type = 'inner') 

#Generate the plot
NTAI_LST_Correlation_Plot <- ggplot(HRD_Score_Combine_df, aes(x = NtAI, y = LST)) +
  geom_point(color = 'blue', size = 3)+
  geom_smooth(method = lm, color = 'red', size = 2, fill = 'red')+
  theme_bw()+
  theme(axis.title = element_text(size = 14, face = 'bold'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none', axis.text = element_text(face = 'bold'))+
  labs(x = 'NtAI Score', y = 'LST Score\n')

#use spearman correlation method to look at both the correlation coefficient (rho) and P value, and add to the graph

NTAI_LST_Correlation_Plot + stat_cor(method = 'spearman', label.x = 30, label.y = 6, size = 6, label.sep = '\n', fontface = 'bold')

ggsave('~/shenz/script/R_plot/BRCA_Genomic_Scar/correlation_plot/NtAI and LST score correlation plot.tiff', width = 8, height = 5)

print('Finish NtAI and LST correlation plot')

#-----------------------------------------------------------------------------------------------------------------------------
#generate spearman correlation plot of NTAI against HRD-LOH score

#Generate the plot
NTAI_HRD_Correlation_Plot <- ggplot(HRD_Score_Combine_df, aes(x = NtAI, y = HRD_LOH)) +
  geom_point(color = 'blue', size = 3)+
  geom_smooth(method = lm, color = 'red', size = 2, fill = 'red')+
  theme_bw()+
  theme(axis.title = element_text(size = 14, face = 'bold'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none', axis.text = element_text(face = 'bold'))+
  labs(x = 'NtAI Score', y = 'HRD-LOH Score\n')

#use spearman correlation method to look at both the correlation coefficient (rho) and P value, and add to the graph

NTAI_HRD_Correlation_Plot + stat_cor(method = 'spearman', label.x = 30, label.y = 6, size = 6, label.sep = '\n', fontface = 'bold')

ggsave('~/shenz/script/R_plot/BRCA_Genomic_Scar/correlation_plot/NtAI and HRD-LOH score correlation plot.tiff', width = 8, height = 5)

print('Finish NtAI and HRD-LOH correlation plot')



