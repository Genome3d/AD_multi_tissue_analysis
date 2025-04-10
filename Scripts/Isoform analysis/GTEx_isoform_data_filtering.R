################################################################################
#integration of dataset to match GTEx median transcription expression
#for 8 tissues (no fetal expression data)
################################################################################

#install.packages("vroom") 
#vroom for fast uploading of datasets
library(vroom)
library(tidyverse)
library(stringr)


#attributes (connect annotations with tissues) - sample ID with tissues
GTex_sample_attributes=vroom("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt") %>%
  filter(SMTSD =="Whole Blood" | 
           SMTSD =="Artery - Coronary" |
           SMTSD =="Artery - Aorta" | 
           SMTSD =="Brain - Cortex" | 
           SMTSD =="Skin - Not Sun Exposed (Suprapubic)" | 
           SMTSD =="Skin - Sun Exposed (Lower leg)" | 
           SMTSD =="Liver" | 
           SMTSD =="Lung" )

#gene names to link to gene ids 
GTEx_gene_info=vroom("gencode.v26.GRCh38.genes.gtf") %>%
  filter(`(GRCh38),` == "AC092066.1;" | 
           `(GRCh38),` == "ACE;" | 
           `(GRCh38),` == "APOC2;" |
           `(GRCh38),` == "BIN1;" |
           `(GRCh38),` == "EARS2;" |
           `(GRCh38),` == "GAL3ST4;" |
           `(GRCh38),` == "ICA1L;" |
           `(GRCh38),` == "LINC02210;" |
           `(GRCh38),` == "PLEKHM1;" |
           `(GRCh38),` == "PRSS36;" |
           `(GRCh38),` == "KAT8;")

#get rid of semi colons and commas from character strings, rename to gene name
GTEx_gene_info=GTEx_gene_info %>%
  mutate(Gene_name=str_sub(GTEx_gene_info$`(GRCh38),`, end = -2),
         gene_id=str_sub(GTEx_gene_info$`evidence-based`, end = -2))

#set system environment to larger to upload TPM info
Sys.setenv(VROOM_CONNECTION_SIZE=500000072)

#TPM transcript info
GTEx_TPMs_old_version=vroom("GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz", skip=2)

GTEx_gene_info %>%
  select(gene_id) %>%
  distinct()
  
GTEx_TPMs_old_version %>%
  select(gene_id) %>%
  distinct()

#filter dataset by linking gene IDs connected to gene names with gene IDs in TPM dataset
GTex_TPM_subset =  filter(GTEx_TPMs_old_version, gene_id %in% GTEx_gene_info$`gene_id`)

write.csv(GTex_TPM_subset, "GTex_TPM_subset.csv")

#pivot longer to make sample ids a single column
GTex_TPM_subset_pivot=pivot_longer(GTex_TPM_subset, cols= c(3:17384), names_to = "SAMPID", values_to = "TPM")

#only keep sample IDs for specified tissues
GTex_TPM_subset_pivot_sub =  filter(GTex_TPM_subset_pivot, SAMPID %in% GTex_sample_attributes$SAMPID)

#merge this dataset with TPM values with dataset with tissues, only keep relevent columns
GTex_TPM_tissue=merge(GTex_TPM_subset_pivot_sub, GTex_sample_attributes, by="SAMPID")%>%
  select(c("SAMPID", "transcript_id", "gene_id" , "TPM" , "SMTSD"))

#merge dataset with gene info to get gene names
GTex_TPM_tissue_gene=merge(GTex_TPM_tissue,GTEx_gene_info, by = 'gene_id') 

#select relevent columns, distinct gets rid of duplicated rows, make transcript id and tissue factors for easy summaries
GTex_TPM_tissue_gene=GTex_TPM_tissue_gene%>%
  select(c("Gene_name", "SAMPID", "transcript_id", "gene_id" , "TPM" , "SMTSD"))%>%
  distinct()%>%
  mutate(transcript_id=as.factor(transcript_id),
         tissue=as.factor(SMTSD))

GTex_TPM_tissue_gene %>%
  select(Gene_name) %>%
  distinct()

#write to csv for ease of access
write.csv(GTex_TPM_tissue_gene, "GTEx_TPM_tissue_gene.csv")
