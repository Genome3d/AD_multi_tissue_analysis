################################################################################
#Creation of master dataframe of MR results from nine tissues
################################################################################

#load libraries
library(tidyverse)
library(TwoSampleMR)

#import MR result datasets

#adult cortex
adult_cortex_sig_MR_results = read_csv("adult_cortex_sig_MR_results.csv")
adult_cortex_sig_MR_results$tissue='Adult Cortex'

#aorta
aorta_sig_MR_results = read.csv("aorta_sig_MR_results.csv")
aorta_sig_MR_results$tissue = 'Artery - Aorta'

#coronary
coronary_sig_MR_results=read_csv("coronary_sig_MR_results.csv")
coronary_sig_MR_results$tissue='Artery - Coronary'

#fetal cortex
fetal_cortex_sig_MR_results = read.csv("fetal_cortex_sig_MR_results.csv")
fetal_cortex_sig_MR_results$tissue  = "Fetal Cortex"

#liver
liver_sig_MR_results = read.csv("liver_sig_MR_results.csv")
liver_sig_MR_results$tissue = 'Liver'

#lung
lung_sig_MR_results = read.csv("lung_sig_MR_results.csv")
lung_sig_MR_results$tissue = 'Lung'

#skin exposed
skin_exposed_sig_MR_results = read.csv("skin_exposed_sig_MR_results.csv")
skin_exposed_sig_MR_results$tissue = 'Skin - Exposed'

#skin unexposed
skin_unexposed_sig_MR_results = read.csv("skin_unexposed_sig_MR_results.csv")
skin_unexposed_sig_MR_results$tissue = 'Skin - Unexposed'

#whole blood
whole_blood_sig_MR_results = read.csv("blood_sig_MR_results.csv")
whole_blood_sig_MR_results$tissue = 'Whole Blood'

#join datasets together
LOAD_2SMR_data=list(adult_cortex_sig_MR_results,
                    aorta_sig_MR_results,
                    coronary_sig_MR_results,
                    fetal_cortex_sig_MR_results,
                    liver_sig_MR_results,
                    lung_sig_MR_results,
                    skin_exposed_sig_MR_results,
                    skin_unexposed_sig_MR_results,
                    whole_blood_sig_MR_results)%>%
  reduce(full_join)

#generate odds ratios from beta values (TwoSampleMR function)
LOAD_2SMR_data=generate_odds_ratios(LOAD_2SMR_data)

#there are 185 unique risk genes across the 9 tissues 
#(further filtering in script below)
length(unique(LOAD_2SMR_data$exposure))

################################################################################
#SNP and regulation (log(aFC)) info dataset

#getting SNP ID and gene regulation info prior to MR 

adult_cortex_harmonised = read.csv("adult_cortex_harmonised.csv")
adult_cortex_harmonised$tissue='Adult Cortex'

aorta_harmonised = read.csv("aorta_harmonised.csv")
aorta_harmonised$tissue='Artery - Aorta'

blood_harmonised = read.csv("blood_harmonised.csv")
blood_harmonised$tissue='Whole Blood'

coronary_harmonised = read.csv("coronary_harmonised.csv")
coronary_harmonised$tissue='Artery - Coronary'

fetal_cortex_harmonised = read.csv("fetal_cortex_harmonised.csv")
fetal_cortex_harmonised$tissue='Fetal Cortex'

liver_harmonised = read.csv("liver_harmonised.csv")
liver_harmonised$tissue='Liver'

lung_harmonised = read.csv("lung_harmonised.csv")
lung_harmonised$tissue='Lung'

skin_exposed_harmonised = read.csv("skin_exposed_harmonised.csv")
skin_exposed_harmonised$tissue='Skin - Exposed'

skin_unexposed_harmonised = read.csv("skin_unexposed_harmonised.csv")
skin_unexposed_harmonised$tissue='Skin - Unexposed'

#join datasets together, keep selected columns
tissue_SNP_list=list(aorta_harmonised,
                     coronary_harmonised,
                     fetal_cortex_harmonised,
                     adult_cortex_harmonised,
                     lung_harmonised,
                     liver_harmonised,
                     skin_unexposed_harmonised,
                     skin_exposed_harmonised,
                     blood_harmonised) %>%
  reduce(full_join) %>%
  dplyr::select(id.exposure,SNP,exposure,tissue,beta.exposure,se.exposure,pval.exposure)

#filter to keep SNPs only significant following 2SMR
#keep rows that match up with significant id exposure
filtered_SNP_df = filter(tissue_SNP_list, id.exposure %in% LOAD_2SMR_data$id.exposure) 


write.csv(filtered_SNP_df, "LOAD_SNP_data.csv")

################################################################################
#Filter MR dataset further to exclude one SNP:multiple exposure in a single tissue
#as this violates independence assumption
filtered_SNP_df = filtered_SNP_df%>%
  mutate(SNP_tissue = paste(SNP, tissue, sep = "-"))

#shows which snps occur more than once in a single tissue 
snp_exposure_freq_table = data.frame(table(filtered_SNP_df$SNP,filtered_SNP_df$tissue)) %>%
  filter(Freq > 1) %>%
  rename("SNP" = "Var1",
         "tissue" = "Var2") %>%
  mutate(SNP_tissue = paste(SNP, tissue, sep = "-"))

#13 eQTLs regulate more than one exposure in each tissue --> exclusion
filtered_SNP_exclusion = filter(filtered_SNP_df, SNP_tissue %in% snp_exposure_freq_table$SNP_tissue) %>%
  mutate(exposure_tissue = paste(exposure, tissue, sep = "-"))

#filter tissue comparison dataset by removing exposure-tissue combos which have multiple 
#exposures for the same SNP --> removing 26 observations

LOAD_2SMR_data = LOAD_2SMR_data%>%
  mutate(exposure_tissue=paste(exposure, tissue, sep = "-"))

LOAD_2SMR_data = filter(LOAD_2SMR_data, !exposure_tissue %in% filtered_SNP_exclusion$exposure_tissue)

length(unique(LOAD_2SMR_data$exposure))

#removes 6 genes entirely from analysis

LOAD_2SMR_data = LOAD_2SMR_data[,-16]

write.csv(LOAD_2SMR_data,"LOAD_2SMR_data.csv")
################################################################################
#joint dataset MR and aFC for switch gene analysis

tissue_comparison_switch.final = LOAD_2SMR_data %>%
  group_by(exposure) %>%
  filter(any(or > 1) & any(or < 1)) %>%
  ungroup() %>%
  dplyr::select(exposure,tissue,or,or_uci95, or_lci95, pval)

#subset of this for switch genes
filtered_SNP_switch_df = filter(filtered_SNP_df, exposure %in% tissue_comparison_switch.final$exposure)
#join with switch data for full info
final_filtered_switch=full_join(filtered_SNP_switch_df,tissue_comparison_switch.final)%>%
  na.omit() %>%
  mutate(aFC.lci95 = beta.exposure - 1.96 * se.exposure,
         aFC.uci95 = beta.exposure + 1.96 * se.exposure,
         exposure.snp=paste(exposure,SNP,sep = " --> "),
         aFC_rank = factor(exposure.snp, levels = unique(exposure.snp[order(-or)])))

write.csv(final_filtered_switch, "switch_genes_or_afc_data.csv")
