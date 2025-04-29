################################################################################
#Transcript isoform analysis script
#Zillah Daysh
################################################################################
#Load libraries
library(tidyverse)
library(rstatix)
library(ggpubr)
library(reshape2)
library(broom)
library(FSA)
library(forcats)

################################################################################
#Read in datasets
################################################################################

#transcript expression info - dataset from GTEx_isoform_data_filtering.R script
GTex_TPM_tissue_gene = read.csv("GTEx_TPM_tissue_gene.csv") 

#extract donor ID information
GTex_TPM_tissue_gene = GTex_TPM_tissue_gene %>%
  mutate(SAMPID = substr(SAMPID,1,10),
         SAMPID = sub("-$", "", SAMPID))

#2SMR data for all nine tissues
AD_2SMR_data = read.csv("AD_2SMR_data.csv")[,-c(1:2,4:5,7:8,18)]

################################################################################
#Data wrangling for analysis
################################################################################

#Subset data for relevant columns
TPM_subset=GTex_TPM_tissue_gene  %>%
  select(c("Gene_name", "SAMPID" , "transcript_id", "tissue", "TPM")) 

#rename tissues
TPM_subset$tissue = fct_recode(TPM_subset$tissue,
                              "Adult Cortex" = "Brain - Cortex",
                              "Skin - Unexposed" = "Skin - Not Sun Exposed (Suprapubic)" , 
                              "Skin - Exposed" = "Skin - Sun Exposed (Lower leg)")

#filter to only keep tissues for each gene which were significant in 2SMR
TPM_tissue_gene = TPM_subset %>%
  inner_join(AD_2SMR_data, by = c("Gene_name" = "exposure", "tissue" = "tissue")) %>%
  mutate(transcript_id = substr(transcript_id, 1, 15),
         AD_risk = ifelse(or > 1 , "Risk", "Protection"))

#remove switch gene transcripts whose medians are absent in all tissues
remove_transcripts = TPM_tissue_gene %>%
  group_by(tissue, transcript_id, Gene_name) %>%
  summarise(`Median TPM` = round (median(TPM), 3)) %>%
  ungroup() %>%
  group_by(transcript_id) %>%
  summarise(absent = all(`Median TPM` ==0))%>%
  filter(absent) %>%
  pull(transcript_id)

#create new columns for total and percent TPM for each transcript for each donor ID 
TPM_tissue_gene_filtered = filter(TPM_tissue_gene, !transcript_id %in% remove_transcripts) %>%
  group_by(SAMPID,tissue,Gene_name) %>%
  mutate(percent_TPM = (`TPM`/sum(`TPM`)),
         total_TPM = sum(`TPM`)) %>%
  ungroup()

#change - to : for tissue codes 
#needed for for loop below to separate pairwise comparison groups
TPM_tissue_gene_filtered$tissue = fct_recode(TPM_tissue_gene_filtered$tissue,
                                       "Skin : Unexposed" = "Skin - Unexposed" , 
                                       "Skin : Exposed" = "Skin - Exposed",
                                       "Artery : Aorta" = "Artery - Aorta",
                                       "Artery : Coronary" = "Artery - Coronary")

#create transcript name column - corresponds to ensembl transcript names linked to transcript IDs
TPM_tissue_gene_filtered$transcript_name = TPM_tissue_gene_filtered$transcript_id %>%
  fct_recode("ACE-201" = "ENST00000290863",
             "ACE-202" = "ENST00000290866",
             "ACE-204" = "ENST00000428043",
             "ACE-206" = "ENST00000578679",
             "ACE-209" = "ENST00000579314",
             "ACE-210" = "ENST00000579462",
             "ACE-211" = "ENST00000579726",
             "ACE-212" = "ENST00000580318",
             "ACE-214" = "ENST00000582244",
             "ACE-216" = "ENST00000582678",
             "ACE-219" = "ENST00000583645",
             "ACE-220" = "ENST00000584529",
             
             "APOC2-201" = "ENST00000252490",
             "APOC2-202" = "ENST00000585786",
             "APOC2-203" = "ENST00000590360",
             "APOC2-204" = "ENST00000591597",
             "APOC2-205" = "ENST00000592257",
             
             "BIN1-201" = "ENST00000259238",
             "BIN1-202" = "ENST00000316724",
             "BIN1-204" = "ENST00000348750",
             "BIN1-205" = "ENST00000351659",
             "BIN1-208" = "ENST00000376113",
             "BIN1-209" = "ENST00000393040",
             "BIN1-210" = "ENST00000393041",
             "BIN1-211" = "ENST00000409400",
             "BIN1-212" = "ENST00000462958",
             "BIN1-213" = "ENST00000466111",
             "BIN1-214" = "ENST00000484253",
             
             "EARS2-201" = "ENST00000449606",
             "EARS2-205" = "ENST00000562799",
             "EARS2-207" = "ENST00000563459",
             "EARS2-210" = "ENST00000564668",
             "EARS2-211" = "ENST00000564759",
             "EARS2-213" = "ENST00000564987",
             "EARS2-215" = "ENST00000565344",
             
             "GAL3ST4-201" = "ENST00000360039",
             "GAL3ST4-203" = "ENST00000413800",
             "GAL3ST4-204" = "ENST00000423751",
             "GAL3ST4-206" = "ENST00000482469",
             "GAL3ST4-207" = "ENST00000495882",
             
             "ICA1L-201" = "ENST00000358299",
             "ICA1L-202" = "ENST00000392237",
             "ICA1L-206" = "ENST00000418208",
             "ICA1L-209" = "ENST00000421334",
             "ICA1L-210" = "ENST00000425178",
             "ICA1L-218" = "ENST00000457524",
             "ICA1L-219" = "ENST00000476602",
             "ICA1L-220" = "ENST00000484561",
             "ICA1L-221" = "ENST00000494560",
             "ICA1L-222" = "ENST00000617388",
             
             "PLEKHM1-201" = "ENST00000430334",
             "PLEKHM1-203" = "ENST00000579131",
             "PLEKHM1-204" = "ENST00000579197",
             "PLEKHM1-205" = "ENST00000580205",
             "PLEKHM1-206" = "ENST00000580404",
             "PLEKHM1-207" = "ENST00000581448",
             "PLEKHM1-210" = "ENST00000584420",
             "PLEKHM1-211" = "ENST00000585506",
             "PLEKHM1-212" = "ENST00000586084",
             "PLEKHM1-213" = "ENST00000586562",
             "PLEKHM1-215" = "ENST00000590991",
             "PLEKHM1-216" = "ENST00000591580",
             "PLEKHM1-217" = "ENST00000581932",
             
             "PRSS36-201" = "ENST00000268281",
             "PRSS36-202" = "ENST00000418068",
             "PRSS36-204" = "ENST00000562368",
             "PRSS36-205" = "ENST00000562390",
             "PRSS36-206" = "ENST00000563693",
             "PRSS36-209" = "ENST00000571878",
             
             "KAT8-201" = "ENST00000219797",
             "KAT8-202" = "ENST00000448516",
             "KAT8-203" = "ENST00000537402",
             "KAT8-204" = "ENST00000538768",
             "KAT8-205" = "ENST00000539683",
             "KAT8-206" = "ENST00000543774",
             "KAT8-207" = "ENST00000573144"
             
  ) 

#there is no current ensembl ID for this PLEKHM1 transcript (now archived transcript)
#so remove from dataset
TPM_tissue_gene_filtered = TPM_tissue_gene_filtered %>%
  filter(!transcript_id == "ENST00000582035")


################################################################################
#Kruskal-Wallis test to compare relative median abundances 
################################################################################
#Loop for each gene

#gene names
genes = unique(TPM_tissue_gene_filtered$Gene_name)

# Create an empty list to store the models and pairwise comparison results
models.kw = list()
merged.data.final_gene =  list()
merged.data.final_gene_risk = list()

# Loop over each gene
for (i in genes) {
  # Filter data for the current gene
  TPM_model_gene = TPM_tissue_gene_filtered %>% 
    filter(Gene_name == i) %>%
    mutate(tissue = factor(tissue, 
                           levels = c(
                             unique(tissue[AD_risk == "Protection"]), 
                             unique(tissue[AD_risk == "Risk"])
                           )))
  
  # Get unique transcripts for gene
  transcripts = unique(TPM_model_gene$transcript_name)
  
  # Create lists to store test results
  pairwise_results_for_gene = list()
  
  # Loop over each transcript
  for (j in transcripts) {
    # Filter data for the current transcript
    TPM_transcript_data = TPM_model_gene %>% 
      filter(transcript_name == j)
    
    unique_tissues = unique(TPM_transcript_data$tissue)
    
    if (length(unique_tissues) > 2) {
      # Fit the Kruskal-Wallis test
      kw.test = kruskal.test(percent_TPM ~ tissue, data = TPM_transcript_data)
      models.kw[[paste(i, j, "kw", sep = "_")]] = kw.test #save in list
      
      # Perform Dunn's test for pairwise comparisons
      dunn_pairwise = dunnTest(data = TPM_transcript_data, percent_TPM ~ tissue, method = "bh")
      #arrange, and get *** to show p val significance
      dunn_pairwise_df = dunn_pairwise$res %>%
        mutate(transcript_name = j) %>%
        rowwise() %>%
        mutate(group1 = str_trim(strsplit(Comparison, "-")[[1]][1]),
               group2 = str_trim(strsplit(Comparison, "-")[[1]][2]),
               p.signif = case_when(P.adj < 0.001 ~ "***",
                                    P.adj < 0.01 ~ "**", 
                                    P.adj < 0.05 ~ "*",
                                    P.adj > 0.05 ~ "ns")) %>%
        as_tibble()
    } else {
      # Perform Wilcoxon rank-sum test (Mann-Whitney U test) for two tissues (PLEKHM1)
      wilcox_test = wilcox.test(percent_TPM ~ tissue, data = TPM_transcript_data)
      #renamed object to dunn test even though it is not not
      #this is to enable the loop to work with PLEKHM1
      #arrange to be same as actual dunn test df above
      dunn_pairwise_df = tibble(
        Comparison = paste(unique_tissues, collapse = " - "),
        Z = wilcox_test$statistic,
        P.adj = wilcox_test$p.value,
        transcript_name = j,
        group1 = unique_tissues[1],
        group2 = unique_tissues[2],
        p.signif = case_when(wilcox_test$p.value < 0.001 ~ "***",
                             wilcox_test$p.value < 0.01 ~ "**",
                             wilcox_test$p.value < 0.05 ~ "*",
                             wilcox_test$p.value > 0.05 ~ "ns")
      )
    }
    #save in list
    pairwise_results_for_gene[[j]] = dunn_pairwise_df
  }
  
  # Combine pairwise transcript results for each gene
  merged_results = bind_rows(pairwise_results_for_gene)
  
  # Max TPM proportion values for plotting significance (just for location on plot)
  max_values = TPM_tissue_gene_filtered %>%
    group_by(Gene_name, transcript_name) %>%
    summarise(max_TPM = max(percent_TPM, na.rm = TRUE), .groups = 'drop')
  
  # Merge max values and calculate y.position for significance stars
  merged_data = merged_results %>%
    left_join(max_values, by = "transcript_name") %>%
    group_by(transcript_name) %>%
    mutate(y.position = max_TPM + 0.05 + (row_number() - 1) * 0.08) %>%
    ungroup()
  

  #for each switch gene, there is a single risk tissue and then numerous protective tissues or vice versa
  #so for plotting, just keeping the single risk tissue vs other tissues for main figure
  #to show sig differences between the risk tissue and ALL protective tissues (or reverse)
  #diff dataset to merged data which has all pairwise comparisons (supp figures)
  
  # Extract tissue and risk information for the current gene
  tissue_risk_info = TPM_tissue_gene_filtered %>%
    filter(Gene_name == i) %>%
    select(tissue, AD_risk) %>%
    distinct()
  
  # Determine tissue_filter based on the presence of "Protection" and "Risk" tissues
  #needed for below plotting to group tissues into "risk" and "protection" and keeping 
  #relevant significance bars
  num_protection = nrow(tissue_risk_info %>% filter(AD_risk == "Protection"))
  num_risk = nrow(tissue_risk_info %>% filter(AD_risk == "Risk"))
  
  if (num_protection == 1 && num_risk > 1) {
    tissue_filter = tissue_risk_info$tissue[tissue_risk_info$AD_risk == "Protection"]
  } else if (num_risk == 1 && num_protection > 1) {
    tissue_filter = tissue_risk_info$tissue[tissue_risk_info$AD_risk == "Risk"]
  } else {
    tissue_filter = tissue_risk_info$tissue
  }
  
  #filter pairwise significance data to be x risk/protective tissues vs the unique risk/protective tissue
  merged_data_risk = merged_data %>%
    filter(grepl(tissue_filter, Comparison)) %>%
    group_by(transcript_name) %>%
    mutate(y.position = max_TPM + 0.05 + (row_number() - 1) * 0.08) %>%
    filter(!any(p.signif == "ns")) %>%
    filter(P.adj == max(P.adj)) %>%
    ungroup() %>%
    mutate(group1 = recode(group1, !!!setNames(tissue_risk_info$AD_risk, tissue_risk_info$tissue)),
           group2 = recode(group2, !!!setNames(tissue_risk_info$AD_risk, tissue_risk_info$tissue)))
  
  # Store the results
  merged.data.final_gene_risk[[i]] = merged_data_risk
  merged.data.final_gene[[i]] = merged_data

}

################################################################################
#Following are plots in paper, change to other gene names for other plots
#(will need to manually change scale fill manual to line up with number of risk/protective tissues)
#number of risk/protective tissues easily found in TPM_tissue_gene_filtered using this code:
#TPM_tissue_gene_filtered %>%
  #select(Gene_name, tissue, AD_risk) %>%
  #distinct()

################################################################################
#BIN1
################################################################################
TPM_model_BIN1 = TPM_tissue_gene_filtered %>%
  filter(Gene_name == "BIN1")

TPM_model_BIN1$tissue = factor(TPM_model_BIN1$tissue, 
                                levels = c(unique(TPM_model_BIN1$tissue[TPM_model_BIN1$AD_risk == "Protection"]), 
                                           unique(TPM_model_BIN1$tissue[TPM_model_BIN1$AD_risk == "Risk"])))

fig_3a = ggplot(TPM_model_BIN1, aes(x = AD_risk, y = percent_TPM)) +  # Tissue risk on x-axis
  geom_boxplot(aes(fill = tissue)) +
  theme_bw(base_size = 16) +
  labs(x = "Tissue risk", y = "Relative TPM abundance", fill = "Tissue") + 
  scale_fill_manual(values = c("#00BFC4", rep("#F8766D", 2))) +
  stat_pvalue_manual(
    merged.data.final_gene_risk$BIN1,
    label = "p.signif",
    tip.length = 0.01,
    hide.ns = TRUE) +
  facet_wrap(~ transcript_name, scales = "free", nrow = 2) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -10, unit = "pt"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

ggsave("fig_3a.tiff", plot = fig_3a, width = 11, height = 6.5, dpi = "print")


#plot with all pairwise results (supp fig) 
ggplot(TPM_model_BIN1, aes(x = tissue, y = percent_TPM)) +  
  geom_boxplot(aes(fill = tissue)) +
  theme_bw(base_size = 16) +
  labs(x = "Tissue", y = "Relative TPM abundance") + 
  scale_fill_manual(values = c(("#00BFC4"), rep("#F8766D",2)))+
  stat_pvalue_manual(
    merged.data.final_gene$BIN1,
    label = "p.signif",
    tip.length = 0.01,
    hide.ns = TRUE) +
  facet_wrap(~ transcript_name, scales = "free", nrow = 2) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

################################################################################
#GAL3ST4
################################################################################

TPM_model_GAL3ST4 = TPM_tissue_gene_filtered %>%
  filter(Gene_name == "GAL3ST4")


TPM_model_GAL3ST4$tissue = factor(TPM_model_GAL3ST4$tissue, 
                                 levels = c(unique(TPM_model_GAL3ST4$tissue[TPM_model_GAL3ST4$AD_risk == "Protection"]), 
                                            unique(TPM_model_GAL3ST4$tissue[TPM_model_GAL3ST4$AD_risk == "Risk"])))

fig_3b = ggplot(TPM_model_GAL3ST4, aes(x = AD_risk, y = percent_TPM)) +  # Tissue risk on x-axis
  geom_boxplot(aes(fill = tissue)) +
  theme_bw(base_size = 16) +
  labs(x = "Tissue risk", y = "Relative TPM abundance", fill = "Tissue")  +
  scale_fill_manual(values = c(rep("#00BFC4",1), rep("#F8766D",2)))+
  stat_pvalue_manual(
    merged.data.final_gene_risk$GAL3ST4,
    label = "p.signif",
    tip.length = 0.01,
    hide.ns = TRUE) +
  facet_wrap(~ transcript_name, scales = "free", nrow = 1) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -10, unit = "pt"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))


ggsave("fig_3b.tiff", plot = fig_3b, width = 11, height = 3.85, dpi = "print")

#plot with all pairwise results (supp fig) 
ggplot(TPM_model_GAL3ST4, aes(x = tissue, y = percent_TPM)) +  
  geom_boxplot(aes(fill = tissue)) +
  theme_bw(base_size = 16) +
  labs(x = "Tissue", y = "Relative TPM abundance") + 
  scale_fill_manual(values = c(("#00BFC4"), rep("#F8766D",2)))+
  stat_pvalue_manual(
    merged.data.final_gene$GAL3ST4,
    label = "p.signif",
    tip.length = 0.01,
    hide.ns = TRUE) +
  facet_wrap(~ transcript_name, scales = "free", nrow = 2) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))


################################################################################
#ICA1L
################################################################################

TPM_model_ICA1L = TPM_tissue_gene_filtered %>%
  filter(Gene_name == "ICA1L")


TPM_model_ICA1L$tissue = factor(TPM_model_ICA1L$tissue, 
                                levels = c(unique(TPM_model_ICA1L$tissue[TPM_model_ICA1L$AD_risk == "Protection"]), 
                                           unique(TPM_model_ICA1L$tissue[TPM_model_ICA1L$AD_risk == "Risk"])))


fig_3c = ggplot(TPM_model_ICA1L, aes(x = AD_risk, y = percent_TPM)) +  # Tissue risk on x-axis
  geom_boxplot(aes(fill = tissue)) +
  theme_bw(base_size = 16) +
  labs(x = "Tissue risk", y = "Relative TPM abundance", fill = "Tissue")  +
  scale_fill_manual(values = c(rep("#00BFC4",2), rep("#F8766D",1)))+
  stat_pvalue_manual(
    merged.data.final_gene_risk$ICA1L,
    label = "p.signif",
    tip.length = 0.01,
    hide.ns = TRUE) +
  facet_wrap(~ transcript_name, scales = "free", nrow = 2) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -10, unit = "pt"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))


ggsave("fig_3c.tiff", plot = fig_3c, width = 11, height = 7, dpi = "print")

#plot with all pairwise results (supp fig) 
ggplot(TPM_model_ICA1L, aes(x = tissue, y = percent_TPM)) +  
  geom_boxplot(aes(fill = tissue)) +
  theme_bw(base_size = 16) +
  labs(x = "Tissue", y = "Relative TPM abundance") + 
  scale_fill_manual(values = c(("#00BFC4"), rep("#F8766D",2)))+
  stat_pvalue_manual(
    merged.data.final_gene$ICA1L,
    label = "p.signif",
    tip.length = 0.01,
    hide.ns = TRUE) +
  facet_wrap(~ transcript_name, scales = "free", nrow = 2) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))





################################################################################
#KAT8
################################################################################

TPM_model_KAT8 = TPM_tissue_gene_filtered %>%
  filter(Gene_name == "KAT8")

#only keep whole blood vs others for plot as whole blood is the anomoly

merged_data_KAT8_blood = merged.data.final_gene$KAT8 %>%
  filter(group1 == "Whole Blood" | group2 == "Whole Blood") %>%
  group_by(transcript_name) %>%
  mutate(y.position = max_TPM + 0.05 + (row_number() - 1) * 0.08) %>%
  ungroup()


#plot of whole blood significance vs other tissues
fig_4b = ggplot(TPM_model_KAT8, aes(x = tissue, y = percent_TPM)) +
  geom_boxplot(aes(fill = "black")) +
  theme_bw(base_size = 16) +
  labs(x = "Tissue", y = "Percent TPM") +
  scale_fill_manual(values = "darkgrey") +
  stat_pvalue_manual(merged_data_KAT8_blood,
                     label = "p.signif",
                     tip.length = 0.01,
                     hide.ns = TRUE) +
  facet_wrap(~as.factor(transcript_name), scales = "free", nrow = 2) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 65, hjust = 1, size = 14))+
  annotate("rect" , xmin = 7.6, xmax = 8.4, ymin= -Inf, ymax = Inf, fill = "#00BFC4",alpha=0.15)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

ggsave("fig_4b.tiff", plot = fig_4b, width = 11, height = 8, dpi = "print")


#all pairwise differences (supp figure)
ggplot(TPM_model_KAT8, aes(x = tissue, y = percent_TPM)) +
  geom_boxplot(aes(fill = "black")) +
  theme_bw(base_size = 16) +
  labs(x = "Tissue", y = "Percent TPM") +
  scale_fill_manual(values = "darkgrey") +
  stat_pvalue_manual(merged.data.final_gene$KAT8,
                     label = "p.signif",
                     tip.length = 0.01,
                     hide.ns = TRUE) +
  facet_wrap(~as.factor(transcript_name), scales = "free", nrow = 2) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 65, hjust = 1, size = 14))+
  annotate("rect" , xmin = 7.6, xmax = 8.4, ymin= -Inf, ymax = Inf, fill = "#00BFC4",alpha=0.15)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))


