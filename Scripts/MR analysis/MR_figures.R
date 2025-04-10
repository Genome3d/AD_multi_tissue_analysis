################################################################################
#MR figure Script 
################################################################################

#load libraries 
library(tidyverse)
library(wesanderson) #colour scheme 
library(pals) #colour scheme
library(RIdeogram) #ideogram
library(BiocManager) # install packages if needed
library(biomaRt) #ensembl gene annotations
library(forestploter) #forest plot

#Read in AD causal risk gene data - created in Multi_tissue_MR_data_filtering.R
LOAD_2SMR_data = read.csv("LOAD_2SMR_data.csv")[,-1]
LOAD_SNP_data = read.csv( "LOAD_SNP_data.csv")[,-1]

################################################################################
#Fig. 2A - Tile Plot

#LOLOAD risk genes present in at least 6 tissues
multi_tissue_LOAD_genes = LOAD_2SMR_data %>%
  group_by(exposure) %>%
  filter(n_distinct(tissue) > 5)

fig_2a = ggplot(multi_tissue_LOAD_genes, aes(y = tissue, x= exposure, fill = or)) +
  geom_tile(color = "black", lwd = 0.75, lty = 1,aes(width = 0.42, height = 1)) +  
  theme_bw(base_size = 18) +                                
  xlab("LOAD risk genes") +   
  ylab("Tissue") +
  labs(fill = "OR (risk)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_blank()) + 
  scale_fill_gradientn(colours = ocean.balance(30),  
                       values = scales::rescale(c(0.5, 1, 3.5)),
                       limits = c(0.5, 3.5)) 

ggsave("fig_2a.tiff", plot = fig_2a, width = 11, height = 4.4, dpi = "print")


################################################################################
#Fig. 2B - switch gene tile plot

LOAD_switch_genes = LOAD_2SMR_data %>%
  group_by(exposure) %>%
  filter(any(or > 1) & any(or < 1)) %>%
  ungroup()


fig_2b = ggplot(LOAD_switch_genes, aes(y = tissue, x= exposure, fill = or)) +
  geom_tile(color = "black", lwd = 0.75, lty = 1,aes(width = 0.38, height = 1)) + 
  theme_bw(base_size = 18) +                                  
  xlab("LOAD risk genes") +   
  ylab("Tissue") +
  labs(fill = "OR (risk)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, ),
        strip.text = element_blank()) + 
  scale_fill_gradientn(colours = ocean.balance(30),  
                       values = scales::rescale(c(0.5, 1, 3.5)),  
                       limits = c(0.5, 3.5)) 


ggsave("fig_2b.tiff", plot = fig_2b, width = 11, height = 4.4, dpi = "print")

################################################################################
#Fig. 3, 4b
# see Isoform_analysis.R script

################################################################################
#Fig. 4a - Forest Plot

KAT8_LOAD_data = full_join(LOAD_SNP_data[,c(2:7)], LOAD_2SMR_data[,c(4,9,10,13:15)]) %>%
  filter(exposure == "KAT8") %>% #select KAT8 gene
  mutate(aFC.lci95 = beta.exposure - 1.96 * se.exposure, #create upper and lower CI for allelic fold change
         aFC.uci95 = beta.exposure + 1.96 * se.exposure,
         'AD OR' = paste(rep(" ", 20), collapse = " "),
         'OR (95% CI)' = ifelse(is.na(or), "", sprintf("%.2f (%.2f to %.2f)",or,or_lci95,or_uci95)),
         'aFC (eQTL effect)' = paste(rep(" ", 20), collapse = " "),
         'log2(aFC) (95% CI)' = ifelse(is.na(beta.exposure), "", sprintf("%.2f (%.2f to %.2f)",beta.exposure,aFC.lci95,aFC.uci95))) %>%
  rename("Gene" = "exposure",
         "Tissue" = "tissue",
         "OR pval" = "pval.exposure",
         "pval.aFC" = "pval") %>%
  arrange(SNP)

#rearrange dataframe for forest plot
KAT8_forest_plot = KAT8_LOAD_data %>%
  group_by(Gene) %>%
  summarise_all(~ NA) %>%
  mutate(`Gene : Tissue` = Gene) %>%
  bind_rows(KAT8_LOAD_data %>% mutate(`Gene : Tissue` = paste0("    ", Tissue))) %>%
  mutate(`Gene : Tissue` = ifelse(is.na(Tissue), Gene, `Gene : Tissue`)) %>%
  dplyr::select(`Gene : Tissue`, everything(), -Gene, -Tissue) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), ifelse(is.numeric(.), NA, ""), .))) %>%
  ungroup() %>%
  mutate_at(c('OR pval','pval.aFC'), ~signif(., 4)) %>%
  mutate_at(c('OR pval','pval.aFC'), ~as.character(.)) %>%
  mutate_at(c('OR pval','pval.aFC'), ~replace_na(.,""))

#set forest plot themes
tm_KAT8 <- forest_theme(
  core=list(bg_params=list(fill = ifelse(KAT8_forest_plot$or < 1 , "#00BFC4" ,"#F8766D" ),
                           alpha = 0.15)),
  colhead=list(fg_params=list(hjust=0.5, x=0.5)),
  ci_Theight = 0.3)



fig_4a = forest(KAT8_forest_plot[c(1,2,12,13,5,14,15)],
       est = list(KAT8_forest_plot$or,
                  KAT8_forest_plot$beta.exposure),
       lower = list(KAT8_forest_plot$or_lci95,
                    KAT8_forest_plot$aFC.lci95),
       upper = list(KAT8_forest_plot$or_uci95,
                    KAT8_forest_plot$aFC.uci95),
       theme = tm_KAT8,
       ci_column = c(3,6),
       ref_line = c(1,0),
       ticks_at = list(c(0.5, 1, 1.5, 2),
                       c(-0.5, 0, 0.5, 1)),
)

ggsave("fig_4a.tiff", plot = fig_4a, width = 11, height = 3.5, dpi = "print")
################################################################################
#Supplementary Figure 1 - Ideogram using ensembl annotations

#Get unique gene list
total_gene_list_LOAD=unique(LOAD_2SMR_data$exposure)

#gene mart for annotating risk genes
gene_mart=useMart("ensembl",dataset = "hsapiens_gene_ensembl" )

#use ensembl to get location of risk genes
LOAD_MR_gene_location = getBM(attributes=c('hgnc_symbol',
                                           'chromosome_name',
                                           'start_position',
                                           'end_position',
                                           'band'),
                              filters = ("hgnc_symbol"),
                              values=list(total_gene_list_LOAD),
                              mart=gene_mart) 

LOAD_MR_gene_location_ideogram = LOAD_MR_gene_location %>%
  filter(nchar(chromosome_name) <= 2) %>%
  dplyr::select(c('chromosome_name', 'start_position', 'end_position'))%>%
  dplyr::rename("Chr" = "chromosome_name", 
                "Start" = "start_position", 
                "End" = "end_position") %>%
  mutate("Shape" = "circle", 
         "color" = "B32357", 
         "Type" = "LOAD risk gene") %>%
  dplyr::select(c('Type', 'Shape', 'Chr', 'Start', 'End', 'color'))

#annotations for ideogram
data("human_karyotype")
data("gene_density")

#create ideogram
ideogram(karyotype = human_karyotype, label = LOAD_MR_gene_location_ideogram,
         overlaid = gene_density,
         colorset1 = pal,
         label_type = "marker")

#convert to png to view
convertSVG("chromosome.svg", device = "png", dpi = 1600)

################################################################################
#Supplementary figure 2
#Forest plot of switch genes 

final_filtered_switch_forest = full_join(LOAD_SNP_data[,c(2:7)], LOAD_2SMR_data[,c(4,9,10,13:15)]) %>%
  filter(exposure %in% LOAD_switch_genes$exposure) %>% 
  mutate(aFC.lci95 = beta.exposure - 1.96 * se.exposure, 
         aFC.uci95 = beta.exposure + 1.96 * se.exposure,
         'AD OR' = paste(rep(" ", 20), collapse = " "),
         'OR (95% CI)' = ifelse(is.na(or), "", sprintf("%.2f (%.2f to %.2f)",or,or_lci95,or_uci95)),
         'aFC (eQTL effect)' = paste(rep(" ", 20), collapse = " "),
         'log2(aFC) (95% CI)' = ifelse(is.na(beta.exposure), "", sprintf("%.2f (%.2f to %.2f)",beta.exposure,aFC.lci95,aFC.uci95))) %>%
  rename("Gene" = "exposure",
         "Tissue" = "tissue",
         "OR pval" = "pval.exposure",
         "pval.aFC" = "pval") %>%
  arrange(Gene)

final_filtered_switch_paper_forest = final_filtered_switch_forest %>% 
  group_by(Gene) %>%
  summarise_all(~ NA) %>%
  mutate(`Gene : Tissue` = Gene) %>%
  bind_rows(final_filtered_switch_forest %>% mutate(`Gene : Tissue` = paste0("    ", Tissue))) %>%
  arrange(Gene, !is.na(Tissue)) %>%
  mutate(`Gene : Tissue` = ifelse(is.na(Tissue), Gene, `Gene : Tissue`)) %>%
  dplyr::select(`Gene : Tissue`, everything(), -Gene, -Tissue) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), ifelse(is.numeric(.), NA, ""), .))) %>%
  ungroup() %>%
  mutate_at(c('OR pval','pval.aFC'), ~signif(., 4)) %>%
  mutate_at(c('OR pval','pval.aFC'), ~as.character(.)) %>%
  mutate_at(c('OR pval','pval.aFC'), ~replace_na(.,""))


tm <- forest_theme(
  core=list(bg_params=list(fill = ifelse(final_filtered_switch_paper_forest$or < 1 , "#00BFC4" ,"#F8766D" ),
                           alpha = 0.15)),
  colhead=list(fg_params=list(hjust=0.5, x=0.5)),
  ci_Theight = 0.3,
  base_size = 15)


supp_fig_2 = forest(final_filtered_switch_paper_forest[c(1,2,12,13,5,14,15)],
       est = list(final_filtered_switch_paper_forest$or,
                  final_filtered_switch_paper_forest$beta.exposure),
       lower = list(final_filtered_switch_paper_forest$or_lci95,
                    final_filtered_switch_paper_forest$aFC.lci95),
       upper = list(final_filtered_switch_paper_forest$or_uci95,
                    final_filtered_switch_paper_forest$aFC.uci95),
       theme = tm,
       ci_column = c(3,6),
       ref_line = c(1,0),
       ticks_at = list(c(0.5, 1, 1.5, 2),
                       c(-0.5, 0, 0.5, 1)),
)

ggsave("supp_fig_2.tiff", plot = supp_fig_2, width = 14, height = 15, dpi = "print")

################################################################################
#Supplementary Figures 3:10
#streamline_transcript_analysis_script.R

################################################################################
#Supplementary Table 2 
LOAD_2SMR_data_table = merge(LOAD_2SMR_data, LOAD_SNP_data, by = c("id.exposure"))
LOAD_2SMR_data_table = LOAD_2SMR_data_table[,c(4,16, 7:10, 11:15)] 

LOAD_2SMR_data_table = LOAD_2SMR_data_table %>%
  rename("Exposure (Gene)" = "exposure.x",
         "beta.MR" = "b",
         "SE.MR" = "se",
         "pval.MR" = "pval",
         "Tissue" = "tissue.x",
         "OR" = "or",
         "OR_lci95" = "or_lci95",
         "OR_uci95" = "or_uci95",
         "lci95" = "lo_ci",
         "uci95" = "up_ci") %>%
  arrange((`Exposure (Gene)`))

write.csv(LOAD_2SMR_data_table, "LOAD_2SMR_data_table.csv")
################################################################################
#Supplementary Table 3
LOAD_2SMR_KAT8_table = LOAD_2SMR_data_table %>%
  filter(`Exposure (Gene)` == "KAT8")

write.csv(LOAD_2SMR_KAT8_table[1:11], "LOAD_2SMR_KAT8_table.csv")


