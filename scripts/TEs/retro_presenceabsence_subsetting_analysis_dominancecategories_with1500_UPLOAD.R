library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(broom)

##R script for analysis of presence/absence of retrotransposons based on homoeolog expression
##bias categories. Only 1.5kb and 5kb subsets were used for further analysis.

##now do analysis for dominance subsets

retros_1kb <- read.csv("Y:/Sophie/transcriptome_paper/Analysis_of_Promoters/Retro/retro_1kb_removed0_bytriad.csv")
retros_1.5kb <- read.csv("Y:/Sophie/transcriptome_paper/Analysis_of_Promoters/Retro/retro_1.5kb_removed0_bytriad.csv")
retros_2kb <- read.csv("Y:/Sophie/transcriptome_paper/Analysis_of_Promoters/Retro/retro_2kb_removed0_bytriad.csv")
retros_5kb <- read.csv("Y:/Sophie/transcriptome_paper/Analysis_of_Promoters/Retro/retro_5kb_removed0_bytriad.csv")


##upload the gene name files for the top and bottom sets of genes (10-80-10, 25-50-25, 5-90-5, 1-98-1)
setwd("Y:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/")
central_genes <- read.csv("HC_CS_no_stress_movement_Central.txt")
A_dom_genes <- read.csv("HC_CS_no_stress_movement_A.dominant.txt")
A_sup_genes <- read.csv("HC_CS_no_stress_movement_A.suppressed.txt")
B_dom_genes <- read.csv("HC_CS_no_stress_movement_B.dominant.txt")
B_sup_genes <- read.csv("HC_CS_no_stress_movement_B.suppressed.txt")
D_dom_genes <- read.csv("HC_CS_no_stress_movement_D.dominant.txt")
D_sup_genes <- read.csv("HC_CS_no_stress_movement_D.suppressed.txt")

##now want to subset into the relevant subsets (i.e. split A_dom into A and B/D samples)

A_dom_genes_A <- A_dom_genes %>% filter(grepl("A01", x))
A_dom_genes_BD <- A_dom_genes %>% filter(!grepl("A01", x))
B_dom_genes_B <- B_dom_genes %>% filter(grepl("B01", x))
B_dom_genes_AD <- B_dom_genes %>% filter(!grepl("B01", x))
D_dom_genes_D <- D_dom_genes %>% filter(grepl("D01", x))
D_dom_genes_AB <- D_dom_genes %>% filter(!grepl("D01", x))

A_sup_genes_A <- A_sup_genes %>% filter(grepl("A01", x))
A_sup_genes_BD <- A_sup_genes %>% filter(!grepl("A01", x))
B_sup_genes_B <- B_sup_genes %>% filter(grepl("B01", x))
B_sup_genes_AD <- B_sup_genes %>% filter(!grepl("B01", x))
D_sup_genes_D <- D_sup_genes %>% filter(grepl("D01", x))
D_sup_genes_AB <- D_sup_genes %>% filter(!grepl("D01",x))

##now merge
dominant_genes <- rbind(A_dom_genes_A, B_dom_genes_B, D_dom_genes_D)
dominant_genes_supp <- rbind(A_dom_genes_BD, B_dom_genes_AD, D_dom_genes_AB)
suppressed_genes <- rbind(A_sup_genes_A, B_sup_genes_B, D_sup_genes_D)
suppressed_genes_dom <- rbind(A_sup_genes_BD, B_sup_genes_AD, D_sup_genes_AB)

##add category labels to the retro table
retros_1kb = within(retros_1kb, {
  category = ifelse(retros_1kb$gene_id %in% dominant_genes$x, "Dom",
                    ifelse(retros_1kb$gene_id %in% dominant_genes_supp$x, "Dom_Low",
                           ifelse(retros_1kb$gene_id %in% suppressed_genes$x, "Supp",
                                  ifelse(retros_1kb$gene_id %in% suppressed_genes_dom$x, "Supp_High",
                                         ifelse(retros_1kb$gene_id %in% central_genes$x, "Central", "NA"))))) 
})

retros_1.5kb = within(retros_1.5kb, {
  category = ifelse(retros_1.5kb$gene_id %in% dominant_genes$x, "Dom",
                    ifelse(retros_1.5kb$gene_id %in% dominant_genes_supp$x, "Dom_Low",
                           ifelse(retros_1.5kb$gene_id %in% suppressed_genes$x, "Supp",
                                  ifelse(retros_1.5kb$gene_id %in% suppressed_genes_dom$x, "Supp_High",
                                         ifelse(retros_1.5kb$gene_id %in% central_genes$x, "Central", "NA"))))) 
})

retros_2kb = within(retros_2kb, {
  category = ifelse(retros_2kb$gene_id %in% dominant_genes$x, "Dom",
                    ifelse(retros_2kb$gene_id %in% dominant_genes_supp$x, "Dom_Low",
                           ifelse(retros_2kb$gene_id %in% suppressed_genes$x, "Supp",
                                  ifelse(retros_2kb$gene_id %in% suppressed_genes_dom$x, "Supp_High",
                                         ifelse(retros_2kb$gene_id %in% central_genes$x, "Central", "NA"))))) 
})

retros_5kb = within(retros_5kb, {
  category = ifelse(retros_5kb$gene_id %in% dominant_genes$x, "Dom",
                    ifelse(retros_5kb$gene_id %in% dominant_genes_supp$x, "Dom_Low",
                           ifelse(retros_5kb$gene_id %in% suppressed_genes$x, "Supp",
                                  ifelse(retros_5kb$gene_id %in% suppressed_genes_dom$x, "Supp_High",
                                         ifelse(retros_5kb$gene_id %in% central_genes$x, "Central", "NA"))))) 
})

##gather the data frames and combine
retros_1kb_long <- gather(retros_1kb, variable, value, -gene_id, -group_id, -category)
retros_1.5kb_long <- gather(retros_1.5kb, variable, value, -gene_id, -group_id, -category)
retros_2kb_long <- gather(retros_2kb, variable, value, -gene_id, -group_id, -category)
retros_5kb_long <- gather(retros_5kb, variable, value, -gene_id, -group_id, -category)

full_data_table <- rbind(retros_1kb_long, retros_1.5kb_long, retros_2kb_long, retros_5kb_long)

#chi square test
full_data_table$category <- as.character(full_data_table$category)

pres_1kb <- table((full_data_table%>%filter(category!="NA")%>%filter(grepl("presence_1kb",variable))%>%filter(!grepl("sum",variable))%>%pull(.,value)),(full_data_table%>%filter(category!="NA")%>%filter(grepl("presence_1kb",variable))%>%filter(!grepl("sum",variable))%>%pull(.,category)))
chisq.test(pres_1kb)
pres_1.5kb <- table((full_data_table%>%filter(category!="NA")%>%filter(grepl("presence_1.5kb",variable))%>%filter(!grepl("sum",variable))%>%pull(.,value)),(full_data_table%>%filter(category!="NA")%>%filter(grepl("presence_1.5kb",variable))%>%filter(!grepl("sum",variable))%>%pull(.,category)))
chisq.test(pres_1.5kb)
pres_2kb <- table((full_data_table%>%filter(category!="NA")%>%filter(grepl("presence_2kb",variable))%>%filter(!grepl("sum",variable))%>%pull(.,value)),(full_data_table%>%filter(category!="NA")%>%filter(grepl("presence_2kb",variable))%>%filter(!grepl("sum",variable))%>%pull(.,category)))
chisq.test(pres_2kb)
pres_5kb <- table((full_data_table%>%filter(category!="NA")%>%filter(grepl("presence_5kb",variable))%>%filter(!grepl("sum",variable))%>%pull(.,value)),(full_data_table%>%filter(category!="NA")%>%filter(grepl("presence_5kb",variable))%>%filter(!grepl("sum",variable))%>%pull(.,category)))
chisq.test(pres_5kb)
##chi-square test is not relevant, as only two comparisons!
anova = full_data_table %>% group_by(variable) %>% do(tidy(aov(value ~ category, data=.)))
##no sig diff 
setwd("Y://Sophie/transcriptome_paper/Analysis_of_Promoters/Retro/")
write.csv(anova, "Retros_anova_between_dominance_categories_0removed_1500inc.csv")

##now get plots
ggplot(full_data_table, aes(x=variable, y=value, fill=category))+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90))

summary_retro <- ddply(full_data_table, c("variable","category"), summarise,
                       N = length(value),
                       mean = mean(value),
                       sd = sd(value),
                       se = sd / sqrt(N))

ggplot(summary_retro, aes(x=variable, y=mean, fill=category)) +
  geom_bar(position=position_dodge(),stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(angle = 90))

##just plot the count data now
ggplot((full_data_table%>%filter(grepl("count",variable))%>%filter(!grepl("sum",variable))%>%filter(category != "NA")), aes(x=variable, y=value, fill=category))+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90))

ggplot((summary_retro%>%filter(grepl("count",variable))%>%filter(!grepl("sum",variable))%>%filter(category != "NA")), aes(x=variable, y=mean, fill=category)) +
  geom_bar(position=position_dodge(),stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(angle = 90))

##just plot the presence absence data now
ggplot((full_data_table%>%filter(grepl("presence",variable))%>%filter(!grepl("sum",variable))%>%filter(category != "NA")), aes(x=variable, y=value, fill=category))+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90))

ggplot((summary_retro%>%filter(grepl("presence",variable))%>%filter(!grepl("sum",variable))%>%filter(category != "NA")), aes(x=variable, y=mean, fill=category)) +
  geom_bar(position=position_dodge(),stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(angle = 90))

count_1kb <- table((full_data_table%>%filter(category!="NA")%>%filter(grepl("count_1kb",variable))%>%pull(.,value)),(full_data_table%>%filter(category!="NA")%>%filter(grepl("count_1kb",variable))%>%pull(.,category)))
chisq.test(count_1kb)
count_1.5kb <- table((full_data_table%>%filter(category!="NA")%>%filter(grepl("count_1.5kb",variable))%>%pull(.,value)),(full_data_table%>%filter(category!="NA")%>%filter(grepl("count_1.5kb",variable))%>%pull(.,category)))
chisq.test(count_1.5kb)
count_2kb <- table((full_data_table%>%filter(category!="NA")%>%filter(grepl("count_2kb",variable))%>%pull(.,value)),(full_data_table%>%filter(category!="NA")%>%filter(grepl("count_2kb",variable))%>%pull(.,category)))
chisq.test(count_2kb)
count_5kb <- table((full_data_table%>%filter(category!="NA")%>%filter(grepl("count_5kb",variable))%>%pull(.,value)),(full_data_table%>%filter(category!="NA")%>%filter(grepl("count_5kb",variable))%>%pull(.,category)))
chisq.test(count_5kb)

full_data_table$category <- as.factor(full_data_table$category)
kruskal = full_data_table %>% filter(category!="NA") %>% filter(grepl("count",variable)) %>% filter(!grepl("sum",variable)) %>%
  group_by(variable) %>% do(tidy(kruskal.test(value ~ category, data=.)))

wilcox = full_data_table %>% filter(category != "NA") %>% filter(grepl("count",variable)) %>% filter(!grepl("sum",variable)) %>% 
  group_by(variable) %>% do(tidy(pairwise.wilcox.test(.$value, .$category,p.adjust.method = "BH")))
wilcox_nopadj = full_data_table %>% filter(category != "NA") %>% filter(grepl("count",variable)) %>% filter(!grepl("sum",variable)) %>% 
  group_by(variable) %>% do(tidy(pairwise.wilcox.test(.$value, .$category,p.adjust.method = "none")))

write.csv(kruskal, "Retros_kruskal_between_dominance_categories_0removed_1500inc.csv")
write.csv(wilcox, "Retros_wilcoxBH_between_dominance_categories_0removed_1500inc.csv")
write.csv(wilcox_nopadj, "Retros_wilcox_between_dominance_categories_0removed_1500inc.csv")
