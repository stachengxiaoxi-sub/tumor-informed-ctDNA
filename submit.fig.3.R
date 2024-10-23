#######Figure.3
rm(list = ls())
gc()

library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)
library(readr)
############## Read in data ########
cli.2<-read.csv("F:/燃石医学/Choice-1-ITH 研究/manuscript/manuscript for Nature cancer/source data 1/source data 2.csv")

# Read the CSV file
plot_mat <- read_csv("F:/燃石医学/Choice-1-ITH 研究/manuscript/manuscript for Nature cancer/source data 1/plot_mat.csv")

# Convert the tibble to a data frame
plot_mat <- as.data.frame(plot_mat)

# Set the first column as rownames
rownames(plot_mat) <- plot_mat[[1]]

# Remove the first column as it is now the rownames
plot_mat <- plot_mat[,-1]

# View the result
head(plot_mat)

pvals <- c()

for (i in 11:18) {
  x <- cli.2[cli.2$ctDNA_positive_C1D1 == "yes", i]
  y <- cli.2[cli.2$ctDNA_positive_C1D1 == "no", i]
  
  pval <- wilcox.test(x, y)$p.value
  pvals <- c(pvals, pval)
}

pvals_df <- data.frame(column = colnames(cli.2)[11:18], pval = pvals)

# Extract the names of the columns
names(cli.2)

cols <- c("ctDNA_positive_C1D1", colnames(cli.2)[11:18])

# Extract data from columns 9 to 16
data<- cli.2[, cols]


# Convert data to long format
data_long <- melt(data, id.vars = "ctDNA_positive_C1D1")

# Plot multiple plots
library(ggprism)

ggplot(data_long, aes(x = ctDNA_positive_C1D1, y = value, fill = ctDNA_positive_C1D1)) +
  geom_violin(alpha=0.5) +
  #geom_dotplot(binaxis='y', stackdir='center', position = "dodge",dotsize = 0.5) +
  stat_summary(aes(colour = ctDNA_positive_C1D1), fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", position=position_dodge(0.9), size=0.5,width = 0.5) +
  stat_summary(aes(colour = ctDNA_positive_C1D1), fun.y=mean, geom="crossbar", position=position_dodge(0.9), size=0.5, width = 0.5) +
  scale_fill_manual(values=c("#008EC8", "#D8433F")) +
  scale_colour_manual(values=c("#008EC8", "#D8433F")) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "",
       x = "ctDNA_positive",
       y = "Value") +
  theme_prism()
########################figure 3H#####################

library(magrittr)
library(dplyr)
library(ComplexHeatmap)
library(ggpubr)
library(Hmisc)
library(circlize)
library(RColorBrewer)
# Define the driver genes to compare
driver_genes <- c('TP53', 'PTEN', 'NOTCH1', 'NFE2L2', 'LRP1B', 'KMT2D', 'KEAP1', 'CSMD3')

# Create an empty vector to store the gene names and p-values
gene_names <- c()
p_values <- c()

# Loop through each gene and perform Fisher's exact test
for (gene in driver_genes) {
  
  # Extract WT and Mut information for the current gene and create a 2x2 contingency table
  gene_table <- table(cli.2[[gene]], cli.2$ctDNA_positive_C1D1)
  
  # Perform Fisher's exact test on the contingency table
  fisher_test <- fisher.test(gene_table)
  
  # Store the gene name and its p-value
  gene_names <- c(gene_names, gene)
  p_values <- c(p_values, round(fisher_test$p.value, 3))
}

# Create a data frame with the results
fisher_results<- data.frame(Gene = gene_names, P_Value = p_values)


############## Differences in other variables between the two groups ########################


## Preparing clinical data for plotting group ----

onco_cli <- cli.2 %>%
  dplyr::select(tumor_sample_uuid2,
                group = ctDNA_positive_C1D1, treatment, age, gender, stage)

{
  # Format the gene names to match the row order with plot.mat
  calc_str <- fisher_results[match(rownames(plot_mat), fisher_results$Gene),] %>%
    dplyr::mutate(rname = paste0(Gene, "\t", P_Value)) 
  }

## plot mutation ----

{
  vc_syn <- c("Missense","Truncating","Inframe","Other")
  # One can use any colors, here in this example color palette from RColorBrewer package is used
  vc_cols = RColorBrewer::brewer.pal(n = length(vc_syn), name = 'Set1')
  names(vc_cols) = c(vc_syn)
  print(vc_cols)
}
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # #8DD3C7
  Truncating = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*1, 
              gp = gpar(fill = vc_cols["Truncating"], col = NA))
  },
  
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*1, 
              gp = gpar(fill = vc_cols["Missense"], col = NA))
  },
  
  Inframe = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*1, 
              gp = gpar(fill = vc_cols["Inframe"], col = NA))
  },
  Other = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, 
              gp = gpar(fill = vc_cols["Other"], col = NA))
  }
)
heatmap_legend_param = list(title = "Alternations", at = vc_syn, 
                            labels = vc_syn,direction = "vertical")

{
  # 左半
  left_clin <- cli.2 %>%
    dplyr::filter(ctDNA_positive_C1D1== 'no') %>%
    dplyr::select(tumor_sample_uuid2,ctDNA_positive_C1D1,treatment,age,gender,stage)
  left_data <- plot_mat %>%
    as.data.frame() %>%
    dplyr::select(left_clin$tumor_sample_uuid2)
  left_data[is.na(left_data)] <- ''
  # heatmap_legend_param = list(title = "Alternations", at = vc_syn, 
  #                             labels = vc_syn,direction = "horizontal",nrow = 1)
  
  ha_l_m = HeatmapAnnotation(ctDNA = factor(left_clin$ctDNA_positive_C1D1,levels = c('No','Yes')),
                             Treatment = factor(left_clin$treatment),
                             Age = anno_barplot(left_clin$age, ylim = c(0, max(left_clin$age, na.rm = TRUE)),border = F,gp=gpar(border =NA,fill="#008001",lty="blank")),
                             Gender = factor(left_clin$gender),
                             Stage = factor(left_clin$stage),
                             col = list(ctDNA = c("No" = ggsci::pal_jama('default')(2)[1], 
                                                  "Yes" = ggsci::pal_jama('default')(2)[2]),
                                        Treatment = c("ICI-chemo" = ggsci::pal_jama('default')(2)[1], 
                                                      "chemo" = ggsci::pal_jama('default')(2)[2]),
                                        Gender = c("Male" = ggsci::pal_jama('default')(2)[1], 
                                                   "Female" = ggsci::pal_jama('default')(2)[2]),
                                        Stage = c("I" = ggsci::pal_jama('default')(4)[1], 
                                                  "II" = ggsci::pal_jama('default')(4)[2],
                                                  "III" = ggsci::pal_jama('default')(4)[3],
                                                  "IV" = ggsci::pal_jama('default')(4)[4])),
                             show_legend = c(F), # 左侧图top 注释不显示图例
                             show_annotation_name = T,annotation_name_side = 'left',
                             gap = unit(2, "mm"),
                             annotation_legend_param = list(
                               ctDNA = list(title = "ctDNA",direction = "vertical"))
  )
  p_left_m <- oncoPrint(left_data,
                        # get_type = function(x) gsub(":.*$", "", strsplit(x, ";")[[1]]), # Take one out of two
                        top_annotation = ha_l_m,
                        # top_annotation = NULL,
                        right_annotation = NULL, # Do not display barplot
                        # Add rank and test P-values on the right side
                        alter_fun = alter_fun, 
                        pct_digits = 2,
                        # column_split = table(onco.cli$rsgroup)[1], # Split position, need to sort onco.cli by rsgroup
                        col = vc_cols, 
                        column_title = NULL, 
                        row_order = c(1,8,5,6,7,2,3,4), # 1
                        # column_order = onco.cli$patient_id,
                        remove_empty_columns = T, # Cannot delete, otherwise it will not work subsequently (actually not useful)
                        remove_empty_rows = TRUE,
                        show_heatmap_legend = F, # Do not display legend for the left plot
                        show_column_names = F,
                        show_row_names = F)
  p_left_m
  # right
  right_clin <- cli.2%>%
    dplyr::filter(ctDNA_positive_C1D1 == 'yes') %>%
    dplyr::select(tumor_sample_uuid2,ctDNA_positive_C1D1,treatment,age,gender,stage)
  right_data <- plot_mat %>%
    as.data.frame() %>%
    dplyr::select(right_clin$tumor_sample_uuid2)
  right_data[is.na(right_data)] <- ''
  ra_r_m <- rowAnnotation(`Chiq test` = anno_text(calc_str$p.res.str,
                                                  location = 1,
                                                  rot = 1,
                                                  just = "right",
                                                  show_name = T,
                                                  gp = gpar(fontsize = 10)),
                          annotation_name_side = 'top',
                          annotation_name_rot = 0)
  
  heatmap_legend_param = list(title = "Alternations", at = vc_syn, 
                              labels = vc_syn,direction = "vertical") # vertical horizontal
  
  ha_r_m = HeatmapAnnotation(ctDNA = factor(right_clin$ctDNA_positive_C1D1,levels = c('No','Yes')),
                             Treatment = factor(right_clin$treatment),
                             Age = anno_barplot(right_clin$age, ylim = c(0, max(right_clin$age, na.rm = TRUE)),border = F,gp=gpar(border =NA,fill="#008001",lty="blank")),
                             Gender = factor(right_clin$gender),
                             Stage = factor(right_clin$stage),
                             col = list(ctDNA = c("No" = ggsci::pal_jama('default')(2)[1], 
                                                  "Yes" = ggsci::pal_jama('default')(2)[2]),
                                        Treatment = c("ICI-chemo" = ggsci::pal_jama('default')(2)[1], 
                                                      "chemo" = ggsci::pal_jama('default')(2)[2]),
                                        Gender = c("Male" = ggsci::pal_jama('default')(2)[1], 
                                                   "Female" = ggsci::pal_jama('default')(2)[2]),
                                        Stage = c( "III" = ggsci::pal_jama('default')(4)[3],
                                                  "IV" = ggsci::pal_jama('default')(4)[4])),
                             show_legend = c(T), # 左侧图top 注释不显示图例
                             show_annotation_name = F,annotation_name_side = 'left',
                             gap = unit(2, "mm"),
                             annotation_legend_param = list(
                               ctDNA = list(title = "ctDNA",direction = "vertical"))
  )
  p_right_m <- oncoPrint(right_data,
                         # get_type = function(x) gsub(":.*$", "", strsplit(x, ";")[[1]]), # 两个取一个
                         top_annotation = ha_r_m,
                         # top_annotation = NULL,
                         right_annotation = ra_r_m, # 不显示barplot
                         # 右侧加秩和检验P值
                         alter_fun = alter_fun, 
                         pct_digits = 2,
                         # column_split = table(onco.cli$rsgroup)[1], # 切割位置,需对onco.cli按rsgroup排序
                         col = vc_cols, 
                         column_title = NULL, 
                         row_order = c(1,8,5,6,7,2,3,4), # 1
                         # column_order = onco.cli$patient_id,
                         remove_empty_columns = T,
                         remove_empty_rows = TRUE,
                         show_heatmap_legend = T, 
                         # show_row_names = F,
                         show_column_names = F,
                         # pct_side = 'right',  # 百分比放在右边
                         # row_names_side = "left",
                         heatmap_legend_param = heatmap_legend_param)
  p_right_m
  
  p_mut <- p_left_m + p_right_m
 
  draw(p_mut, 
       merge_legend = TRUE,
       heatmap_legend_side = "right", 
       annotation_legend_side = "right")
 
}



#########################figure 2I####################################

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(gridExtra)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Define the function to perform Fisher's exact test
fisher_test_func <- function(data, group_var, mut_var) {
  # Create a contingency table
  contingency_table <- data %>%
    group_by(!!group_var, !!mut_var) %>%
    summarise(count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = !!mut_var, values_from = count, values_fill = list(count = 0)) %>%
    mutate(across(c(`WT`, `Mut`), ~ .x + (is.na(.) * 0)))  # Replace NA with 0
  
  # Perform Fisher's exact test if the contingency table has enough data
  test_result <- contingency_table %>%
    filter(`WT` > 0 & `Mut` > 0) %>%
    mutate(p_value = fisher.test(.[[paste(group_var, "WT", sep = "_")]], .[[paste(group_var, "Mut", sep = "_")]])$p.value)
  
  # Return the results with p-values
  return(test_result)
}

# List of genes to perform Fisher's exact test
genes <- c("PI3K_Akt", "ErbB", "Ras", "mTOR", "MAPK", "HR", "Base.excision.repair", "p53",
           "Wnt", "Nucleotide.excision.repair", "Notch", "Mismatch.repair")

# Store results in a list
results_list <- list()

# Perform Fisher's exact test for each gene and store results
for (gene in genes) {
  if (gene %in% names(cli.2)) {  # Check if the gene column exists in the data frame
    results_list[[gene]] <- fisher_test_func(cli.2, group_var = "ctDNA_positive_C1D1", mut_var = gene)
  } else {
    message(paste("Column", gene, "does not exist in the data frame."))
  }
}

# Combine all results into a single data frame
results_df <- bind_rows(results_list) %>%
  mutate(p_str = ifelse(p_value < 0.01, "<0.01", sprintf("%.3f", p_value)))

# Calculate the mutation proportion for each gene and ctDNA status
calc_str <- results_df %>%
  mutate(p_str = ifelse(mut_p_value < 0.01,"<0.01",sprintf("%.3f",mut_p_value))) %>%
  select(gene,ctDNA_positive_C1D1,p_str) %>%
  mutate(gene = case_when(gene =='Base.excision.repair'~"BER",
                          gene =='Mismatch.repair'~"MMR",
                          gene =='Nucleotide.excision.repair'~"NER",
                          TRUE ~ gene))

# Pivot long data for plotting
tmp_data <- cli.2 %>%
  pivot_longer(cols = -c('Tumor_Sample_UUID','ctDNA_positive_C1D1'),
               names_to = "gene", values_to = "type") %>%
  group_by(gene,type,ctDNA_positive_C1D1) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(gene,ctDNA_positive_C1D1) %>%
  mutate(freq = round(n/sum(n),4)*100) %>%
  ungroup() %>%
  filter(type == 'Mut') %>%
  mutate(gene = case_when(gene =='Base.excision.repair'~"BER",
                          gene =='Mismatch.repair'~"MMR",
                          gene =='Nucleotide.excision.repair'~"NER",
                          TRUE ~ gene))

# Plot barplot
p_out <- ggplot(data=tmp_data, aes(x=gene, y=freq, fill=ctDNA_positive_C1D1)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=freq, group=ctDNA_positive_C1D1),
            fontface="bold", vjust = 1.2,
            position=position_dodge(width=0.9), size=2.5) + 
  scale_fill_manual(name = "",
                    values = c("#377EB8","#E41A1C"),
                    labels = c("Yes" = "ctDNA+", "No" = "ctDNA-")) +
  theme_classic() + 
  labs(title = "Mutation Proportion by ctDNA Positive C1D1 Status", x = "Gene", y = "Proportion of Mutations") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(calc_str, label = "{p_str}",
                     x = "gene",
                     y.position = max(tmp_data$freq)+1) +
  facet_grid( ~ gene, space  = "free", scales = "free", drop = TRUE, shrink = TRUE)

# Save the plot
ggsave("./results/ctDNA_WES_pathway_barplot.pdf", p_out, width = 8, height = 4)

# Visualization of the Results of Genetic Mutation Pathway Classification----

library(magrittr)
library(dplyr)
library(ggstatsplot)
library(ggpubr)
library(ggplot2)
library(ggstats)
library(Hmisc)
library(circlize)
library(RColorBrewer)
library(tidyverse)
library(ggprism)
names(cli.2)
# plot
plot_data <- cli.2 %>%
 dplyr::select("tumor_sample_uuid2", "ctDNA_positive_C1D1","PI3K_Akt","ErbB","Ras","mTOR","MAPK", "HR","Base.excision.repair","p53",                       
               "Wnt","Nucleotide.excision.repair","Notch","Mismatch.repair" ) 

# barplot
statistic_res <- data.frame()
# plot_list <- list()
for (p_name in colnames(plot_data)[-c(1,2)]) {
  # p_name <- 'ErbB'
  statistic_data <- plot_data %>%
    dplyr::select(ctDNA_positive_C1D1,pathway = p_name)
  p_value <- chisq.test(with(statistic_data,table(ctDNA_positive_C1D1, pathway)),correct = F)%>%
    .$p.value %>%
    signif(3)
  pval_pos <- data.frame(
    pathway = p_name,
    group1 = "ctDNA-",
    group2 = "ctDNA+",
    label = p_value,
    y.position = 1.05)
  statistic_res <- dplyr::bind_rows(statistic_res,pval_pos)
  
  # p_out <- ggplot(statistic_data) +
  #   geom_bar(aes(x = pathway, fill = ctDNA),position = "fill")+
  #   labs(y = "Percentage") +
  #   scale_y_continuous(labels = scales::percent) +
  #   scale_x_discrete() +
  #   scale_fill_discrete(name = "", labels = c("Yes" = "cdDNA+", "No" = "ctDNA-")) +
  #   theme_classic() +
  #   add_pvalue(pval_pos, label = "Chiq test's P = {label}")+
  #   theme(aspect.ratio = 2/1.5, 
  #         plot.title = element_text(hjust = 0.5, size = 12),
  #         axis.text.y = element_text(size = 12, vjust = -0.2), 
  #         axis.text.x = element_text(size = 10),
  #         axis.title = element_text(size = 12), 
  #         axis.title.x.bottom = element_text(size = 12, vjust = 0.5)) +
  #   ggtitle(p_name)
  # plot_list[[p_name]] <- p_out
}
# p_out <- ggarrange(plotlist = plot_list,ncol = 4,nrow = 3,common.legend = T)
# ggsave("./results/ctDNA_WES_pathway1.pdf",p_out,width = 10,height = 14)

# barplot 
{
  calc_str <- statistic_res %>%
    dplyr::mutate(p_str = ifelse(label < 0.01,"<0.01",sprintf("%.3f",label))) %>%
    dplyr::select(pathway,group1,group2,p_str)%>%
    dplyr::mutate(pathway = case_when(pathway =='Base excision repair'~"BER",
                                      pathway =='Mismatch repair'~"MMR",
                                      pathway =='Nucleotide excision repair'~"NER",
                                      TRUE ~ pathway))
 
  tmp_data <- plot_data %>%
    tidyr::pivot_longer(cols = -c('tumor_sample_uuid2','ctDNA_positive_C1D1'),
                        names_to = "pathway",values_to = "type") %>%
    dplyr::group_by(pathway,type,ctDNA_positive_C1D1) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(pathway,ctDNA_positive_C1D1) %>%
    dplyr::mutate(freq = round(n/sum(n),4)*100)%>%
    dplyr::ungroup() %>%
    dplyr::filter(type == 'Mut') %>%
    dplyr::mutate(pathway = case_when(pathway =='Base excision repair'~"BER",
                                      pathway =='Mismatch repair'~"MMR",
                                      pathway =='Nucleotide excision repair'~"NER",
                                      TRUE ~ pathway))
  
  p_out <- ggplot(data=tmp_data, aes(x=pathway, y=freq)) +
    geom_bar(aes(fill=ctDNA_positive_C1D1),stat="identity", position=position_dodge(width=0.9))+
    geom_text(aes(label=freq,group=ctDNA_positive_C1D1),
              fontface="bold",vjust = 1.2,
              position=position_dodge(width=0.9), size=2.5) + 
    scale_fill_manual(name = "",
                      values = c("#377EB8","#E41A1C"),
                      labels = c("Yes" = "ctDNA+", "No" = "ctDNA-"))+
    theme_classic() + 
    xlab("")+
    ylab("Alterations detect rate (%)")+
    # add_pvalue(calc_str, x = "pathway")
    stat_pvalue_manual(calc_str, label = "{p_str}",
                       x = "pathway",
                       y.position = max(tmp_data$freq)+1)+
    theme(legend.position = "top",
          axis.text.x = element_blank())+
    facet_grid( ~ pathway,  space  = "free", scales = "free", drop = TRUE, shrink = TRUE)
  p_out
  }

