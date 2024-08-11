library(DESeq2)
library(ComplexHeatmap)
library(circlize)

####################
##### Heatmaps #####
####################

#INPUT are results from (1) PcKD v. SfmbtKD and PcKD v. ScmKD, or (2) PcKD v. toySfmbtKD and PcKD v. toyScmKD DESeq 
# found in "~/Brown_etal_2024/Figure_6/INPUT/"

# Get Gene IDs for significantly changed genes in the PcKD vs. EAD and WD comparisons, remove duplicates, 
# and use them to extract their corresponding log2FCs

sfmbt_for_heatmap <- res_PcKD_v_SfmbtKD_sig %>%
  dplyr::select("geneid")

scm_for_heatmap <- res_PcKD_v_ScmKD_sig %>%
  dplyr::select("geneid")



toy.sfmbt_for_heatmap <- res_PcKD_v_toySfmbtKD_sig %>%
  dplyr::select("geneid")

toy.scm_for_heatmap <- res_PcKD_v_toyScmKD_sig %>%
  dplyr::select("geneid")




heatmap_genes_DE <- bind_rows(sfmbt_for_heatmap, scm_for_heatmap) %>% 
  distinct

heatmap_genes_toy = bind_rows(toy.sfmbt_for_heatmap, toy.scm_for_heatmap) %>%
  distinct





sfmbt_log2fc <- filter(res_PcKD_v_SfmbtKD, geneid %in% heatmap_genes_DE$geneid) %>%
  dplyr::select("geneid", "log2FoldChange")

scm_log2fc <- filter(res_PcKD_v_ScmKD, geneid %in% heatmap_genes_DE$geneid) %>%
  dplyr::select("geneid", "log2FoldChange")


toy.sfmbt_log2fc <- filter(res_PcKD_v_toySfmbtKD, geneid %in% heatmap_genes_toy$geneid) %>%
  dplyr::select("geneid", "log2FoldChange")

toy.scm_log2fc <- filter(res_PcKD_v_toyScmKD, geneid %in% heatmap_genes_toy$geneid) %>%
  dplyr::select("geneid", "log2FoldChange")



sig_matrix_DE <- left_join(sfmbt_log2fc, scm_log2fc, by = "geneid", suffix = c("Sfmbt", "Scm")) %>%
  column_to_rownames("geneid") %>%
  data.matrix

sig_matrix_toy <- left_join(toy.sfmbt_log2fc, toy.scm_log2fc, by = "geneid", suffix = c("toy-Sfmbt", "toy-Scm")) %>%
  column_to_rownames("geneid") %>%
  data.matrix

# Replace all rownames except those corresponding to genes of interest
# Note: if you want to change the list of genes to be labeled, first rerun the above pipe
rownames(sig_matrix_DE)[!rownames(sig_matrix_DE) %in% c("eya", "toy", "dac", "vg", "Ubx", "sd", "Scr", "Abd-B", "Antp", "ey")] <- ""

Heatmap(sig_matrix_DE, col = colorRamp2(c(-10, -5, 0, 5, 10), c("midnightblue", "cadetblue3", "snow1", "chocolate2", "firebrick")),
        cluster_rows = T, cluster_columns = T, row_gap = unit(2, "mm"), 
        heatmap_legend_param = list(title = "log2(FC)"),
        row_dend_reorder = F, show_row_dend = T, column_gap = unit(1, "mm"),
        row_names_gp = gpar(fontsize = 8, fontface = 3), border = T, row_km = 4,
        column_names_rot = 360, column_names_centered = T, show_column_names = T,
        column_labels = c("SfmbtKD EAD", "ScmKD EAD"), column_names_side = "top", 
        column_title_side = "top", column_title = "PcKD eye-antennal disc expression vs.")

rownames(sig_matrix_toy)[!rownames(sig_matrix_toy) %in% c("eya", "toy", "dac", "vg", "Ubx", "sd", "Scr", "Abd-B", "Antp", "ey")] <- ""

Heatmap(sig_matrix_toy, col = colorRamp2(c(-10, -5, 0, 5, 10), c("midnightblue", "cadetblue3", "snow1", "chocolate2", "firebrick")),
        cluster_rows = T, cluster_columns = T, row_gap = unit(2, "mm"), 
        heatmap_legend_param = list(title = "log2(FC)"),
        row_dend_reorder = T, show_row_dend = T, row_km = 3,
        row_names_gp = gpar(fontsize = 8, fontface = 3), border = T, column_names_rot = 360, 
        column_names_centered = T, show_column_names = T,
        column_labels = c("toy-SfmbtKD EAD", "toy-ScmKD EAD"), column_names_side = "top", 
        column_title_side = "top", column_title = "PcKD eye-antennal disc expression vs.")

#column_split = 2,
#row_km = # -> separate rows based on clusered km (mean)
#row_split = # -> split by categorical variables