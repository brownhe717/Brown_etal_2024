library(ggpubr)
library(viridis)

#INPUT are results from WT v. PcKD, WT v. SfmbtKD, and WT v. ScmKD DESeq 
  # found in "~/Brown_etal_2024/Figure_4/INPUT/"

# Generate MA plot (WT v. PcKD)
ggmaplot(res_WT_vs_PcKD, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_WT_vs_PcKD$Geneid), top = 2, 
         label.select = c("Antp", "vg", "nub", "ap", "Sox15", "Dr", "Ubx", "ey", "toy", "eya", "Scr"), 
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", 
         palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))

# Generate MA plot (WT v. SfmbtKD)
ggmaplot(res_wt_v_SfmbtKD, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_wt_v_SfmbtKD$geneid), top = 2, 
         label.select = c("Antp", "vg", "nub", "ap", "Sox15", "Dr", "Ubx", "ey", "toy", "eya", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", 
         palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))

# Generate MA plot (WT v. ScmKD)
ggmaplot(res_wt_v_ScmKD, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_wt_v_ScmKD$geneid), top = 2, 
         label.select = c("Antp", "vg", "nub", "ap", "Sox15", "Dr", "Ubx", "ey", "toy", "eya", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", 
         palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))
