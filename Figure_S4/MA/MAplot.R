library(ggpubr)
library(viridis)

#INPUT are results from PcKD v. toyPcKD, SfmbtKD v. toySfmbtKD, and ScmKD v. toyScmKD DESeq 
  # found in "~/Brown_etal_2024/Figure_S4/INPUT/"

# Generate MA plot (PcKD v. toyPcKD)
ggmaplot(res_PcKD_v_toyPcKD, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_PcKD_v_toyPcKD$geneid), 
         top = 20, #label.select = c("Antp", "vg", "nub", "ap", "Sox15", "Dr", "Ubx", "ey", "toy", "eya", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", 
         palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))

# Generate MA plot (SfmbtKD v. toySfmbtKD)
ggmaplot(res_SfmbtKD_v_toySfmbtKD, fdr = 0.05, fc = 2, size = 0.4, 
         genenames = as.vector(res_SfmbtKD_v_toySfmbtKD$geneid), top = 2, 
         label.select = c("Antp", "vg", "nub", "ap", "Sox15", "Dr", "Ubx", "ey", "toy", "eya", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", 
         palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))

# Generate MA plot (ScmKD v. toyScmKD)
ggmaplot(res_ScmKD_v_toyScmKD, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_ScmKD_v_toyScmKD$geneid), 
         top = 20, #label.select = c("Antp", "vg", "nub", "ap", "Sox15", "Dr", "Ubx", "ey", "toy", "eya", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", 
         palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))