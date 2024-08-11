library(ggpubr)
library(viridis)

#INPUT are results from (1) PcKD v. SfmbtKD and PcKD v. ScmKD, or (2) PcKD v. toySfmbtKD and PcKD v. toyScmKD DESeq 
# found in "~/Brown_etal_2024/Figure_6/INPUT/"

# Generate MA plot (PcKD v. SfmbtKD)
ggmaplot(res_PcKD_v_SfmbtKD, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_PcKD_v_SfmbtKD$geneid), 
         top = 10, label.select = c("Antp", "Abd-B", "vg","Ubx", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", 
         font.main = "bold", palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))

# Generate MA plot (PcKD v. ScmKD)
ggmaplot(res_PcKD_v_ScmKD, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_PcKD_v_ScmKD$geneid), 
         top = 10, label.select = c("Antp", "Abd-B", "vg","Ubx", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", 
         font.main = "bold", palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))

# Generate MA plot (PcKD v. toySfmbtKD)
ggmaplot(res_PcKD_v_toySfmbtKD, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_PcKD_v_toySfmbtKD$geneid), 
         top = 10, label.select = c("Antp", "Abd-B", "vg","Ubx", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", 
         font.main = "bold", palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))

# Generate MA plot (PcKD v. toyScmKD)
ggmaplot(res_PcKD_v_toyScmKD, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_PcKD_v_toyScmKD$geneid), 
         top = 10, label.select = c("Antp", "Abd-B", "vg","Ubx", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", 
         palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))

