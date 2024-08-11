library(ggpubr)
library(viridis)

#INPUT are results from toyKD v. toyPcKD, toyKD v. toySfmbtKD, and toyKD v. toyScmKD DESeq 
  # found in "~/Brown_etal_2024/Figure_5/INPUT/"

# Generate MA plot (toyKD v. toyPcKD)
ggmaplot(res_toy_v_toyPcKD, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_toy_v_toyPcKD$geneid), 
         top = 10, label.select = c("Antp", "Abd-B", "vg","Ubx", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", 
         font.main = "bold", palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))

# Generate MA plot (toyKD v. toySfmbtKD)
ggmaplot(res_toy_v_toySfmbtKD, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_toy_v_toySfmbtKD$geneid), 
         top = 10, label.select = c("Antp", "Abd-B", "vg","Ubx", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold", 
         palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10), max.overlaps = FALSE)

# Generate MA plot (toyKD v. toyScmKD)
ggmaplot(res_toy_v_toyScmKD, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_toy_v_toyScmKD$geneid), 
         top = 10, label.select = c("Antp", "Abd-B", "vg","Ubx", "Scr"),
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", 
         font.main = "bold", palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))
