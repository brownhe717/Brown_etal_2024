library(ggpubr)
library(viridis)

#INPUT are results from WT v. toyKD DESeq -- can be found in "~/Brown_etal_2024/Figure_3/INPUT/WT.v.toyKD.diff.txt"

# Generate MA plot
ggmaplot(res_wt_vs_toyKD, fdr = 0.05, fc = 2, size = 0.4, genenames = as.vector(res_wt_vs_toyKD$geneid), top = 10, 
         font.label = c("bold", 11), label.rectangle = TRUE, font.legend = "bold", font.main = "bold",
         palette = c("steelblue3", "yellowgreen","slategrey"),
         ggtheme = ggplot2::theme_classic(), select.top.method = "padj", ylim = c(-10,10))
#bgcolor("black") changes the background color to black
#Where up = WD enriched (EAD depleted); and down = WD depleted (EAD enirched)