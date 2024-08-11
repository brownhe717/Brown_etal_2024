library(clusterProfiler)
library(org.Dm.eg.db)
library(GOplot)

#INPUT are results from WT v. toyKD DESeq -- can be found in "~/Brown_etal_2024/Figure_3/INPUT/WT.v.toyKD.diff.txt"

geneList = res_wt_vs_toyKD_sig[,1]
head(geneList)

ego = enrichGO(geneList, OrgDb = "org.Dm.eg.db", ont = "BP", pAdjustMethod = "fdr", keyType = "SYMBOL")
head(ego)

dotplot(ego, showCategory=10)