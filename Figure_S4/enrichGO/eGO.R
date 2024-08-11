library(clusterProfiler)
library(org.Dm.eg.db)
library(GOplot)

#INPUT are results from PcKD v. toyPcKD, SfmbtKD v. toySfmbtKD, and ScmKD v. toyScmKD DESeq 
# found in "~/Brown_etal_2024/Figure_S4/INPUT/"

#PcKD v. toyPcKD
geneList = res_PcKD_v_toyPcKD_sig[,1]
head(geneList)

ego = enrichGO(geneList, OrgDb = "org.Dm.eg.db", ont = "BP", pAdjustMethod = "fdr", keyType = "SYMBOL")
head(ego)

dotplot(ego, showCategory=12)

#SfmbtKD v. toySfmbtKD
geneList = res_SfmbtKD_v_toySfmbtKD_sig[,1]
head(geneList)

ego = enrichGO(geneList, OrgDb = "org.Dm.eg.db", ont = "BP", pAdjustMethod = "fdr", keyType = "SYMBOL")
head(ego)

dotplot(ego, showCategory=15)

#ScmKD v. toyScmKD
geneList = res_ScmKD_v_toyScmKD_sig[,1]
head(geneList)

ego = enrichGO(geneList, OrgDb = "org.Dm.eg.db", ont = "BP", pAdjustMethod = "fdr", keyType = "SYMBOL")
head(ego)

dotplot(ego, showCategory=15)