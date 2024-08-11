library(DESeq2)
library(ggplot2)

#INPUT is dds object from DESeq2 -- script found in "Brown_etal_2024/RNAseq/DESeq2/DeSeq2.R"

#########################
###### Plot Counts ######
#########################

dds$tissue <- factor(dds$tissue, levels=c("WT","toyKD","PcKD","toy-PcKD",
                                          "SfmbtKD", "toy-SfmbtKD","ScmKD","toy-ScmKD"))

#Pc counts
d <- plotCounts(dds, gene="Pc", intgroup = "tissue", returnData = TRUE)

ggplot(d, aes(x=tissue,y=count, color=tissue)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  theme_bw() +
  scale_color_viridis_d() +
  theme(axis.text.x = element_text(angle = 90))

#Sfmbt counts
d <- plotCounts(dds, gene="Sfmbt", intgroup = "tissue", returnData = TRUE)

ggplot(d, aes(x=tissue,y=count, color=tissue)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  theme_bw() +
  scale_color_viridis_d() +
  theme(axis.text.x = element_text(angle = 90))

#Scm counts
d <- plotCounts(dds, gene="Scm", intgroup = "tissue", returnData = TRUE)

ggplot(d, aes(x=tissue,y=count, color=tissue)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  theme_bw() +
  scale_color_viridis_d() +
  theme(axis.text.x = element_text(angle = 90))

#upd1 counts
d <- plotCounts(dds, gene="upd1", intgroup = "tissue", returnData = TRUE)

ggplot(d, aes(x=tissue,y=count, color=tissue)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  theme_bw() +
  scale_color_viridis_d() +
  theme(axis.text.x = element_text(angle = 90))

#upd2 counts
d <- plotCounts(dds, gene="upd2", intgroup = "tissue", returnData = TRUE)

ggplot(d, aes(x=tissue,y=count, color=tissue)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  theme_bw() +
  scale_color_viridis_d() +
  theme(axis.text.x = element_text(angle = 90))

#upd3 counts
d <- plotCounts(dds, gene="upd3", intgroup = "tissue", returnData = TRUE)

ggplot(d, aes(x=tissue,y=count, color=tissue)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  theme_bw() +
  scale_color_viridis_d() +
  theme(axis.text.x = element_text(angle = 90))
