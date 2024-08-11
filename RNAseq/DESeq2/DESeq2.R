#!/usr/env/bin Rscript

library(DESeq2)
library(tidyverse)
library(org.Dm.eg.db)

# Read count matrix
counts <- read.delim("~/Brown_etal_2024/RNAseq/DESeq2/INPUT/counts.txt",sep="\t",row.names="geneid")

# Generate experimental design formula (rewritten to generate rownames in the same command)
coldata <- data.frame(row.names = colnames(counts), "tissue" = c(rep("WT",3),
                                                                 rep("toyKD",3),
                                                                 rep("PcKD",3),
                                                                 rep("toy-PcKD",3),
                                                                 rep("SfmbtKD",3),
                                                                 rep("toy-SfmbtKD",3),
                                                                 rep("ScmKD",3),
                                                                 rep("toy-ScmKD",3)))

# Input count matrix
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~tissue)

# Filter out rows with less than 10 counts total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Perform Wald test
dds <- DESeq(dds)

#########################################################
## Group 1: What changes occur when toy is lost alone? ##
################## WT EAD vs. toyKD EAD #################
#########################################################

# Perform DE testing ("treated", "untreated")
res_wt_vs_toyKD <- results(dds, contrast=c("tissue","toyKD","WT"))

# Add fold change classification to complete results
res_wt_vs_toyKD <- as.data.frame(res_wt_vs_toyKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_wt_vs_toyKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in toyKD",
    log2FoldChange <= -1 ~ "depleted in toyKD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_wt_vs_toyKD_sig <- subset(res_wt_vs_toyKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_wt_vs_toyKD_sig, "WT.v.toyKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

####################################################################
### Group 2: What changes occur when a PcG member is lost alone? ###
########################### WT vs. PcKD ############################
####################################################################

# Perform DE testing 
res_WT_vs_PcKD <- results(dds, contrast=c("tissue","PcKD", "WT"))

# Add fold change classification to complete results
res_WT_vs_PcKD <- as.data.frame(res_WT_vs_PcKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_WT_vs_PcKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in absence of Pc",
    log2FoldChange <= -1 ~ "depleted in absence of Pc",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_WT_vs_PcKD_sig <- subset(res_WT_vs_PcKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_WT_vs_PcKD_sig, "WT.v.PcKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

########################
#### WT vs. SfmbtKD ####
########################

# Perform DE testing 
res_wt_v_SfmbtKD <- results(dds, contrast=c("tissue","SfmbtKD", "WT"))

# Add fold change classification to complete results
res_wt_v_SfmbtKD <- as.data.frame(res_wt_v_SfmbtKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_wt_v_SfmbtKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in SfmbtKD",
    log2FoldChange <= -1 ~ "depleted in SfmbtKD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_wt_v_SfmbtKD_sig <- subset(res_wt_v_SfmbtKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_wt_v_SfmbtKD_sig, "WT.v.SfmbtKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

########################
#### WT vs. ScmKD ####
########################

# Perform DE testing 
res_wt_v_ScmKD <- results(dds, contrast=c("tissue","ScmKD", "WT"))

# Add fold change classification to complete results
res_wt_v_ScmKD <- as.data.frame(res_wt_v_ScmKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_wt_v_ScmKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in ScmKD",
    log2FoldChange <= -1 ~ "depleted in ScmKD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_wt_v_ScmKD_sig <- subset(res_wt_v_ScmKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_wt_v_ScmKD_sig, "WT.v.ScmKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

############################################################################
### Group 3: What changes occur when both toy and a PcG member are lost? ### 
############################ toyKD vs. toy-PcKD ############################
############################################################################

# Perform DE testing 
res_toy_v_toyPcKD <- results(dds, contrast=c("tissue","toy-PcKD", "toyKD"))

# Add fold change classification to complete results
res_toy_v_toyPcKD <- as.data.frame(res_toy_v_toyPcKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_toy_v_toyPcKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in toy-PcKD",
    log2FoldChange <= -1 ~ "depleted in toy-PcKD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_toy_v_toyPcKD_sig <- subset(res_toy_v_toyPcKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_toy_v_toyPcKD_sig, "toyKD.v.toy-PcKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

###############################
#### toyKD vs. toy-SfmbtKD ####
###############################

# Perform DE testing 
res_toy_v_toySfmbtKD <- results(dds, contrast=c("tissue","toy-SfmbtKD", "toyKD"))

# Add fold change classification to complete results
res_toy_v_toySfmbtKD <- as.data.frame(res_toy_v_toySfmbtKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_toy_v_toySfmbtKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in toy-SfmbtKD",
    log2FoldChange <= -1 ~ "depleted in toy-SfmbtKD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_toy_v_toySfmbtKD_sig <- subset(res_toy_v_toySfmbtKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_toy_v_toySfmbtKD_sig, "toyKD.v.toy-SfmbtKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

#############################
#### toyKD vs. toy-ScmKD ####
#############################

# Perform DE testing 
res_toy_v_toyScmKD <- results(dds, contrast=c("tissue","toy-ScmKD", "toyKD"))

# Add fold change classification to complete results
res_toy_v_toyScmKD <- as.data.frame(res_toy_v_toyScmKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_toy_v_toyScmKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in toy-ScmKD",
    log2FoldChange <= -1 ~ "depleted in toy-ScmKD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_toy_v_toyScmKD_sig <- subset(res_toy_v_toyScmKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_toy_v_toyScmKD_sig, "toyKD.v.toy-ScmKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

##############################################################################
### Group 4: What toy-dependent changes occur to allow for an eye-to-wing? ### 
########################## SfmbtKD vs. toy-SfmbtKD ###########################
##############################################################################

# Perform DE testing 
res_SfmbtKD_v_toySfmbtKD <- results(dds, contrast=c("tissue","toy-SfmbtKD", "SfmbtKD"))

# Add fold change classification to complete results
res_SfmbtKD_v_toySfmbtKD <- as.data.frame(res_SfmbtKD_v_toySfmbtKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_SfmbtKD_v_toySfmbtKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in toy-SfmbtKD",
    log2FoldChange <= -1 ~ "depleted in toy-SfmbtKD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_SfmbtKD_v_toySfmbtKD_sig <- subset(res_SfmbtKD_v_toySfmbtKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_SfmbtKD_v_toySfmbtKD_sig, "SfmbtKD.v.toy-SfmbtKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

#############################
#### ScmKD vs. toy-ScmKD ####
#############################

# Perform DE testing 
res_ScmKD_v_toyScmKD <- results(dds, contrast=c("tissue","toy-ScmKD", "ScmKD"))

# Add fold change classification to complete results
res_ScmKD_v_toyScmKD <- as.data.frame(res_ScmKD_v_toyScmKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_ScmKD_v_toyScmKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in toy-ScmKD",
    log2FoldChange <= -1 ~ "depleted in toy-ScmKD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_ScmKD_v_toyScmKD_sig <- subset(res_ScmKD_v_toyScmKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_ScmKD_v_toyScmKD_sig, "ScmKD.v.toy-ScmKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

##############################################################################
#### Group 5: What toy-dependent changes occur to allow for extra growth? ####
############################# PcKD vs. toy-PcKD ##############################
##############################################################################

# Perform DE testing 
res_PcKD_v_toyPcKD <- results(dds, contrast=c("tissue","toy-PcKD", "PcKD"))

# Add fold change classification to complete results
res_PcKD_v_toyPcKD <- as.data.frame(res_PcKD_v_toyPcKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_PcKD_v_toyPcKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in toy-PcKD",
    log2FoldChange <= -1 ~ "depleted in toy-PcKD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_PcKD_v_toyPcKD_sig <- subset(res_PcKD_v_toyPcKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_PcKD_v_toyPcKD_sig, "PcKD.v.toy-PcKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

############################################################
#### Group 6: Is there more than one way to skin a cat? ####
################## PcKD vs. toy-SfmbtKD ####################
############################################################

# Perform DE testing 
res_PcKD_v_toySfmbtKD <- results(dds, contrast=c("tissue","toy-SfmbtKD", "PcKD"))

# Add fold change classification to complete results
res_PcKD_v_toySfmbtKD <- as.data.frame(res_PcKD_v_toySfmbtKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_PcKD_v_toySfmbtKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in toy-SfmbtKD",
    log2FoldChange <= -1 ~ "depleted in toy-SfmbtKD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_PcKD_v_toySfmbtKD_sig <- subset(res_PcKD_v_toySfmbtKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_PcKD_v_toySfmbtKD_sig, "PcKD.v.toy-SfmbtKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

##########################
### PcKD vs. toy-ScmKD ###
##########################

# Perform DE testing 
res_PcKD_v_toyScmKD <- results(dds, contrast=c("tissue","toy-ScmKD", "PcKD"))

# Add fold change classification to complete results
res_PcKD_v_toyScmKD <- as.data.frame(res_PcKD_v_toyScmKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_PcKD_v_toyScmKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in toy-ScmKD",
    log2FoldChange <= -1 ~ "depleted in toy-ScmKD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_PcKD_v_toyScmKD_sig <- subset(res_PcKD_v_toyScmKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_PcKD_v_toyScmKD_sig, "PcKD.v.toy-ScmKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

##############################################################
#### Group 7: How does Pc eye-to-wing differ from others? ####
##################### PcKD vs. SfmbtKD #######################
##############################################################

# Perform DE testing 
res_PcKD_v_SfmbtKD <- results(dds, contrast=c("tissue","SfmbtKD", "PcKD"))

# Add fold change classification to complete results
res_PcKD_v_SfmbtKD <- as.data.frame(res_PcKD_v_SfmbtKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_PcKD_v_SfmbtKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in SfmbtKD",
    log2FoldChange <= -1 ~ "depleted in SfmbtKD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_PcKD_v_SfmbtKD_sig <- subset(res_PcKD_v_SfmbtKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_PcKD_v_SfmbtKD_sig, "PcKD.v.SfmbtKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)

##########################
### PcKD vs. ScmKD ###
##########################

# Perform DE testing 
res_PcKD_v_ScmKD <- results(dds, contrast=c("tissue","ScmKD", "PcKD"))

# Add fold change classification to complete results
res_PcKD_v_ScmKD <- as.data.frame(res_PcKD_v_ScmKD) %>%
  rownames_to_column("geneid") %>%
  mutate(res_PcKD_v_ScmKD_change = case_when(
    log2FoldChange >= 1 ~ "enriched in ScmKD",
    log2FoldChange <= -1 ~ "depleted in ScmKD",
    TRUE ~ "ns"
  ))

# Pull significantly changed genes
res_PcKD_v_ScmKD_sig <- subset(res_PcKD_v_ScmKD, padj < 0.05 & abs(log2FoldChange) >= 1)

# Write to table
write.table(res_PcKD_v_ScmKD_sig, "PcKD.v.ScmKD.diff.txt", sep="\t", col.names=T, row.names=F, quote=F)