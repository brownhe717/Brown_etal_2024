#library(tsrexplorer)
library(GenomicFeatures)
library(Gviz)
library(viridis)
library(scales)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

# This script imports signal tracks and displays them at the specified regions of the Drosophila 

gviz_dir = file.path("~/Desktop/Brown et al 2024/scripts/Figure 7")
bigwig_dir = file.path("~/Desktop/Brown et al 2024/scripts/Figure 7/INPUT/bigs")

# Create genomic axis track
axis.track <- GenomeAxisTrack(col = "black", scale = 0.1, col.range = "black")

options(ucscChromosomeNames = FALSE)

# Create gene annotation track
txdb <- makeTxDbFromGFF("~/Desktop/Brown et al 2024/scripts/Figure 7/INPUT/Dm_genome/genes.gtf")

genome.track <- GeneRegionTrack(txdb, genome = "dm6", shape = "arrow", names = "Genes", col = "black",
                                showId = FALSE, fill = "black", trancriptAnnotation = "gene_symbol", collapseTranscripts = "meta")


# Create data tracks

# Get colors from the viridis palette for the number of tracks to be plotted
show_col(viridis_pal()(40))


#############################
### Set y-axis limits EAD ###
#############################


#H3K27me3 tracks
WT_pos_lim = c(0,25000)
toy_pos_lim = c(0,25000)
Pc_pos_lim = c(0,25000)
toyPc_pos_lim = c(0,25000)
Sfmbt_pos_lim = c(0,25000)
toySfmbt_pos_lim = c(0,25000)

#IgG tracks
IgG_WT_pos_lim = c(0,25000)
IgG_toy_pos_lim = c(0,25000)
IgG_Pc_pos_lim = c(0,25000)
IgG_toyPc_pos_lim = c(0,25000)
IgG_Sfmbt_pos_lim = c(0,25000)
IgG_toySfmbt_pos_lim = c(0,25000)


# H3K27me3
WT_EAD <- DataTrack(range = file.path(bigwig_dir, "WT_H3K27me3_1_SF.bigwig"), genome = "dm6", 
                          name = "DE-GAL4", col.histogram = "#440154FF", fill.histogram = "#440154FF", ylim = WT_pos_lim)
#try changing 120hr K27me3 to #404588FF

toy_EAD <- DataTrack(range = file.path(bigwig_dir, "toyKD_H3K27me3_2_SF.bigwig"), genome = "dm6", 
                    name = "toyKD", col.histogram = "#453581FF", fill.histogram = "#453581FF", ylim = toy_pos_lim)
#try changing 120hr K27me3 to #404588FF

Pc_EAD <- DataTrack(range = file.path(bigwig_dir, "PcKD_H3K27me3_1_SF.bigwig"), genome = "dm6", 
                    name = "PcKD", col.histogram = "#3A538BFF", fill.histogram = "#3A538BFF", ylim = Pc_pos_lim)
#try changing 120hr K27me3 to #404588FF

toyPc_EAD <- DataTrack(range = file.path(bigwig_dir, "toyPcKD_H3K27me3_1_SF.bigwig"), genome = "dm6", 
                    name = "toyPcKD", col.histogram = "#26818EFF", fill.histogram = "#26818EFF", ylim = toyPc_pos_lim)
#try changing 120hr K27me3 to #404588FF

Sfmbt_EAD <- DataTrack(range = file.path(bigwig_dir, "SfmbtKD_H3K27me3_2_SF.bigwig"), genome = "dm6", 
                    name = "SfmbtKD", col.histogram = "#40BC72FF", fill.histogram = "#40BC72FF", ylim = Sfmbt_pos_lim)
#try changing 120hr K27me3 to #404588FF

toySfmbt_EAD <- DataTrack(range = file.path(bigwig_dir, "toySfmbtKD_H3K27me3_1_SF.bigwig"), genome = "dm6", 
                    name = "toySfmbtKD", col.histogram = "#76D153FF", fill.histogram = "#76D153FF", ylim = toySfmbt_pos_lim)
#try changing 120hr K27me3 to #404588FF



# IgG
#IgG_EAD <- DataTrack(range = file.path(bigwig_dir4, "10WT_3rd_instar_IgG_Ub_2_normalized.bigwig"), genome = "dm6", 
#name = "IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_EAD_pos_lim)
IgG_WT <- DataTrack(range = file.path(bigwig_dir, "WT_IgG_H3K27me3_2_SF.bigwig"), genome = "dm6", 
                         name = "WT 3rd Instar IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_WT_pos_lim)

IgG_toy <- DataTrack(range = file.path(bigwig_dir, "toyKD_IgG_H3K27me3_3_SF.bigwig"), genome = "dm6", 
                    name = "WT 3rd Instar IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_toy_pos_lim)

IgG_Pc <- DataTrack(range = file.path(bigwig_dir, "PcKD_IgG_H3K27me3_1_SF.bigwig"), genome = "dm6", 
                    name = "WT 3rd Instar IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_Pc_pos_lim)

IgG_toyPc <- DataTrack(range = file.path(bigwig_dir, "toyPcKD_IgG_H3K27me3_1_SF.bigwig"), genome = "dm6", 
                    name = "WT 3rd Instar IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_toyPc_pos_lim)

IgG_Sfmbt <- DataTrack(range = file.path(bigwig_dir, "SfmbtKD_IgG_H3K27me3_2_SF.bigwig"), genome = "dm6", 
                    name = "WT 3rd Instar IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_Sfmbt_pos_lim)

IgG_toySfmbt <- DataTrack(range = file.path(bigwig_dir, "toySfmbtKD_IgG_H3K27me3_3_SF.bigwig"), genome = "dm6", 
                    name = "WT 3rd Instar IgG", col.histogram = "#FDE725FF", fill.histogram = "#FDE725FF", ylim = IgG_toySfmbt_pos_lim)

#####################
### Overlay Track ###
#####################

displayPars(IgG_WT) = list(alpha.title=1,alpha=1)
displayPars(WT_EAD) = list(alpha.title=1, alpha=1)
ot_WT_EAD = OverlayTrack(trackList = list(WT_EAD, IgG_WT))

displayPars(IgG_toy) = list(alpha.title=1,alpha=1)
displayPars(toy_EAD) = list(alpha.title=1, alpha=1)
ot_toy_EAD = OverlayTrack(trackList = list(toy_EAD, IgG_toy))

displayPars(IgG_Pc) = list(alpha.title=1,alpha=1)
displayPars(Pc_EAD) = list(alpha.title=1, alpha=1)
ot_Pc_EAD = OverlayTrack(trackList = list(Pc_EAD, IgG_Pc))

displayPars(IgG_toyPc) = list(alpha.title=1,alpha=1)
displayPars(toyPc_EAD) = list(alpha.title=1, alpha=1)
ot_toyPc_EAD = OverlayTrack(trackList = list(toyPc_EAD, IgG_toyPc))

displayPars(IgG_Sfmbt) = list(alpha.title=1,alpha=1)
displayPars(Sfmbt_EAD) = list(alpha.title=1, alpha=1)
ot_Sfmbt_EAD = OverlayTrack(trackList = list(Sfmbt_EAD, IgG_Sfmbt))

displayPars(IgG_toySfmbt) = list(alpha.title=1,alpha=1)
displayPars(toySfmbt_EAD) = list(alpha.title=1, alpha=1)
ot_toySfmbt_EAD = OverlayTrack(trackList = list(toySfmbt_EAD, IgG_toySfmbt))


###################
### Plot Tracks ###
###################

#all vg locus
plotTracks(list(axis.track, 
                ot_WT_EAD, #6,000
                ot_toy_EAD, #12,000
                ot_Pc_EAD, #4,000
                ot_toyPc_EAD, #12,000
                ot_Sfmbt_EAD, #25,000
                ot_toySfmbt_EAD, #8,000
                genome.track), 
           chromosome = "chr2R",from = 12884201, to = 12900600, 
           background.title = "white", 
           col.title = "black", 
           col.axis= "black", 
           type="histogram", 
           baseline = 0, 
           col.baseline = "black"
)

#DE>PcG KD vg locus
plotTracks(list(axis.track, 
                ot_WT_EAD, #6,000
                ot_Pc_EAD, #4,000
                ot_Sfmbt_EAD, #25,000
                genome.track), 
           chromosome = "chr2R",from = 12884201, to = 12900600, 
           background.title = "white", 
           col.title = "black", 
           col.axis= "black", 
           type="histogram", 
           baseline = 0, 
           col.baseline = "black"
)

#DE-toy>PcG KD vg locus
plotTracks(list(axis.track, 
                ot_toy_EAD, #6,000
                ot_toyPc_EAD, #4,000
                ot_toySfmbt_EAD, #25,000
                genome.track), 
           chromosome = "chr2R",from = 12884201, to = 12900600, 
           background.title = "white", 
           col.title = "black", 
           col.axis= "black", 
           type="histogram", 
           baseline = 0, 
           col.baseline = "black"
)
