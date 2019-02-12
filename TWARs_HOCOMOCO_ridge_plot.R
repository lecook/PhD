## Laura E Cook, University of Melbourne
## 17 Oct 2018
### Makes a ridge plot for all craniofacial PWMs which were found to have matched binding sites in TWAR sequences from MEME FIMO output

library(plyr)
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(ggridges)
library(reshape2)
library(RColorBrewer)
library(scales) # for ggplot
library(gridExtra) # for ggplot 


setwd("~/phd/OneDrive - The University of Melbourne/tfbs_project/results/meme")

options(width=150)


## Need to plot the distribution of scores for each TF PWM for each species

## Need column with scores, column with PWM name, column with species


tcyn <- read.table("tcyn_CFmotifs_HOCOMOCO_meme_fimo_0.01_80percent_filtered_output.tsv", stringsAsFactors=F, sep="\t", header=T)
tcyn["species"] <- "tcyn"

clup <- read.table("Clup_CFmotifs_HOCOMOCO_meme_fimo_0.01_80percent_filtered_output.tsv", stringsAsFactors=F, sep="\t", header=T)
clup["species"] <- "clup"

ailMel1 <- read.table("ailMel1_CFmotifs_HOCOMOCO_meme_fimo_0.01_80percent_filtered_output.tsv", stringsAsFactors=F, sep="\t", header=T)
ailMel1["species"] <- "ailMel1"

sarHar1 <- read.table("sarHar1_CFmotifs_HOCOMOCO_meme_fimo_0.01_80percent_filtered_output.tsv", stringsAsFactors=F, sep="\t", header=T)
sarHar1["species"] <- "sarHar1"


# Combine species tables
data <- rbind(tcyn, clup, ailMel1, sarHar1)

# Clean up motif column
data$motif_id <- gsub(pattern="_H.*", replacement="", x=data$motif_id)

pdf("ridge_plots_all2.pdf", height=6)
    ggplot(data, aes(x=score, y=motif_id, fill=species)) +
        geom_density_ridges(scale=4, alpha=0.3) +
        #geom_vline(xintercept=log(0.99/0.01), linetype="dashed", color="grey30") +
        scale_fill_brewer(palette="Dark2", labels=c("thylacine", "wolf", "panda", "devil"), guide="legend") +
        theme_ridges() +
        theme(axis.text.y = element_text(size=rel(0.5))) +
        xlim(-20,20)
dev.off()


