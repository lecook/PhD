## 06 Feb 2019
## Laura E Cook, University of Melbourne
## This script is a work in progress for analysing the loss and/or gain of TFBS in the TWAR sequences

###### ------------------------------- LOAD IN PACKAGES AND DATA ----------------------------------- ######

library(plyr)
library(data.table) 
library(dplyr)
library(GenomicRanges)
library(stringr)
library(tidyverse)


# Set working directory
setwd("~/phd/OneDrive - The University of Melbourne/tfbs_project/results/meme/")

# Load in table of craniofacial motif matches from MEME FIMO for all species: Thylacine, Wolf, Panda and Devil
m <- read.table("all_CFmotifs_HOCOMOCO_core_80%_overlap_0.001.tsv", header=TRUE, sep="")



###### ------------------------------- FINDING UNIQUE MATCHES BASED ON INTERVALS ----------------------------------- ######

## Try and do this for every TWAR

for(twar in unique(m$sequence_name)){
	twar <- subset(m, sequence_name=="chr1")


	# Create an empty list
	my.list <- list()

	# Loop over every unique TF id within the list of matches for this TWAR
	for(tf in unique(twar$motif_id)){

	# Generate a list of the tf names
	list.name <- as.character(tf)

		# Convert dataframe into a GRanges object to find overlaps and unique sites to compare
		GR <- makeGRangesFromDataFrame(subset(twar, motif_id==tf), keep.extra.columns=TRUE, seqnames.field="sequence_name", start.field="start", end.field="stop", strand.field="strand")
		
		# Extract the interval data according to the metadata "species"
		ailMel1 <- GR[which(elementMetadata(GR)[, "species"] %in% "ailMel1")]
		tcyn <- GR[which(elementMetadata(GR)[, "species"] %in% "tcyn")]
		sarHar1 <- GR[which(elementMetadata(GR)[, "species"] %in% "sarHar1")]
		clup <- GR[which(elementMetadata(GR)[, "species"] %in% "clup")]
		

		# Extract the unique matches for each species
		# Can't work out a way to do this in one go so have to look at thylacine no devil, thylacine no panda and thylacine no wolf and then overlap those lists
		# subsetByOverlaps() looks for the overlapping sequences between two GRanges objects, invert=TRUE means that it does the opposite and outputs the sequence found in thylacine but not in devil
		# I've ignored strand (-/+) in this because it kept mucking up and I don't know that I'm necessarily interested in strand?? Check with Irene about this.


		## Unique to thylacine
		thylacine_not_devil <- subsetByOverlaps(tcyn, sarHar1, invert=TRUE, ignore.strand=TRUE)
		thylacine_not_panda <- subsetByOverlaps(tcyn, ailMel1, invert=TRUE, ignore.strand=TRUE)
		thylacine_not_wolf <- subsetByOverlaps(tcyn, clup, invert=TRUE, ignore.strand=TRUE)

		# Overlap the list above to find those matches unique to thylacine
		uT <- as(Reduce(subsetByOverlaps, list(thylacine_not_wolf, thylacine_not_panda, thylacine_not_devil)), "data.frame")

		# If there is data present (there isn't always a unique match for every TF) add a column to the data.frame so I can keep track of which match is unique in which species
		# If the number of rows in the dataframe is greater the zero add the column with species name
		# Repeat this for the other species - perhaps I can do this in a loop if I had a list of the species names?? ASK IRENE
		if(nrow(uT)>0){
			uT$unique <- "thylacine"
		}	

		## Unique to wolf
		wolf_not_devil <- subsetByOverlaps(clup, sarHar1, invert=TRUE, ignore.strand=TRUE)
		wolf_not_panda <- subsetByOverlaps(clup, ailMel1, invert=TRUE, ignore.strand=TRUE)
		wolf_not_thylacine <- subsetByOverlaps(clup, tcyn, invert=TRUE, ignore.strand=TRUE)
		uW <- as(Reduce(subsetByOverlaps, list(wolf_not_thylacine, wolf_not_panda, wolf_not_devil)), "data.frame")
		if(nrow(uW)>0){
			uW$unique <- "wolf"
		}	

		## Unique to panda
		panda_not_devil <- subsetByOverlaps(ailMel1, sarHar1, invert=TRUE, ignore.strand=TRUE)
		panda_not_thylacine <- subsetByOverlaps(ailMel1, tcyn, invert=TRUE, ignore.strand=TRUE)
		panda_not_wolf <- subsetByOverlaps(ailMel1, clup, invert=TRUE, ignore.strand=TRUE)
		uP <- as(Reduce(subsetByOverlaps, list(panda_not_wolf, panda_not_thylacine, panda_not_devil)), "data.frame")
		if(nrow(uP)>0){
			uP$unique <- "panda"
		}	

		## Unique to devil
		devil_not_thylacine <- subsetByOverlaps(sarHar1, tcyn, invert=TRUE, ignore.strand=TRUE)
		devil_not_panda <- subsetByOverlaps(sarHar1, ailMel1, invert=TRUE, ignore.strand=TRUE)
		devil_not_wolf <- subsetByOverlaps(sarHar1, clup, invert=TRUE, ignore.strand=TRUE)
		uD <- as(Reduce(subsetByOverlaps, list(devil_not_wolf, devil_not_panda, devil_not_thylacine)), "data.frame")
		if(nrow(uD)>0){
			uD$unique <- "devil"
		}

		# Look at some other combinations

		## Unique to thylacine and wolf
		uTW <-as(Reduce(subsetByOverlaps, list(thylacine_not_panda, thylacine_not_devil, clup)), "data.frame")
		if(nrow(uTW)>0){
			uTW$unique <- "thylacine_wolf"
		}

		## Unique to thylacine and panda
		uTP <-as(Reduce(subsetByOverlaps, list(thylacine_not_wolf, thylacine_not_devil, ailMel1)), "data.frame")
		if(nrow(uTP)>0){
			uTP$unique <- "thylacine_panda"
		}

		## Unique to thylacine and devil
		uTD <-as(Reduce(subsetByOverlaps, list(sarHar1, thylacine_not_panda, thylacine_not_wolf)), "data.frame")
		if(nrow(uTD)>0){
			uTD$unique <- "thylacine_devil"
		}

		## Unique to wolf and panda
		uWP <-as(Reduce(subsetByOverlaps, list(wolf_not_thylacine, wolf_not_devil, ailMel1)), "data.frame")
		if(nrow(uWP)>0){
			uWP$unique <- "wolf_panda"
		}

		## Unique to wolf and panda
		uWD <-as(Reduce(subsetByOverlaps, list(sarHar1, wolf_not_thylacine, wolf_not_panda)), "data.frame")
		if(nrow(uWD)>0){
			uWD$unique <- "wolf_devil"
		}

		## Unique to panda and devil
		uPD <-as(Reduce(subsetByOverlaps, list(sarHar1, panda_not_thylacine, panda_not_wolf)), "data.frame")
		if(nrow(uPD)>0){
			uPD$unique <- "panda_devil"
		}


		## Match is present in all four species
		uTWPD <- as(Reduce(subsetByOverlaps, list(tcyn, ailMel1, sarHar1, clup)), "data.frame")
		if(nrow(uTWPD)>0){
			uTWPD$unique <- "thylacine_wolf_panda_devil"
		}

		# For each TF bind all the combinations into 1 dataframe
		df_temp <- rbind(uT, uW, uD, uP, uTW, uTP, uTD, uWP, uWD, uPD, uTWPD)
		
		# Add the TF that the loop is on to the dataframe so I know which TF this data is for
		if(nrow(df_temp)>0){
			df_temp$TF <- tf
		}


			# For each TF in my list add the information from the dataframe (dt_temp) 
			my.list[[list.name]] <- df_temp

	}

	# Combine all the data for every TF in the list
	final.data <- do.call(rbind, my.list)
	write.table(final.data, file="final_data.txt")


	### Need to work out a way to keep both regions that overlap
	### So that I can later compare the scores for the same regions


	###### ------------------------------- COUNTING MATCHES ----------------------------------- ######

	# Create an empty list
	my.list2 <- list()

	# Loop over every unique TF id within the list of matches for this TWAR
	for(tf in unique(final.data$TF)){

	# Generate a list of the tf names and a list of the unique species
	list.name2 <- as.character(tf)
		
		# Subset the data by tf
		tf_subset <- subset(final.data, TF == tf)
		tf_twar_subset <- subset(twar, motif_id == tf)

		# Count occurances in unique matches dataset (final.data)
		uT_count <- sum(str_count(tf_subset$unique, "\\bthylacine\\b"))
		uD_count <-  sum(str_count(tf_subset$unique, "\\bdevil\\b"))
		uP_count <- sum(str_count(tf_subset$unique, "\\bpanda\\b"))
		uW_count <- sum(str_count(tf_subset$unique, "\\bwolf\\b"))
		uTW_count <- sum(str_count(tf_subset$unique, "\\bthylacine_wolf\\b"))
		uTP_count <- sum(str_count(tf_subset$unique, "\\bthylacine_panda\\b"))
		uTD_count <- sum(str_count(tf_subset$unique, "\\bthylacine_devil\\b"))
		uWD_count <- sum(str_count(tf_subset$unique, "\\bwolf_devil\\b"))
		uPD_count <- sum(str_count(tf_subset$unique, "\\bpanda_devil\\b"))
		uTWPD_count <- sum(str_count(tf_subset$unique, "\\bthylacine_wolf_panda_devil\\b"))

		# Count all occurances from original dataset
		T <- sum(str_count(tf_twar_subset$species, "\\btcyn\\b"))
		W <- sum(str_count(tf_twar_subset$species, "clup"))
		P <- sum(str_count(tf_twar_subset$species, "ailMel1"))
		D <- sum(str_count(tf_twar_subset$species, "sarHar1"))


		# Save these into a new df
		df_count <- data.frame(T, W, P, D, uT_count, uW_count, uD_count, uP_count, uTW_count, uTP_count, uTD_count, uWD_count, uPD_count, uTWPD_count)
			my.list2[[list.name2]] <- df_count
	}


	# Combine all the count data for every TF in the list
	count.data <- do.call(rbind, my.list2)
	write.table(count.data, file="count_data.txt")


## DONE! For now.. urghh

### Work out how to do this for ALL the twars.. maybe need to write it as a function? 
### Then have a list of TWARs and list of TFs and loop through both lists