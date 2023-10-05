#data processing for trawl 2018 CO1 full dataset

####Libraries####
library(dada2)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(ape)
library(qualpalr)
library(viridis)
library(ShortRead)
library(Biostrings)
library(seqinr)

####Environment Setup####
theme_set(theme_bw())
setwd("~/projects/CO1_runs/trawl-2018/dada2_all")

####File Path Setup####
####NB: IN THIS DATASET WE ARE USING PREVIOUSLY PREPARED INPUT FILES FOR THE FILTERANDTRIM() COMMAND, THEY ARE ALREADY IN THE INPUT FOLDERS, SO WE CAN JUST CREATE THE FILE PATH INPUTS AND SKIP THE FILE CREATION STEPS (N FILTERING, PRIMER TRIMMING, ETC.)
####IMPORTANT NOTE: this involves combining processed fastq files (we started out combining raw fastqs and processing them as one batch, but some lanes require different trimming parameters and the best read-retention can be achieved by processing each lane individually and combining reads post-processing steps)
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- "~/projects/CO1_runs/trawl-2018/demultiplexed/merged"
list.files(path)
fnFs <- sort(list.files(path, pattern="R1", full.names = TRUE)) #change the pattern to match all your R1 files
fnRs <- sort(list.files(path, pattern="R2", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 2) #change the delimiter in quotes and the number at the end of this command to decide how to split up the file name, and which element to extract for a unique sample name

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # creating fn*s.filtN objects so they refer to existing filtered files
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE)) # creating cut*s objects so they refer to existing filtered files
cutRs <- sort(list.files(path.cut, pattern = "R2", full.names = TRUE))

####filter and trim reads####
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

####trim & filter####
#filter and trim command. dada2 can canonically handle lots of errors, I am typically permissive in the maxEE parameter set here, in order to retain the maximum number of reads possible. error correction steps built into the dada2 pipeline have no trouble handling data with this many expected errors.
#it is best, after primer removal, to not truncate with 18s data, or with data from any region in which the length is broadly variable. you may exclude organisms that have a shorter insert than the truncation length (definitely possible, good example is giardia). defining a minimum sequence length is best.
#trimright 10 bp to improve merging rates downstream
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, trimRight = c(10,10), truncLen=c(0,0), minLen = c(150,150),
                     maxN=c(0,0), maxEE=c(4,6), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE,
                     compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100

####learn error rates####
#the next three sections (learn error rates, dereplication, sample inference) are the core of dada2's sequence processing pipeline. read the dada2 paper and their online documentation (linked at top of this guide) for more information on how these steps work
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

pdf("error_rates.dada2.R1s.pdf", width = 10, height = 10) # define plot width and height. completely up to user.
  plotErrors(errF, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered before moving on.
dev.off()
pdf("error_rates.dada2.R2s.pdf", width = 10, height = 10) # define plot width and height. completely up to user.
  plotErrors(errR, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered before moving on.
dev.off()

####dereplication####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names #this is just to ensure that all your R objects have the same sample names in them
names(derepFs) <- sample.names
names(derepRs) <- sample.names

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

####OPTIONAL: remove low-sequence samples before merging####
#a "subscript out of bounds" error at the next step (merging) may indicate that you aren't merging any reads in one or more samples.
#NB, NOT getting this error doesn't necessarily mean that all of your samples end up with more than 0 merged reads, as i found out while processing a large 18s dataset. your guess is as good as mine as to why this error does or does not appear, but filtering out the samples that cause it is necessary for completion of the pipeline.
#samples_to_keep <- as.numeric(out[,"reads.out"]) > 100 #example of simple method used above after the filter and trim step. if you already did this but still got an error when merging, try the steps below
getN <- function(x) sum(getUniques(x)) #keeping track of read retention, number of unique sequences after ASV inference
track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
samples_to_keep <- track[,4] > 50 #your threshold. try different ones to get the lowest one that will work. #this method accounts for dereplication/ASVs left after inference
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)] #record names of samples you have the option of removing

####merge paired reads####
mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

####construct sequence table####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

####View Sequence Length Distribution Post-Merging####
#most useful with merged data. this plot will not show you much for forward reads only, which should have a uniform length distribution.
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
pdf("length_histogram.merged_reads.pdf", width = 10, height = 8) # define plot width and height. completely up to user.
  plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot
dev.off()

####remove low-count singleton ASVs####
#create phyloseq otu_table
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

#some metrics from the sequence table
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV
length(otu_pres_abs_rowsums) #how many ASVs
#23558
length(which(otu_pres_abs_rowsums == 1)) #how many ASVs only present in one sample
#34170

#what are the counts of each ASV
otu_rowsums <- rowSums(otus) #raw counts per ASV
otu_singleton_rowsums <- as.data.frame(otu_rowsums[which(otu_pres_abs_rowsums == 1)]) #raw read counts in ASVs only presesnt in one sample
hist(otu_singleton_rowsums[,1], breaks=500, xlim = c(0,200), xlab="# Reads in ASV") #histogram plot of above
length(which(otu_singleton_rowsums <= 1)) #how many are there with N reads or fewer? (N=1 in example)
#[1] 227

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "transform_sample_counts" to create a relative abundance table, and then filter your ASVs by choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(seqtab) # sanity check
dim(otus) # (this should be the same as last command, but the dimensions reversed)
#[1] 34170   190
otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step (a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need to be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those only present in one sample, this is a value we can use to filter them for relative abundance on a per-sample basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.001) #here is where you set your relative abundance threshold #which ASVs pass our filter for relative abundance
length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
#[1] 19397
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abundance value than our threshold)
otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
dim(otus_filt) #how many ASVs did you retain?
#[1] 13320   190
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=TRUE, verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
dim(seqtab.nosingletons.nochim)
#190 13320
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons) #proportion of nonchimeras #it should be relatively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low
#[1] 0.9765551 data looks to be good quality w/r/t chimeras

####track read retention through steps####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab.nosingletons), rowSums(seqtab.nosingletons.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100, track[,7]/track[,1]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nosingletons", "nochimeras", "percent_singletons", "percent_chimeras", "percent_retained_of_total")

#overall, Run 1 (sept) has higher proportions of both singleton ASVs, and chimeric ASVs than Run 2 (Oct), but read retention is higher for Run 1, likely meaning that read quality was higher in Run1 than Run2.
#overall, read depth appears greater in Run 1 than Run 2

####save output from sequnce table construction steps####
write.table(data.frame("row_names"=rownames(track),track),"read_retention.CO1_merged.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"sequence_table.CO1.merged.txt", row.names=FALSE, quote=F, sep="\t")

#if you must save your sequence table and load it back in before doing taxonomy assignments, here is how to read and format the object so that dada2 will accept it again
seqtab.nosingletons.nochim <- fread("sequence_table.CO1.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.nosingletons.nochim) <- seqtab.nosingletons.nochim[,1] #set row names
seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,-1] #remove column with row names in it
seqtab.nosingletons.nochim <- as.matrix(seqtab.nosingletons.nochim) #cast the object as a matrix
mode(seqtab.nosingletons.nochim) <- "numeric"


##### now replace the long ASV names (the actual sequences) with human-readable names####
#save the new names and sequences as a .fasta file in your project working directory, and save a table that shows the mapping of sequences to new ASV names
my_otu_table <- t(as.data.frame(seqtab.nosingletons.nochim)) #transposed (OTUs are rows) data frame. unclassing the otu_table() output avoids type/class errors later on
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') #create new names
write.table(cbind(ASV.num, ASV.seq), "sequence_ASVname_mapping.txt", sep="\t", quote=F, row.names=F, col.names=F)
write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "CO1_ASV_sequences.fasta") #save sequences with new names in fasta format
#IMPORTANT: sanity checks
colnames(seqtab.nosingletons.nochim) == ASV.seq #only proceed if this tests as true for all elements

#rename your ASVs in the taxonomy table and sequence table objects
seqtab.withnames <- seqtab.nosingletons.nochim
colnames(seqtab.withnames) <- ASV.num

#re-save sequence table with updated names
write.table(data.frame("row_names"=rownames(seqtab.withnames),seqtab.withnames),"sequence_table.CO1.merged.corrected_ASV_names.txt", row.names=FALSE, quote=F, sep="\t")


####in this section we are collapsing ASVs with vsearch and a custom script, because some are in reverse-complement orientation. alternative method rather than using swarm to process this data ####

#cluster ASVs with vsearch at 100% (this collapses ASVs with opposite orientations and outputs all sequences in the same orientation, if --strand both is specified), then collapse OTU table based on the clustering map produced
vsearch --cluster_fast CO1_ASV_sequences.fasta --iddef 1 --id 1 --strand both --centroids CO1_ASV_sequences.collapsed.fasta --uc alignment_orientation_correction/clusters.uc
#use mafft to make an alignment of the orientation-corrected sequences. sanity check to see if anything is still in the opposite orientation
mafft --thread 38 --reorder --auto CO1_ASV_sequences.collapsed.fasta > CO1_ASV_sequences.collapsed.aligned.fasta #view with AliView, re-orient sequences that appear to be in the wrong orientation, remove gaps, save over CO1_ASV_sequences.oriented.fasta

#reoriented file is saved as:
CO1_ASV_sequences.oriented.fasta
#and R-compatible phylip format saved as:
CO1_ASV_sequences.oriented.phy

#next step is to collapse ASVs based on the clusters generated. i think i can do this with my old python code from the MED/QIIME pipeline. let's check.
#inputs for this script are: 
#--matrix_count (a tab-delimited OTU table with OTUs in rows and samples in columns)
#--otu_map (tab delimited map of OTUs that should be collapsed, no headers. can be a modified version of the the clusters.uc file from usearch output above)
#--taxonomy_assignments (tab delimited taxonomy assignments for the OTUs you want to keep, no headers. maybe rather than modifying this script too heavily, i will dummy up a fake file for this, we will see, since i think i depend on it for which OTUs to actually keep in the output)
#--outputfile just the name of your output file, a collapsed OTU table based on the map you provide, with taxonomy assignments for remaining OTUs

#need to transpose our sequence table for this
awk '{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}' sequence_table.CO1.merged.corrected_ASV_names.txt > sequence_table.CO1.merged.corrected_ASV_names.transposed.txt

#need to create a map file from the .uc file by deleting the first 8 columns and retaining only those with something in column 10
cut -f 9- alignment_orientation_correction/clusters.uc | grep -v "*" | awk '{print $2"\t"$1}' | sort | awk -F'\t' -v OFS='\t' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' | tr -s '\t' > alignment_orientation_correction/otu_map.txt

python2 ~/programs/hakai_genomics_repo/pipeline_scripts/collapse_OTU_table.py --matrix_count sequence_table.CO1.merged.corrected_ASV_names.transposed.txt \
--otu_map alignment_orientation_correction/otu_map.txt \
--outputfile sequence_table.CO1.merged.corrected_ASV_names.collapsed.txt


##### now replace the long ASV names (the actual sequences) with human-readable names####
#save the new names and sequences as a .fasta file in your project working directory, and save a table that shows the mapping of sequences to new ASV names
my_otu_table <- t(as.data.frame(seqtab.nosingletons.nochim.collapsed)) #transposed (OTUs are rows) data frame. unclassing the otu_table() output avoids type/class errors later on
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') #create new names
write.table(cbind(ASV.num, ASV.seq), "sequence_ASVname_mapping.final.txt", sep="\t", quote=F, row.names=F, col.names=F)
write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "CO1_ASV_sequences.final.fasta") #save sequences with new names in fasta format
#IMPORTANT: sanity checks
colnames(seqtab.nosingletons.nochim.collapsed) == ASV.seq #only proceed if this tests as true for all elements

#rename your ASVs in the taxonomy table and sequence table objects
colnames(seqtab.nosingletons.nochim.collapsed) <- ASV.num


#re-save sequence table with updated names
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim.collapsed),seqtab.nosingletons.nochim.collapsed),"sequence_table.CO1.merged.collapsed.w_ASV_names.txt", row.names=FALSE, quote=F, sep="\t")


#### taxonomy assignment with blast and LCA, Nov 2022 ####
#assign taxonomy for dada2-processed CO1 amplicon data with blast using NCBI NT and custom CO1 databases together
mkdir blast_96_sim
# blast against custom blast DB
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/CO1_database/blast_DB/CO1.BOLD_genbank_combined.rep_set.blast_DB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query CO1_ASV_sequences.final.fasta  -out blast_96_sim/CO1_ASV_sequences.customDB.blast.out
# blast against genbank NT blast DB
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/NCBI_NT/2021-11-05/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query CO1_ASV_sequences.final.fasta  -out blast_96_sim/CO1_ASV_sequences.blast.out


#filter input fasta to only contain sequences with no hits in the above two blast runs
cat blast_96_sim/CO1_ASV_sequences*blast.out | cut -f1,1 | sort | uniq > blast_96_sim/blast_hit_ASVs.txt
grep -wv -f blast_96_sim/blast_hit_ASVs.txt sequence_ASVname_mapping.final.txt | cut -f1,1 | sort > blast_96_sim/no_blast_hit_ASVs.txt
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' blast_96_sim/no_blast_hit_ASVs.txt CO1_ASV_sequences.final.fasta > blast_96_sim/CO1_ASV_sequences.no_blast_hit.fasta

#blast this output with lower thresholds for similarity
mkdir blast_90_sim
# blast against custom blast DB
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 90 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/CO1_database/blast_DB/CO1.BOLD_genbank_combined.rep_set.blast_DB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query blast_96_sim/CO1_ASV_sequences.no_blast_hit.fasta  -out blast_90_sim/CO1_ASV_sequences.customDB.blast.out
# blast against genbank NT blast DB
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 90 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/NCBI_NT/2021-11-05/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query blast_96_sim/CO1_ASV_sequences.no_blast_hit.fasta  -out blast_90_sim/CO1_ASV_sequences.blast.out


#IMPORTANT: next steps are done in R, for simplicity's sake. a custom script here would also work. this is simpler.
#we need to add taxonIDs for the customDB (adding them directly to the blast DB has not worked in the past, they don't get returned in the blast output). Using the blast output and a map of accessions to taxonIDs, we add the taxonIDs for each blast result.
library(tidyverse)
library(data.table)
library(ShortRead)
library(Biostrings)
library(seqinr)
taxonidmap <- read.delim("~/projects/taxonomyDBs/CO1_database/taxonID_map/CO1.BOLD_genbank_combined.taxonID_map.w_hakai_barcodes.txt", sep="\t", header=F)
blastfile <- read.delim("blast_96_sim/CO1_ASV_sequences.customDB.blast.out", sep="\t", header=F) #read in blast results for 96% similarity iteration
blastfile2 <- read.delim("blast_90_sim/CO1_ASV_sequences.customDB.blast.out", sep="\t", header=F) #read in blast results for 90% similarity iteration
blastfile <- rbind(blastfile, blastfile2) #join iterations for customDB blast
colnames(taxonidmap) <- c("accession", "taxonID")
colnames(blastfile) <- c("asv", "col2", "accession", "blasttaxid", "col5", "col6", "col7", "col8")
taxonidmap$accession <- trimws(taxonidmap$accession, which = c("both"))
blastfile_wtaxids <- merge(blastfile,taxonidmap, by="accession", all.x=TRUE)
blastfile_output <- subset(blastfile_wtaxids, select=c("asv", "col2", "accession", "taxonID", "col5", "col6", "col7", "col8")) #it's okay here that "col2" is just all NA values, we need it to conform to the blast output from the NT database, but it doesn't get used by the taxonomy assignment scripts
blastfile_output <- blastfile_output[order(blastfile_output$asv),]
write.table(blastfile_output, "blast_96_sim/CO1_ASV_sequences.customDB_96_90.blast.out", row.names=F, col.names=F, quote=F, sep="\t") #overwriting input

#combine blast results in a way that the add taxonomy and LCA scripts can handle
blastout_customDB <- read.delim("blast_96_sim/CO1_ASV_sequences.customDB_96_90.blast.out", sep="\t", header=F)
blastout_NCBINT <- read.delim("blast_96_sim/CO1_ASV_sequences.blast.out", sep="\t", header=F)
blastout_NCBINT_2 <- read.delim("blast_90_sim/CO1_ASV_sequences.blast.out", sep="\t", header=F)
blastout_NCBINT <- rbind(blastout_NCBINT, blastout_NCBINT_2) #join iterations for NT blast
blastout_combined <- rbind(blastout_customDB, blastout_NCBINT) #combine tables
tmp <- blastout_combined[order(-blastout_combined$V5),] #order descending order for percent identity
tmp <- tmp[order(tmp$V1),] #order by ASV name
blastout_combined <- tmp 
#write to file
write.table(blastout_combined, "blast_96_sim/CO1_ASV_sequences.combined_all.blast.out", row.names=F, col.names=F, quote=F, sep="\t")
#now quit R and continue with the remaining code in your bash shell

#avoid " keyerror: 'NA' " with python script by filtering on "NA" as a wholeword string
#explanation: occasionally, an output line from blast has an NA, which causes an error with the Simple-LCA script below. quick fix is to remove these lines (they're quite rare anyway)
grep -v -w "NA" blast_96_sim/CO1_ASV_sequences.combined_all.blast.out > blast_96_sim/tmp

#execute first step for the LCA program (adding taxonomy strings based on taxonIDs in blast output)
python2 ~/programs/galaxy-tool-BLAST/blastn_add_taxonomy_lite.py -i blast_96_sim/tmp -t ~/projects/taxonomyDBs/NCBI_taxonomy/2021-11-05/rankedlineage.dmp -m ~/projects/taxonomyDBs/NCBI_taxonomy/2021-11-05/merged.dmp -o taxonomy
cat <(head -n 1 ~/programs/galaxy-tool-lca/example/example.tabular) taxonomy_tmp > tmp

#in this section, a series of taxonomy string modifications are made. This makes the final output more readable/easier to understand, but is optional.
#IMPORTANT: please note that the filtering criteria in the final step depend on some of this filtering (e.g. blast hits with the word "phylum" will be removed. see -fh parameter in final step) 
# because i'm replacing "unknown phylum" in these sed commands, these sequences are retained.
# if you choose not to do the replacement, the blast hits with "unknown plylum" will not be used in LCA determination unless you also change the filtering criteria set in the -fh parameter during the final step.
# also note, "unknown phylum" is present in any taxonomy where the clade's phylum is uncertain in the NCBI taxonomy system, it doesn't indicate any other kind of uncertainty about the provenance of the sequence.

#label fix for clades missing "kingdom" label
sed -i 's/unknown kingdom \/ Bacillariophyta/Eukaryota \/ Bacillariophyta/g' tmp #Bacillariophyta
sed -i 's/unknown kingdom \/ Ciliophora/Eukaryota \/ Ciliophora/g' tmp #Ciliophora
sed -i 's/unknown kingdom \/ Discosea/Eukaryota \/ Discosea/g' tmp #Discosea
sed -i 's/unknown kingdom \/ Evosea/Eukaryota \/ Evosea/g' tmp #Evosea
sed -i 's/unknown kingdom \/ Haptista/Eukaryota \/ Haptista/g' tmp #Haptista
sed -i 's/unknown kingdom \/ Rhodophyta/Eukaryota \/ Rhodophyta/g' tmp #Rhodophyta

#and for those missing kingdom + phylum labels
sed -i 's/Eukaryota \/ unknown phylum \/ Chrysophyceae/Eukaryota \/ Chrysophyceae \/ Chrysophyceae/g' tmp #Chrysophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Cryptophyceae/Eukaryota \/ Cryptophyceae \/ Cryptophyceae/g' tmp #Cryptophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Oomycota/Eukaryota \/ Oomycota \/ Oomycota/g' tmp #Oomycota
sed -i 's/Eukaryota \/ unknown phylum \/ Phaeophyceae/Eukaryota \/ Phaeophyceae \/ Phaeophyceae/g' tmp #Phaeophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Phaeophyceae/Eukaryota \/ Phaeophyceae \/ Phaeophyceae/g' tmp #Phaeophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Bigyra/Eukaryota \/ Bigyra \/ Bigyra/g' tmp #Bigyra
sed -i 's/Eukaryota \/ unknown phylum \/ Dictyochophyceae/Eukaryota \/ Dictyochophyceae \/ Dictyochophyceae/g' tmp #Dictyochophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Dinophyceae/Eukaryota \/ Dinophyceae \/ Dinophyceae/g' tmp #Dinophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Pelagophyceae/Eukaryota \/ Pelagophyceae \/ Pelagophyceae/g' tmp #Pelagophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Raphidophyceae/Eukaryota \/ Raphidophyceae \/ Raphidophyceae/g' tmp #Raphidophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Synurophyceae/Eukaryota \/ Synurophyceae \/ Synurophyceae/g' tmp #Synurophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Bolidophyceae/Eukaryota \/ Bolidophyceae \/ Bolidophyceae/g' tmp #Bolidophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Polycystinea/Eukaryota \/ Polycystinea \/ Polycystinea/g' tmp #Polycystinea
sed -i 's/Eukaryota \/ unknown phylum \/ Choanoflagellata/Eukaryota \/ Choanoflagellata \/ Choanoflagellata/g' tmp #Choanoflagellata
sed -i 's/Eukaryota \/ unknown phylum \/ Filasterea/Eukaryota \/ Filasterea \/ Filasterea/g' tmp #Filasterea

#label for those missing phylum + class labels
sed -i 's/Eukaryota \/ unknown phylum \/ unknown class \/ Telonemida/Eukaryota \/ Telonemida \/ Telonemida \/ Telonemida/g' tmp #Telonemida
sed -i 's/Eukaryota \/ unknown phylum \/ unknown class \/ Jakobida/Eukaryota \/ Jakobida \/ Jakobida \/ Jakobida/g' tmp #Jakobida
sed -i 's/Eukaryota \/ unknown phylum \/ unknown class \/ Pirsoniales/Eukaryota \/ Pirsoniales \/ Pirsoniales \/ Pirsoniales/g' tmp #Pirsoniales

#execute final step for the LCA program (forming consensus taxonomies for each ASV)
mv -f tmp blast_96_sim/CO1_ASV_sequences.combined_all.blast.out #put tmp file where it belongs, add label (this is an overwrite, careful!!)
python2 ~/programs/galaxy-tool-lca/lca.species.py -i blast_96_sim/CO1_ASV_sequences.combined_all.blast.out -o taxonomy_table.CO1.NCBI_NT+customDB.iterative_blast.txt -b 100 -id 90 -cov 50 -t best_hit -tid 98 -tcov 50 -fh environmental,unidentified,kingdom -flh unclassified

####same blast hit resutls processed with LCA only####
python2 ~/programs/galaxy-tool-lca/lca.species.py -i blast_96_sim/CO1_ASV_sequences.combined_all.blast.out -o taxonomy_table.CO1.NCBI_NT+customDB.iterative_blast.LCA_ONLY.txt -b 100 -id 90 -cov 50 -t only_lca -fh environmental,unidentified,kingdom -flh unclassified

#cleanup
rm taxonomy_tmp blast_96_sim/tmp #remove redundant files
