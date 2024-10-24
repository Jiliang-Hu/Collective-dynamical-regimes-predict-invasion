#library("knitr")
#library("gridExtra")
library("devtools")
library("dada2")
packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

miseq_path <- "./16S data"
list.files(miseq_path)

# Sort ensures forward\reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq.gz"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq.gz"))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
#just show that our vectors point to the actual files (just showing the first 3)
fnFs[1:3]
fnRs[1:3]

#show quality of the first two sequences
#look here for quality score: https://www.illumina.com/science/education/sequencing-quality-scores.html
plotQualityProfile(fnFs[1:5])
plotQualityProfile(fnRs[1:5])


# FILTERING
#We are going to filter the fastq, cutting at a lower quality threshold 
#we determine: 230 bases for forward reads, 180 for reverse

#define new (filtered) filenames
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

#we also apply standard filtering techniques
# that aim at maximum 2 expected errors per read
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(145,140), trimLeft = 0,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

head(out)


# identifies identical sequences, groups them into 'unique sequences', and assigns 'abundance'

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames
names(derepFs)
names(derepRs)

#distinguish sequencing errors from real biological variation
#uses unsupervised learning to find the model parameters
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#inspect the forward unique sequences (from sample 1) in the dada file after errors have been processed
dadaFs[[1]]
# and the reverse sequences from sample 1
dadaRs[[1]]


# CONSTRUCT SEQUENCE TABLE AND REMOVE CHIMERAS
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


seqtabAllM <- makeSequenceTable(mergers,orderBy = "abundance")
# maybe try this to dereplicate?
#mergeSequenceTables(mergers, orderBy= "abundance")

dim(seqtabAllM)
table(nchar(getSequences(seqtabAllM)))


# remove chimeras from 'unique sequences'
seqtabMNoC <- removeBimeraDenovo(seqtabAllM)
#seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtabMNoC)

#fraction of actual (non-chimeric) species
sum(seqtabMNoC)/sum(seqtabAllM)

# Taxa based on the GreenGenes data, species assigned from the Ribosome Database Project
#fastaRefTaxa <- "./gg_13_8_train_set_97.fa.gz"
GreenGenesTaxM <- assignTaxonomy(seqtabMNoC,"silva_nr_v138_train_set.fa")


#Adding Species from the RdP on the top of the GreenGenes classification did not 
#improved the resolution. Maybe it's a way of not giving false Species-level assignments
GreenGenes.print <- GreenGenesTaxM # Removing sequence rownames for display only
rownames(GreenGenes.print) <- NULL
head(GreenGenes.print)
write.csv(GreenGenes.print, "M_GreenGenesTaxa.csv")
#RdpSpeciesTab.print <- RdpSpeciesTab # Removing sequence rownames for display only
#rownames(RdpSpeciesTab.print) <- NULL
#head(RdpSpeciesTab.print)
#write.csv(RdpSpeciesTab.print, "RdpSpeciesTaxa.csv")

#get OTU table
ps <- phyloseq(otu_table(seqtabMNoC, taxa_are_rows=FALSE),
               sampleNames,
               tax_table(GreenGenesTaxM))

# assign new names to species for a more readable output
new.names <- paste0("ASV", seq(ntaxa(ps))) # Define new names ASV1, ASV2, ...
seqs <- taxa_names(ps) # Store sequences
names(seqs) <- new.names # Make map from ASV1 to full sequence
taxa_names(ps) <- new.names # Rename to human-friendly format

write.csv(seqs, "M_UniqueSequences.csv")
write.csv(otu_table(ps), "M_OTUtableGreenGenes.csv")
write.csv(tax_table(ps), "M_TAXAtableGreenGenes.csv")

