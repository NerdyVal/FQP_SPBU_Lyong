library(ggplot2)
library(dada2)
library(phyloseq)
library(Biostrings)

# list amplicons into forward/reverse reads
amplicons_F <- sort(list.files(path = "/path_to_reads", 
                               pattern = "R1", # forward reads
                               full.names = TRUE))

amplicons_R <- sort(list.files(path = "/path_to_reads",
                               pattern = "R2", # reverse reads
                               full.names = TRUE))


sample.names <- sapply(strsplit(basename(amplicons_F), "_"), `[`, 1)
head (sample.names)
plotQualityProfile(amplicons_F[1:3])

# trimming & filtering
root.qc = './qc'
forward.qc <- file.path(root.qc, paste0(sample.names, "_R1.fastq.gz"))
reverse.qc <- file.path(root.qc, paste0(sample.names, "_R2.fastq.gz"))

qc.out <- filterAndTrim(amplicons_F, forward.qc, amplicons_R, reverse.qc, 
                        truncLen=c(225, 225),maxN=0, maxEE=c(2,2), truncQ = 2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
head(qc.out)

# error model for forward reads
errF <- learnErrors(forward.qc, multithread=TRUE)

# error model for reverse reads
errR <- learnErrors(reverse.qc, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(forward.qc, err=errF, multithread=TRUE)
dadaRs <- dada(reverse.qc, err=errR, multithread=TRUE)
dadaFs[1]

# merge Forward & Reverse Reads
merged <- mergePairs(dadaF = dadaFs, # dada result
                     derepF = forward.qc, # path of filtered reads
                     dadaR = dadaRs,
                     derepR = reverse.qc,
                     minOverlap = 8,
                     maxMismatch = 2)
head(merged[[1]])

# construct Sequence Table
seqtab <- makeSequenceTable(merged)
dim(seqtab)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# track the number of reads or sequences that made it through each step of the processes
getN <- function(x) sum(getUniques(x))
track <- cbind(qc.out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merged, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

saveRDS(seqtab.nochim, "seqtab")

# assign Taxonomy
memory.size()
taxa <- assignTaxonomy(seqtab.nochim, "/path_to_assignmentr", multithread=TRUE)
taxa.print = taxa
rownames(taxa.print) <- NULL
head(taxa.print)

bad_taxa = as.character(rownames(taxa))[(taxa[, 1] == 'Eukaryota' & 
                                           is.na(taxa[, 3])) | is.na(taxa[, 2])]
print(length(bad_taxa))
print(ncol(seqtab.nochim))
seqtab.nochim = seqtab.nochim[, sapply(as.character(colnames(seqtab.nochim)), 
                                       function(x) !(x %in% bad_taxa))]
print(nrow(seqtab.nochim))
taxa.print <- taxa # removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa_sp <- addSpecies(taxa, "/path_to_species_assignment", tryRC=TRUE)

# phyloseq object
samples.out <- rownames(seqtab.nochim)
metadata <- data.frame(SampleID=sample.names)
rownames(metadata) <- samples.out
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))

# get sequences
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps) # connect dna data to taxa names
ps <- merge_phyloseq(ps, dna) # merge data
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps))) # change taxa names into shorter id (ASVn)
ps

#save the phyloseq object:
saveRDS(ps, file = "name.RDS")
dir.create("MicrobiomeAnalyst_data")

# For the representative sequences:
df_refseq <- as.data.frame(refseq(ps))
df_refseq <- tibble::rownames_to_column(df_refseq, var = "#NAME")
write.table(df_refseq, "MicrobiomeAnalyst_data/rep_seqs.csv", row.names=FALSE, col.names=FALSE)

# For the ASV table:
df_asv_table <- as.data.frame(t(otu_table(ps)))
df_asv_table <- tibble::rownames_to_column(df_asv_table, var = "#NAME")
write.csv(df_asv_table, "MicrobiomeAnalyst_data/asv_table.csv", row.names=FALSE)

# For the taxonomy table:
df_tax_table <- as.data.frame(tax_table(ps))
df_tax_table <- tibble::rownames_to_column(df_tax_table, var = "#TAXONOMY")
write.csv(df_tax_table, "MicrobiomeAnalyst_data/tax_table.csv", row.names=FALSE)

# For the metadata:
resulting_table <- tibble::rownames_to_column(metadata, var = "#NAME")
write.csv(metadata, "MicrobiomeAnalyst_data/metadata.csv", row.names=FALSE)
