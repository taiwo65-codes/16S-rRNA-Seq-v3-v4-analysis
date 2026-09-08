#

## 4/9/26
#  CFI 16SrRNA data analysis (v3-v4) * Taiwo Bankole*


# install DADA2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.22")


# ----- load package -----
library(dada2); packageVersion("dada2")



#----- set path ____________
path <- "/Users/banko/Desktop/Li_lab/trimmed"
list.files(path) # view files in path 



# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))



# Extract the sample names from the file names for later use
# This strips away "trim_" and everything after "_R1" to give you clean sample names 
# (e.g., "Li_pos_control_S414_L001")
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)
sample.names <- gsub("trim_", "", sample.names)

# Print the sample names to make sure they look correct
sample.names


# plot the quality of the first two forward reads
plotQualityProfile(fnFs[c(11:20)])

# plot the quality of the first two reverse reads (might be worse)
plotQualityProfile(fnRs[c(6:10)])

# plotQualityProfile(fnFs[c(1,41)])
# plotQualityProfile(fnRs[c(1,41)])



# Filter and Trim

# Create file paths for the new filtered files
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Assign your sample names to these new files
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Run the actual filtering command
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(240, 230),
                     maxN=0, 
                     maxEE=c(2, 2), 
                     truncQ=2, 
                     rm.phix=TRUE,
                     compress=TRUE, 
                     multithread=FALSE)

# Print the summary to see how many reads survived
head(out)




# Error rates

# 1. Learn the error rates for the Forward reads
errF <- learnErrors(filtFs, multithread=FALSE)

# 2. Learn the error rates for the Reverse reads
errR <- learnErrors(filtRs, multithread=FALSE)

# 3. Plot the Forward error rates to visualize the machine learning model
plotErrors(errF, nominalQ=TRUE)

# 3. Plot the Forward error rates to visualize the machine learning model
plotErrors(errR, nominalQ=TRUE)




# sample inference
# apply sample inference to forward reads
dadaFs <- dada(filtFs, err=errF, multithread=F)

# apply sample inference to forward reads
dadaRs <- dada(filtRs, err=errR, multithread=F)

# inspecting the returned dada-class object
dadaFs[[1]]


# merge forward and revrese reads together
# This uses those identical "name tags" we set up earlier to zipper them up.
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])


# Construct your first Amplicon Sequence Variant (ASV) table
seqtab <- makeSequenceTable(mergers)

# 5. Check the dimensions of your new table
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


## ************************************* 


# Keep only sequences that fall within the normal biological V3-V4 length range
seqtab.clean <- seqtab[, nchar(colnames(seqtab)) %in% 399:430]

# Check the new distribution to confirm it worked
table(nchar(getSequences(seqtab.clean)))

# Check the dimensions to see how many ASVs survived
dim(seqtab.clean)


### *******************************************

# 1. Remove the chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab.clean, method="consensus", multithread=FALSE, verbose=TRUE)

# 2. Check how many ASVs are left
dim(seqtab.nochim)

# 3. Check what percentage of your total reads survived the whole pipeline
sum(seqtab.nochim)/sum(seqtab.clean)


### **********************************

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

## ********************************************

# Taxonomy assignments

# 1. Assign Kingdom down to Genus
taxa <- assignTaxonomy(seqtab.nochim, file.path(path, "silva_nr99_v138.2_toGenus_trainset.fa.gz"), multithread=FALSE)

# 2. Add the exact Species names where possible
taxa <- addSpecies(taxa, file.path(path, "silva_v138.2_assignSpecies.fa.gz"))

# 3. Print the top 5 most abundant bacteria to see your results!
taxa.print <- taxa # Removing sequence names for display only
rownames(taxa.print) <- NULL
head(taxa.print)



# save output.
saveRDS(seqtab.nochim, file.path(path, "seqtab_dada2_raw.rds"))
saveRDS(taxa, file.path(path, "taxa_dada2_raw.rds"))






#####################################################

# *****************************
# 1. Load the required libraries
library(phyloseq)
library(ggplot2) # We use this to make the plot look professional

# 2. Load your newly minted metadata spreadsheet
meta <- read.csv(file.path(path, "CFI metadata.csv"), row.names=1)

# 3. Glue everything together into a phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa))

# 4. Filter out Mitochondria, Chloroplasts, and unclassified artifacts
ps.clean <- subset_taxa(ps, Family != "Mitochondria" | is.na(Family))
ps.clean <- subset_taxa(ps.clean, Order != "Chloroplast" | is.na(Order))
ps.clean <- subset_taxa(ps.clean, !is.na(Phylum))

# 5. Remove the sequencing controls from the plot (we don't want them in the biological graph)
ps.mice <- subset_samples(ps.clean, Group != "Negative" & Group != "Positive")



# **********************************
# ==========================================
# PATH A: ALPHA DIVERSITY (No Filtering)
# ==========================================
# You will use ps.clean directly for this.
# Example: plot_richness(ps.clean, x="Treatment", measures=c("Shannon", "Chao1"))


# Generate the Alpha Diversity Boxplot!
plot_richness(ps.mice, x="Group", measures=c("Observed", "Chao1")) + 
  geom_boxplot(aes(fill=Group), alpha=0.7) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold")) +
  labs(title="Alpha Diversity: Gut Microbiome Richness", x="Treatment Group")




# ************** alpha diversity estimation **************
# 1. Extract all alpha diversity metrics
alpha_div <- estimate_richness(ps.mice)

# 2. Save to CSV (using row.names=TRUE to keep your Sample IDs!)
write.csv(alpha_div, file.path(path, "alpha_diversities_complete.csv"), row.names=TRUE)

# 3. Extract the metadata from your phyloseq object as a standard dataframe
meta_mice <- as(sample_data(ps.mice), "data.frame")

# 4. Merge the alpha diversity table with your metadata so R knows which mouse is which
# We use merge() by row.names to ensure the DNA data perfectly matches the treatment group
alpha_meta <- merge(alpha_div, meta_mice, by="row.names")



# " R reread my hyphen as periods for rownames, change to underscore nexttime


# recalculate
# 1. Add the math columns directly into your metadata (bypassing the merge completely)
meta_mice$Shannon <- alpha_div$Shannon
meta_mice$Simpson <- alpha_div$Simpson
meta_mice$Chao1 <- alpha_div$Chao1
meta_mice$Observed <- alpha_div$Observed

# 2. Re-run the Global Statistical Test to ensure it sees all 4 groups
cat("\n=== Kruskal-Wallis Test for Simpson Index ===\n")
kruskal.test(Simpson ~ Group, data = meta_mice)

# 3. Re-run the Pairwise Post-Hoc Test
cat("\n=== Pairwise Wilcoxon Test for Simpson Index (FDR Corrected) ===\n")
pairwise.wilcox.test(meta_mice$Simpson, meta_mice$Group, p.adjust.method = "BH")



###############################################################




# ==========================================
# PATH B: BETA DIVERSITY (Filtering Noise)
# ==========================================

# Step 1: Filter by Total Abundance
# Remove any ASV that does not have at least 10 total reads across the ENTIRE experiment.
# (This clears out the extreme sequencing artifacts that survived DADA2).
ps.beta <- prune_taxa(taxa_sums(ps.mice) >= 10, ps.mice)

# Step 2: Filter by Prevalence
# An ASV might have 500 reads, but if all 500 reads are in just ONE single mouse, 
# it's likely a fluke or an individual infection, not a true group-level feature.
# This command forces an ASV to be present in at least 2 different samples to survive.
ps.beta <- filter_taxa(ps.beta, function(x) sum(x > 0) >= 2, TRUE)

# Step 3: Compare your dimensions
cat("Original ASVs for Alpha Diversity:", ntaxa(ps.mice), "\n")
cat("Filtered ASVs for Beta Diversity:", ntaxa(ps.beta), "\n")

###########################################################################

# Beta diversity

# 1. Calculate the Bray-Curtis distance matrix and run the PCoA algorithm
ord.bray <- ordinate(ps.beta, method="PCoA", distance="bray")

# 2. Generate a professional PCoA plot with clustering ellipses
plot_ordination(ps.beta, ord.bray, color="Group") +
  geom_point(size=4, alpha=0.8) +
  stat_ellipse(aes(group=Group), type="t", linetype=2, alpha=0.5) + 
  theme_bw() +
  theme(legend.position="right", text=element_text(size=14, face="bold")) +
  labs(title="Beta Diversity: PCoA (Bray-Curtis)", 
       subtitle="Do the microbiomes cluster by treatment?",
       color="Treatment Group")



# Permanova calculations

# 1. Load the vegan package for community ecology statistics
library(vegan)

# 2. Extract the metadata from your stringently filtered beta object
meta_beta <- as(sample_data(ps.beta), "data.frame")

# 3. Calculate the exact Bray-Curtis distance matrix
bray_dist <- phyloseq::distance(ps.beta, method="bray")

# 4. Run the global PERMANOVA test
# set.seed ensures your p-values are identical every time you re-run the code
set.seed(123) 
adonis2(bray_dist ~ Group, data = meta_beta, permutations = 999)



# pairwise permanova
# ==========================================
# TEST 1: Did CFI alter the healthy controls? pvalue= 0.032
# ==========================================
# 1. Keep only the Control and Control-CFI mice
ps.control_pair <- subset_samples(ps.beta, Group %in% c("Control", "Control-CFI"))

# 2. Re-calculate the distance and metadata for just these mice
bray_control <- phyloseq::distance(ps.control_pair, method="bray")
meta_control <- as(sample_data(ps.control_pair), "data.frame")

# 3. Run the PERMANOVA
cat("\n=== PERMANOVA: Control vs Control-CFI ===\n")
set.seed(123)
adonis2(bray_control ~ Group, data = meta_control, permutations = 999)


# ==========================================
# TEST 2: Did CFI rescue the HFD mice? p value = 0.04
# ==========================================
# 1. Keep only the HFD and HFD-CFI mice
ps.hfd_pair <- subset_samples(ps.beta, Group %in% c("HFD", "HFD-CFI"))

# 2. Re-calculate the distance and metadata for just these mice
bray_hfd <- phyloseq::distance(ps.hfd_pair, method="bray")
meta_hfd <- as(sample_data(ps.hfd_pair), "data.frame")

# 3. Run the PERMANOVA
cat("\n=== PERMANOVA: HFD vs HFD-CFI ===\n")
set.seed(123)
adonis2(bray_hfd ~ Group, data = meta_hfd, permutations = 999)





#######################################################################

#Taxonomy check for 'Phylum'

# 1. Transform counts to Relative Abundance (Percentages)
ps.rel <- transform_sample_counts(ps.mice, function(x) x / sum(x))

# 2. Agglomerate the data to the Phylum level
ps.phylum <- tax_glom(ps.rel, taxrank = "Phylum")

# 3. Build the professional Bar Plot
plot_bar(ps.phylum, x="Name", fill="Phylum") + 
  # facet_grid automatically groups your mice by their treatment group!
  facet_grid(~Group, scales="free_x", space="free_x") +
  # This line removes the ugly black borders from the boxes to make it look clean
  geom_bar(aes(fill=Phylum), color=NA, stat="identity", position="stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face="bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "right") +
  labs(title="Gut Microbiome Composition (Phylum Level)",
       y="Relative Abundance (1.0 = 100%)", x="Mouse ID")



# for genus
# 2. Agglomerate the data to the Phylum level
ps.genus <- tax_glom(ps.rel, taxrank = "Genus")

# 3. Build the professional Bar Plot
plot_bar(ps.genus, x="Name", fill="Genus") + 
  # facet_grid automatically groups your mice by their treatment group!
  facet_grid(~Group, scales="free_x", space="free_x") +
  # This line removes the ugly black borders from the boxes to make it look clean
  geom_bar(aes(fill=Genus), color=NA, stat="identity", position="stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face="bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "right") +
  labs(title="Gut Microbiome Composition (Genus Level)",
       y="Relative Abundance", x="Mouse ID")



# ######## extract relative abundance table ##############


################################################# ignore ############

# 1. Transform to Relative Abundance (Percentages)
ps.prop <- transform_sample_counts(ps.mice, function(otu) otu / sum(otu) * 100)

# 2. Aggregate to Genus level
ps.genus <- tax_glom(ps.prop, taxrank = "Genus")

# 3. Melt everything together (This saves your metadata!)
melted_data <- psmelt(ps.genus)

# 4. Use tidyr to make it "wide" just like your script
library(tidyr)
library(dplyr)

wide_genus_df <- melted_data %>%
  select(Sample, Group, Genus, Abundance) %>%
  pivot_wider(names_from = Sample, values_from = Abundance)

# 5. Save as a text file
write.table(wide_genus_df,
            file = "genus_rel_abundance.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


##########################################################################################


# 1. Transform to Proportions (Decimals only, NO multiplying by 100)
ps.prop_1 <- transform_sample_counts(ps.mice, function(otu) otu / sum(otu))

# 2. Aggregate to the Phylum level
ps.phylum_1 <- tax_glom(ps.prop, taxrank = "Phylum")

# 3. Melt everything together to save metadata
melted_phylum <- psmelt(ps.phylum)

# 4. Use tidyr to make it "wide" (just like your screenshot)
library(tidyr)
library(dplyr)

wide_phylum_df <- melted_phylum %>%
  select(Sample, Phylum, Abundance) %>%
  pivot_wider(
    names_from = Sample, 
    values_from = Abundance,
    # Fills in empty spaces with clean 0s
    values_fill = list(Abundance = 0) 
  )

# 5. Save as a text file
write.csv(wide_phylum_df,
            file = "phylum_rel_abundance_wide.csv",
            row.names = T, col.names = T)







# Now doing the above phylum codes for other taxa.

# 1. Transform to Proportions (Decimals)
ps.prop_2 <- transform_sample_counts(ps.mice, function(otu) otu / sum(otu))

# 3. Create the custom function with the updated 'Name' variable
export_wide_taxonomy <- function(physeq_prop, tax_level, filename) {
  
  # Aggregate and melt
  ps_glom <- tax_glom(physeq_prop, taxrank = tax_level)
  melted_data <- psmelt(ps_glom)
  
  # Pivot using your metadata 'Name' column instead of the raw Sample ID
  wide_df <- melted_data %>%
    # We select Name here so it is available for the pivot
    select(Name, all_of(tax_level), Abundance) %>%
    pivot_wider(
      names_from = Name, # This turns your clean labels into the column headers
      values_from = Abundance,
      values_fill = list(Abundance = 0)
    ) %>%
    # This forcefully renames the first column (e.g., "Phylum") to "Name"
    rename(Name = all_of(tax_level)) 
  
  # Save the file
  write.table(wide_df,
              file = file.path(path, filename),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  cat(paste("Successfully saved:", filename, "\n"))
}

# 4. Execute the function for every taxonomic level
export_wide_taxonomy(ps.prop_2, "Class",   "class_rel_abundance_wide.txt")
export_wide_taxonomy(ps.prop_2, "Order",   "order_rel_abundance_wide.txt")
export_wide_taxonomy(ps.prop_2, "Family",  "family_rel_abundance_wide.txt")
export_wide_taxonomy(ps.prop_2, "Genus",   "genus_rel_abundance_wide.txt")
export_wide_taxonomy(ps.prop_2, "Species", "species_rel_abundance_wide.txt")
export_wide_taxonomy(ps.prop_2, "Phylum",   "phylum_rel_abundance_wide.txt")



# 1. Transform to Proportions (Decimals) ONLY ONCE to save processing time
ps.prop <- transform_sample_counts(ps.mice, function(otu) otu / sum(otu))

# 2. Load required packages
library(tidyr)
library(dplyr)

# 3. Create a custom function to handle the melting, pivoting, and saving
export_wide_taxonomy <- function(physeq_prop, tax_level, filename) {
  
  # Aggregate to the specific taxonomic level
  ps_glom <- tax_glom(physeq_prop, taxrank = tax_level)
  
  # Melt to extract data and metadata
  melted_data <- psmelt(ps_glom)
  
  # Pivot to the exact wide format from your screenshot
  wide_df <- melted_data %>%
    # all_of() tells R to dynamically select the column matching your tax_level
    select(Sample, all_of(tax_level), Abundance) %>%
    pivot_wider(
      names_from = Sample, 
      values_from = Abundance,
      values_fill = list(Abundance = 0) # Fills empty spaces with clean 0s
    )
  
  # Save as a clean text file
  write.table(wide_df,
              file = file.path(path, filename),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  # Print a confirmation message to the console
  cat(paste("Successfully saved:", filename, "\n"))
}

# 4. Execute the function for every remaining taxonomic level!
export_wide_taxonomy(ps.prop, "Class",   "class_rel_abundance_wide.txt")
export_wide_taxonomy(ps.prop, "Order",   "order_rel_abundance_wide.txt")
export_wide_taxonomy(ps.prop, "Family",  "family_rel_abundance_wide.txt")
export_wide_taxonomy(ps.prop, "Genus",   "genus_rel_abundance_wide.txt")
export_wide_taxonomy(ps.prop, "Species", "species_rel_abundance_wide.txt")
















# ********* Genus rel abundance ******

# ps.prop <- transform_sample_counts(ps_filtered, function(otu) otu / sum(otu)) # if desired in %, just multiply by 100

# Aggregrate to desired rank/ taxon
# ps.genus <- tax_glom(ps.prop, taxrank = "Genus")

# Extract OTU table as data frame
otu_genus_df <- as.data.frame(otu_table(ps.genus))
if(!taxa_are_rows(ps.genus)){
  otu_genus_df <- t(otu_genus_df)
}

# Extract taxonomy table as data frame
tax_genus_df <- as.data.frame(tax_table(ps.genus))

# Merge taxonomy and abundance
merged_genus_df <- cbind(otu_genus_df, tax_genus_df)

# 5. Summarize to class level (no ASV IDs)
genus_abundance <- merged_genus_df %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum))

# Save as TXT
write.table(genus_abundance,
            file = " genus rel_abundance.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)





#########################################################
# Deseq2

# 1. Load the package (If it asks to update other packages, you can usually type 'n' for no)
library(DESeq2)

# 2. Convert your phyloseq object into a DESeq2 object
# We use ps.mice because DESeq2 requires RAW counts, not percentages!
dds <- phyloseq_to_deseq2(ps.mice, ~ Group)

# 3. Handle the "sparse" microbiome zeroes (The modern fix for 16S data)
dds <- estimateSizeFactors(dds, type="poscounts")

# 4. Run the differential abundance algorithm
dds <- DESeq(dds, fitType="local")

# 5. Extract the specific comparison: HFD-CFI vs HFD
# The order is: c("ColumnName", "Numerator (Treatment)", "Denominator (Baseline)")
res <- results(dds, contrast=c("Group", "HFD-CFI", "HFD"))

# 6. Filter for ONLY the statistically significant bacteria (adjusted p-value < 0.05)
sigtab <- res[which(res$padj < 0.05), ]

# 7. Attach the SILVA bacterial taxonomy names to the math results
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps.mice)[rownames(sigtab), ], "matrix"))

# 8. Sort the results so the most highly enriched bacteria are at the top
sigtab <- sigtab[order(sigtab$log2FoldChange, decreasing=TRUE), ]

# 9. Save the full results to your folder for your manuscript!
write.csv(sigtab, file.path(path, "DESeq2_HFD_vs_HFD-CFI.csv"), row.names=TRUE)

# 10. Print the top results to the screen
head(sigtab[, c("log2FoldChange", "padj", "Phylum", "Family", "Genus", "Species")], 10)


########################################################################

# ==========================================================
# Comparison 1: Control vs. Control-CFI
# Question: What exactly did Caffeine do to the healthy mice?
# ==========================================================
# Note: "Control-CFI" is the numerator (positive fold change), "Control" is the baseline
res_ctrl <- results(dds, contrast=c("Group", "Control-CFI", "Control"))

# Filter for significance and add taxonomy names
sigtab_ctrl <- res_ctrl[which(res_ctrl$padj < 0.05), ]
sigtab_ctrl <- cbind(as(sigtab_ctrl, "data.frame"), as(tax_table(ps.mice)[rownames(sigtab_ctrl), ], "matrix"))
sigtab_ctrl <- sigtab_ctrl[order(sigtab_ctrl$log2FoldChange, decreasing=TRUE), ]

# Save to your folder and print the top 10 results
write.csv(sigtab_ctrl, file.path(path, "DESeq2_Control_vs_Control-CFI.csv"), row.names=TRUE)
cat("\n=== Top 10 Bacteria Altered by CFI in Healthy Controls ===\n")
head(sigtab_ctrl[, c("log2FoldChange", "padj", "Phylum", "Family", "Genus")], 10)


# ==========================================================
# Comparison 2: Control-CFI vs. HFD
# Question: How different is a healthy-treated gut from a severe obese gut?
# ==========================================================
# Note: "HFD" is the numerator, "Control-CFI" is the baseline
res_hfd_c <- results(dds, contrast=c("Group", "HFD", "Control-CFI"))

# Filter for significance and add taxonomy names
sigtab_hfd_c <- res_hfd_c[which(res_hfd_c$padj < 0.05), ]
sigtab_hfd_c <- cbind(as(sigtab_hfd_c, "data.frame"), as(tax_table(ps.mice)[rownames(sigtab_hfd_c), ], "matrix"))
sigtab_hfd_c <- sigtab_hfd_c[order(sigtab_hfd_c$log2FoldChange, decreasing=TRUE), ]

# Save to your folder and print the top 10 results
write.csv(sigtab_hfd_c, file.path(path, "DESeq2_Control-CFI_vs_HFD.csv"), row.names=TRUE)
cat("\n=== Top 10 Bacteria Different Between HFD and Control-CFI ===\n")
head(sigtab_hfd_c[, c("log2FoldChange", "padj", "Phylum", "Family", "Genus")], 10)








































































































































