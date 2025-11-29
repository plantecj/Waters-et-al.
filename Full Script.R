##Sequenced Read Processing

#Load "dada2"

#Load data into R environment
#Sample names
samples <- scan("samples", what = "character")
#Forward and reverse reads
forward_reads<- paste0(samples, "_R1_001_trimmed.fastq.gz")
reverse_reads<- paste0(samples, "_R2_001_trimmed.fastq.gz")
#Filtered forward and reverse reads
filtered_forward_reads<-paste0(samples, "_R1_filtered.fastq.gz")
filtered_reverse_reads<-paste0(samples, "_R2_filtered.fastq.gz")

#Trim reads
#Check quality plots
plotQualityProfile(forward_reads)
plotQualityProfile(reverse_reads)
#Trim
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxEE=c(2,2), rm.phix=TRUE, minLen=175, truncLen=c(250,200))
#Check what was filtered
class(filtered_out)
dim (filtered_out)
filtered_out
#Check quality plots after trimming
plotQualityProfile(filtered_forward_reads)
plotQualityProfile(filtered_reverse_reads)

#Create error model with data
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)

#Plot error model
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)

#Dereplicate reads
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples 
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples

##Infer ASVs
dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread = TRUE)
#Merge forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse, derep_reverse, trimOverhang=TRUE)
#Generate OTU matrix (Count Table)
seqtab <- makeSequenceTable(merged_amplicons)
#Save as RDS
saveRDS(seqtab, file = "seqtab.rds")

#Several samples were sequenced at a separate time.
#Pre-prossessing was the same up to this point.
#Merge the previous count table with the current table.
seqtab_combo <- mergeSequenceTables(seqtab_1, seqtab_2, repeats = "sum")
#Save as RDS
saveRDS(seqtab_combo, file= "seqtab_combo.rds")

#Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab_combo, verbose=T)

#Assign Taxonomy using SILVA v 132
download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r132_March2018.RData", destfile="SILVA_SSU_r132_2018.RData")
#load "DECIPHER"
library(DECIPHER)
#Create DNAStringSet using ASVs
dna_silva <- DNAStringSet(getSequences(seqtab.nochim))
#Classify ASVs
tax_info_silva <- IdTaxa(test=dna_silva, trainingSet=trainingSet, strand="both", processors=NULL)

#Extract taxa tables, ASV tables, and asv_fasta data
ranks_silva <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax_silva <- t(sapply(tax_info_silva, function(x) {
  m_silva <- match(ranks_silva, x$rank)
  taxa_silva <- x$taxon[m_silva]
  taxa_silva[startsWith(taxa_silva, "unclassified_")] <- NA
  taxa_silva
}))
colnames(asv_tax_silva) <- ranks_silva
rownames(asv_tax_silva) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax_silva, "ASVs_taxonomy_silva.tsv", sep = "\t", quote=F, col.names=NA)

#Remove contaminants 
#Install and load "decontam"
library(decontam)
vector_for_decontam <- c(rep(TRUE, 1), rep(FALSE, 26), rep(TRUE, 2), rep(FALSE, 54))
contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)
table(contam_df$contaminant)

#Initial analysis
library(tidyverse)  
library(phyloseq)
library(vegan)
library(DESeq2)
library(dendextend)
library(viridis)
library(ggplot2)
#Create Count Table
count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1, check.names=F, sep="\t")[ , -c(1, 28, 29)]
#Create Taxonomy Table
tax_tab_silva <- as.matrix(read.table("ASVs_taxonomy_silva.tsv", header=T, row.names=1, check.names=F, sep="\t"))
#Add in sample metadata



##Make phylogenetic tree using "phangorn"
library(ape)
library(phangorn)
library(phyloseq)
library(dada2)
library(DECIPHER)

#Extract sequences from dada2 output
class(asv_seqs)
sequences<- getSequences(asv_seqs)
names(sequences)= sequences


#Sequence alignment with DECIPHER
alignment <- AlignSeqs(DNAStringSet(sequences), anchor = NA)

#Change the sequence alignment output to a phyDat structure (data stored as integers in
# list with as many vectors as sequences)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

#Create distance matrix
dm <- dist.ml(phang.align)

#Create neigbor joining tree
treeNJ <- NJ(dm) 

#Set up internal maximum likelihood
fit = pml(treeNJ, data=phang.align) #negative edges length changed to 0!
#pml computes the likelihood of a phylogenetic tree given a sequence alignment and a model

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,rearrangement = "stochastic", control = pml.control(trace = 0))
# optim.pml optimizes the different model parameters

saveRDS(fitGTR, file="fitGTR.rds")


#Variance Stabilizing Transformation
# Create DESeq2 object 
deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sampledatv1, design = ~Season)
#Season should represent the season that the samples where taken from
#Use a Variance Stabilizing Transformation (McMurdie and Hlomes, 2014: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531) to normalize sampling depth
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
#Extract transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)

#Create Phyloseq Objects
#Normalized to Sampling Depth
#Make a count table without ASV labels (ASV_1, ASV_2, etc.) in order to add phylogenetic tree to the phyloseq object
asv_tab_phylo<-t(seqtab.nochim) 

write.table(asv_tab_phylo, "asv_tab_phylo.tsv",
            sep="\t", quote=F, col.names=NA)
count_tab_phylo <- read.table("asv_tab_phylo.tsv", header=T, row.names=1,
                              check.names=F, sep="\t")[ , -c(1, 28, 29)]
deseq_counts_phylo <- DESeqDataSetFromMatrix(count_tab_phylo, colData = sampledatv1, design = ~Season.base)
deseq_counts_vst_phylo <- varianceStabilizingTransformation(deseq_counts_phylo)
vst_trans_count_tab_phylo <- assay(deseq_counts_vst_phylo)

#Use transformed table to create a Phyloseq Object
#OTU (count) table
vst_count_phy <- otu_table(vst_trans_count_tab_phylo, taxa_are_rows=T)
#Then the sample table
sample_info_tab_phy <- sample_data(sampledata)
#Then the phylogenetic tree (fitGTR)
tree<-phy_tree(fitGTR$tree)
#Combine 
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy, tree)



#Not normalized to sampling depth
#asv_tax_silva
ranks_silva <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax_silva_fasta <- t(sapply(tax_info, function(x) {
  m_silva <- match(ranks_silva, x$rank)
  taxa_silva <- x$taxon[m_silva]
  taxa_silva[startsWith(taxa_silva, "unclassified_")] <- NA
  taxa_silva
}))
colnames(asv_tax_silva_fasta) <- ranks_silva

#Changing the asv_tax_silva_fasta headings from numbers into the fasta sequences using "count_tab_phylo"
rownames(asv_tax_silva_fasta)<-count_tab_phylo$X

#Save files
saveRDS(asv_tax_silva_fasta, file="tax_tab_silva_fasta.rds")

#ASV Physeq
count_tab_phy <- otu_table(as.matrix(count_tab_phylo), taxa_are_rows=T)
tax_tab_phy_silva_fasta <- tax_table(as.matrix(tax_tab_silva_fasta))
ASV_physeq_silva <- phyloseq(count_tab_phy, tax_tab_phy_silva_fasta, sample_info_tab_phy, tree)

saveRDS (ASV_physeq_silva, file="silva_physeq")

#Filter out everything but diatoms
Silva_physeq_Diatom<-subset_taxa(ASV_physeq_silva, class=="Diatomea")
saveRDS (Silva_physeq_Diatom, file="silva_physeq_diatom")

#Prepare table for future analyses
silva_diatable<-otu_table(silva_physeq_diatom)
ASV_Table<-t(silva_diatable)








## Remove Singletons
#Load "dplyr" package
install.packages("dplyr")
library(dplyr)
#Load in ASV table with SampleID and ASVs as columns

#remove the singletons
ASV_Table_no_singletons <- ASV_Table
%>% select(where(~ sum(.x)>1))

#Save as csv
write.csv(ASV_Table_no_singletons, "Diatom_Table.csv", row.names = FALSE)

#Remove singletons from ASV_Table with additional samples collected from the water.
water_table_no_singletons <- ASV_Table_Water
%>% select(where(~ sum(.x)>1))

#Save as csv
write.csv(water_table_no_singletons, "Water_Diatom_Table.csv")





##Square-Root Transform Diatom Table and Water Diatom Table
Diatom_Table_sqrt <-sqrt(Diatom_Table)
Water_Diatom_Table_sqrt <-sqrt(Water_Diatom_Table)


##Extract Rare and Abundant ASVs for both transformed and non-transformed tables
#Tranformed
total_abund_sqrt <- colSums (Diatom.Table.sqrt)
#Relative abundance
rel_abund.sqrt <- total_abund.sqrt / sum(total_abund.sqrt)
#Rare and Abundant ASVs
rare_asvs.sqrt <- Diatom.Table.sqrt[, rel_abund.sqrt < 0.001]
abundant_asvs.sqrt <- Diatom.Table.sqrt[, rel_abund.sqrt >= 0.001]

#Not Transformed
total_abund <- colSums (Diatom_Table)
#Relative abundance
rel_abund <- total_abund / sum(total_abund)
#Rare and Abundant ASVs
rare_asvs <- Diatom_Table[, rel_abund < 0.001]
abundant_asvs <- Diatom_Table[, rel_abund >= 0.001]


##Rarefaction Curves
#Total (using "Diatom_Table")

library(vegan)

rarecurve(Diatom_Table, step=1000, col=Metadata$Color, lwd=2, ylab="ASV Richness", label=T)
abline(v=(min(rowSums(Diatom_Table))))
legend("bottomright", legend = Legend_dat$Season.base, col = Legend_dat$Color, lty=1)

#Rare rarefaction curve
rarecurve(rare_asvs, step=100, col=Metadata$Color, lwd=2, ylab="ASV Richness", label=T)
abline(v=(min(rowSums(rare_asvs))))
legend("bottomright", legend = Legend_dat$Season.base, col = Legend_dat$Color, lty=1)





##Distance Decay Model
##Distance Decay Model
#Load "vegan"
library(vegan)
#Use square-root transformed diatom table Diatom.Table.sqrt
#Install time matrix "time_dat"


#Create euclidean distance matrix for time
euc_time <- dist(time_dat$hrs, method = "euclidean")



#Total ASVs

#Create a bray-curtis matrix
bray_dist<-(vegdist(Diatom.Table.sqrt, method = "bray"))

#Convert to vectors 
time_vec <- as.vector(euc_time)
bc_vec <- as.vector(bray_dist)

dist_decay <- cor.test(time_vec, bc_vec, method = "pearson")

#data:  time_vec and bc_vec
#t = 12.161, df = 2554, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.1969724 0.2702758
#sample estimates:
#  cor 
#0.2339566 

plot(bc_vec ~ time_vec, ylab = "Community Dissimilarity", xlab = "Temporal Distance (Hours)")

model <- lm(bc_vec ~ time_vec)
abline(model)



#Rare ASVs
#Extract rare and abundant ASVs
total_abund <- colSums (Diatom.Table.sqrt)
#Relative abundance
rel_abund <- total_abund / sum(total_abund)
#Rare and Abundant ASVs
rare_asvs <- Diatom.Table.sqrt[, rel_abund < 0.001]
abundant_asvs <- Diatom.Table.sqrt[, rel_abund >= 0.001]

#Create a bray-curtis matrix
bray_dist_rare<-(vegdist(rare_asvs, method = "bray"))

#Convert to vectors 
bc_vec_rare <- as.vector(bray_dist_rare)

dist_decay_rare <- cor.test(time_vec, bc_vec_rare, method = "pearson")
#data:  time_vec and bc_vec_rare
#t = 16.35, df = 2554, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.2722923 0.3424969
#sample estimates:
#  cor 
#0.3078135

plot(bc_vec_rare ~ time_vec, ylab = "Community Dissimilarity", xlab = "Temporal Distance (Hours)")

model_rare <- lm(bc_vec_rare ~ time_vec)
abline(model_rare)


#Abundant ASVs
#Create a bray-curtis matrix
bray_dist_abundant<-(vegdist(abundant_asvs, method = "bray"))

#Convert to vectors 
bc_vec_abundant <- as.vector(bray_dist_abundant)

dist_decay_abundant <- cor.test(time_vec, bc_vec_abundant, method = "pearson")
#data:  time_vec and bc_vec_abundant
#t = 10.43, df = 2554, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.1646466 0.2390248
#sample estimates:
#  cor 
#0.2021272 

plot(bc_vec_abundant ~ time_vec, ylab = "Community Dissimilarity", xlab = "Temporal Distance (Hours)")

model_abundant <- lm(bc_vec_abundant ~ time_vec)
abline(model_abundant)




##NMDS 
library(vegan)
#Total ASVs
#Calculate NMDS using square-root transformed table
nmds_bray <- metaMDS(Diatom_Table_sqrt, distance="bray")


#Extract NMDS scores:
data.scores.nmds = as.data.frame(scores(nmds_bray)$sites)
#Save as csv to add in season, phase, tide, and zone metadata manually
write.csv(data.scores.nmds, "data.scores.nmds.csv")


#Plot season and phase 
library(ggplot2)
ggplot(data.scores.nmds, aes(x=NMDS1, y=NMDS2)) +geom_point(size=2, aes(color=Season, shape=Phase))

#Order seasons and separate out Phases by Season.
data.scores.nmds$Season <- factor(data.scores.nmds$Season, levels = c("SUMMER", "FALL","WINTER","SPRING"))

#Plot
ggplot(data.scores.nmds, aes(x=NMDS1, y=NMDS2)) +geom_point(size=2, aes(color=Phase, shape=Phase)) + coord_equal() + facet_wrap(~Season)

#Incorporate the Tide factor
ggplot(data.scores.nmds, aes(x=NMDS1, y=NMDS2)) +geom_point(size=2, aes(color=Tide, shape=Tide)) + coord_equal() + facet_grid(Phase~Season)


#Repeat NMDS for ASV table containing water samples
#Square root transformation
Water_Diatom_Table_sqrt <-sqrt(Water_Diatom_Table)

#NMDS
water_nmds <- metaMDS(Water_Diatom_Table_sqrt, distance="bray")

#Extract NMDS scores:
data.scores.nmds.water = as.data.frame(scores(water_nmds)$sites)

#Save as a CSV and add in metadata manually
write.csv(data.scores.nmds.water, "Water.NMDS.Scores.csv")

#Now to plot
#Identify the water samples so they can be labeled
to_label<-c("H20FALN","H20FALS", "H20SUMN", "H20SUMS", "H20WINS", "H20WINN", "H20SPRN", "H20SPRS")

#Plot
library(ggrepel)
ggplot(Water_NMDS_Scores, aes(x=NMDS1, y=NMDS2)) +geom_point(size=2, aes(color=Season, shape=Phase)) + geom_text_repel(data = subset(Water_NMDS_Scores,SampleID %in% to_label), aes(label=SampleID))



##Rare
nmds_bray_rare <- metaMDS(rare_asvs.sqrt, distance="bray")
#Extract the scores:
data.scores.nmds_rare = as.data.frame(scores(nmds_bray_rare)$sites)
#Save as csv to add in season, phase, tide and zone metadata manually
write.csv(data.scores.nmds_rare, "data.scores.nmds.rare.csv")

#Plot Season and Phase
ggplot(data.scores.nmds.rare, aes(x=NMDS1, y=NMDS2)) +geom_point(size=2, aes(color=Season, shape=Phase))
#Order seasons and separate out Phases by Season.
data.scores.nmds.rare$Season <- factor(data.scores.nmds.rare$Season, levels = c("SUMMER", "FALL","WINTER","SPRING"))

ggplot(data.scores.nmds.rare, aes(x=NMDS1, y=NMDS2)) +geom_point(size=2, aes(color=Phase, shape=Phase)) + coord_equal() + facet_wrap(~Season)

#Add Tide
ggplot(data.scores.nmds.rare, aes(x=NMDS1, y=NMDS2)) +geom_point(size=2, aes(color=Tide, shape=Tide)) + coord_equal() + facet_grid(Phase~Season)



##Abundant
nmds_bray_abundant <- metaMDS(abundant_asvs.sqrt, distance="bray")
#Extract the scores:
data.scores.nmds_abundant = as.data.frame(scores(nmds_bray_abundant)$sites)
#Save as csv to add in season, phase, tide and zone metadata manually
write.csv(data.scores.nmds_abundant, "data.scores.nmds.abundant.csv")


#Plot Season and Phase 
ggplot(data.scores.nmds.abundant, aes(x=NMDS1, y=NMDS2)) +geom_point(size=2, aes(color=Season, shape=Phase))
#Order seasons and separate out Phases by Season.
data.scores.nmds.abundant$Season <- factor(data.scores.nmds.abundant$Season, levels = c("SUMMER", "FALL","WINTER","SPRING"))

ggplot(data.scores.nmds.abundant, aes(x=NMDS1, y=NMDS2)) +geom_point(size=2, aes(color=Phase, shape=Phase)) + coord_equal() + facet_wrap(~Season)

#Add Tide
ggplot(data.scores.nmds.abundant, aes(x=NMDS1, y=NMDS2)) +geom_point(size=2, aes(color=Tide, shape=Tide)) + coord_equal() + facet_grid(Phase~Season)






##PERMANOVA and PERMDISP
#PERMANOVA for Total ASVs
library(vegan)
#Using Square root transformation table
#Run the addonis function for Season
permanovafit_Season<-adonis2(Diatom_Table_sqrt ~ Season.base, data=sampledata, permutations=999, method="bray")
#Season and phase
permanovafit_Season_Phase<-adonis2(Diatom_Table_sqrt ~ Season, data=sampledata, permutations=999, method="bray")
#Tide
permanovafit_Tide<-adonis2(Diatom_Table_sqrt ~ Tide, data=sampledata, permutations=999, method="bray")
#Phase
permanovafit_Phase<-adonis2(Diatom_Table_sqrt ~ Phase, data=sampledata, permutations=999, method="bray")
#Zone
permanovafit_Zone<-adonis2(Diatom_Table_sqrt ~ Zone, data=sampledata, permutations=999, method="bray")

## PERMDSIP for Total ASVs
#Check for homogeneity of multivariate dispersion using betadisper (PERMDISP)
Distance_Data <-vegdist(Diatom_Table_sqrt)

#Season
anova(betadisper(Distance_Data, sampledata$Season.base))
#Season and Phase
anova(betadisper(Distance_Data, sampledata$Season))
#Tide
anova(betadisper(Distance_Data, sampledata$Tide))
#Phase
anova(betadisper(Distance_Data, sampledata$Phase))
#Zone
anova(betadisper(Distance_Data, sampledata$Zone))


##Compare groups within Phase
#Make individual tables for each season
#one for FALL
FAL_ASV<-Diatom_Table_sqrt[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),]
#one for SPRING
SPR_ASV<-Diatom_Table_sqrt[c(19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36),]
#one for SUMMER
SUM_ASV<-Diatom_Table_sqrt[c(37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54),]
#one for WINTER
WIN_ASV<-Diatom_Table_sqrt[c(55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72),]

#Do the same for metadata
SPR_sampledat<-sampledata[c(19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36),]
FAL_sampledat<-sampledata[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),]
WIN_sampledat<-sampledata[c(55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72),]
SUM_sampledat<-sampledata[c(37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54),]

#PERMANOVA for Phases within each Season
permanovafit_Phase_FAL<-adonis2(FAL_ASV ~ Phase, data=FAL_sampledat, permutations=999, method="bray")
permanovafit_Phase_SPR<-adonis2(SPR_ASV ~ Phase, data=SPR_sampledat, permutations=999, method="bray")
permanovafit_Phase_WIN<-adonis2(WIN_ASV ~ Phase, data=WIN_sampledat, permutations=999, method="bray")
permanovafit_Phase_SUM<-adonis2(SUM_ASV ~ Phase, data=SUM_sampledat, permutations=999, method="bray")

#Check with PERMDISP
FAL_distdata<-vegdist(FAL_ASV)
WIN_distdata<-vegdist(WIN_ASV)
SPR_distdata<-vegdist(SPR_ASV)
SUM_distdata<-vegdist(SUM_ASV)
anova(betadisper(FAL_distdata, FAL_sampledat$Phase)) 
anova(betadisper(WIN_distdata, WIN_sampledat$Phase)) 
anova(betadisper(SPR_distdata, SPR_sampledat$Phase)) 
anova(betadisper(SUM_distdata, SUM_sampledat$Phase)) 


#PERMANOVA for Tides within Season and Phase
#Create sub-tables for each phase per season
FAL_ASV_Neap<-Diatom_Table_sqrt[c(1,2,3,4,5,6,7,8,9),]
FAL_ASV_Spring<-Diatom_Table_sqrt[c(10,11,12,13,14,15,16,17,18),]
SPR_ASV_Neap<-Diatom_Table_sqrt[c(19,20,21,22,23,24,25,26,27),]
SPR_ASV_Spring<-Diatom_Table_sqrt[c(28,29,30,31,32,33,34,35,36),]
SUM_ASV_Neap<-Diatom_Table_sqrt[c(37,38,39,40,41,42,43,44,45),]
SUM_ASV_Spring<-Diatom_Table_sqrt[c(46,47,48,49,50,51,52,53,54),]
WIN_ASV_Neap<-Diatom_Table_sqrt[c(55,56,57,58,59,60,61,62,63),]
WIN_ASV_Spring<-Diatom_Table_sqrt[c(64,65,66,67,68,69,70,71,72),]

#Metadata
FAL_sampledat_Neap<-sampledata[c(1,2,3,4,5,6,7,8,9),]
FAL_sampledat_Spring<-sampledata[c(10,11,12,13,14,15,16,17,18),]
SPR_sampledat_Neap<-sampledata[c(19,20,21,22,23,24,25,26,27),]
SPR_sampledat_Spring<-sampledata[c(28,29,30,31,32,33,34,35,36),]
SUM_sampledat_Neap<-sampledata[c(37,38,39,40,41,42,43,44,45),]
SUM_sampledat_Spring<-sampledata[c(46,47,48,49,50,51,52,53,54),]
WIN_sampledat_Neap<-sampledata[c(55,56,57,58,59,60,61,62,63),]
WIN_sampledat_Spring<-sampledata[c(64,65,66,67,68,69,70,71,72),]

permanovafit_Tide_FAL_Neap<-adonis2(FAL_ASV_Neap ~ Tide, data=FAL_sampledat_Neap, permutations=999, method="bray")
permanovafit_Tide_SPR_Neap<-adonis2(SPR_ASV_Neap ~ Tide, data=SPR_sampledat_Neap, permutations=999, method="bray")
permanovafit_Tide_WIN_Neap<-adonis2(WIN_ASV_Neap ~ Tide, data=WIN_sampledat_Neap, permutations=999, method="bray")
permanovafit_Tide_SUM_Neap<-adonis2(SUM_ASV_Neap ~ Tide, data=SUM_sampledat_Neap, permutations=999, method="bray")
permanovafit_Tide_FAL_Spring<-adonis2(FAL_ASV_Spring ~ Tide, data=FAL_sampledat_Spring, permutations=999, method="bray")
permanovafit_Tide_SPR_Spring<-adonis2(SPR_ASV_Spring ~ Tide, data=SPR_sampledat_Spring, permutations=999, method="bray")
permanovafit_Tide_WIN_Spring<-adonis2(WIN_ASV_Spring ~ Tide, data=WIN_sampledat_Spring, permutations=999, method="bray")
permanovafit_Tide_SUM_Spring<-adonis2(SUM_ASV_Spring ~ Tide, data=SUM_sampledat_Spring, permutations=999, method="bray")

#PERMDIST
FAL_distdata_Neap<-vegdist(FAL_ASV_Neap)
WIN_distdata_Neap<-vegdist(WIN_ASV_Neap)
SPR_distdata_Neap<-vegdist(SPR_ASV_Neap)
SUM_distdata_Neap<-vegdist(SUM_ASV_Neap)
FAL_distdata_Spring<-vegdist(FAL_ASV_Spring)
WIN_distdata_Spring<-vegdist(WIN_ASV_Spring)
SPR_distdata_Spring<-vegdist(SPR_ASV_Spring)
SUM_distdata_Spring<-vegdist(SUM_ASV_Spring)

anova(betadisper(FAL_distdata_Neap, FAL_sampledat_Neap$Tide)) 
anova(betadisper(WIN_distdata_Neap, WIN_sampledat_Neap$Tide)) 
anova(betadisper(SPR_distdata_Neap, SPR_sampledat_Neap$Tide)) 
anova(betadisper(SUM_distdata_Neap, SUM_sampledat_Neap$Tide)) 
anova(betadisper(FAL_distdata_Spring, FAL_sampledat_Spring$Tide)) 
anova(betadisper(WIN_distdata_Spring, WIN_sampledat_Spring$Tide)) 
anova(betadisper(SPR_distdata_Spring, SPR_sampledat_Spring$Tide)) 
anova(betadisper(SUM_distdata_Spring, SUM_sampledat_Spring$Tide)) 




##PERMANOVA for Rare ASVs
#Season - rare
adonis2(rare_asvs.sqrt ~ Season.base, data=sampledata, permutations=999, method="bray")
#Season and phase - rare
adonis2(rare_asvs.sqrt ~ Season, data=sampledata, permutations=999, method="bray")
#Tide - rare
adonis2(rare_asvs.sqrt ~ Tide, data=sampledata, permutations=999, method="bray")
#Phase - rare
adonis2(rare_asvs.sqrt ~ Phase, data=sampledata, permutations=999, method="bray")
#Zone - rare
adonis2(rare_asvs.sqrt ~ Zone, data=sampledata, permutations=999, method="bray")

# PERMDISP - rare
#Check for homogeneity of multivariate dispersion using betadisper (PERMDISP)
Dis_Dat_rare <-vegdist(rare_asvs.sqrt)

#Season - rare
anova(betadisper(Dis_Dat_rare, sampledata$Season.base))
#Season and Phase - rare
anova(betadisper(Dis_Dat_rare, sampledata$Season))
#Tide - rare
anova(betadisper(Dis_Dat_rare, sampledata$Tide))
#Phase - rare
anova(betadisper(Dis_Dat_rare, sampledata$Phase))
#Zone - rare
anova(betadisper(Dis_Dat_rare, sampledata$Zone))

#Compare groups within Phase
#Make individual tables for each season
#FALL
FAL_ASV_rare<-rare_asvs.sqrt[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),]
#SPRING
SPR_ASV_rare<-rare_asvs.sqrt[c(19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36),]
#SUMMER
SUM_ASV_rare<-rare_asvs.sqrt[c(37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54),]
#WINTER
WIN_ASV_rare<-rare_asvs.sqrt[c(55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72),]

#Do the same form the metadata
SPR_sampledat<-sampledata[c(19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36),]
FAL_sampledat<-sampledata[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),]
WIN_sampledat<-sampledata[c(55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72),]
SUM_sampledat<-sampledata[c(37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54),]

#PERMANOVA for Phases within each Season
permanovafit_Phase_FAL_rare<-adonis2(FAL_ASV_rare ~ Phase, data=FAL_sampledat, permutations=999, method="bray")
permanovafit_Phase_FAL_rare  

permanovafit_Phase_SPR_rare<-adonis2(SPR_ASV_rare ~ Phase, data=SPR_sampledat, permutations=999, method="bray")
permanovafit_Phase_SPR_rare

permanovafit_Phase_WIN_rare<-adonis2(WIN_ASV_rare ~ Phase, data=WIN_sampledat, permutations=999, method="bray")
permanovafit_Phase_WIN_rare

permanovafit_Phase_SUM_rare<-adonis2(SUM_ASV_rare ~ Phase, data=SUM_sampledat, permutations=999, method="bray")
permanovafit_Phase_SUM_rare

#Check with PERMDISP
FAL_distdata_rare<-vegdist(FAL_ASV_rare)
WIN_distdata_rare<-vegdist(WIN_ASV_rare)
SPR_distdata_rare<-vegdist(SPR_ASV_rare)
SUM_distdata_rare<-vegdist(SUM_ASV_rare)
anova(betadisper(FAL_distdata_rare, FAL_sampledat$Phase)) 
anova(betadisper(WIN_distdata_rare, WIN_sampledat$Phase)) 
anova(betadisper(SPR_distdata_rare, SPR_sampledat$Phase)) 
anova(betadisper(SUM_distdata_rare, SUM_sampledat$Phase)) 


#PERMANOVA for Tides within Season and Phase
#Create sub-tables for each phase per season
FAL_ASV_Neap_rare<-rare_asvs.sqrt[c(1,2,3,4,5,6,7,8,9),]
FAL_ASV_Spring_rare<-rare_asvs.sqrt[c(10,11,12,13,14,15,16,17,18),]
SPR_ASV_Neap_rare<-rare_asvs.sqrt[c(19,20,21,22,23,24,25,26,27),]
SPR_ASV_Spring_rare<-rare_asvs.sqrt[c(28,29,30,31,32,33,34,35,36),]
SUM_ASV_Neap_rare<-rare_asvs.sqrt[c(37,38,39,40,41,42,43,44,45),]
SUM_ASV_Spring_rare<-rare_asvs.sqrt[c(46,47,48,49,50,51,52,53,54),]
WIN_ASV_Neap_rare<-rare_asvs.sqrt[c(55,56,57,58,59,60,61,62,63),]
WIN_ASV_Spring_rare<-rare_asvs.sqrt[c(64,65,66,67,68,69,70,71,72),]

#Metadata
FAL_sampledat_Neap<-sampledata[c(1,2,3,4,5,6,7,8,9),]
FAL_sampledat_Spring<-sampledata[c(10,11,12,13,14,15,16,17,18),]
SPR_sampledat_Neap<-sampledata[c(19,20,21,22,23,24,25,26,27),]
SPR_sampledat_Spring<-sampledata[c(28,29,30,31,32,33,34,35,36),]
SUM_sampledat_Neap<-sampledata[c(37,38,39,40,41,42,43,44,45),]
SUM_sampledat_Spring<-sampledata[c(46,47,48,49,50,51,52,53,54),]
WIN_sampledat_Neap<-sampledata[c(55,56,57,58,59,60,61,62,63),]
WIN_sampledat_Spring<-sampledata[c(64,65,66,67,68,69,70,71,72),]

permanovafit_Tide_FAL_Neap_rare<-adonis2(FAL_ASV_Neap_rare ~ Tide, data=FAL_sampledat_Neap, permutations=999, method="bray")
permanovafit_Tide_SPR_Neap_rare<-adonis2(SPR_ASV_Neap_rare ~ Tide, data=SPR_sampledat_Neap, permutations=999, method="bray")
permanovafit_Tide_WIN_Neap_rare<-adonis2(WIN_ASV_Neap_rare ~ Tide, data=WIN_sampledat_Neap, permutations=999, method="bray")
permanovafit_Tide_SUM_Neap_rare<-adonis2(SUM_ASV_Neap_rare ~ Tide, data=SUM_sampledat_Neap, permutations=999, method="bray")
permanovafit_Tide_FAL_Spring_rare<-adonis2(FAL_ASV_Spring_rare ~ Tide, data=FAL_sampledat_Spring, permutations=999, method="bray")
permanovafit_Tide_SPR_Spring_rare<-adonis2(SPR_ASV_Spring_rare ~ Tide, data=SPR_sampledat_Spring, permutations=999, method="bray")
permanovafit_Tide_WIN_Spring_rare<-adonis2(WIN_ASV_Spring_rare ~ Tide, data=WIN_sampledat_Spring, permutations=999, method="bray")
permanovafit_Tide_SUM_Spring_rare<-adonis2(SUM_ASV_Spring_rare ~ Tide, data=SUM_sampledat_Spring, permutations=999, method="bray")


#PERMDIST
FAL_distdata_Neap_rare<-vegdist(FAL_ASV_Neap_rare)
WIN_distdata_Neap_rare<-vegdist(WIN_ASV_Neap_rare)
SPR_distdata_Neap_rare<-vegdist(SPR_ASV_Neap_rare)
SUM_distdata_Neap_rare<-vegdist(SUM_ASV_Neap_rare)
FAL_distdata_Spring_rare<-vegdist(FAL_ASV_Spring_rare)
WIN_distdata_Spring_rare<-vegdist(WIN_ASV_Spring_rare)
SPR_distdata_Spring_rare<-vegdist(SPR_ASV_Spring_rare)
SUM_distdata_Spring_rare<-vegdist(SUM_ASV_Spring_rare)

anova(betadisper(FAL_distdata_Neap_rare, FAL_sampledat_Neap$Tide)) 
anova(betadisper(WIN_distdata_Neap_rare, WIN_sampledat_Neap$Tide))
anova(betadisper(SPR_distdata_Neap_rare, SPR_sampledat_Neap$Tide)) 
anova(betadisper(SUM_distdata_Neap_rare, SUM_sampledat_Neap$Tide)) 
anova(betadisper(FAL_distdata_Spring_rare, FAL_sampledat_Spring$Tide)) 
anova(betadisper(WIN_distdata_Spring_rare, WIN_sampledat_Spring$Tide)) 
anova(betadisper(SPR_distdata_Spring_rare, SPR_sampledat_Spring$Tide)) 
anova(betadisper(SUM_distdata_Spring_rare, SUM_sampledat_Spring$Tide)) 



#PERMANOVA for abundant
#Season - abundant
adonis2(abundant_asvs.sqrt ~ Season.base, data=sampledata, permutations=999, method="bray")
#Season and phase - abundant
adonis2(abundant_asvs.sqrt ~ Season, data=sampledata, permutations=999, method="bray")
#Now for tide - abundant
adonis2(abundant_asvs.sqrt ~ Tide, data=sampledata, permutations=999, method="bray")
#Phase - abundant
adonis2(abundant_asvs.sqrt ~ Phase, data=sampledata, permutations=999, method="bray")
#Zone - abundant
adonis2(abundant_asvs.sqrt ~ Zone, data=sampledata, permutations=999, method="bray")

# PERMDISP - abundant
#Check for homogeneity of multivariate dispersion using betadisper (PERMDISP)
Dis_Dat_abundant <-vegdist(abundant_asvs.sqrt)

#Season - abundant
anova(betadisper(Dis_Dat_abundant, sampledata$Season.base))
#Season and Phase - abundant
anova(betadisper(Dis_Dat_abundant, sampledata$Season))
#Tide - abundant
anova(betadisper(Dis_Dat_abundant, sampledata$Tide))
#Phase - abundant
anova(betadisper(Dis_Dat_abundant, sampledata$Phase))
#Zone - abundant
anova(betadisper(Dis_Dat_abundant, sampledata$Zone))






#BIOENV

#Load vegan
library(vegan)
library(usdm)

#Upload the Square-root transformed Diatom_Table and remove the SUMNTE50 Row because it had some NAs

#Scale the environmental data
scaled_dat<-as.data.frame(scale(env_dat_noNA))
#Check the variance inflation factor for each variable
vif(env_dat_full)

#Percent sand, percent mud, and the meiofauna parameters have high 
#"infinite" variance inflation factor values (i.e. perfect correlation). 
#The eight meiofauna parameters are highly correlated. Remove
#ostrocod number parameter because ostracods are not nearly as prevalent
No_ostracod_number <- scaled_dat[ , !(names(scaled_dat)) %in% "Ostracod.."]
vif(No_ostracod_number)

#This reduced VIF values for each of the meiofauna number parameters (no longer INF)
#Repeat with Ostracod mass (Ost.Mass)
No_ostracod_mass <- No_ostracod_number[ , !(names(No_ostracod_number)) %in% "Ost.Mass"]
vif(No_ostracod_mass)

#Percent sand and mud are the only parameters with infinite VIFs now.
#The samples are over 90 percent mud that value is likely more useful
#Remove percent sand parameter
No_sand <- No_ostracod_mass[ , !(names(No_ostracod_mass)) %in% "X..SAND."]
vif(No_sand)

#No parameters have infinite VIFs now.
#Continue removing the highest VIFs until all are below a threshold of 10
No_total_meio <- No_sand[ , !(names(No_sand)) %in% "Total.."]
vif(No_total_meio)
No_mud <- No_total_meio[ , !(names(No_total_meio)) %in% "X..MUD."]
vif(No_mud)
No_sed <- No_mud[ , !(names(No_mud)) %in% "Sed.temp"]
vif(No_sed)
No_air <- No_sed[ , !(names(No_sed)) %in% "air.temp"]
vif(No_air)

#All VIFs below threshold of 10

#BIOENV Up to 10 combinations
BIOENV <- "bioenv"(Diatom_Table, No_air, method = "spearman", index = "bray", upto = 10, trace = T, partial = NULL)
#Best model has 4 parameters (max. 10 allowed):
#  Water.temp Spartina....alive. Spartine.Height..Alive. Spartine.Height..Dead.
#with correlation  0.4285206 

#Use output of BIOENV in RDA
bray_dist<-vegdist(Diatom_Table, method="bray")
bray_dist<-as.matrix(bray_dist)
bray_dist<-as.data.frame(bray_dist)

euc_time <- vegdist(time_dat$hrs, method = "euc")
euc_time<- as.matrix(euc_time)
euc_time<-as.data.frame(euc_time)
scaled_dat<- as.data.frame(scaled_dat)
time_pcnm<- pcnm(euc_time)

dbRDA.bray = capscale(bray_dist ~ Water.temp + Spartina....alive. + Spartine.Height..Alive. + Spartine.Height..Dead. + Condition(time_pcnm$vectors), scaled_dat)

#Water contribution to RDA
dbRDA.water = capscale(bray_dist ~ Water.temp + Condition(time_pcnm$vectors), scaled_dat)
RsquareAdj(dbRDA.water)$adj.r.squared * 100
#0.2270481 percent

#Spartina....alive.   contribution to RDA
dbRDA.sparalive = capscale(bray_dist ~ Spartina....alive. + Condition(time_pcnm$vectors), scaled_dat)
RsquareAdj(dbRDA.sparalive)$adj.r.squared * 100
#1.638168 percent

#Spartina.height.alive.   contribution to RDA
dbRDA.sparheightalive = capscale(bray_dist ~ Spartine.Height..Alive. + Condition(time_pcnm$vectors), scaled_dat)
RsquareAdj(dbRDA.sparheightalive)$adj.r.squared * 100
#2.201032 percent

#Spartina.Height..Dead.   contribution to RDA
dbRDA.sparheightdead = capscale(bray_dist ~ Spartine.Height..Dead. + Condition(time_pcnm$vectors), scaled_dat)
RsquareAdj(dbRDA.sparheightdead)$adj.r.squared * 100
#0.5047229 percent


#View Results
anova(dbRDA.bray)
anova(dbRDA.bray, by="terms", permu=200)

#Extract the scores
site.scores <- scores(dbRDA.bray, display = "sites", scaling =1)
env.scores <- scores(dbRDA.bray, display = "bp", scaling =1)
#Convert to dataframes and add in season and phase meta data to site scores
sites_df<-as.data.frame(site.scores)
sites_df$Season <-sampledata$Season.base
sites_df$Phase <- sampledata$Phase


env_df <- as.data.frame(env.scores)
#Re-label variables for plot
env_df$Variable <- c("WT", "LSN", "LSH", "DSH" )


#Plot RDA
library (ggplot2)
library(ggrepel)
g.bray<- ggplot()
RDS_plot<- g.bray + geom_point(data=sites_df, aes(x=CAP1, y=CAP2, color=Season, shape=Phase), size = 3) + geom_segment(data = env_df, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), arrow = arrow(length = unit(0.25, "cm")), colour='red')+geom_text_repel (data = env_df, aes(x=CAP1, y=CAP2, label = Variable), color = "black", max.overlaps = 20)

#Variation partitioning analysis
euc_time<- as.matrix(euc_time)
bray_varpart<-varpart(bray_dist, ~scaled_dat$Water.temp + scaled_dat$Spartina....alive. + scaled_dat$Spartine.Height..Alive. + scaled_dat$Spartine.Height..Dead., ~euc_time)

#Plot
VarPlot<-plot(bray_varpart, Xnames = c("ENV", "TEMP"))


#Rare BIOENV
#Use the environmental data culled in the BIOENV analysis for Total ASVs
#BIOENV analysis (BIOENV_Parameters.csv)
#Should be 24 variables

#Use rare_asvs_noNA.sqrt table. Removed SUMNTE50 because of NAs

#Compute BIOENV with up to 10 combinations
BIOENV_Rare <- "bioenv"(rare_asvs_noNA.sqrt, BIOENV_Parameters, method = "spearman", index = "bray", upto = 10, trace = T, partial = NULL)

#Use output of BIOENV in RDA
scaled_data <-scale(env_dat)
bray_dist_rare<-vegdist(rare_asvs.sqrt, method="bray")
bray_dist_rare<-as.matrix(bray_dist_rare)
bray_dist_rare<-as.data.frame(bray_dist_rare)

euc_time <- vegdist(time_dat$hrs, method = "euc")
euc_time<- as.matrix(euc_time)
euc_time<-as.data.frame(euc_time)
scaled_data<- as.data.frame(scaled_data)
time_pcnm<- pcnm(euc_time)

dbRDA.bray.rare = capscale(bray_dist_rare ~ Water.temp + Spartina....alive. + Spartine.Height..Alive. + Spartine.Height..Dead. + NO23..mg.N.L. + MEAN..um. + Condition(time_pcnm$vectors), scaled_data)

#View results
anova(dbRDA.bray.rare)
anova(dbRDA.bray.rare, by="terms", permu=200)

#Extract the scores
site.scores.rare <- scores(dbRDA.bray.rare, display = "sites", scaling =1)
env.scores.rare <- scores(dbRDA.bray.rare, display = "bp", scaling =1)
#Convert to dataframes and add in season and phase meta data to site scores
sites_df_rare<-as.data.frame(site.scores.rare)
sites_df_rare$Season <-sampledata$Season
sites_df_rare$Phase <- sampledata$Phase


env_df_rare <- as.data.frame(env.scores.rare)
#Update variable names for plotting
env_df_rare$Variable <- c("WT", "LSN", "LSH", "DSH", "NO23", "MGS")


#Plot RDA
g.bray.rare<- ggplot()
RDS_plot_rare<- g.bray.rare + geom_point(data=sites_df_rare, aes(x=CAP1, y=CAP2, color=Season, shape=Phase), size = 3) + geom_segment(data = env_df_rare, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), arrow = arrow(length = unit(0.25, "cm")), colour='red')+geom_text_repel (data = env_df_rare, aes(x=CAP1, y=CAP2, label = Variable), color = "black", max.overlaps = 20)

#Variation partitioning analysis.
euc_time<- as.matrix(euc_time)
bray_varpart_rare<-varpart(bray_dist_rare, ~scaled_data$Water.temp + scaled_data$Spartina....alive. + scaled_data$Spartine.Height..Alive. + scaled_data$Spartine.Height..Dead. + scaled_data$NO23..mg.N.L.+scaled_data$MEAN..um., ~euc_time)

#Plot
VarPlot_rare<-plot(bray_varpart_rare, Xnames = c("ENV", "TEMP"))


#Abundant ASVs
#BIOENV
#Use the environmental data culled in the BIOENV analysis for Total ASVs
#BIOENV analysis (BIOENV_Parameters.csv)
#Should be 24 variables

#Use abundant_asvs.noNA.sqrt table for BIOENV. Remove SUMNTE50 because of NAs


#BIOENV with up to 10 combinations
BIOENV_abundant <- "bioenv"(abundant_asvs.noNA.sqrt, BIOENV_Parameters, method = "spearman", index = "bray", upto = 10, trace = T, partial = NULL)

#Use output of BIOENV in RDA
scaled_data <-scale(env_dat)
bray_dist_abundant<-vegdist(abundant_asvs.sqrt, method="bray")
bray_dist_abundant<-as.matrix(bray_dist_abundant)
bray_dist_abundant<-as.data.frame(bray_dist_abundant)

euc_time <- vegdist(time_dat$hrs, method = "euc")
euc_time<- as.matrix(euc_time)
euc_time<-as.data.frame(euc_time)
scaled_data<- as.data.frame(scaled_data)
time_pcnm<- pcnm(euc_time)

dbRDA.bray.abundant = capscale(bray_dist_abundant ~ Water.temp + Spartina....alive. + Spartine.Height..Alive. + Spartine.Height..Dead. + Condition(time_pcnm$vectors), scaled_data)

#View results
anova(dbRDA.bray.abundant)
anova(dbRDA.bray.abundant, by="terms", permu=200)


#Extract RDA scores
site.scores.abundant <- scores(dbRDA.bray.abundant, display = "sites", scaling =1)
env.scores.abundant <- scores(dbRDA.bray.abundant, display = "bp", scaling =1)
#Convert to dataframes and add in season and phase meta data to site scores
sites_df_abundant<-as.data.frame(site.scores.abundant)
sites_df_abundant$Season <-sampledata$Season.base
sites_df_abundant$Phase <- sampledata$Phase


env_df_abundant <- as.data.frame(env.scores.abundant)
#Update variable names for plotting
env_df_abundant$Variable <- c("WT", "LSN", "LSH", "DSH")

#Plot
g.bray.abundant<- ggplot()
RDS_plot_abundant<- g.bray.abundant + geom_point(data=sites_df_abundant, aes(x=CAP1, y=CAP2, color=Season, shape=Phase), size = 3) + geom_segment(data = env_df_abundant, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), arrow = arrow(length = unit(0.25, "cm")), colour='red')+geom_text_repel (data = env_df_abundant, aes(x=CAP1, y=CAP2, label = Variable), color = "black", max.overlaps = 20)


#Variation partitioning analysis.
euc_time<- as.matrix(euc_time)
bray_varpart_abundant<-varpart(bray_dist_abundant, ~scaled_data$Water.temp + scaled_data$Spartina....alive. + scaled_data$Spartine.Height..Alive. + scaled_data$Spartine.Height..Dead., ~euc_time)
#Plot
VarPlot_abundant<-plot(bray_varpart_abundant, Xnames = c("ENV", "TEMP"))

              
              




##Alpha Diversity for Total ASVs
library(tidyverse)
library(vegan)
library(ggplot2)
library(FSA)

#Calculate Diversity Indices
#Richness
richness<-specnumber(Diatom_Table)
#shannon-weiner
shannon <- diversity(Diatom_Table, index = "shannon")
#Pielou's Evenness
pielou <- shannon / log(richness)
#Chao1
chao1 <- estimateR(Diatom_Table)["S.chao1",]

#Combine results
alpha_div <- data.frame(SampleID = metadat$SampleID, Richness = richness, Chao1 = chao1, Shannon = shannon, Pielou = pielou)%>% left_join(metadat, by = "SampleID")
#Transform to plot
alpha_long <- alpha_div %>% pivot_longer(cols = c("Richness", "Chao1", "Shannon", "Pielou"), names_to = "Metric", values_to = "Value")
#Set the order I want the seasons to appear
alpha_long$Season <- factor(alpha_long$Season, levels = c("SUMMER", "FALL", "WINTER","SPRING"))
#Combine season and phase factors
alpha_long <- alpha_long %>% mutate (SeasonPhase = interaction(Season, Phase, sep = "_"))
#Reorder neap and spring
alpha_long$SeasonPhase <- factor(alpha_long$SeasonPhase, levels = c("SUMMER_N", "SUMMER_S", "FALL_N", "FALL_S", "WINTER_N", "WINTER_S", "SPRING_N", "SPRING_S"))

#Plot results
ggplot(alpha_long, aes(x=SeasonPhase, y = Value, fill = Season)) + geom_boxplot(alpha = 0.85) + facet_wrap(~ Metric, scales = "free_y") + scale_x_discrete(labels = function(x) gsub("_","\n", x)) + labs(x = "Season and Tidal Phase", y = "Diversity Value", fill = "Season")



#Add tide
#Combine results
alpha_div <- data.frame(SampleID = sampledata$SampleID, Richness = richness, Chao1 = chao1, Shannon = shannon, Pielou = pielou)%>% left_join(sampledata, by = "SampleID")
#Transform to plot
alpha_long_tide <- alpha_div_tide %>% pivot_longer(cols = c("Richness", "Chao1", "Shannon", "Pielou"), names_to = "Metric", values_to = "Value")
#Set the order for season
alpha_long_tide$Season <- factor(alpha_long_tide$Season, levels = c("SUMMER", "FALL", "WINTER","SPRING"))
#Combine season and phase factors
alpha_long_tide <- alpha_long_tide %>% mutate (SeasonPhase = interaction(Season, Phase, Tide, sep = "_"))
#Reorder neap and spring
alpha_long_tide$SeasonPhaseTide <- factor(alpha_long_tide$SeasonPhaseTide, levels = c("SUMMER_N_TE","SUMMER_N_TL", "SUMMER_N_TF", "SUMMER_S_TE", "SUMMER_S_TL", "SUMMER_S_TF","FALL_N_TE","FALL_N_TL", "FALL_N_TF", "FALL_S_TE", "FALL_S_TL", "FALL_S_TF", "WINTER_N_TE", "WINTER_N_TL", "WINTER_N_TF", "WINTER_S_TE", "WINTER_S_TL", "WINTER_S_TF", "SPRING_N_TE", "SPRING_N_TL", "SPRING_N_TF", "SPRING_S_TE", "SPRING_S_TL", "SPRING_S_TF"))

#Plot results
ggplot(alpha_long_tide, aes(x=SeasonPhase, y = Value, fill = Season)) + geom_boxplot(alpha = 0.85) + facet_wrap(~ Metric, scales = "free_y") + scale_x_discrete(labels = function(x) gsub("_","\n", x)) + labs(x = "Season, Lunar Phase, and Tidal Phase", y = "Diversity Value", fill = "Season.base")


#Compare with ANOVA or Krusakal-Wallis
library(dplyr)
metrics <- c("Richness", "Chao1", "Shannon", "Pielou")
results <- list()

for (m in metrics) {
  for (s in unique(alpha_div$Season)) {
    data_sub <- alpha_div %>% filter(Season == s)
    
    # Check normality for each lunar phase group
    shapiro_p <- data_sub %>%
      group_by(Phase) %>%
      summarise(p_value = shapiro.test(get(m))$p.value)
    
    # If *all* groups are roughly normal (p > 0.05), use ANOVA
    if (all(shapiro_p$p_value > 0.05)) {
      test <- aov(get(m) ~ Phase, data = data_sub)
      pval <- summary(test)[[1]][["Pr(>F)"]][1]
      method <- "ANOVA"
    } else {
      test <- kruskal.test(get(m) ~ Phase, data = data_sub)
      pval <- test$p.value
      method <- "Kruskal-Wallis"
    }
    
    results[[paste(m, s, sep = "_")]] <- data.frame(
      Metric = m,
      Season = s,
      Method = method,
      P_value = pval
    )
  }
}

results_df <- do.call(rbind, results)
results_df


write.csv(results_df, "Anova_results_Season_Phase.csv")



#Do the same for Tide
results_tide <- list()

for (m in metrics) {
  for (s in unique(Alpha_Diversity_Values$Season)) {
    for (l in unique(Alpha_Diversity_Values$Phase)) {
      
      data_sub <- Alpha_Diversity_Values %>% 
        filter(Season == s, Phase == l)
      
      # Skip if there’s not enough data (e.g., only one tide category)
      if (length(unique(data_sub$Tide)) < 2 || nrow(data_sub) < 5) next
      
      # Check normality per tidal cycle group
      shapiro_p <- data_sub %>%
        group_by(Tide) %>%
        summarise(p_value = tryCatch(shapiro.test(get(m))$p.value, 
                                     error = function(e) NA))
      
      # Decide test based on normality
      if (all(shapiro_p$p_value > 0.05, na.rm = TRUE)) {
        test <- aov(get(m) ~ Tide, data = data_sub)
        pval <- summary(test)[[1]][["Pr(>F)"]][1]
        method <- "ANOVA"
      } else {
        test <- kruskal.test(get(m) ~ Tide, data = data_sub)
        pval <- test$p.value
        method <- "Kruskal-Wallis"
      }
      
      results_tide[[paste(m, s, l, sep = "_")]] <- data.frame(
        Metric = m,
        Season = s,
        LunarPhase = l,
        Method = method,
        P_value = pval
      )
    }
  }
}

results_tide_df <- do.call(rbind, results_tide)
results_tide_df

write.csv (results_tide_df, "Anova_Tide.csv")




#Season Only
#Richness
richness_anova <- aov(Alpha_Diversity_Values$Richness ~ Alpha_Diversity_Values$Season)
shapiro.test(residuals(richness_anova))
#p-value = 0.1323 normal
summary(richness_anova)
TukeyHSD(richness_anova)

#Shannon
shannon_anova <- aov(Alpha_Diversity_Values$Shannon ~ Alpha_Diversity_Values$Season)
shapiro.test(residuals(shannon_anova))
#p-value = 0.0003948 Not normal switch to Kruskal Wallis
kruskal.test(Alpha_Diversity_Values$Shannon ~ Alpha_Diversity_Values$Season)
dunnTest(Alpha_Diversity_Values$Shannon ~ Alpha_Diversity_Values$Season, method = "bh")


#Chao1
Chao1_anova <- aov(Alpha_Diversity_Values$Chao1 ~ Alpha_Diversity_Values$Season)
shapiro.test(residuals(Chao1_anova))
#p-value = 0.4088 normal
summary(Chao1_anova)
TukeyHSD(Chao1_anova)


#Pielou
Pielou_anova <- aov(Alpha_Diversity_Values$Pielou ~ Alpha_Diversity_Values$Season)
shapiro.test(residuals(Pielou_anova))
#p-value = 2.38e-05 not normal switch to Krusakl-Wallis
kruskal.test(Alpha_Diversity_Values$Pielou ~ Alpha_Diversity_Values$Season)
dunnTest(Alpha_Diversity_Values$Pielou ~ Alpha_Diversity_Values$Season, method = "bh")




    

##Alpha Diversity for Rare ASVs
library(tidyverse)
library(vegan)
library(ggplot2)
library(FSA)

#Extract rare ASVs
total_abund <- colSums (Diatom_Table)
#Relative abundance
rel_abund <- total_abund / sum(total_abund)
#Rare ASVs
rare_asvs <- Diatom_Table[, rel_abund < 0.001]


#Calculate Diversity Indices
#Richness
richness<-specnumber(rare_asvs)
#shannon-weiner
shannon <- diversity(rare_asvs, index = "shannon")
#Pielou's Evenness
pielou <- shannon / log(richness)
#Chao1
chao1 <- estimateR(rare_asvs)["S.chao1",]

#Combine results
alpha_div <- data.frame(SampleID = metadat$SampleID, Richness = richness, Chao1 = chao1, Shannon = shannon, Pielou = pielou)%>% left_join(metadat, by = "SampleID")
#Add tide
alpha_div <- data.frame(SampleID = sampledata$SampleID, Richness = richness, Chao1 = chao1, Shannon = shannon, Pielou = pielou)%>% left_join(sampledata, by = "SampleID")
#save as CSV
write.csv(alpha_div, "Rare_alpha.csv")


#Compare with ANOVA or Krusakal-Wallis
library(dplyr)
metrics <- c("Richness", "Chao1", "Shannon", "Pielou")
results <- list()

for (m in metrics) {
  for (s in unique(alpha_div$Season)) {
    data_sub <- alpha_div %>% filter(Season == s)
    
    # Check normality for each lunar phase group
    shapiro_p <- data_sub %>%
      group_by(Phase) %>%
      summarise(p_value = shapiro.test(get(m))$p.value)
    
    # If *all* groups are roughly normal (p > 0.05), use ANOVA
    if (all(shapiro_p$p_value > 0.05)) {
      test <- aov(get(m) ~ Phase, data = data_sub)
      pval <- summary(test)[[1]][["Pr(>F)"]][1]
      method <- "ANOVA"
    } else {
      test <- kruskal.test(get(m) ~ Phase, data = data_sub)
      pval <- test$p.value
      method <- "Kruskal-Wallis"
    }
    
    results[[paste(m, s, sep = "_")]] <- data.frame(
      Metric = m,
      Season = s,
      Method = method,
      P_value = pval
    )
  }
}

results_df <- do.call(rbind, results)
results_df


write.csv(results_df, "Rare_Anova_results_Season_Phase.csv")



#Do the same for Tide
results_tide <- list()

for (m in metrics) {
  for (s in unique(alpha_div_tide$Season.base)) {
    for (l in unique(alpha_div_tide$Phase)) {
      
      data_sub <- alpha_div_tide %>% 
        filter(Season.base == s, Phase == l)
      
      # Skip if there’s not enough data (e.g., only one tide category)
      if (length(unique(data_sub$Tide)) < 2 || nrow(data_sub) < 5) next
      
      # Check normality per tidal cycle group
      shapiro_p <- data_sub %>%
        group_by(Tide) %>%
        summarise(p_value = tryCatch(shapiro.test(get(m))$p.value, 
                                     error = function(e) NA))
      
      # Decide test based on normality
      if (all(shapiro_p$p_value > 0.05, na.rm = TRUE)) {
        test <- aov(get(m) ~ Tide, data = data_sub)
        pval <- summary(test)[[1]][["Pr(>F)"]][1]
        method <- "ANOVA"
      } else {
        test <- kruskal.test(get(m) ~ Tide, data = data_sub)
        pval <- test$p.value
        method <- "Kruskal-Wallis"
      }
      
      results_tide[[paste(m, s, l, sep = "_")]] <- data.frame(
        Metric = m,
        Season = s,
        Phase = l,
        Method = method,
        P_value = pval
      )
    }
  }
}

results_tide_df <- do.call(rbind, results_tide)
results_tide_df

write.csv (results_tide_df, "Rare_Anova_Results_Tide.csv")




#For Season
#Richness
richness_anova <- aov(alpha_div$Richness ~ alpha_div$Season)
shapiro.test(residuals(richness_anova))
#p-value = 0.04474 not normal
kruskal.test(alpha_div$Richness ~ alpha_div$Season)
dunnTest(alpha_div$Richness ~ alpha_div$Season, method = "bh")

#Shannon
shannon_anova <- aov(alpha_div$Shannon ~ alpha_div$Season)
shapiro.test(residuals(shannon_anova))
#p-value = 0.1811 normal 
TukeyHSD(shannon_anova)

#Chao1
Chao1_anova <- aov(alpha_div$Chao1 ~ alpha_div$Season)
shapiro.test(residuals(Chao1_anova))
#p-value = 0.2157 normal
summary(Chao1_anova)
TukeyHSD(Chao1_anova)

#Pielou
Pielou_anova <- aov(alpha_div$Pielou ~ alpha_div$Season)
shapiro.test(residuals(Pielou_anova))
#p-value = 0.03692 not normal switch to Krusakl-Wallis
kruskal.test(alpha_div$Pielou ~ alpha_div$Season)
dunnTest(alpha_div$Pielou ~ alpha_div$Season, method = "bh")






##Alpha Diversity for Abundant ASVs
library(tidyverse)
library(vegan)
library(ggplot2)
library(FSA)

#Extract rare and abundant ASVs
total_abund <- colSums (Diatom_Table)
#Relative abundance
rel_abund <- total_abund / sum(total_abund)
#Rare and Abundant ASVs
rare_asvs <- Diatom_Table[, rel_abund < 0.001]
abundant_asvs <- Diatom_Table[, rel_abund >= 0.001]

#Calculate Diversity Indices
#Richness
richness<-specnumber(abundant_asvs)
#shannon-weiner
shannon <- diversity(abundant_asvs, index = "shannon")
#Pielou's Evenness
pielou <- shannon / log(richness)
#Chao1
chao1 <- estimateR(abundant_asvs)["S.chao1",]

#Combine results
alpha_div <- data.frame(SampleID = metadat$SampleID, Richness = richness, Chao1 = chao1, Shannon = shannon, Pielou = pielou)%>% left_join(metadat, by = "SampleID")
#Add tide
alpha_div <- data.frame(SampleID = metadat$SampleID, Richness = richness, Chao1 = chao1, Shannon = shannon, Pielou = pielou)%>% left_join(metadat, by = "SampleID")
#save as CSV
write.csv(alpha_div, "abundant_alpha.csv")


#Compare with ANOVA or Krusakal-Wallis
library(dplyr)
metrics <- c("Richness", "Chao1", "Shannon", "Pielou")
results <- list()

for (m in metrics) {
  for (s in unique(alpha_div$Season)) {
    data_sub <- alpha_div %>% filter(Season == s)
    
    # Check normality for each lunar phase group
    shapiro_p <- data_sub %>%
      group_by(Phase) %>%
      summarise(p_value = shapiro.test(get(m))$p.value)
    
    # If *all* groups are roughly normal (p > 0.05), use ANOVA
    if (all(shapiro_p$p_value > 0.05)) {
      test <- aov(get(m) ~ Phase, data = data_sub)
      pval <- summary(test)[[1]][["Pr(>F)"]][1]
      method <- "ANOVA"
    } else {
      test <- kruskal.test(get(m) ~ Phase, data = data_sub)
      pval <- test$p.value
      method <- "Kruskal-Wallis"
    }
    
    results[[paste(m, s, sep = "_")]] <- data.frame(
      Metric = m,
      Season = s,
      Method = method,
      P_value = pval
    )
  }
}

results_df <- do.call(rbind, results)
results_df


write.csv(results_df, "Abundant_Anova_results_Season_Phase.csv")



#Do the same for Tide
results_tide <- list()

for (m in metrics) {
  for (s in unique(alpha_div$Season)) {
    for (l in unique(alpha_div$Phase)) {
      
      data_sub <- alpha_div %>% 
        filter(Season == s, Phase == l)
      
      # Skip if there’s not enough data (e.g., only one tide category)
      if (length(unique(data_sub$Tide)) < 2 || nrow(data_sub) < 5) next
      
      # Check normality per tidal cycle group
      shapiro_p <- data_sub %>%
        group_by(Tide) %>%
        summarise(p_value = tryCatch(shapiro.test(get(m))$p.value, 
                                     error = function(e) NA))
      
      # Decide test based on normality
      if (all(shapiro_p$p_value > 0.05, na.rm = TRUE)) {
        test <- aov(get(m) ~ Tide, data = data_sub)
        pval <- summary(test)[[1]][["Pr(>F)"]][1]
        method <- "ANOVA"
      } else {
        test <- kruskal.test(get(m) ~ Tide, data = data_sub)
        pval <- test$p.value
        method <- "Kruskal-Wallis"
      }
      
      results_tide[[paste(m, s, l, sep = "_")]] <- data.frame(
        Metric = m,
        Season = s,
        Phase = l,
        Method = method,
        P_value = pval
      )
    }
  }
}

results_tide_df <- do.call(rbind, results_tide)
results_tide_df

write.csv (results_tide_df, "Abundant_Anova_Results_Tide.csv")




#Just for Season
#Richness
richness_anova <- aov(alpha_div$Richness ~ alpha_div$Season)
shapiro.test(residuals(richness_anova))
#p-value = 0.001162 not normal
kruskal.test(alpha_div$Richness ~ alpha_div$Season)
#p-value = 1.507e-06 
dunnTest(alpha_div$Richness ~ alpha_div$Season, method = "bh")
#Comparison           Z     P.unadj       P.adj
#   FALL - SPRING  0.7776437 4.367791e-01 4.367791e-01
#   FALL - SUMMER  3.6170404 2.979908e-04 8.939725e-04
# SPRING - SUMMER  2.8393966 4.519894e-03 9.039788e-03
#   FALL - WINTER -1.7307558 8.349533e-02 1.001944e-01
# SPRING - WINTER -2.5083995 1.212794e-02 1.819192e-02
# SUMMER - WINTER -5.3477962 8.903165e-08 5.341899e-07


#Shannon
shannon_anova <- aov(alpha_div$Shannon ~ alpha_div$Season)
shapiro.test(residuals(shannon_anova))
#p-value = 0.0001579 not normal 
kruskal.test(alpha_div$Shannon ~ alpha_div$Season)
#p-value =  0.001528
dunnTest(alpha_div$Shannon ~ alpha_div$Season, method = "bh")
# Comparison          Z      P.unadj       P.adj
#   FALL - SPRING  3.4562204 0.0005478073 0.003286844
#   FALL - SUMMER  2.5881834 0.0096483582 0.019296716
# SPRING - SUMMER -0.8680369 0.3853741270 0.462448952
#   FALL - WINTER  0.7406186 0.4589246979 0.458924698
# SPRING - WINTER -2.7156017 0.0066155452 0.019846636
# SUMMER - WINTER -1.8475648 0.0646653259 0.096997989

#Chao1
Chao1_anova <- aov(alpha_div$Chao1 ~ alpha_div$Season)
shapiro.test(residuals(Chao1_anova))
#p-value = 0.00102 not normal
kruskal.test(alpha_div$Chao1 ~ alpha_div$Season)
#p-value = 1.712e-06
dunnTest(alpha_div$Chao1 ~ alpha_div$Season, method = "bh")
# Comparison          Z      P.unadj        P.adj
#   FALL - SPRING  0.7735997 4.391676e-01 4.391676e-01
#   FALL - SUMMER  3.6207655 2.937326e-04 8.811979e-04
# SPRING - SUMMER  2.8471658 4.411038e-03 8.822075e-03
#   FALL - WINTER -1.6987292 8.937022e-02 1.072443e-01
# SPRING - WINTER -2.4723289 1.342360e-02 2.013539e-02
# SUMMER - WINTER -5.3194947 1.040558e-07 6.243350e-07

#Pielou
Pielou_anova <- aov(alpha_div$Pielou ~ alpha_div$Season)
shapiro.test(residuals(Pielou_anova))
#p-value = 3.298e-05 not normal 
kruskal.test(alpha_div$Pielou ~ alpha_div$Season)
#p-value = 0.008254 significant
dunnTest(alpha_div$Pielou ~ alpha_div$Season, method = "bh")
#        Comparison          Z      P.unadj       P.adj
#   FALL - SPRING  3.3765839 0.0007339197 0.004403518
#   FALL - SUMMER  1.6723647 0.0944524307 0.141678646
# SPRING - SUMMER -1.7042193 0.0883401333 0.176680267
#   FALL - WINTER  1.1945462 0.2322643774 0.278717253
# SPRING - WINTER -2.1820377 0.0291067477 0.087320243
# SUMMER - WINTER -0.4778185 0.6327794051 0.632779405




##Alpha Diversity plots for Rare and Abundant
library(tidyverse)
library(vegan)
library(ggplot2)

#Extract Rare and Abundant ASVS 
#Total abundance across samples
total_abund <- colSums (Diatom_Table)
#Relative abundance
rel_abund <- total_abund / sum(total_abund)
#Rare and Abundant ASVs
rare_asvs <- Diatom_Table[, rel_abund < 0.001]
abundant_asvs <- Diatom_Table[, rel_abund >= 0.001]


#Rare
#Calculate Diversity Indices
#Richness
richness_rare<-specnumber(rare_asvs)
#shannon-weiner
shannon_rare <- diversity(rare_asvs, index = "shannon")
#Pielou's Evenness
pielou_rare <- shannon_rare / log(richness_rare)
#Chao1
chao1 <- estimateR(rare_asvs)["S.chao1",]

#Combine results
alpha_div_rare <- data.frame(SampleID = metadat$SampleID, Richness = richness_rare, Chao1 = chao1, Shannon = shannon_rare, Pielou = pielou_rare)%>% left_join(metadat, by = "SampleID")
#Transform to plot
alpha_long_rare <- alpha_div_rare %>% pivot_longer(cols = c("Richness", "Chao1", "Shannon", "Pielou"), names_to = "Metric", values_to = "Value")
#Set the order for Season
alpha_long_rare$Season <- factor(alpha_long_rare$Season, levels = c("SUMMER", "FALL", "WINTER","SPRING"))
#Combine season and phase factors
alpha_long_rare <- alpha_long_rare %>% mutate (SeasonPhase = interaction(Season, Phase, sep = "_"))
#Reorder neap and spring
alpha_long_rare$SeasonPhase <- factor(alpha_long_rare$SeasonPhase, levels = c("SUMMER_N", "SUMMER_S", "FALL_N", "FALL_S", "WINTER_N", "WINTER_S", "SPRING_N", "SPRING_S"))

#Plot results
ggplot(alpha_long_rare, aes(x=SeasonPhase, y = Value, fill = Season)) + geom_boxplot(alpha = 0.85) + facet_wrap(~ Metric, scales = "free_y") + scale_x_discrete(labels = function(x) gsub("_","\n", x)) + labs(x = "Season and Tidal Phase", y = "Diversity Value", fill = "Season")


#abundant
#Calculate Diversity Indices
#Richness
richness_abundant<-specnumber(abundant_asvs)
#shannon-weiner
shannon_abundant <- diversity(abundant_asvs, index = "shannon")
#Pielou's Evenness
pielou_abundant <- shannon_abundant / log(richness_abundant)
#Chao1
chao1 <- estimateR(abundant_asvs)["S.chao1",]

#Combine results
alpha_div_abundant <- data.frame(SampleID = metadat$SampleID, Richness = richness_abundant, Chao1 = chao1, Shannon = shannon_abundant, Pielou = pielou_abundant)%>% left_join(metadat, by = "SampleID")
#Transform to plot
alpha_long_abundant <- alpha_div_abundant %>% pivot_longer(cols = c("Richness", "Chao1", "Shannon", "Pielou"), names_to = "Metric", values_to = "Value")
#Set the order for Season
alpha_long_abundant$Season <- factor(alpha_long_abundant$Season, levels = c("SUMMER", "FALL", "WINTER","SPRING"))
#Combine season and phase factors
alpha_long_abundant <- alpha_long_abundant %>% mutate (SeasonPhase = interaction(Season, Phase, sep = "_"))
#Reorder neap and spring
alpha_long_abundant$SeasonPhase <- factor(alpha_long_abundant$SeasonPhase, levels = c("SUMMER_N", "SUMMER_S", "FALL_N", "FALL_S", "WINTER_N", "WINTER_S", "SPRING_N", "SPRING_S"))

#Plot results
ggplot(alpha_long_abundant, aes(x=SeasonPhase, y = Value, fill = Season)) + geom_boxplot(alpha = 0.85) + facet_wrap(~ Metric, scales = "free_y") + scale_x_discrete(labels = function(x) gsub("_","\n", x)) + labs(x = "Season and Tidal Phase", y = "Diversity Value", fill = "Season")




##Rare vs. Total ALpha Diversity Values
Diatom_Table$Dataset <- "Total"
Rare_asvs$Dataset <- "Rare"

library(dplyr)
library(tidyr)

#Combine into long format
alpha_compared_long <- bind_rows(Diatom_Table, Rare_asvs) %>% pivot_longer(cols = c(Richness,Chao1, Shannon, Pielou), names_to = "Metric", values_to = "Value")

library(ggplot2)
library(patchwork)

ggplot(alpha_compared_long, aes(x = Dataset, y = Value, fill = Dataset)) + geom_boxplot(alpha = 0.8) + facet_wrap(~ Metric, scales = "free_y") + labs(x = "", y = "Value")




##Table of ASV Abundance by Genus
library(tidyr)
library(dplyr)
library(tibble)

#Long format for merge
long <- Diatom_Table %>% rownames_to_column("SampleID") %>% pivot_longer(-SampleID, names_to = "ASVNumber", values_to = "Abundance")

#Merge with Taxaonomy data and Metadata
long_tax <- merge(long, ID, by = "ASVNumber")
tax_meta <- long_tax %>% left_join(metadata, by = "SampleID")

#Summarize by genus within each Phase nested in Season
Genus_Table <- tax_meta %>% group_by(Season, Phase, genus) %>% summarise(TotalAbundance = sum(Abundance, na.rm = T)) %>% group_by(Season, Phase) %>% mutate(RelAbund = TotalAbundance / sum(TotalAbundance)) %>% ungroup()

write.csv(Genus_Table, "Genus_Table.csv")



###Rare Taxa Relative Abundance
#Extracting Rare and Abundant ASVs

#Starting with the ASV table (Diatom_Table)

#Total abundance across samples
total_abund <- colSums (Diatom_Table)
#Relative abundance
rel_abund <- total_abund / sum(total_abund)
#Rare and Abundant ASVs
rare_asvs <- Diatom_Table[, rel_abund < 0.001]
abundant_asvs <- Diatom_Table[, rel_abund >= 0.001]

#Ensure ASV column in "ID" matches with the rare_asvs rows
colnames(rare_asvs) <- gsub("_", "", colnames(rare_asvs))

#Long format for merge
library(tidyr)
library(dplyr)
library(tibble)
rare_long <- rare_asvs %>% rownames_to_column("SampleID") %>% pivot_longer(-SampleID, names_to = "ASVNumber", values_to = "Abundance")

#Merge with Taxaonomy data and Metadata
rare_long_tax <- merge(rare_long, ID_Other, by = "ASVNumber")
rare_tax_meta <- rare_long_tax %>% left_join(metadata, by = "SampleID")

#Summarize by genus within each Phase nested in Season
rare_genus_phase <- rare_tax_meta %>% group_by(Season, Phase, genus) %>% summarise(TotalAbundance = sum(Abundance, na.rm = T)) %>% group_by(Season, Phase) %>% mutate(RelAbund = TotalAbundance / sum(TotalAbundance)) %>% ungroup()

#Plot stacked bar chart
library(ggplot2)
ggplot(rare_genus_phase,aes(x = Phase, y = RelAbund, fill = genus)) + geom_bar(stat = "identity") + facet_wrap (~ Season, scales = "free_x")

#A table may be more helpful. Save as csv for creation of Excel table.
write.csv(rare_genus_phase, "Rare_Genus_Abundance.csv")





##Comparing ASV percentages with less than 1% abundance and less than 0.1% abundance

#Extracting Rare (<0.1%) and 1% ASVs

#Total abundance across samples
total_abund <- colSums (Diatom_Table)
#Relative abundance
rel_abund <- total_abund / sum(total_abund)
#Rare and 1% ASVs
rare_asvs <- Diatom_Table[, rel_abund < 0.001]
onepercent_asvs <- Diatom_Table[, rel_abund < 0.01]
#Save
write.csv (rare_asvs, "Rare_ASV.csv")
write.csv (onepercent_asvs, "onepercent_ASV.csv")
#Compare each of these tables to Diatom_Table to find percent of ASVs in rare and 1 percent groups. 

#Percentages are much higher than expected. Try with water samples included.
total_abund_water <- colSums (silva_diatable_water)
#Relative abundance
rel_abund_water <- total_abund_water / sum(total_abund_water)
#Rare and 1% ASVs
rare_asvs_water <- silva_diatable_water[, rel_abund_water < 0.001]
onepercent_asvs_water <- silva_diatable_water[, rel_abund_water < 0.01]
#Save
write.csv (rare_asvs_water, "Rare_ASV_water.csv")
write.csv (onepercent_asvs_water, "onepercent_ASV_water.csv")


#Total abundance across all samples
asv_sums <- colSums (Diatom_Table)
total_reads <- sum(asv_sums)
asv_rel_abund <- (asv_sums / total_reads) * 100

#Extract low abundance ASVs and 1% ASVs
rare_abundance_asvs <- Diatom_Table[, asv_rel_abund < 0.001]
low_abundance_asvs <- Diatom_Table[, asv_rel_abund < 0.01]

#Save
write.csv(low_abundance_asvs, "Low_Abundance_ASVs.csv")
write.csv(rare_abundance_asvs, "Rare_Abundance_ASVs.csv")

#Total Abundances across all samples with water
asv_sums_water <- colSums (silva_diatable_water)
total_reads_water <- sum(asv_sums_water)
asv_rel_abund_water <- (asv_sums_water / total_reads_water) * 100

#Extract low abundance ASVs and 1% ASVs
rare_abundance_asvs_water <- silva_diatable_water[, asv_rel_abund_water < 0.001]
low_abundance_asvs_water <- silva_diatable_water[, asv_rel_abund_water < 0.01]

#Save
write.csv(low_abundance_asvs_water, "Low_Abundance_ASVs_water.csv")
write.csv(rare_abundance_asvs_water, "Rare_Abundance_ASVs_water.csv")







##Presence absence tables for Venn Diagrams
asv_pa <- Diatom_Table
asv_pa[asv_pa > 0] <- 1
write.csv(asv_pa, "asv_pa.csv")


#Add metadata in excel

#Add up asvs present
asv_by_season<- asv_pa %>% group_by(Season) %>% summarise(across(starts_with("ASV"), sum), .groups = "drop")
asv_by_season
#Do the same for lunar phase
asv_by_phase<- asv_pa %>% group_by(Phase) %>% summarise(across(starts_with("ASV"), sum), .groups = "drop")
asv_by_phase
#Tide
asv_by_tide<- asv_pa %>% group_by(Tide) %>% summarise(across(starts_with("ASV"), sum), .groups = "drop")
asv_by_tide


#Transform vertical
library(janitor)
asv_season <- t(asv_by_season)
asv_season <- row_to_names(asv_season, row_number = 1)
asv_phase <- t(asv_by_phase)
asv_phase <- row_to_names(asv_phase, row_number = 1)
asv_tide <- t(asv_by_tide)
asv_tide <- row_to_names(asv_tide, row_number = 1)

#Save these to work with in excel
write.csv(asv_season, "asv_season.csv")
write.csv(asv_phase, "asv_phase.csv")
write.csv(asv_tide, "asv_tide.csv")


#Water vs. sediment
asv_pa_water <- water_diatom_table
asv_pa_water[asv_pa_water > 0] <- 1


#add metadata in excel
write.csv(asv_pa_water, "asv_water.csv")

#Add up asvs present
asv_by_sedwat<- asv_water %>% group_by(SedWat) %>% summarise(across(starts_with("X.ASV"), sum), .groups = "drop")
asv_by_sedwat





#Community_Stability Metric (ASV Table without Transformations)
library(codyn)
library (tidyr)
library (dplyr)

#There are a few Zones within Sampling times that were collected at the exact same time.
#This does not work well with variance ratio and synchrony
#Average the time variable for each zone within the sampling events
#Denoted as E1 through E24 (Event)
df_avg <- Diatom_Time %>% group_by(Event) %>% summarise(mean_time = mean(Time, na.rm = TRUE), across(starts_with("ASV_"), ~mean(.x, na.rm = T)) %>% ungroup())
df_avg_sorted <- df_avg %>% arrange(mean_time)

#Long format
Long_Tab <- df_avg_sorted %>% pivot_longer (cols = starts_with("ASV_"), names_to = "ASV", values_to = "abundance")

#Community Stability
Comm_stability <- community_stability(Long_Tab, time.var = "mean_time", abundance.var = "abundance")
Comm_stability
# 2.424508

#Variance ratio 
var_ratio <- variance_ratio(Long_Tab, time.var = "mean_time", abundance.var = "abundance", species.var = "ASV", bootnumber = 1000)
#    lowerCI  upperCI nullmean       VR
#    0.4124164 2.035575 1.000444 4.308251

#Synchrony Loreau
sync <- synchrony(Long_Tab, time.var = "mean_time", abundance.var = "abundance", species.var = "ASV", metric ="Loreau")
sync
#0.1424745

#Synchrony Gross
sync_gross <- synchrony(Long_Tab, time.var = "mean_time", abundance.var = "abundance", species.var = "ASV", metric ="Gross")
sync_gross
#0.259726

#Total Turnover
turn_total<- turnover(Long_Tab, time.var = "mean_time", species.var = "ASV", abundance.var = "abundance", metric = "total")
turn_total
#Gain Turnover
turn_gain<- turnover(Long_Tab, time.var = "mean_time", species.var = "ASV", abundance.var = "abundance", metric = "appearance")
turn_gain
#Loss Turnover
turn_loss<- turnover(Long_Tab, time.var = "mean_time", species.var = "ASV", abundance.var = "abundance", metric = "disappearance")
turn_loss


#Combine total, gain, and loss
turn_all <- turn_total %>% left_join((turn_gain), by="mean_time") %>% left_join((turn_loss), by = "mean_time")
turn_all
write.csv(turn_all, "turnover_notransform.csv")

#Plot
library(ggplot2)
library(tidyr)

#In hours
turn_long <- pivot_longer(turn_all, cols = c("appearance", "disappearance", "total"), names_to = "type", values_to = "value")
turnover_hours<-ggplot(turn_long, aes(x = mean_time, y = value, color = type)) + geom_line(linewidth = 1) +labs(x = "Time in Hours", y = "Proportion of Species")

#By sampling event
turn_long_event <- pivot_longer(turnover_notransform, cols = c("appearance", "disappearance", "total"), names_to = "type", values_to = "value")
turnover_event<-ggplot(turn_long_event, aes(x = Event, y = value, color = type)) + geom_line(linewidth = 1) +labs(x = "Sampling Event", y = "Proportion of Species")




#Break into rare and common ASVs
#Add time data to rare and common ASV tables in Excel

df_avg_rare <- Rare_Diatom_Time %>% group_by(Event) %>% summarise(mean_time = mean(Time, na.rm = TRUE), across(starts_with("ASV_"), ~mean(.x, na.rm = T)) %>% ungroup())
df_avg_sorted_rare <- df_avg_rare %>% arrange(mean_time)
#Perform community stability metrics with the rare asvs.
#Long format
Long_Tab_rare <- df_avg_sorted_rare %>% pivot_longer (cols = starts_with("ASV_"), names_to = "ASV", values_to = "abundance")

#Community Stability
Comm_stability_rare <- community_stability(Long_Tab_rare, time.var = "mean_time", abundance.var = "abundance")
Comm_stability_rare
# 1.955587

#ariVance Ratio 
var_ratio_rare <- variance_ratio(Long_Tab_rare, time.var = "mean_time", abundance.var = "abundance", species.var = "ASV", bootnumber = 1000)
#    lowerCI  upperCI  nullmean       VR
#1 0.4290932 1.954778 0.9828788 28.86773

#Synchrony Loreau
sync_rare <- synchrony(Long_Tab_rare, time.var = "mean_time", abundance.var = "abundance", species.var = "ASV", metric ="Loreau")
sync_rare
#0.1235129

#Synchrony Gross
sync_rare_Gross <- synchrony(Long_Tab_rare, time.var = "mean_time", abundance.var = "abundance", species.var = "ASV", metric ="Gross")
sync_rare_Gross
#0.2375695


#Total Turnover
turn_total_rare<- turnover(Long_Tab_rare, time.var = "mean_time", species.var = "ASV", abundance.var = "abundance", metric = "total")
turn_total_rare
#Gain Turnover
turn_gain_rare<- turnover(Long_Tab_rare, time.var = "mean_time", species.var = "ASV", abundance.var = "abundance", metric = "appearance")
turn_gain_rare
#Loss Turnover
turn_loss_rare<- turnover(Long_Tab_rare, time.var = "mean_time", species.var = "ASV", abundance.var = "abundance", metric = "disappearance")
turn_loss_rare


#Combine total, gain, and loss
turn_all_rare <- turn_total_rare %>% left_join((turn_gain_rare), by="mean_time") %>% left_join((turn_loss_rare), by = "mean_time")
turn_all_rare
write.csv(turn_all_rare, "turnover.rare.notransform.csv")

#Plot
library(ggplot2)
library(tidyr)

#By sampling event
turn_long_rare <- pivot_longer(turnover.rare.notransform, cols = c("appearance", "disappearance", "total"), names_to = "type", values_to = "value")
turnover_event_rare<-ggplot(turn_long_rare, aes(x = Event, y = value, color = type)) + geom_line(linewidth = 1) +labs(x = "Sampling Event", y = "Proportion of Species")




#Perform community stability metrics with abundant ASVs.
df_avg_abundant <- Abundant_Diatom_Time %>% group_by(Event) %>% summarise(mean_time = mean(Time, na.rm = TRUE), across(starts_with("ASV_"), ~mean(.x, na.rm = T)) %>% ungroup())
df_avg_sorted_abundant <- df_avg_abundant %>% arrange(mean_time)
#Long format
Long_Tab_abundant <- df_avg_sorted_abundant %>% pivot_longer (cols = starts_with("ASV_"), names_to = "ASV", values_to = "abundance")

#Community Stability
Comm_stability_abundant <- community_stability(Long_Tab_abundant, time.var = "mean_time", abundance.var = "abundance")
Comm_stability_abundant
# 2.450375

#Variance Ratio 
var_ratio_abundant <- variance_ratio(Long_Tab_abundant, time.var = "mean_time", abundance.var = "abundance", species.var = "ASV", bootnumber = 1000)
var_ratio_abundant
#    lowerCI  upperCI nullmean      VR
#  1 0.4458857 2.105821 1.028032 3.33125

#Synchrony Loreau
sync_abundant <- synchrony(Long_Tab_abundant, time.var = "mean_time", abundance.var = "abundance", species.var = "ASV", metric ="Loreau")
sync_abundant
#0.1519999

#Synchrony Gross
sync_abundant_gross <- synchrony(Long_Tab_abundant, time.var = "mean_time", abundance.var = "abundance", species.var = "ASV", metric ="Gross")
sync_abundant_gross
#0.3865044

#Total Turnover
turn_total_abundant<- turnover(Long_Tab_abundant, time.var = "mean_time", species.var = "ASV", abundance.var = "abundance", metric = "total")
turn_total_abundant
#Gain Turnover
turn_gain_abundant<- turnover(Long_Tab_abundant, time.var = "mean_time", species.var = "ASV", abundance.var = "abundance", metric = "appearance")
turn_gain_abundant
#Loss Turnover
turn_loss_abundant<- turnover(Long_Tab_abundant, time.var = "mean_time", species.var = "ASV", abundance.var = "abundance", metric = "disappearance")
turn_loss_abundant


#Combine total, gain, and loss
turn_all_abundant <- turn_total_abundant %>% left_join((turn_gain_abundant), by="mean_time") %>% left_join((turn_loss_abundant), by = "mean_time")
turn_all_abundant
write.csv(turn_all_abundant, "turnover.abund.notransform.csv")

#Plot
library(ggplot2)
library(tidyr)

#By sampling event
turn_long_abundant <- pivot_longer(turnover.abund.notransform, cols = c("appearance", "disappearance", "total"), names_to = "type", values_to = "value")
turnover_event_abundant<-ggplot(turn_long_abundant, aes(x = Event, y = value, color = type)) + geom_line(linewidth = 1) +labs(x = "Sampling Event", y = "Proportion of Species")





##SIMPER analysis
library(vegan)
#Run SIMPER
Simper_season<-simper(Diatom_Table, MET$Season, permutations = 999, ord = TRUE )
summary(Simper_season)


Simper_season<-summary(Simper_season, ordered = T)

#Extract Season Combinations
Simper_Summer_Fall<-as.data.frame(Simper_season$Summer_Fall)
Simper_Summer_Winter<-as.data.frame(Simper_season$Summer_Winter)
Simper_Summer_Spring<-as.data.frame(Simper_season$Summer_Spring)
Simper_Fall_Winter<-as.data.frame(Simper_season$Fall_Winter)
Simper_Fall_Spring<-as.data.frame(Simper_season$Fall_Spring)
Simper_Winter_Spring<-as.data.frame(Simper_season$Winter_Spring)

#Compute percent contribution for each ASV
Simper_Summer_Fall$percent <- 100 * Simper_Summer_Fall$average / sum(Simper_Summer_Fall$average)
Simper_Summer_Winter$percent <- 100 * Simper_Summer_Winter$average / sum(Simper_Summer_Winter$average)
Simper_Summer_Spring$percent <- 100 * Simper_Summer_Spring$average / sum(Simper_Summer_Spring$average)
Simper_Fall_Winter$percent <- 100 * Simper_Fall_Winter$average / sum(Simper_Fall_Winter$average)
Simper_Fall_Spring$percent <- 100 * Simper_Fall_Spring$average / sum(Simper_Fall_Spring$average)
Simper_Winter_Spring$percent <- 100 * Simper_Winter_Spring$average / sum(Simper_Winter_Spring$average)

#Order by descending contribution
Simper_Summer_Fall <- Simper_Summer_Fall[order(Simper_Summer_Fall$percent, decreasing = T), ]
Simper_Summer_Winter <- Simper_Summer_Winter[order(Simper_Summer_Winter$percent, decreasing = T), ]
Simper_Summer_Spring <- Simper_Summer_Spring[order(Simper_Summer_Spring$percent, decreasing = T), ]
Simper_Fall_Winter <- Simper_Fall_Winter[order(Simper_Fall_Winter$percent, decreasing = T), ]
Simper_Fall_Spring <- Simper_Fall_Spring[order(Simper_Fall_Spring$percent, decreasing = T), ]
Simper_Winter_Spring <- Simper_Winter_Spring[order(Simper_Winter_Spring$percent, decreasing = T), ]

#Compute cumulative sum
Simper_Summer_Fall$cumsum_percent <- cumsum(Simper_Summer_Fall$percent)
Simper_Summer_Winter$cumsum_percent <- cumsum(Simper_Summer_Winter$percent)
Simper_Summer_Spring $cumsum_percent <- cumsum(Simper_Summer_Spring$percent)
Simper_Fall_Winter$cumsum_percent <- cumsum(Simper_Fall_Winter$percent)
Simper_Fall_Spring$cumsum_percent <- cumsum(Simper_Fall_Spring$percent)
Simper_Winter_Spring$cumsum_percent <- cumsum(Simper_Winter_Spring$percent)

#Filter top contributors (70%)
top70_Summer_fall <- subset(Simper_Summer_Fall, cumsum_percent <=70)
top70_Summer_Winter <- subset(Simper_Summer_Winter, cumsum_percent <=70)
top70_Summer_Spring <- subset(Simper_Summer_Spring, cumsum_percent <=70)
top70_Fall_Winter <- subset(Simper_Fall_Winter, cumsum_percent <=70)
top70_Fall_Spring <- subset(Simper_Fall_Spring, cumsum_percent <=70)
top70_Winter_Spring <- subset(Simper_Winter_Spring, cumsum_percent <=70)


#Save as csvs for organization into chart
write.csv(top70_Summer_fall, "top70_Summer_fall.csv")
write.csv(top70_Summer_Winter, "top70_Summer_Winter.csv")
write.csv(top70_Summer_Spring, "top70_Summer_Spring.csv")
write.csv(top70_Fall_Winter, "top70_Fall_Winter.csv")
write.csv(top70_Fall_Spring, "top70_Fall_Spring.csv")
write.csv(top70_Winter_Spring, "top70_Winter_Spring.csv")