# BiocManager::install("apeglm")
# BiocManager::install("DESeq2")
# BiocManager::install("DEGreport")
# browseVignettes(package = "DESeq2")
# install.packages(c("ashr", "ggbeeswarm", "pheatmap"))
# devtools::install_github("teunbrand/ggh4x")

library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(DEGreport)
library(vegan)

# set colors for plots 
zero <- "#18f27f" #green
ten <- "#099c4e"
thirty <- "#022412"

srur <- "#18f27f"
scan <- "#FF00FF"

BR.zero <- "#18f27f"
BR.ten <- "#088c47"
BR.thirty <- "#033018"
SC.zero <- "#FF00FF"
SC.ten <-  "#9d009d"
SC.thirty <- "#3b003b"

# for sig cateogory plots
both <- "purple3"
ten.only <- "darkorange"
thirty.only <- "skyblue"

  
setwd("/Users/jennaekwealor/Documents/dissertation_repositories/UV_syntrichia_acute")

directory <- ("data/htseq-count_10142020/")

# load sample table
sampleTable <- read.csv("sample_key.csv", stringsAsFactors = TRUE)

# filter for just the species i want to look at here
sampleTable_SC <- sampleTable %>% filter(species == "SC")
sampleTable_BR <- sampleTable %>% filter(species == "BR")

# load annotations
func_annot <- read.csv(file = "2021_func_annot_GO_terms_SC.csv")
func_annot <- subset(func_annot, select = -X)

#### UV effect in SC ####

dds_SC <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_SC,
                                      directory = directory,
                                      design = ~ time) # does the UV  have an effect
                                      

# relevel treatments so that one is the 'reference' in differential abundance analyses
dds_SC$time <- relevel(dds_SC$time, ref = "zero") # level treatments so one is marked as 'zero' as the baseline in terms of diff abundance

# run the DE analysis
dds_SC <- DESeq(dds_SC)

resultsNames(dds_SC)
res_SC <- results(dds_SC)
summary(res_SC)
# 6868 total genes

# filter out rows with less than 10 reads
keep <- rowSums(counts(dds_SC)) >= 10 #  rows with at least 10 reads
dds_SC <- dds_SC[keep,] # keep rows from above


#### plot dispersion ####
plotDispEsts(dds_SC)

#### plot PCA ####
vsd_SC <- vst(dds_SC, blind=T)
plotPCA(vsd_SC, intgroup="time") + theme_minimal()
# ntop = number of top genes to use for principal components, selected by highest row variance

# use returnData = TRUE should the function only return the data.frame of PC1 and PC2 with intgroup covariates for custom plotting
pca_df_SC <- plotPCA(vsd_SC, intgroup=c("time"), returnData = TRUE) 


p_pca_all_SC <- pca_df_SC %>% ggplot(.,aes(x=PC1,y=PC2)) +
  geom_point(aes(color = time), size = 3) + 
  ggforce::geom_mark_ellipse(aes(color = time)) +  
  scale_color_manual(name = "UV Treatment", labels = c(expression(T[0]), expression(T[10]), expression(T[30])),values = c(SC.zero, SC.ten, SC.thirty)) +
  theme_minimal() +
  theme(legend.text.align = 0,
        text = element_text(size = 14)) +
  xlab("PC1 (45%)") + ylab("PC2 (20%)") +
  scale_x_continuous(limits = c(-25,20)) +
  scale_y_continuous(limits = c(-10,15)) 

p_pca_all_SC

pdf("plots/p_pca_all_SC.pdf") 
p_pca_all_SC
dev.off()


#### permanova to test if groups are different? ####


#### ten vs. zero ####

resultsNames(dds_SC) # gives a list of "name" comparisons you can use as shortcut instead of "contrasts"
results(dds_SC, name = "time_ten_vs_zero", alpha = 0.05)  -> time_ten_vs_zero_SC
# make df of UV genes with pvalue of less than 0.05 and shrunken LFC of at least 2
time_ten_vs_zero_SC %>%
  as.data.frame() %>%
  mutate(transcript_id = rownames(.) ) %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  filter(padj < 0.05) -> time_ten_vs_zero_df_SC

# join annotations to results file
inner_join(x = func_annot, y = time_ten_vs_zero_df_SC, by = "transcript_id") %>%
  distinct(transcript_id,  .keep_all = T) -> time_ten_vs_zero_ann_SC

write.csv(time_ten_vs_zero_ann_SC, "tables/SC/time_ten_vs_zero_ann.csv", row.names=FALSE)

length(time_ten_vs_zero_ann_SC$transcript_id)
# 10 genes
time_ten_vs_zero_ann_SC %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
  arrange(abs_log2FoldChange) %>%
  tail(n = 10) -> time_ten_vs_zero_ann_strict_SC
write.csv(time_ten_vs_zero_ann_strict_SC, "tables/SC/time_ten_vs_zero_ann_strict.csv", row.names=FALSE)


#### Plot results  ####

# MAplot 
plotMA(time_ten_vs_zero_SC, ylim=c(-4,4), cex=0.5, alpha = 0.05)
abline(h=c(-1,1), col="yellow", lwd=2) 
while (!is.null(dev.list()))  dev.off()

pdf("plots/maplot_time_ten_vs_zero_SC.pdf") 
plotMA(time_ten_vs_zero_SC, ylim=c(-4,4), cex=0.5, alpha = 0.05,
       xlab = "Mean of Normalized Counts",
       ylab = "Log2 Fold Change")
abline(h=c(-1,1), col="yellow", lwd=2) 
while (!is.null(dev.list()))  dev.off()


#### thirty vs. zero ####

resultsNames(dds_SC) # gives a list of "name" comparisons you can use as shortcut instead of "contrasts"
results(dds_SC, name = "time_thirty_vs_zero", alpha = 0.05)  -> time_thirty_vs_zero_SC
# make df of UV genes with svalue of less than 0.05 and shrunken LFC of at least 2
time_thirty_vs_zero_SC %>%
  as.data.frame() %>%
  mutate(transcript_id = rownames(.) ) %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  filter(padj < 0.05) -> time_thirty_vs_zero_df_SC

# join annotations to results file
inner_join(x = func_annot, y = time_thirty_vs_zero_df_SC, by = "transcript_id") %>%
  distinct(transcript_id,  .keep_all = T) -> time_thirty_vs_zero_ann_SC

write.csv(time_thirty_vs_zero_ann_SC, "tables/SC/time_thirty_vs_zero_ann.csv", row.names=FALSE)

length(time_thirty_vs_zero_ann_SC$transcript_id)
# 126 genes
time_thirty_vs_zero_ann_SC %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
  arrange(abs_log2FoldChange) %>%
  tail(n = 10) -> time_thirty_vs_zero_ann_strict_SC
write.csv(time_thirty_vs_zero_ann_strict_SC, "tables/SC/time_thirty_vs_zero_ann_strict.csv", row.names=FALSE)


#### Plot results  ####

# MAplot 
plotMA(time_thirty_vs_zero_SC, ylim=c(-4,4), cex=0.5, alpha = 0.05)
abline(h=c(-1,1), col="yellow", lwd=2) 
while (!is.null(dev.list()))  dev.off()

pdf("plots/maplot_time_thirty_vs_zero_SC.pdf") 
plotMA(time_thirty_vs_zero_SC, ylim=c(-4,4), cex=0.5, alpha = 0.05,
       xlab = "Mean of Normalized Counts",
       ylab = "Log2 Fold Change")
abline(h=c(-1,1), col="yellow", lwd=2) 
while (!is.null(dev.list()))  dev.off()



#### UV effect in BR ####

dds_BR <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_BR,
                                     directory = directory,
                                     design = ~ time) # does the UV  have an effect


# relevel treatments so that one is the 'reference' in differential abundance analyses
dds_BR$time <- relevel(dds_BR$time, ref = "zero") # level treatments so one is marked as 'zero' as the baseline in terms of diff abundance

# run the DE analysis
dds_BR <- DESeq(dds_BR)

resultsNames(dds_BR)
res_BR <- results(dds_BR)
summary(res_BR)
# 6851 total genes

# filter out rows with less than 10 reads
keep <- rowSums(counts(dds_BR)) >= 10 #  rows with at least 10 reads
dds_BR <- dds_BR[keep,] # keep rows from above


# plot dispersion
plotDispEsts(dds_BR)

#### plot PCA ####
vsd_BR <- vst(dds_BR, blind=T)
plotPCA(vsd_BR, intgroup="time") + theme_minimal()
# ntop = number of top genes to use for principal components, selected by highest row variance

# use returnData = TRUE should the function only return the data.frame of PC1 and PC2 with intgroup covariates for custom plotting
pca_df_BR <- plotPCA(vsd_BR, intgroup=c("time"), returnData = TRUE) 


p_pca_all_BR <- pca_df_BR %>% ggplot(.,aes(x=PC1,y=PC2)) +
  geom_point(aes(color = time), size = 3) + 
  ggforce::geom_mark_ellipse(aes(color = time)) +  
  scale_color_manual(name = "UV Treatment", labels = c(expression(T[0]), expression(T[10]), expression(T[30])),values = c(BR.zero, BR.ten, BR.thirty)) +
  theme_minimal() +
  theme(legend.text.align = 0,
        text = element_text(size = 14)) +
  xlab("PC1 (43%)") + ylab("PC2 (23%)") +
  scale_x_continuous(limits = c(-25,20)) +
  scale_y_continuous(limits = c(-10,15)) 

p_pca_all_BR

pdf("plots/p_pca_all_BR.pdf") 
p_pca_all_BR
dev.off()



#### ten vs. zero ####

resultsNames(dds_BR) # gives a list of "name" comparisons you can use as shortcut instead of "contrasts"
results(dds_BR, name = "time_ten_vs_zero")  -> time_ten_vs_zero_BR
# make df of UV genes with svalue of less than 0.05 and shrunken LFC of at least 2
time_ten_vs_zero_BR %>%
  as.data.frame() %>%
  mutate(transcript_id = rownames(.) ) %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  filter(padj < 0.05) -> time_ten_vs_zero_df_BR

# join annotations to results file
inner_join(x = func_annot, y = time_ten_vs_zero_df_BR, by = "transcript_id") %>%
  distinct(transcript_id,  .keep_all = T) -> time_ten_vs_zero_ann_BR

write.csv(time_ten_vs_zero_ann_BR, "tables/BR/time_ten_vs_zero_ann.csv", row.names=FALSE)

length(time_ten_vs_zero_ann_BR$transcript_id)
# 18 genes
time_ten_vs_zero_ann_BR %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
  arrange(abs_log2FoldChange) %>%
  tail(n = 10) -> time_ten_vs_zero_ann_strict_BR
write.csv(time_ten_vs_zero_ann_strict_BR, "tables/BR/time_ten_vs_zero_ann_strict.csv", row.names=FALSE)


#### Plot results  ####

# MAplot 
plotMA(time_ten_vs_zero_BR, ylim=c(-4,4), cex=0.5, alpha = 0.05)
abline(h=c(-1,1), col="yellow", lwd=2) 
while (!is.null(dev.list()))  dev.off()

pdf("plots/maplot_time_ten_vs_zero_BR.pdf") 
plotMA(time_ten_vs_zero_BR, ylim=c(-4,4), cex=0.5, alpha = 0.05,
       xlab = "Mean of Normalized Counts",
       ylab = "Log2 Fold Change")
abline(h=c(-1,1), col="yellow", lwd=2) 
while (!is.null(dev.list()))  dev.off()




#### thirty vs. zero ####

resultsNames(dds_BR) # gives a list of "name" comparisons you can use as shortcut instead of "contrasts"
results(dds_BR, name = "time_thirty_vs_zero")  -> time_thirty_vs_zero_BR
# make df of UV genes with svalue of less than 0.05 and shrunken LFC of at least 2
time_thirty_vs_zero_BR %>%
  as.data.frame() %>%
  mutate(transcript_id = rownames(.) ) %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  filter(padj < 0.05) -> time_thirty_vs_zero_df_BR

# join annotations to results file
inner_join(x = func_annot, y = time_thirty_vs_zero_df_BR, by = "transcript_id") %>%
  distinct(transcript_id,  .keep_all = T) -> time_thirty_vs_zero_ann_BR

write.csv(time_thirty_vs_zero_ann_BR, "tables/BR/time_thirty_vs_zero_ann.csv", row.names=FALSE)

length(time_thirty_vs_zero_ann_BR$transcript_id)
# 38 genes
time_thirty_vs_zero_ann_BR %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
  arrange(abs_log2FoldChange) %>%
  tail(n = 10) -> time_thirty_vs_zero_ann_strict_BR
write.csv(time_thirty_vs_zero_ann_strict_BR, "tables/BR/time_thirty_vs_zero_ann_strict.csv", row.names=FALSE)


#### Plot results  ####

# MAplot 
plotMA(time_thirty_vs_zero_BR, ylim=c(-4,4), cex=0.5, alpha = 0.05)
abline(h=c(-1,1), col="yellow", lwd=2) 
while (!is.null(dev.list()))  dev.off()

pdf("plots/maplot_time_thirty_vs_zero_BR.pdf") 
plotMA(time_thirty_vs_zero_BR, ylim=c(-4,4), cex=0.5, alpha = 0.05,
       xlab = "Mean of Normalized Counts",
       ylab = "Log2 Fold Change")
abline(h=c(-1,1), col="yellow", lwd=2) 
while (!is.null(dev.list()))  dev.off()
















#### interaction: do the species differ in their response to UV? #####
# design formula that models the species difference at time 0, the difference over time, and any species-specific differences over time (the interaction term species:time).
# convert HTSeqCount to DESeqDataSet
dds_int <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ species + # does the species  have an effect in the reference time 0
                                    time + # does the time (UV treatment) have an effect in the reference species?
                                    species:time) # does the ltime (UV treatment)  have different effects based on species

# relevel treatments so that one is the 'reference' in differential abundance analyses
dds_int$time <- relevel(dds_int$time, ref = "zero") # level treatments so one is marked as 'zero' as the baseline in terms of diff abundance
dds_int$species <- relevel(dds_int$species, ref = "BR")

# run the DE analysis
dds_int <- DESeq(dds_int)
resultsNames(dds_int)
res_int <- results(dds_int)
summary(res_int)
# 6883 total genes

# filter out rows with less than 10 reads
keep <- rowSums(counts(dds_int)) >= 10 #  rows with at least 10 reads
dds_int <- dds_int[keep,] # keep rows from above


# plot dispersion
plotDispEsts(dds_int)

# plot PCA
vsd_int <- vst(dds_int, blind=T)
plotPCA(vsd_int, intgroup=c("species", "time")) + theme_minimal()
# ntop = number of top genes to use for principal components, selected by highest row variance

# use returnData = TRUE should the function only return the data.frame of PC1 and PC2 with intgroup covariates for custom plotting
pca_df <- plotPCA(vsd_int, intgroup=c("species", "time"), returnData = TRUE) 


p_pca_all <- pca_df %>% ggplot(.,aes(x=PC1,y=PC2)) +
  geom_point(aes(color = interaction(species, time)), size = 2) + 
  scale_color_manual(name = "Species", values = c(BR.zero=BR.zero, BR.ten=BR.ten, BR.thirty=BR.thirty, SC.zero=SC.zero, SC.ten=SC.ten, SC.thirty=SC.thirty, BR=srur, SC=scan)) +
  theme_minimal() +
  theme(legend.text.align = 0,
        text = element_text(size = 14)) +
  xlab("PC1 (93%)") + ylab("PC2 (3%)")  +
  stat_ellipse(aes(x=PC1,y=PC2,group=species, color = species),type = "norm", level = 0.95) 

p_pca_all

pdf("plots/p_pca_all.pdf") 
p_pca_all
dev.off()


#### LRT to compare species-specific effects ####
# Genes with small p values from this test are those which at one or more time points after time 0 showed a species-specific effect. 
# Note therefore that this will not give small p values to genes that moved up or down over time in the same way in both species
# dds_int = full model

dds_intLRT <- DESeq(dds_int, test="LRT", reduced = ~ species + time)
resLRT <- results(dds_intLRT)
resLRT$symbol <- mcols(dds_intLRT)$symbol
head(resLRT[order(resLRT$padj),], 4)
summary(resLRT)
# 6859 total genes

# list of genes diff in the two species
results(dds_intLRT) -> LRT

# relevel treatments so that one is the 'reference' in differential abundance analyses
dds_intLRT$time <- relevel(dds_int$time, ref = "zero") # level treatments so one is marked as 'zero' as the baseline in terms of diff abundance
dds_intLRT$species <- relevel(dds_int$species, ref = "BR")

# make df of UV genes with svalue of less than 0.05 
LRT %>%
  as.data.frame() %>%
  mutate(transcript_id = rownames(.) ) %>%
  filter(padj < 0.05) -> LRT_df

# join annotations to results file
inner_join(x = func_annot, y = LRT_df, by = "transcript_id") %>%
  distinct(transcript_id,  .keep_all = T) %>% select(-c(log2FoldChange, lfcSE))-> LRT_ann

# baseMean: mean of normalized counts for all samples
# stat: the difference in deviance between the reduced model and the full model
# pvalue: the stat value is compared to a chi-squared distribution to generate a pvalue
# padj: BH adjusted p-values

write.csv(LRT_ann, "tables/LRT_ann.csv", row.names=FALSE)

length(LRT_ann$transcript_id)
# 69 genes
LRT_ann %>%
  arrange(padj) %>%
  head(n = 10) -> LRT_ann_strict
write.csv(LRT_ann_strict, "tables/LRT_ann_strict.csv", row.names=FALSE)

#### cluster analysis for LRT sig genes ####
# identify groups of genes that share a pattern of expression change across the sample groups (levels). 
# To do this we will be using a clustering tool called degPatterns from the ‘DEGreport’ package. 
# The degPatterns tool uses a hierarchical clustering approach based on pair-wise correlations between genes, then cuts the hierarchical tree to generate groups of genes with similar expression profiles.
# The tool cuts the tree in a way to optimize the diversity of the clusters, such that the variability inter-cluster > the variability intra-cluster.

# vsd transform counts from LRT
vsd_intLRT <- vst(dds_intLRT, blind=T)
vsd_intLRT_mat <- assay(vsd_intLRT)


# Obtain rlog values for those significant genes
cluster_rlog <- vsd_intLRT_mat[LRT_df$transcript_id, ]

# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
# metadata: the metadata dataframe that corresponds to samples
# time: character column name in metadata that will be used as variable that changes
# col: character column name in metadata to separate samples
# metadata file needs samplenames as rownames
sampleTable_meta <- sampleTable[,-1]
rownames(sampleTable_meta) <- sampleTable[,1]
sampleTable_meta$time <- factor(sampleTable_meta$time, levels=c("zero", "ten", "thirty"))

clusters <- DEGreport::degPatterns(cluster_rlog, metadata = sampleTable_meta, time = "time", col="species", summarize = "condition", minc = 5, reduce = F)

# write csvs of clusters, transcript_ids and annotations
clusters_genes <- clusters[["normalized"]] %>% select(genes, cluster)
colnames(clusters_genes) <- c("transcript_id", "cluster")
# there are two of each row, one for each species. don't need; remove
clusters_genes <- unique(clusters_genes)
clusters_ann <- left_join(clusters_genes, func_annot, by = "transcript_id")
# sort by cluster
clusters_ann <- clusters_ann[order(clusters_ann$cluster),] 
#gsub to replace cluster names
clusters_ann$cluster <- clusters_ann$cluster %>% gsub("1", "I", .) %>% gsub("2", "II", .) %>% gsub("3", "III", .) %>% gsub("4", "IV", .) %>% gsub("5", "V", .) %>% gsub("6", "VI", .)

# write clusters to csv
write.csv(clusters_ann, "tables/clusters_ann.csv", row.names=FALSE)


# extract cluster data for custom figures
cluster_labs <- c("Cluster I: 13 transcripts", "Cluster II: 13 transcripts", "Cluster III: 6 transcripts", "Cluster IV: 19 transcripts", "Cluster V: 10 transcripts", "Cluster VI: 8 transcripts")
names(cluster_labs) <- c("1", "2", "3", "4", "5", "6")


ggplot(clusters[["normalized"]],
       aes(time, value, species), guide = 'none') +
  geom_boxplot(aes(color = interaction(species, time)), outlier.shape=NA, alpha = 0.2, position = position_dodge(width = 0.9)) +
  facet_wrap(~cluster, labeller = labeller(cluster = cluster_labs)) +
  scale_color_manual(name = "Species", values = c(BR.zero=BR.zero, BR.ten=BR.ten, BR.thirty=BR.thirty, SC.zero=SC.zero, SC.ten=SC.ten, SC.thirty=SC.thirty, BR=srur, SC=scan)) +
  scale_fill_manual(name = "Species", values = c(BR.zero=BR.zero, BR.ten=BR.ten, BR.thirty=BR.thirty, SC.zero=SC.zero, SC.ten=SC.ten, SC.thirty=SC.thirty, BR=srur, SC=scan)) +
  stat_summary(fun = median, geom = 'line', size=0.75, aes(color = species, group = interaction(cluster, species)), position = position_dodge(width = 0.9)) +
  geom_point(aes(fill=interaction(species, time)), size=2, shape=21, stroke=0, position = position_jitterdodge(dodge.width = 0.9)) +
  # geom_line(data = clusters[["normalized"]],aes(color=species, group = interaction(species, genes)), alpha = 0.1) +
  theme_light() +
  theme(strip.background =element_rect(fill="black")) +
  ylab("Transcript abundance Z-score") +
  xlab("Minutes of UV Radiation Exposure") +
  scale_x_discrete(labels=c("zero" = "0", "ten" = "10",
                            "thirty" = "30")) -> cluster_plot

cluster_plot 


ggsave("plots/clusters.pdf", cluster_plot, width = 4, height = 3, units = "in", scale = 2.5)




#### individual count plots of LRT sig genes (different in two species) ####
  

# We can plot the counts for the groups over time using ggplot2, 
# for the gene with the smallest adjusted p value, 
# testing for condition-dependent time profile and accounting for differences at time 0 
# Keep in mind that the interaction terms are the difference 
# between the two groups at a given time after accounting for the difference at time 0

# will need to get the list of genes then plotCounts(), looping over the the list
# creating a data frame for ggplot to facet plot them all at once 
num_transcripts <- length(LRT_ann_strict$transcript_id)
count_data <- NULL
final_counts_LRT_ann_strict <- NULL

for(g in 1:num_transcripts) {
  
  one_transcript <- LRT_ann_strict$transcript_id[g]
  
  # returnData = true returns data frame so it can be plotted with ggplot2, etc.
  count_data <- plotCounts(dds_intLRT, gene=one_transcript, intgroup=c("time","species"), 
                           returnData=TRUE)
  count_data$transcript_id <- rep(one_transcript)
  
  # convert rownames to column and group by treatments
  count_data <- count_data %>% tibble::rownames_to_column(., var = "sample") 
  final_counts_LRT_ann_strict <- rbind(final_counts_LRT_ann_strict, count_data)
}
final_counts_LRT_ann_strict

# log transform the counts
final_counts_LRT_ann_strict <- dplyr::mutate(final_counts_LRT_ann_strict, log2_count = log2(count))

# # grep to replace species names
# final_counts_LRT_ann_strict$species <- final_counts_LRT_ann_strict$species %>% gsub("BR", "S. ruralis", .) %>% gsub("SC", "S. caninervis", .)

# add gene names
final_counts_LRT_ann_strict <- left_join(final_counts_LRT_ann_strict, func_annot, by = "transcript_id", keep.both)
 
final_counts_LRT_ann_strict$transcript_id <- factor(final_counts_LRT_ann_strict$transcript_id, levels = c("Sc_g08728", "Sc_g00568", "Sc_g00609", "Sc_g13403", "Sc_g02848", "Sc_g10759", "Sc_g14586", "Sc_g14304", "Sc_g12657", "Sc_g05211"))




ggline <- ggplot(final_counts_LRT_ann_strict,
       aes(x = time, y = log2_count, color = species, fill = interaction(species, time), group = species)) + 
      facet_wrap("transcript_id", ncol = 5) +
  geom_point(aes(fill=interaction(species, time)), size=2, shape=21, stroke=0) + 
  stat_summary(fun=mean, geom="line")  +
  scale_color_manual(name = "Species", values = c(BR=srur, SC=scan), labels = c(expression(italic("S. ruralis")), expression(italic("S. caninervis")))) +
  scale_fill_manual(values = c(BR.zero=BR.zero, BR.ten=BR.ten, BR.thirty=BR.thirty, SC.zero=SC.zero, SC.ten=SC.ten, SC.thirty=SC.thirty)) +
  theme_light() +
  theme(text = element_text(size = 10),
        legend.position="bottom",
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.y=element_text(size = 10),
        strip.background = element_rect(fill = "black")) +
        scale_x_discrete(labels=c("zero" = "0", "ten" = "10",
                              "thirty" = "30")) +
  labs(x="\nMinutes of UV Radiation Exposure", y = "Log2 Transformed Normalized Transcript Count")

ggline 

ggsave("plots/ggline.pdf", ggline, width = 3, height = 2, units = "in", scale = 3)





#### individual LFC plots of SC sig genes that overlap and are distinct in 10 and 30 minutes (different in two species) ####

# first get lists of uniqueness and overlapping in the sig at 10 and 30 in SC
overlap_SC <- time_ten_vs_zero_df_SC %>% filter(transcript_id %in% time_thirty_vs_zero_df_SC$transcript_id) %>% select(transcript_id)
ten_only_SC <- time_ten_vs_zero_df_SC %>% filter(!transcript_id %in% overlap_SC$transcript_id) %>% select(transcript_id)
thirty_only_SC <- time_thirty_vs_zero_df_SC %>% filter(!transcript_id %in% overlap_SC$transcript_id) %>% select(transcript_id)

# take the transcripts that are significant from ten and calculate their absolute value LFC at ten
time_ten_vs_zero_ann_SC_10_absLFC <- time_ten_vs_zero_ann_SC %>% mutate(absLFC = abs(log2FoldChange)) %>% select(transcript_id, absLFC, lfcSE)
time_ten_vs_zero_ann_SC_10_absLFC$time <- rep("ten")

# take these same transcripts and make an identical data frame that has LFC at 0 set to 0 (normalization so they can all be plotted together) 
time_zero_ann_SC_10_absLFC <- time_ten_vs_zero_ann_SC_10_absLFC %>% mutate(time = rep("zero"), absLFC = rep(0), lfcSE = rep(0))

# take these same transcripts and pull their LFC at time 30
# first convert results to df and make rownames into a column
SC_thirty_10_res <- (as.data.frame(results(dds_SC, name = "time_thirty_vs_zero")))
SC_thirty_10_res$transcript_id <- rownames(SC_thirty_10_res)

time_thirty_vs_zero_ann_SC_10_absLFC <- subset(SC_thirty_10_res, transcript_id %in% time_zero_ann_SC_10_absLFC$transcript_id) %>% 
                                                mutate(absLFC = abs(log2FoldChange)) %>% 
                                                select(transcript_id, absLFC, lfcSE) %>%
                                                mutate(time = rep("thirty"))
  
# remove rownames 
rownames(time_thirty_vs_zero_ann_SC_10_absLFC) <- NULL

# put these together for plotting
time_all_ann_SC_10_absLFC <- rbind(time_zero_ann_SC_10_absLFC, time_ten_vs_zero_ann_SC_10_absLFC, time_thirty_vs_zero_ann_SC_10_absLFC) %>% mutate(species = rep("SC"))
# order the times as a factor
time_all_ann_SC_10_absLFC$time <- factor(time_all_ann_SC_10_absLFC$time, levels = c("zero", "ten", "thirty"))
# convert abs(LFC) to numeric

# filter only those that are uniquely sig in t = 10! then repeat for those sig at 30
# overlap between 10 and 30 minutes
# first need to turn rownames into a column for each
time_ten_vs_zero_df_SC$transcript_id <- rownames(time_ten_vs_zero_df_SC)
# remove rownames 
rownames(time_ten_vs_zero_df_SC) <- NULL

time_thirty_vs_zero_df_SC$transcript_id <- rownames(time_thirty_vs_zero_df_SC)
# remove rownames 
rownames(time_thirty_vs_zero_df_SC) <- NULL

# filter time_all_ann_SC_10_absLFC to select for unique only
time_all_ann_SC_10_only_absLFC <- time_all_ann_SC_10_absLFC %>% filter(transcript_id %in% ten_only_SC$transcript_id) %>% mutate(category = rep("ten_only"))

# take the overlap genes and make into an absLFC table
time_all_ann_SC_overlap_absLFC <- time_all_ann_SC_10_absLFC %>% filter(transcript_id %in% overlap_SC$transcript_id) %>% mutate(category = rep("both"))

# take the overlap genes and make into an absLFC table, including LFCs from thirty and from ten as above
###
###
###
###
# take the transcripts that are significant from thirty and calculate their absolute value LFC at ten
time_thirty_vs_zero_ann_SC_30_absLFC <- time_thirty_vs_zero_ann_SC %>% mutate(absLFC = abs(log2FoldChange)) %>% select(transcript_id, absLFC, lfcSE)
time_thirty_vs_zero_ann_SC_30_absLFC$time <- rep("thirty")

# take these same transcripts and make an identical data frame that has LFC at 0 set to 0 (normalization so they can all be plotted together) 
time_zero_ann_SC_30_absLFC <- time_thirty_vs_zero_ann_SC_30_absLFC %>% mutate(time = rep("zero"), absLFC = rep(0), lfcSE = rep(0))

# take these same transcripts and pull their LFC at time 10
# first convert results to df and make rownames into a column
SC_ten_30_res <- (as.data.frame(results(dds_SC, name = "time_ten_vs_zero")))
SC_ten_30_res$transcript_id <- rownames(SC_ten_30_res)

time_ten_vs_zero_ann_SC_30_absLFC <- subset(SC_ten_30_res, transcript_id %in% time_zero_ann_SC_30_absLFC$transcript_id) %>% 
  mutate(absLFC = abs(log2FoldChange)) %>% 
  select(transcript_id, absLFC, lfcSE) %>%
  mutate(time = rep("ten"))

# remove rownames 
rownames(time_ten_vs_zero_ann_SC_30_absLFC) <- NULL

# put these together for plotting
time_all_ann_SC_30_absLFC <- rbind(time_zero_ann_SC_30_absLFC, time_ten_vs_zero_ann_SC_30_absLFC, time_thirty_vs_zero_ann_SC_30_absLFC) %>% mutate(species = rep("SC"))
# order the times as a factor
time_all_ann_SC_30_absLFC$time <- factor(time_all_ann_SC_30_absLFC$time, levels = c("zero", "ten", "thirty"))

# filter time_all_ann_SC_30_absLFC to select for unique only
time_all_ann_SC_30_only_absLFC <- time_all_ann_SC_30_absLFC %>% filter(transcript_id %in% thirty_only_SC$transcript_id) %>% mutate(category = rep("thirty_only"))

# bind ten only, thirty only, and unique together for plotting
time_SC_absLFC <- rbind(time_all_ann_SC_10_only_absLFC, time_all_ann_SC_30_only_absLFC, time_all_ann_SC_overlap_absLFC)

# change order of lines plotted; 10, 30, then both
time_SC_absLFC$category <- factor(time_SC_absLFC$category, levels = c("ten_only", "thirty_only", "both"))


# change time to numbers and numeric
time_SC_absLFC$time <- as.numeric(time_SC_absLFC$time %>% gsub("ten", "10", .)  %>% gsub("thirty", "30", .)  %>% gsub("zero", "0", .))

# prepare dataset for labels; filtering for top LFC in each category
label_30_SC <- time_SC_absLFC %>% filter(category == "thirty_only", time == "30") %>% slice_max(order_by = absLFC, n = 2)
label_10_SC <- time_SC_absLFC %>% filter(category == "ten_only", time == "10") %>% slice_max(order_by = absLFC, n = 2)
label_both_SC <- time_SC_absLFC %>% filter(category == "both", time == "30")%>% slice_max(order_by = absLFC, n = 2)
# add other time points of these transcript_ids so the line can be plotted thicker, don't use left_join need to subset 
thickline_30_SC <- subset(time_SC_absLFC, transcript_id %in% label_30_SC$transcript_id)
thickline_10_SC <- subset(time_SC_absLFC, transcript_id %in% label_10_SC$transcript_id)
thickline_both_SC <- subset(time_SC_absLFC, transcript_id %in% label_both_SC$transcript_id)

# PLOT
ggline_time_SC <- ggplot(data = time_SC_absLFC, aes(x = time, y = absLFC, group = transcript_id, color = category)) + 
  geom_line(alpha = 0) +
  geom_hline(aes(yintercept=1), linetype = "dashed", color = "gray40") +
  guides(color = guide_legend(override.aes = list(alpha = 1) )) +
  geom_line(data = (time_SC_absLFC %>% filter(category == "thirty_only")), aes(x = time, y = absLFC, group = transcript_id), size = 0.4, color = thirty.only, alpha = 0.4) +
  geom_line(data = thickline_30_SC, aes(x = time, y = absLFC, group = transcript_id), size = 0.5, color = thirty.only) +
  geom_line(data = (time_SC_absLFC %>% filter(category == "ten_only")), aes(x = time, y = absLFC, group = transcript_id), size = 0.5, color = ten.only, alpha = 0.6) +
  geom_line(data = thickline_10_SC, aes(x = time, y = absLFC, group = transcript_id), size = 0.7, color = ten.only) +
  geom_line(data = (time_SC_absLFC %>% filter(category == "both")), aes(x = time, y = absLFC, group = transcript_id), size = 0.5, color = both, alpha = 0.8) +
  geom_line(data = thickline_both_SC, aes(x = time, y = absLFC, group = transcript_id), size = 0.7, color = both) +
  scale_color_manual(name = "Significance Category", values = c(both = both, ten_only = ten.only, thirty_only = thirty.only), labels = c(both = "Both", ten_only = expression(T[10]~only), thirty_only= expression(T[30]~only))) +
  scale_x_continuous(labels = c("0" = "0", "10" = "10", "20" = "", "30" = "30")) +
  theme_light() +
  theme(text = element_text(size = 10),
        legend.position=c(0.13,0.79),
        legend.background = element_rect(colour = "gray"),
        legend.text.align = 0,
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.y=element_text(size = 10)) +
  geom_text(data = label_30_SC, aes(label = transcript_id), nudge_x = 1.3, nudge_y = 0.07, size = 2.6, show.legend = F, alpha = 1, fontface = "bold") +
  geom_text(data = label_10_SC, aes(label = transcript_id), nudge_x = 1.3, nudge_y = 0.07, size = 2.6, show.legend = F, alpha = 1, fontface = "bold") +
  geom_text(data = label_both_SC, aes(label = transcript_id), nudge_x = 1.3, nudge_y = 0.07, size = 2.6, show.legend = F, alpha = 1, fontface = "bold") +
  labs(x="\nMinutes of UV Radiation Exposure", y = expression(Absolute~Value~of~Log2~Fold~Change~from~T[0])) 

ggline_time_SC  


ggsave("plots/ggline_time_SC.pdf", ggline_time_SC, width = 2.5, height = 2, units = "in", scale = 3)
  

#### individual LFC plots of BR sig genes that overlap and are distinct in 10 and 30 minutes (different in two species) ####

# first get lists of uniqueness and overlapping in the sig at 10 and 30 in BR
overlap_BR <- time_ten_vs_zero_df_BR %>% filter(transcript_id %in% time_thirty_vs_zero_df_BR$transcript_id) %>% select(transcript_id)
ten_only_BR <- time_ten_vs_zero_df_BR %>% filter(!transcript_id %in% overlap_BR$transcript_id) %>% select(transcript_id)
thirty_only_BR <- time_thirty_vs_zero_df_BR %>% filter(!transcript_id %in% overlap_BR$transcript_id) %>% select(transcript_id)

# take the transcripts that are significant from ten and calculate their absolute value LFC at ten
time_ten_vs_zero_ann_BR_10_absLFC <- time_ten_vs_zero_ann_BR %>% mutate(absLFC = abs(log2FoldChange)) %>% select(transcript_id, absLFC, lfcSE)
time_ten_vs_zero_ann_BR_10_absLFC$time <- rep("ten")

# take these same transcripts and make an identical data frame that has LFC at 0 set to 0 (normalization so they can all be plotted together) 
time_zero_ann_BR_10_absLFC <- time_ten_vs_zero_ann_BR_10_absLFC %>% mutate(time = rep("zero"), absLFC = rep(0), lfcSE = rep(0))

# take these same transcripts and pull their LFC at time 30
# first convert results to df and make rownames into a column
BR_thirty_10_res <- (as.data.frame(results(dds_BR, name = "time_thirty_vs_zero")))
BR_thirty_10_res$transcript_id <- rownames(BR_thirty_10_res)

time_thirty_vs_zero_ann_BR_10_absLFC <- subset(BR_thirty_10_res, transcript_id %in% time_zero_ann_BR_10_absLFC$transcript_id) %>% 
  mutate(absLFC = abs(log2FoldChange)) %>% 
  select(transcript_id, absLFC, lfcSE) %>%
  mutate(time = rep("thirty"))

# remove rownames 
rownames(time_thirty_vs_zero_ann_BR_10_absLFC) <- NULL

# put these together for plotting
time_all_ann_BR_10_absLFC <- rbind(time_zero_ann_BR_10_absLFC, time_ten_vs_zero_ann_BR_10_absLFC, time_thirty_vs_zero_ann_BR_10_absLFC) %>% mutate(species = rep("BR"))
# order the times as a factor
time_all_ann_BR_10_absLFC$time <- factor(time_all_ann_BR_10_absLFC$time, levels = c("zero", "ten", "thirty"))
# convert abs(LFC) to numeric

# filter only those that are uniquely sig in t = 10! then repeat for those sig at 30
# overlap between 10 and 30 minutes
# first need to turn rownames into a column for each
time_ten_vs_zero_df_BR$transcript_id <- rownames(time_ten_vs_zero_df_BR)
# remove rownames 
rownames(time_ten_vs_zero_df_BR) <- NULL

time_thirty_vs_zero_df_BR$transcript_id <- rownames(time_thirty_vs_zero_df_BR)
# remove rownames 
rownames(time_thirty_vs_zero_df_BR) <- NULL

# filter time_all_ann_BR_10_absLFC to select for unique only
time_all_ann_BR_10_only_absLFC <- time_all_ann_BR_10_absLFC %>% filter(transcript_id %in% ten_only_BR$transcript_id) %>% mutate(category = rep("ten_only"))

# take the overlap genes and make into an absLFC table
time_all_ann_BR_overlap_absLFC <- time_all_ann_BR_10_absLFC %>% filter(transcript_id %in% overlap_BR$transcript_id) %>% mutate(category = rep("both"))

# take the overlap genes and make into an absLFC table, including LFCs from thirty and from ten as above
###
###
###
###
# take the transcripts that are significant from thirty and calculate their absolute value LFC at ten
time_thirty_vs_zero_ann_BR_30_absLFC <- time_thirty_vs_zero_ann_BR %>% mutate(absLFC = abs(log2FoldChange)) %>% select(transcript_id, absLFC, lfcSE)
time_thirty_vs_zero_ann_BR_30_absLFC$time <- rep("thirty")

# take these same transcripts and make an identical data frame that has LFC at 0 set to 0 (normalization so they can all be plotted together) 
time_zero_ann_BR_30_absLFC <- time_thirty_vs_zero_ann_BR_30_absLFC %>% mutate(time = rep("zero"), absLFC = rep(0), lfcSE = rep(0))

# take these same transcripts and pull their LFC at time 10
# first convert results to df and make rownames into a column
BR_ten_30_res <- (as.data.frame(results(dds_BR, name = "time_ten_vs_zero")))
BR_ten_30_res$transcript_id <- rownames(BR_ten_30_res)

time_ten_vs_zero_ann_BR_30_absLFC <- subset(BR_ten_30_res, transcript_id %in% time_zero_ann_BR_30_absLFC$transcript_id) %>% 
  mutate(absLFC = abs(log2FoldChange)) %>% 
  select(transcript_id, absLFC, lfcSE) %>%
  mutate(time = rep("ten"))

# remove rownames 
rownames(time_ten_vs_zero_ann_BR_30_absLFC) <- NULL

# put these together for plotting
time_all_ann_BR_30_absLFC <- rbind(time_zero_ann_BR_30_absLFC, time_ten_vs_zero_ann_BR_30_absLFC, time_thirty_vs_zero_ann_BR_30_absLFC) %>% mutate(species = rep("BR"))
# order the times as a factor
time_all_ann_BR_30_absLFC$time <- factor(time_all_ann_BR_30_absLFC$time, levels = c("zero", "ten", "thirty"))

# filter time_all_ann_BR_30_absLFC to select for unique only
time_all_ann_BR_30_only_absLFC <- time_all_ann_BR_30_absLFC %>% filter(transcript_id %in% thirty_only_BR$transcript_id) %>% mutate(category = rep("thirty_only"))

# bind ten only, thirty only, and unique together for plotting
time_BR_absLFC <- rbind(time_all_ann_BR_10_only_absLFC, time_all_ann_BR_30_only_absLFC, time_all_ann_BR_overlap_absLFC)

# change order of lines plotted; 10, 30, then both
time_BR_absLFC$category <- factor(time_BR_absLFC$category, levels = c("ten_only", "thirty_only", "both"))


# change time to numbers and numeric
time_BR_absLFC$time <- as.numeric(time_BR_absLFC$time %>% gsub("ten", "10", .)  %>% gsub("thirty", "30", .)  %>% gsub("zero", "0", .))

# prepare dataset for labels; filtering for top LFC in each category
label_30_BR <- time_BR_absLFC %>% filter(category == "thirty_only", time == "30") %>% slice_max(order_by = absLFC, n = 2)
label_10_BR <- time_BR_absLFC %>% filter(category == "ten_only", time == "10") %>% slice_max(order_by = absLFC, n = 2)
label_both_BR <- time_BR_absLFC %>% filter(category == "both", time == "30")%>% slice_max(order_by = absLFC, n = 2)
# add other time points of these transcript_ids so the line can be plotted thicker, don't use left_join need to subset 
thickline_30_BR <- subset(time_BR_absLFC, transcript_id %in% label_30_BR$transcript_id)
thickline_10_BR <- subset(time_BR_absLFC, transcript_id %in% label_10_BR$transcript_id)
thickline_both_BR <- subset(time_BR_absLFC, transcript_id %in% label_both_BR$transcript_id)

# PLOT
ggline_time_BR <- ggplot(data = time_BR_absLFC, aes(x = time, y = absLFC, group = transcript_id, color = category)) + 
  geom_line(alpha = 0) +
  geom_hline(aes(yintercept=1), linetype = "dashed", color = "gray40") +
  guides(color = guide_legend(override.aes = list(alpha = 1) )) +
  geom_line(data = (time_BR_absLFC %>% filter(category == "thirty_only")), aes(x = time, y = absLFC, group = transcript_id), size = 0.4, color = thirty.only, alpha = 0.4) +
  geom_line(data = thickline_30_BR, aes(x = time, y = absLFC, group = transcript_id), size = 0.5, color = thirty.only) +
  geom_line(data = (time_BR_absLFC %>% filter(category == "ten_only")), aes(x = time, y = absLFC, group = transcript_id), size = 0.5, color = ten.only, alpha = 0.6) +
  geom_line(data = thickline_10_BR, aes(x = time, y = absLFC, group = transcript_id), size = 0.7, color = ten.only) +
  geom_line(data = (time_BR_absLFC %>% filter(category == "both")), aes(x = time, y = absLFC, group = transcript_id), size = 0.5, color = both, alpha = 0.8) +
  geom_line(data = thickline_both_BR, aes(x = time, y = absLFC, group = transcript_id), size = 0.7, color = both) +
  scale_color_manual(name = "Significance Category", values = c(both = both, ten_only = ten.only, thirty_only = thirty.only), labels = c(both = "Both", ten_only = expression(T[10]~only), thirty_only= expression(T[30]~only))) +
  scale_x_continuous(labels = c("0" = "0", "10" = "10", "20" = "", "30" = "30")) +
  theme_light() +
  theme(text = element_text(size = 10),
        legend.position=c(0.13,0.79),
        legend.background = element_rect(colour = "gray"),
        legend.text.align = 0,
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.y=element_text(size = 10)) +
  geom_text(data = label_30_BR, aes(label = transcript_id), nudge_x = 1.3, nudge_y = 0.1, size = 2.6, show.legend = F, alpha = 1, fontface = "bold") +
  geom_text(data = label_10_BR, aes(label = transcript_id), nudge_x = 1.3, nudge_y = 0.1, size = 2.6, show.legend = F, alpha = 1, fontface = "bold") +
  geom_text(data = label_both_BR, aes(label = transcript_id), nudge_x = 1.3, nudge_y = 0.1, size = 2.6, show.legend = F, alpha = 1, fontface = "bold") +
  labs(x="\nMinutes of UV Radiation Exposure", y = expression(Absolute~Value~of~Log2~Fold~Change~from~T[0])) 

ggline_time_BR  


ggsave("plots/ggline_time_BR.pdf", ggline_time_BR, width = 2.5, height = 2, units = "in", scale = 3)


