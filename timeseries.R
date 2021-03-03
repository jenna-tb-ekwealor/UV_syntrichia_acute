# BiocManager::install("GOsummaries")
# BiocManager::install("apeglm")
# BiocManager::install("DESeq2")
# browseVignettes(package = "DESeq2")
# browseVignettes(package = "GOsummaries")
# install.packages(c("ashr", "ggbeeswarm", "pheatmap"))
# devtools::install_github("teunbrand/ggh4x")

library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)

# set colors for plots 
zero <- "#18f27f" #green
ten <- "#099c4e"
thirty <- "#022412"

srur <- "#aa0a3c"
scan <- "#f0f032"
  
setwd("/Users/jennaekwealor/Documents/dissertation_repositories/UV_syntrichia_acute")

directory <- ("data/htseq-count_10142020/")

# load sample table
sampleTable <- read.csv("sample_key.csv", stringsAsFactors = TRUE)

# filter for just the species i want to look at here
sampleTable_SC <- sampleTable %>% filter(species == "SC")
sampleTable_BR <- sampleTable %>% filter(species == "BR")

# load annotations
func_annot <- read.csv(file = "NEW_func_annot_GO_terms_SC.csv")


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


p_pca_all_SC <- pca_df_SC %>% group_by(time) %>% ggplot(.,aes(x=PC1,y=PC2))
p_pca_all_SC <- p_pca_all_SC +
  geom_point(aes(color = time), size = 3) + 
  ggforce::geom_mark_ellipse(aes(color = time)) +  
  scale_color_manual(name = "UV Treatment", labels = c(expression(T[0]), expression(T[10]), expression(T[30])),values = c(zero, ten, thirty)) +
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

# 
# #### plot heat map top varying 1000 GENES ####
# topVarGenes_SC <- head(order(rowVars(assay(vsd_SC)), decreasing = TRUE), 1000)
# mat_SC  <- assay(vsd_SC)[ topVarGenes_SC, ]
# mat_SC  <- mat_SC - rowMeans(mat_SC)
# anno_SC <- as.data.frame(colData(vsd_SC)[, c("time")])
# ann_colors = list(time = c(zero = zero, ten = ten, thirty = thirty))
# 
# heatmap_SC <- pheatmap(mat = mat_SC, 
#                     annotation_col = anno_SC, 
#                     annotation_colors = ann_colors, 
#                     show_rownames = FALSE,
#                     show_colnames = F, 
#                     drop_levels = TRUE,
#                     legend = TRUE, 
#                     annotation_legend = TRUE,
#                     fontsize = 16,
#                     main = "TOP 1000 Heatmap")
# heatmap_SC
# 
# pdf("heatmap_1000.pdf") 
# heatmap
# dev.off


#### ten vs. zero ####

resultsNames(dds_SC) # gives a list of "name" comparisons you can use as shortcut instead of "contrasts"
results(dds_SC, name = "time_ten_vs_zero")  -> time_ten_vs_zero_SC
# make df of UV genes with svalue of less than 0.005 and shrunken LFC of at least 2
time_ten_vs_zero_SC %>%
  as.data.frame() %>%
  mutate(transcript_id = rownames(.) ) %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  filter(padj < 0.005) -> time_ten_vs_zero_df_SC

# join annotations to results file
inner_join(x = func_annot, y = time_ten_vs_zero_df_SC, by = "transcript_id") %>%
  distinct(transcript_id,  .keep_all = T) -> time_ten_vs_zero_ann_SC

write.csv(time_ten_vs_zero_ann_SC, "tables/SC/time_ten_vs_zero_ann.csv")

length(time_ten_vs_zero_ann_SC$transcript_id_SC)
# 0 genes
time_ten_vs_zero_ann_SC %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
  arrange(abs_log2FoldChange) %>%
  tail(n = 10) -> time_ten_vs_zero_ann_strict_SC
write.csv(time_ten_vs_zero_ann_strict_SC, "tables/SC/time_ten_vs_zero_ann_strict.csv")


#### Plot LFC  ####

# lfcshrink 
lfcshrink_time_ten_vs_zero_SC <- lfcShrink(dds_SC, coef="time_ten_vs_zero", lfcThreshold=2, type="apeglm")

# MAplot 
plotMA(lfcshrink_time_ten_vs_zero_SC, ylim=c(-4,4), cex=0.5, alpha = 0.005)
abline(h=c(-2,2), col="lightgray", lwd=2) 
while (!is.null(dev.list()))  dev.off()

pdf("plots/maplot_lfcshrink_time_ten_vs_zero_SC.pdf") 
plotMA(lfcshrink_time_ten_vs_zero_SC, ylim=c(-4,4), cex=0.8, alpha = 0.005,
       xlab = "Mean of Normalized Counts",
       ylab = "Log Fold Change")
abline(h=c(-2,2), col="lightgray", lwd=2) 
while (!is.null(dev.list()))  dev.off()


#### thirty vs. zero ####

resultsNames(dds_SC) # gives a list of "name" comparisons you can use as shortcut instead of "contrasts"
results(dds_SC, name = "time_thirty_vs_zero")  -> time_thirty_vs_zero_SC
# make df of UV genes with svalue of less than 0.005 and shrunken LFC of at least 2
time_thirty_vs_zero_SC %>%
  as.data.frame() %>%
  mutate(transcript_id = rownames(.) ) %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  filter(padj < 0.005) -> time_thirty_vs_zero_df_SC

# join annotations to results file
inner_join(x = func_annot, y = time_thirty_vs_zero_df_SC, by = "transcript_id") %>%
  distinct(transcript_id,  .keep_all = T) -> time_thirty_vs_zero_ann_SC

write.csv(time_thirty_vs_zero_ann_SC, "tables/SC/time_thirty_vs_zero_ann.csv")

length(time_thirty_vs_zero_ann_SC$transcript_id)
# 56 genes
time_thirty_vs_zero_ann_SC %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
  arrange(abs_log2FoldChange) %>%
  tail(n = 10) -> time_thirty_vs_zero_ann_strict_SC
write.csv(time_thirty_vs_zero_ann_strict_SC, "tables/SC/time_thirty_vs_zero_ann_strict.csv")


#### Plot LFC  ####

# lfcshrink 
lfcshrink_time_thirty_vs_zero_SC <- lfcShrink(dds_SC, coef="time_thirty_vs_zero", lfcThreshold=2, type="apeglm")

# MAplot 
plotMA(lfcshrink_time_thirty_vs_zero_SC, ylim=c(-4,4), cex=0.5, alpha = 0.005)
abline(h=c(-2,2), col="lightgray", lwd=2) 
while (!is.null(dev.list()))  dev.off()

pdf("plots/maplot_lfcshrink_time_thirty_vs_zero_SC.pdf") 
plotMA(lfcshrink_time_thirty_vs_zero_SC, ylim=c(-4,4), cex=0.8, alpha = 0.005,
       xlab = "Mean of Normalized Counts",
       ylab = "Log Fold Change")
abline(h=c(-2,2), col="lightgray", lwd=2) 
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

# plot PCA
vsd_BR <- vst(dds_BR, blind=T)
plotPCA(vsd_BR, intgroup="time") + theme_minimal()
# ntop = number of top genes to use for principal components, selected by highest row variance

# use returnData = TRUE should the function only return the data.frame of PC1 and PC2 with intgroup covariates for custom plotting
pca_df_BR <- plotPCA(vsd_BR, intgroup=c("time"), returnData = TRUE) 


p_pca_all_BR <- pca_df_BR %>% group_by(time) %>% ggplot(.,aes(x=PC1,y=PC2))
p_pca_all_BR <- p_pca_all_BR +
  geom_point(aes(color = time), size = 3) + 
  # ADD ggforce's ellipses
  ggforce::geom_mark_ellipse(aes(color = time)) +  
  scale_color_manual(name = "UV Treatment", labels = c(expression(T[0]), expression(T[10]), expression(T[30])),values = c(zero, ten, thirty)) +
  theme_minimal() +
  theme(legend.text.align = 0,
        text = element_text(size = 14)) +
  xlab("PC1 (43%)") + ylab("PC2 (23%)") +
  xlim(-6, 15) +
  ylim(-7, 8)
p_pca_all_BR

pdf("plots/p_pca_all_BR.pdf") 
p_pca_all_BR
dev.off()

# #### plot heat map top varying 1000 GENES ####
# topVarGenes_BR <- head(order(rowVars(assay(vsd_BR)), decreasing = TRUE), 1000)
# mat_BR  <- assay(vsd_BR)[ topVarGenes_BR, ]
# mat_BR  <- mat_BR - rowMeans(mat_BR)
# anno_BR <- as.data.frame(colData(vsd_BR)[, c("time")])
# ann_colors = list(time = c(zero = zero, ten = ten, thirty = thirty))
# 
# heatmap_BR <- pheatmap(mat = mat_BR, 
#                        annotation_col = anno_BR, 
#                        annotation_colors = ann_colors, 
#                        show_rownames = FALSE,
#                        show_colnames = F, 
#                        drop_levels = TRUE,
#                        legend = TRUE, 
#                        annotation_legend = TRUE,
#                        fontsize = 16,
#                        main = "TOP 1000 Heatmap")
# heatmap_BR
# 
# pdf("heatmap_1000.pdf") 
# heatmap
# dev.off()


### ten vs. zero 

resultsNames(dds_BR) # gives a list of "name" comparisons you can use as shortcut instead of "contrasts"
results(dds_BR, name = "time_ten_vs_zero")  -> time_ten_vs_zero_BR
# make df of UV genes with svalue of less than 0.005 and shrunken LFC of at least 2
time_ten_vs_zero_BR %>%
  as.data.frame() %>%
  mutate(transcript_id = rownames(.) ) %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  filter(padj < 0.005) -> time_ten_vs_zero_df_BR

# join annotations to results file
inner_join(x = func_annot, y = time_ten_vs_zero_df_BR, by = "transcript_id") %>%
  distinct(transcript_id,  .keep_all = T) -> time_ten_vs_zero_ann_BR

write.csv(time_ten_vs_zero_ann_BR, "tables/BR/time_ten_vs_zero_ann.csv")

length(time_ten_vs_zero_ann_BR$transcript_id_BR)
# 0 genes
time_ten_vs_zero_ann_BR %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
  arrange(abs_log2FoldChange) %>%
  tail(n = 10) -> time_ten_vs_zero_ann_strict_BR
write.csv(time_ten_vs_zero_ann_strict_BR, "tables/BR/time_ten_vs_zero_ann_strict.csv")


#### Plot LFC #### 

# lfcshrink 
lfcshrink_time_ten_vs_zero_BR <- lfcShrink(dds_BR, coef="time_ten_vs_zero", lfcThreshold=2, type="apeglm")

# MAplot 
plotMA(lfcshrink_time_ten_vs_zero_BR, ylim=c(-4,4), cex=0.5, alpha = 0.005)
abline(h=c(-2,2), col="lightgray", lwd=2) 
while (!is.null(dev.list()))  dev.off()

pdf("plots/maplot_lfcshrink_time_ten_vs_zero_BR.pdf") 
plotMA(lfcshrink_time_ten_vs_zero_BR, ylim=c(-4,4), cex=0.8, alpha = 0.005,
       xlab = "Mean of Normalized Counts",
       ylab = "Log Fold Change")
abline(h=c(-2,2), col="lightgray", lwd=2) 
while (!is.null(dev.list()))  dev.off()


#### thirty vs. zero ####

resultsNames(dds_BR) # gives a list of "name" comparisons you can use as shortcut instead of "contrasts"
results(dds_BR, name = "time_thirty_vs_zero")  -> time_thirty_vs_zero_BR
# make df of UV genes with svalue of less than 0.005 and shrunken LFC of at least 2
time_thirty_vs_zero_BR %>%
  as.data.frame() %>%
  mutate(transcript_id = rownames(.) ) %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
  filter(padj < 0.005) -> time_thirty_vs_zero_df_BR

# join annotations to results file
inner_join(x = func_annot, y = time_thirty_vs_zero_df_BR, by = "transcript_id") %>%
  distinct(transcript_id,  .keep_all = T) -> time_thirty_vs_zero_ann_BR

write.csv(time_thirty_vs_zero_ann_BR, "tables/BR/time_thirty_vs_zero_ann.csv")

length(time_thirty_vs_zero_ann_BR$transcript_id)
# 32 genes
time_thirty_vs_zero_ann_BR %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
  arrange(abs_log2FoldChange) %>%
  tail(n = 10) -> time_thirty_vs_zero_ann_strict_BR
write.csv(time_thirty_vs_zero_ann_strict_BR, "tables/BR/time_thirty_vs_zero_ann_strict.csv")


#### Plot LFC  ####

# lfcshrink 
lfcshrink_time_thirty_vs_zero_BR <- lfcShrink(dds_BR, coef="time_thirty_vs_zero", lfcThreshold=2, type="apeglm")

# MAplot 
plotMA(lfcshrink_time_thirty_vs_zero_BR, ylim=c(-4,4), cex=0.5, alpha = 0.005)
abline(h=c(-2,2), col="lightgray", lwd=2) 
while (!is.null(dev.list()))  dev.off()

pdf("plots/maplot_lfcshrink_time_thirty_vs_zero_BR.pdf") 
plotMA(lfcshrink_time_thirty_vs_zero_BR, ylim=c(-4,4), cex=0.8, alpha = 0.005,
       xlab = "Mean of Normalized Counts",
       ylab = "Log Fold Change")
abline(h=c(-2,2), col="lightgray", lwd=2) 
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


p_pca_all <- pca_df %>% group_by(species, time) %>% ggplot(.,aes(x=PC1,y=PC2))
p_pca_all <- p_pca_all +
  geom_point(aes(color = time), size = 2) + 
  stat_ellipse(aes(x=PC1,y=PC2,group=species, linetype = species),type = "norm", level = 0.95) +
  scale_linetype_manual(values=c("solid", "longdash"), name = "Species", labels = c(expression(italic("S. ruralis")), expression(italic("S. caninervis")))) +
  scale_color_manual(name = "UV Treatment", labels = c(expression(T[0]), expression(T[10]), expression(T[30])),values = c(zero, ten, thirty)) +
  theme_minimal() +
  theme(legend.text.align = 0,
        text = element_text(size = 14)) +
  guides(linetype = guide_legend(order = 1),
         shape = guide_legend(order = 2),
         color = guide_legend(order = 3)) +
  xlab("PC1 (93%)") + ylab("PC2 (3%)") 
p_pca_all

pdf("plots/p_pca_all.pdf") 
p_pca_all
dev.off()

# #### plot heat map top varying 1000 GENES ####
# topVarGenes <- head(order(rowVars(assay(vsd_BR)), decreasing = TRUE), 1000)
# mat  <- assay(vsd)[ topVarGenes, ]
# mat  <- mat - rowMeans(mat)
# anno <- as.data.frame(colData(vsd)[, c("species", "time")])
# ann_colors = list(
#   time = c(zero = zero, ten = ten, thirty = thirty),
#   species = c(BR = "black", SC = "darkgray")
# )
# 
# heatmap <- pheatmap(mat = mat, 
#                     annotation_col = anno, 
#                     annotation_colors = ann_colors, 
#                     show_rownames = FALSE,
#                     show_colnames = F, 
#                     drop_levels = TRUE,
#                     legend = TRUE, 
#                     annotation_legend = TRUE,
#                     fontsize = 16,
#                     main = "TOP 1000 Heatmap")
# heatmap
# 
# pdf("heatmap_1000.pdf") 
# heatmap
# dev.off()




#### LRT to remove species-specific effects ####
# Genes with small p values from this test are those which at one or more time points after time 0 showed a species-specific effect. 
# Note therefore that this will not give small p values to genes that moved up or down over time in the same way in both species

dds_intLRT <- DESeq(dds_int, test="LRT", reduced = ~ species + time)
resLRT <- results(dds_intLRT)
resLRT$symbol <- mcols(dds_intLRT)$symbol
head(resLRT[order(resLRT$padj),], 4)
summary(resLRT)
# 6859 total genes


# list of genes diff in the two species
resultsNames(dds_intLRT) # gives a list of "name" comparisons you can use as shortcut instead of "contrasts"
results(dds_intLRT, name = "speciesSC.timethirty")  -> speciesSC.timethirty

# relevel treatments so that one is the 'reference' in differential abundance analyses
dds_intLRT$time <- relevel(dds_int$time, ref = "zero") # level treatments so one is marked as 'zero' as the baseline in terms of diff abundance
dds_intLRT$species <- relevel(dds_int$species, ref = "BR")

# make df of UV genes with svalue of less than 0.005 
speciesSC.timethirty %>%
  as.data.frame() %>%
  mutate(transcript_id = rownames(.) ) %>%
  filter(padj < 0.005) -> speciesSC.timethirty_df

# join annotations to results file
inner_join(x = func_annot, y = speciesSC.timethirty_df, by = "transcript_id") %>%
  distinct(transcript_id,  .keep_all = T) -> speciesSC.timethirty_ann

write.csv(speciesSC.timethirty_ann, "tables/speciesSC.timethirty_ann.csv")

length(speciesSC.timethirty_ann$transcript_id)
# 24 genes
speciesSC.timethirty_ann %>%
  arrange(padj) %>%
  head(n = 10) -> speciesSC.timethirty_ann_strict
write.csv(speciesSC.timethirty_ann_strict, "tables/speciesSC.timethirty_ann_strict.csv")




# We can plot the counts for the groups over time using ggplot2, 
# for the gene with the smallest adjusted p value, 
# testing for condition-dependent time profile and accounting for differences at time 0 
# Keep in mind that the interaction terms are the difference 
# between the two groups at a given time after accounting for the difference at time 0




# will need to get the list of genes then plotCounts(), looping over the the list
# creating a data frame for ggplot to facet plot them all at once 
num_transcripts <- length(speciesSC.timethirty_ann_strict$transcript_id)
count_data <- NULL
final_counts_speciesSC.timethirty_ann_strict <- NULL

for(g in 1:num_transcripts) {
  
  one_transcript <- speciesSC.timethirty_ann_strict$transcript_id[g]
  
  # returnData = true returns data frame so it can be plotted with ggplot2, etc.
  count_data <- plotCounts(dds_intLRT, gene=one_transcript, intgroup=c("time","species"), 
                           returnData=TRUE)
  count_data$transcript_id <- rep(one_transcript)
  
  # convert rownames to column and group by treatments
  count_data <- count_data %>% tibble::rownames_to_column(., var = "sample") 
  final_counts_speciesSC.timethirty_ann_strict <- rbind(final_counts_speciesSC.timethirty_ann_strict, count_data)
}
final_counts_speciesSC.timethirty_ann_strict

# log transform the counts
final_counts_speciesSC.timethirty_ann_strict <- dplyr::mutate(final_counts_speciesSC.timethirty_ann_strict, log2_count = log2(count))

# grep to replace species names
final_counts_speciesSC.timethirty_ann_strict$species <- final_counts_speciesSC.timethirty_ann_strict$species %>% gsub("BR", "S. ruralis", .) %>% gsub("SC", "S. caninervis", .)

# add gene names
final_counts_speciesSC.timethirty_ann_strict <- left_join(final_counts_speciesSC.timethirty_ann_strict, func_annot, by = "transcript_id", keep.both)
 
final_counts_speciesSC.timethirty_ann_strict$transcript_id <- factor(final_counts_speciesSC.timethirty_ann_strict$transcript_id, levels = c("Sc_g08728", "Sc_g00568", "Sc_g00609", "Sc_g13403", "Sc_g02848", "Sc_g10759", "Sc_g14586", "Sc_g14304", "Sc_g12657", "Sc_g05211"))

ggline <- ggplot(final_counts_speciesSC.timethirty_ann_strict,
       aes(x = time, y = log2_count, color = species, group = species)) + 
      facet_wrap("transcript_id", ncol = 5) +
  geom_point(aes(fill=time), size=2, shape=21, stroke=1) + 
  stat_summary(fun=mean, geom="line")  +
  scale_color_manual(values = c("gray", "black")) +
  scale_fill_manual(values = c(zero=zero, ten=ten, thirty=thirty), guide = 'none') +
  theme_light() +
  theme(text = element_text(size = 12),
        legend.position="bottom",
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.y=element_text(size = 14),
        strip.background = element_rect(fill = "black")) +
        scale_x_discrete(labels=c("zero" = "0", "ten" = "10",
                              "thirty" = "30")) +
  # scale_linetype_manual(values=c("longdash", "solid"), name="Species") +
  labs(x="\nMinutes of UV Exposure", y = "Log2 Transformed Normalized Transcript Count")
  # ggtitle("Top 10 Transcripts \n") 

ggline 

ggsave("plots/ggline.pdf", ggline, width = 6, height = 6, units = "in", scale = 1.3)






