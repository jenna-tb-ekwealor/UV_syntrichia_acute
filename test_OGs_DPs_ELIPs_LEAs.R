library("dplyr")

setwd("/Users/jennaekwealor/Documents/dissertation_repositories/UV_syntrichia_acute")


srur <- "#18f27f"
scan <- "#FF00FF"


species_effect <- read.csv(file = "tables/LRT_ann.csv") %>% dplyr::select(transcript_id:padj)

UV_SC_10 <- read.csv(file = "tables/SC/time_ten_vs_zero_ann.csv") %>% dplyr::select(names:padj)

UV_BR_10 <- read.csv(file = "tables/BR/time_ten_vs_zero_ann.csv") %>% dplyr::select(names:padj)

UV_SC_30 <- read.csv(file = "tables/SC/time_thirty_vs_zero_ann.csv") %>% dplyr::select(names:padj)

UV_BR_30 <- read.csv(file = "tables/BR/time_thirty_vs_zero_ann.csv") %>% dplyr::select(names:padj)


elip <- read.csv(file = "caninervis_elips.csv")
lea <- read.csv(file = "caninervis_leas.csv")


#### check if elips in each of above ####
species_effect_elip <- species_effect %>% filter(transcript_id %in% elip$transcript_id)
length(unique(species_effect$transcript_id))
# [1] 69
length(unique(species_effect_elip$transcript_id))
# [1] 1



UV_SC_10_elip <- UV_SC_10 %>% filter(transcript_id %in% elip$transcript_id)
length(unique(UV_SC_10$transcript_id))
# [1] 10
length(unique(UV_SC_10_elip$transcript_id))
# [1] 1

UV_BR_10_elip <- UV_BR_10 %>% filter(transcript_id %in% elip$transcript_id)
length(unique(UV_BR_10$transcript_id))
# [1] 18
length(unique(UV_BR_10_elip$transcript_id))
# [1] 0


UV_SC_30_elip <- UV_SC_30 %>% filter(transcript_id %in% elip$transcript_id)
length(unique(UV_SC_30$transcript_id))
# [1] 126
length(unique(UV_SC_30_elip$transcript_id))
# [1] 2

UV_BR_30_elip <- UV_BR_30 %>% filter(transcript_id %in% elip$transcript_id)
length(unique(UV_BR_30$transcript_id))
# [1] 38
length(unique(UV_BR_30_elip$transcript_id))
# [1] 0


#### check if leas in each of above ####
species_effect_lea <- species_effect %>% filter(transcript_id %in% lea$transcript_id)
length(unique(species_effect$transcript_id))
# [1] 69
length(unique(species_effect_lea$transcript_id))
# [1] 1


UV_SC_10_lea <- UV_SC_10 %>% filter(transcript_id %in% lea$transcript_id)
length(unique(UV_SC_10$transcript_id))
# [1] 10
length(unique(UV_SC_10_lea$transcript_id))
# [1] 0

UV_BR_10_lea <- UV_BR_10 %>% filter(transcript_id %in% lea$transcript_id)
length(unique(UV_BR_10$transcript_id))
# [1] 18
length(unique(UV_BR_10_lea$transcript_id))
# [1] 0



UV_SC_30_lea <- UV_SC_30 %>% filter(transcript_id %in% lea$transcript_id)
length(unique(UV_SC_30$transcript_id))
# [1] 126
length(unique(UV_SC_30_lea$transcript_id))
# [1] 4

UV_BR_30_lea <- UV_BR_30 %>% filter(transcript_id %in% lea$transcript_id)
length(unique(UV_BR_30$transcript_id))
# [1] 38
length(unique(UV_BR_30_lea$transcript_id))
# [1] 2


### venn diagram of UV results ####

# euler diagram of number sig mets in each treatment 10 mins
VennDiagram::venn.diagram(
  x = list(unique(UV_BR_10$transcript_id),unique(UV_SC_10$transcript_id)),
  euler.d =T,
  scaled = T,
  category.names = c(expression(italic("S. ruralis")) , expression(italic("S. caninervis"))),
  filename  = "venn_acute_10.tiff",
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  height = 800 , 
  width = 800 , 
  resolution = 600,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c(srur, scan),
  
  # Numbers
  cex = 1.25,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0), 
  cat.dist = c(0.03, 0.03),
  cat.fontfamily = "sans"
) 

# euler diagram of number sig mets in each treatment 30 mins
VennDiagram::venn.diagram(
  x = list(unique(UV_BR_30$transcript_id),unique(UV_SC_30$transcript_id)),
  euler.d =T,
  scaled = T,
  category.names = c(expression(italic("S. ruralis")) , expression(italic("S. caninervis"))),
  filename  = "venn_acute_30.tiff",
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  height = 800 , 
  width = 800 , 
  resolution = 600,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c(srur, scan),
  
  # Numbers
  cex = .65,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0), 
  cat.dist = c(0.04, 0.04),
  cat.fontfamily = "sans"
) 



# what's in the venn overlap for the two species? 

overlap_10 <- UV_SC_10 %>% filter(transcript_id %in% UV_BR_10$transcript_id)
overlap_10$transcript_id
# none

overlap_30 <- UV_SC_30 %>% filter(transcript_id %in% UV_BR_30$transcript_id)
overlap_30$transcript_id
# [1] "Sc_g02351" "Sc_g05313" "Sc_g05631" "Sc_g05714" "Sc_g06791" "Sc_g08899" "Sc_g10016" "Sc_g10057" "Sc_g10629" "Sc_g10637"
# [11] "Sc_g11968" "Sc_g13670" "Sc_g13874"

write.csv(overlap_30, "tables/overlap_30_ann.csv", row.names = F)


#### ELIPs and LEAs in clusters ####

clusters_ann <- read.csv(file = "tables/clusters_ann.csv")

#### check if elips in each clusters I - VI above ####
cluster_i_elip <- clusters_ann %>% filter(cluster == "I") %>% filter(transcript_id %in% elip$transcript_id)
length(unique(cluster_i_elip$transcript_id))
# [1] 0

cluster_ii_elip <- clusters_ann %>% filter(cluster == "II") %>% filter(transcript_id %in% elip$transcript_id)
length(unique(cluster_i_elip$transcript_id))
# [1] 0

cluster_iii_elip <- clusters_ann %>% filter(cluster == "III") %>% filter(transcript_id %in% elip$transcript_id)
length(unique(cluster_i_elip$transcript_id))
# [1] 0

cluster_iv_elip <- clusters_ann %>% filter(cluster == "IV") %>% filter(transcript_id %in% elip$transcript_id)
length(unique(cluster_i_elip$transcript_id))
# [1] 0

cluster_v_elip <- clusters_ann %>% filter(cluster == "V") %>% filter(transcript_id %in% elip$transcript_id)
length(unique(cluster_i_elip$transcript_id))
# [1] 0

cluster_vi_elip <- clusters_ann %>% filter(cluster == "VI") %>% filter(transcript_id %in% elip$transcript_id)
length(unique(cluster_i_elip$transcript_id))
# [1] 0

#### check if leas in each clusters I - VI above ####
cluster_i_lea <- clusters_ann %>% filter(cluster == "I") %>% filter(transcript_id %in% lea$transcript_id)
length(unique(cluster_i_lea$transcript_id))
# [1] 0

cluster_ii_lea <- clusters_ann %>% filter(cluster == "II") %>% filter(transcript_id %in% lea$transcript_id)
length(unique(cluster_i_lea$transcript_id))
# [1] 0

cluster_iii_lea <- clusters_ann %>% filter(cluster == "III") %>% filter(transcript_id %in% lea$transcript_id)
length(unique(cluster_i_lea$transcript_id))
# [1] 0

cluster_iv_lea <- clusters_ann %>% filter(cluster == "IV") %>% filter(transcript_id %in% lea$transcript_id)
length(unique(cluster_i_lea$transcript_id))
# [1] 0

cluster_v_lea <- clusters_ann %>% filter(cluster == "V") %>% filter(transcript_id %in% lea$transcript_id)
length(unique(cluster_i_lea$transcript_id))
# [1] 0

cluster_vi_lea <- clusters_ann %>% filter(cluster == "VI") %>% filter(transcript_id %in% lea$transcript_id)
length(unique(cluster_i_lea$transcript_id))
# [1] 0

