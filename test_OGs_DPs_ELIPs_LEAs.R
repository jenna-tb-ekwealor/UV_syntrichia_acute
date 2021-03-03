library("dplyr")

setwd("/Users/jennaekwealor/Documents/dissertation_repositories/UV_syntrichia_acute")

species_effect <- read.csv(file = "tables/speciesSC.timethirty_ann.csv") %>% select(transcript_id:padj)

UV_SC <- read.csv(file = "tables/SC/time_thirty_vs_zero_ann.csv") %>% select(names:padj)

UV_BR <- read.csv(file = "tables/BR/time_thirty_vs_zero_ann.csv") %>% select(names:padj)


elip <- read.csv(file = "caninervis_elips.csv")
lea <- read.csv(file = "caninervis_leas.csv")


#### check if elips in each of above ####
species_effect_elip <- species_effect %>% filter(transcript_id %in% elip$transcript_id)
length(unique(species_effect$transcript_id))
# [1] 24
length(unique(species_effect_elip$transcript_id))
# [1] 1

UV_SC_elip <- UV_SC %>% filter(transcript_id %in% elip$transcript_id)
length(unique(UV_SC$transcript_id))
# [1] 56
length(unique(UV_SC_elip$transcript_id))
# [1] 2

UV_BR_elip <- UV_BR %>% filter(transcript_id %in% elip$transcript_id)
length(unique(UV_BR$transcript_id))
# [1] 32
length(unique(UV_BR_elip$transcript_id))
# [1] 0


#### check if leas in each of above ####
species_effect_lea <- species_effect %>% filter(transcript_id %in% lea$transcript_id)
length(unique(species_effect$transcript_id))
# [1] 24
length(unique(species_effect_lea$transcript_id))
# [1] 0

UV_SC_lea <- UV_SC %>% filter(transcript_id %in% lea$transcript_id)
length(unique(UV_SC$transcript_id))
# [1] 56
length(unique(UV_SC_lea$transcript_id))
# [1] 2

UV_BR_lea <- UV_BR %>% filter(transcript_id %in% lea$transcript_id)
length(unique(UV_BR$transcript_id))
# [1] 32
length(unique(UV_BR_lea$transcript_id))
# [1] 2


### venn diagram of UV results ####

# euler diagram of number sig mets in each treatment 
VennDiagram::venn.diagram(
  x = list(unique(UV_BR$transcript_id),unique(UV_SC$transcript_id)),
  euler.d =T,
  scaled = T,
  category.names = c(expression(italic("S. ruralis")) , expression(italic("S. caninervis"))),
  filename  = "venn_acute.tiff",
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
  fill = '#099c4e',
  
  # Numbers
  cex = .65,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0), 
  cat.dist = c(0.02, 0.02),
  cat.fontfamily = "sans"
) 



# what's in the venn overlap?

overlap <- UV_SC %>% filter(transcript_id %in% UV_BR$transcript_id)
overlap$transcript_id

