# yeast-oaktree

# Bioinformatic Pipeline associated to the Metabarcoding Project on Oak Tree

## 1 - Download the data

Raw reads were deposited on NCBI's Sequence Read Archive (SRA) database under BioProject PRJNA1114957.

## 2 - Download the software(s)

You will need: nf-core, R (Nextflow)

## 3 - Run nf-core/ampliseq

The data needs to be detailed in a samplesheet.

Our samples are named with numbers from 101 to 194, so we need to create the file Samplesheet.tsv:

> sampleID	forwardReads	reverseReads
> S101	/proj/yeast-genomics/MetaOak/DataDelivery_2024-02-13_15-39-23_ngisthlm00569/files/P30009/P30009_101/02-FASTQ/240209_VH00203_376_AAFFLMLM5/P30009_101_S1_L001_R1_001.fastq.gz	/proj/yeast-genomics/MetaOak/DataDelivery_2024-02-13_15-39-23_ngisthlm00569/files/P30009/P30009_101/02-FASTQ/240209_VH00203_376_AAFFLMLM5/P30009_101_S1_L001_R2_001.fastq.gz
> S102	/proj/yeast-genomics/MetaOak/DataDelivery_2024-02-13_15-39-23_ngisthlm00569/files/P30009/P30009_102/02-FASTQ/240209_VH00203_376_AAFFLMLM5/P30009_102_S2_L001_R1_001.fastq.gz	/proj/yeast-genomics/MetaOak/DataDelivery_2024-02-13_15-39-23_ngisthlm00569/files/P30009/P30009_102/02-FASTQ/240209_VH00203_376_AAFFLMLM5/P30009_102_S2_L001_R2_001.fastq.gz
(...)

Now we can run this script on our data:
```{bash}
nextflow run nf-core/ampliseq\
 -r dev -profile uppmax --project naiss2023-22-1116\
 --input "/your-path/MetaOak/Samplesheet.tsv"\
 --illumina_pe_its\ 
 --FW_primer TCCGTAGGTGAACCTGCGG\
 --RV_primer GCTGCGTTCTTCATCGATGC\
 --dada_ref_taxonomy unite-fungi\
 --outdir "/your-path/MetaOak/results"
```

Default setting for taxonomic classification is DADA2 with the SILVA reference taxonomy database. Instead, we use 'unite-fungi' (eukaryotic nuclear ribosomal ITS region).

Running with nf-core/ampliseq v2.9.0dev-g2a23b1b

## 4 - Analysis of the output files produced:

```{r}
ASV_tax.unite.fungi <- read.delim("/your-path/ASV_tax.unite-fungi.tsv") 
ASV_table <- read.delim("/your-path/ASV_table.tsv")
samples_infos <- read.delim("/your-path/R.Stelkens_23_01_sample_info.txt")

# Trimming with 97% confidence threshold
ASV_tax.unite.fungi_trim <- subset(ASV_tax.unite.fungi, confidence>=0.97)
# from 353 to 149
rm(ASV_tax.unite.fungi)

ASV_table <- ASV_table[ASV_table$ASV_ID %in% ASV_tax.unite.fungi_trim$ASV_ID,] 

library(tidyr)
samples_infos <- separate(samples_infos, col=NGI.ID, into=c('lib', 'ID'), sep='_')
samples_infos$ID <- paste0("S", samples_infos$ID)
samples_infos <- samples_infos[,c(2,3)]
```

We can first check where are the hits located.

The resulting abundance and taxonomy tables can be imported into R studio for further analysis.

```{r}
library(vegan)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(stringr)

setwd("/your-path/")
ASVtab = ASV_table

rownames(ASVtab) <- ASVtab[,1]
        
ASVtab <- lapply(ASVtab[,2:95], FUN = as.numeric)
ASVtab <- as.data.frame(ASVtab)

TAXtab = ASV_tax.unite.fungi_trim
metadata = samples_infos
metadata$loc <- str_split(metadata$User.ID, "_") %>% sapply(function(x) x[1])

# Remove pools
metadata_clean <- metadata[!grepl("POOL", metadata$User.ID), ]
rm(metadata)
rm(samples_infos)

ASVtab = t(ASVtab)
rownames(ASVtab) <- samples_infos$User.ID

# number of reads (max and min), number of ASVs per sample, and a quick overview of the taxonomy of the reads (at phylum level).
ASVtab <- as.data.frame(ASVtab)
sort(rowSums(ASVtab)) # Number of reads - between 0 and 489,676 per sample
sort(rowSums(ASVtab != 0)) # Number of ASVs - between 0 and 26 per sample
# Between 1 and 26 ASVs detected per tree
# 3 trees with no ASVs detected = BLA_1,  BJO_2,  VAR1_2 

# remove the pools
ASVtab <- ASVtab[!grepl("POOL", rownames(ASVtab)), ]
```

### Diversity indexes calculation

Our reads have different sequencing depth. So, to avoid problems connected to differences in sequencing depth, we need to rarefy the table. For a couple of samples, the sequencing depth is way lower that 120k:

S136_1.filt.fastq.gz S183_1.filt.fastq.gz
               20351                36977   

We are therefore filtering out those samples before choosing a sequencing depth of 121000, since it is a very similar depth to the samples with the lowest amount of reads.

```{r}
ASVtab1 = readRDS("/your-path/DADA2_table.rds")

sort(rowSums(ASVtab1))
ASVtab1 <- as.data.frame(ASVtab1)

# Remove pools
list <- c("S188_1.filt.fastq.gz","S187_1.filt.fastq.gz","S180_1.filt.fastq.gz","S179_1.filt.fastq.gz",
          "S172_1.filt.fastq.gz","S164_1.filt.fastq.gz","S156_1.filt.fastq.gz","S148_1.filt.fastq.gz",
          "S140_1.filt.fastq.gz","S132_1.filt.fastq.gz","S124_1.filt.fastq.gz","S116_1.filt.fastq.gz",
          "S108_1.filt.fastq.gz")

ASVtab2 <- ASVtab1[!(rownames(ASVtab1) %in% list),]
ASVtab2 <- ASVtab2[(colnames(ASVtab2) %in%
                      ASV_tax.unite.fungi_trim$sequence)] # 137

# Select columns where the sum is not zero
col_sums <- colSums(ASVtab1, na.rm = TRUE)
ASVtab2_non_empty <- ASVtab2[rowSums(ASVtab2) > 0, ] # 78 / 137

sort(rowSums(ASVtab2_non_empty)) # All >20k

metadata_clean$richness = vegan::specnumber(ASVtab2) # calculate ASV richness
metadata_clean$shannon = vegan::diversity(ASVtab2, index = "shannon")
metadata_clean$evenness = metadata_clean$shannon / log(metadata_clean$richness)
metadata_clean$simpson <- vegan::diversity(ASVtab2, "simpson")

# Akso add the true number of species per sample
species_tab = aggregate(t(ASVtab) ~ TAXtab$Species, FUN = "sum") # Sum all the ASV together in the same sample
species_tab_melted = melt(species_tab, value = species_tab$`TAXtab$Species`) ## This melt the table

species_tab_clean <- subset(species_tab_melted, value != 0)
species_tab_clean$value_unique <- 1
species_tab_clean <- species_tab_clean %>%
  group_by(variable) %>%
  summarise(sp_nb = sum(value_unique)) 

metadata_clean <- merge(metadata_clean, species_tab_clean, by.x="User.ID", by.y="variable")
```


### Add environmental and tree-related data

```{r}
library(xlsx)
tree <- read.xlsx("/your-path/bark_2024_yeasts.xlsx", 
                  sheetName = "tree-sp")
tree$Yeast_Individuals <- paste0("S1",tree$Yeast_Individuals)
list_QuPe <- tree$Yeast_Individuals[tree$Tree_species=="QUPE"]
list_QuRo <- tree$Yeast_Individuals[tree$Tree_species=="QURO"]

metadata_clean$tree_sp <- "NA"
metadata_clean$tree_sp[metadata_clean$ID %in% list_QuPe] <- "QuPe"
metadata_clean$tree_sp[metadata_clean$ID %in% list_QuRo] <- "QuRo"

# Add other ecological factors 
tree_eco_fac <- read.xlsx("/your-path/bark_2024_yeasts.xlsx", 
                  sheetName = "eco-fac")
tree_eco_fac$Yeast_Individuals <- paste0("S1",tree_eco_fac$Yeast_Individuals)

metadata_clean <- merge(metadata_clean, tree_eco_fac, by.x="ID", by.y="Yeast_Individuals", all.x=T)

metadata_clean$island <- "no"
metadata_clean$island[metadata_clean$loc %in% c("FH","GAR","BLA","RAL")] <- "yes"

# Add longitude and latitude
envir_metrics <- read.xlsx("/your-path/bark_2024_yeasts.xlsx", 
                  sheetName = "unique_loc")
envir_metrics <- envir_metrics[,c(1,3:4)]

metadata_clean <- merge(metadata_clean, envir_metrics, by.x="loc", by.y="SiteCode", all.x=T)
```


## 5 - Plotting the data - main text

### Figure 3.A. Percentage of ASVs mapping to each genus across all samples.

```{r}
genus_tab = aggregate(t(ASVtab) ~ TAXtab$Genus, FUN = "sum") # Sum all the ASV together in the same sample
genus_tab_melted = melt(genus_tab, value = genus_tab$`TAXtab$Genus`) 
rownames(genus_tab) = c(genus_tab$`TAXtab$Genus`)
genus_tab = genus_tab[,2:82]
total_ab_genus = rowSums(genus_tab)/sum(rowSums(genus_tab))
total_ab_genus = total_ab_genus*100
names(total_ab_genus) = rownames(genus_tab)
total_ab_genus <- as.data.frame(total_ab_genus)
total_ab_genus$names <- row.names(total_ab_genus)
total_ab_genus <- subset(total_ab_genus, names!="NA")

library(dplyr)
genus_identified <- genus_tab_melted %>%
  filter(value > 0) %>%          # Keep only rows where value is greater than 0
  group_by(`TAXtab$Genus`) %>%          # Group by Species
  summarise(places_identified = n_distinct(variable))  # Count distinct places (variable) where species is identified

num_colors <- 20
custom_palette <- colorRampPalette(brewer.pal(11, "Spectral"))(num_colors)

total_ab_genus$formatted_percentage <- ifelse(total_ab_genus$total_ab_genus < 0.1, "<0.01",
                                              round(total_ab_genus$total_ab_genus, 2))
lbls <- paste(total_ab_genus$names, " (",total_ab_genus$formatted_percentage, "%)", sep = "")

lbls <- sapply(1:nrow(total_ab_genus), function(i) {
  bquote(italic(.(total_ab_genus$names[i])) ~ "(" ~ .(total_ab_genus$formatted_percentage[i]) ~ "%)")
})

ggplot(total_ab_genus, aes(x="", y=total_ab_genus, fill=names)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)+
  scale_fill_manual(values = scales::hue_pal()(length(unique(total_ab_genus$names))), 
                    labels = lbls) +
  theme_void()+
  theme(legend.text = element_text()) +  # Adjust text size if needed
  labs(fill = "Genus")
```


### Figure 3.B. Projection of samples using non-metric multidimensional scaling (NMDS).

Ordination is based on ASVs richness with colours indicating the dominant genus. Outlier samples LANY-1, TAN-1 and VASV-1 were removed for display purposes. 

```{r}
# Find the most represented genus in number of ASV reads for each samples
library(dplyr)
colnames(genus_tab_melted) <- c("Genus","Sample","value")
max_value_rows <- genus_tab_melted %>%
  group_by(Sample) %>%
  top_n(1, value) %>%
  ungroup()
max_value_rows <- subset(max_value_rows, value != 0)

metadata_clean <- merge(metadata_clean, max_value_rows[,c(1,2)], 
                        by.x=c("User.ID"), by.y=c("Sample"), all.x=T)

head(ASVtab)
asv_matrix <- as.matrix(ASVtab)
asv_matrix <- na.omit(asv_matrix)
asv_matrix <- asv_matrix[rowSums(asv_matrix) > 0, ]
asv_matrix <- subset(asv_matrix, !(rownames(asv_matrix) %in% c("LANY_1","TAN_1","VASV_1"))) 

filtered_data <- subset(metadata_clean, sp_nb > 0)
filtered_data <- subset(filtered_data, sp_nb!='NA')

custom_palette <- c("#e68511ff","#cc9400ff","#7aae00ff","#0ab500ff", "#00bd66ff", "#00bf99ff",
                    "#00bdc2ff","#00b8e6ff","#8294ffff", "#c77affff","#ff61ccff","#ff66a1ff")

filtered_data$Genus <- as.factor(filtered_data$Genus) # 10 levels

genus_colors <- setNames(custom_palette, levels(filtered_data$Genus))
filtered_data <- filtered_data %>%
  mutate(Genus_Color = genus_colors[Genus])

colnames(filtered_data)[1] <- "Sample"
filtered_data <- subset(filtered_data,!(Sample %in% c("LANY_1","TAN_1","VASV_1")))

nmds <- metaMDS(asv_matrix, distance = "bray", k = 2, trymax = 1000)
# Run 357 stress 0.1707845 
# ... Procrustes: rmse 0.001342651  max resid 0.006971928 

nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$Sample <- rownames(nmds_scores)

plot_data <- left_join(nmds_scores, filtered_data, by = "Sample")

p <- ggplot(plot_data, aes(x = NMDS1, y = NMDS2, color = Genus_Color, size = sp_nb)) +
  geom_point(alpha = 0.7) +
  scale_color_identity() +
  theme_minimal() +
  labs(title = '2D NMDS Ordination of ASVs', x = 'NMDS1', y = 'NMDS2')
p
```


### Figure 4. Chord diagram depicting yeast species co-occurrence, based on shared presence across sampled trees. 
The connections (arcs) between species indicate how frequently they co-occur, with the thickness of the arcs representing the strength of co-occurrence ( thicker lines corresponding to more frequent co-occurrence). 

```{r}
library(vegan)
library(cooccur)
library(corrplot)
library(dplyr)
library(igraph)

merged_data <- merge(TAXtab, ASV_table, by="ASV_ID")
merged_data <- merged_data[, c(8,11:104)]

species_counts <- merged_data %>%
  group_by(Species) %>%
  summarise(across(starts_with("S"), sum, na.rm = TRUE))

id_to_userid <- setNames(metadata_clean$User.ID, metadata_clean$ID)
colnames(species_counts)[-1] <- id_to_userid[colnames(species_counts)[-1]]

species_counts <- species_counts[, !is.na(colnames(species_counts)) & colnames(species_counts) != ""]
species_counts <- species_counts[-1,] # remove unknown species
species_counts <- as.data.frame(species_counts)
rownames(species_counts) <- species_counts$Species
species_counts <- species_counts[,-1] # remove species column

presence_absence <- species_counts
presence_absence[presence_absence > 0] <- 1
presence_absence_clean <- presence_absence[rowSums(presence_absence) != 0, ]

presence_absence_binary <- as.data.frame(sapply(presence_absence_clean, function(x) ifelse(x > 0, 1, 0)))

# Create a species x species co-occurrence matrix
species_names <- rownames(presence_absence_clean)
co_occurrence_matrix <- matrix(0, nrow = length(species_names), ncol = length(species_names),
                               dimnames = list(species_names, species_names))

# Fill the co-occurrence matrix with counts
for (i in 1:nrow(presence_absence_clean)) {
  for (j in 1:nrow(presence_absence_clean)) {
    if (i != j) {
      # Count the number of times both species are present and sum their counts
      co_occurrence_matrix[i, j] <- sum(presence_absence_clean[i, ] * presence_absence_clean[j, ])
    }
  }
}

# Convert the matrix to a data frame for better readability
co_occurrence_df <- as.data.frame(co_occurrence_matrix)

library(reshape2)

# Convert the co-occurrence matrix to a long format
co_occurrence_long <- melt(co_occurrence_matrix)
co_occurrence_long$Var1 <- gsub("_", " ", co_occurrence_long$Var1)
co_occurrence_long$Var2 <- gsub("_", " ", co_occurrence_long$Var2)

library(network)
library(sna)

co_occurrence_pairs <- subset(co_occurrence_long, value >= 1)

species_list <- unique(c(co_occurrence_pairs$Var1, co_occurrence_pairs$Var2))
n_species <- length(species_list)

# Create an adjacency matrix
adj_matrix <- matrix(0, nrow = n_species, ncol = n_species,
                     dimnames = list(species_list, species_list))

# Fill the adjacency matrix with values
for (i in 1:nrow(co_occurrence_pairs)) {
  adj_matrix[co_occurrence_pairs$Var1[i], co_occurrence_pairs$Var2[i]] <- co_occurrence_pairs$value[i]
}

# Create a network object
net <- network(adj_matrix, directed = FALSE)

# Assign edge attributes (co-occurrence values)
for (i in 1:n_species) {
  for (j in 1:n_species) {
    if (adj_matrix[i, j] > 0) {
      net[i, j] <- adj_matrix[i, j]
    }
  }
}

library(circlize)

species_colors <- c(
  "Candida castellii" = "#ff6b66",
  "Citeromyces matritensis" = "#f47f00",
  "Debaryomyces delbrueckii" = "#d59300", #ok
  "Hanseniaspora osmophila" = "#6ab100",
  "Hyphopichia burtonii" = "#00ba00", #ok
  "Issatchenkia orientalis" = "#00c25b", #ok
  "Kazachstania servazzii" = "#00c397", #ok
  "Kluyveromyces dobzhanskii" = "#00bdc2ff", #ok
  "Kluyveromyces lactis" = "#00c0c4",
  "Kregervanrija fluxuum" = "#00baeb", #ok
  "Lachancea kluyveri" = "#00a8ff",
  "Ogataea dorogensis" = "#8091ff", #ok
  "Ogataea nonfermentans" = "#6C7AD8",
  "Pichia mandshurica" = "#9545B7",
  "Pichia membranifaciens" = "#d570ff", #ok
  "Pichia sp" = "#AF61D0",
  "Saccharomyces paradoxus" = "#ff4cd1",
  "Saccharomycodes ludwigii" = "#ff56a2"  #ok
)

hist(adj_matrix[adj_matrix > 0], breaks = 20, 
     main = "Distribution of Co-occurrence Values", xlab = "Co-occurrence Value")
threshold <- quantile(adj_matrix[adj_matrix > 0], 0.75)

circos.clear()  # Clear the previous plot

circos.par(
  gap.after = rep(2, length(species_names)),  # Adds a small gap between each species
  track.margin = c(0, 0),               # Increase the margin around the circle
  canvas.xlim = c(-1.5, 1.5),                 # Adjust canvas limits for more room
  canvas.ylim = c(-1.5, 1.5)                  # Same for Y axis, giving more space
)

# Re-draw the chord diagram
chordDiagram(
  adj_matrix, 
  transparency = 0.5, 
  annotationTrack = "grid", 
  preAllocateTracks = 1, 
  grid.col = species_colors  # Pass custom colors
)

# Add species names, avoiding overlap
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    circos.text(
      CELL_META$xcenter, 
      CELL_META$ylim[1] + 0.1,  # Move labels slightly outwards
      sector.name, 
      facing = "clockwise", 
      niceFacing = TRUE, 
      adj = c(0, 0.5), 
      cex = 0.5  # Adjust the text size to fit
    )
  }
)

# Run statistics on co-occurrence 
cooccurrence_results <- cooccur(mat = presence_absence_clean, type = "spp_site")
summary(cooccurrence_results)
results <- cooccurrence_results$results

avoiding_pairs <- subset(results, p_gt >= 0.05 & p_lt < 0.05) 
# 16/17: "Pichia_sp" and "Saccharomyces_paradoxus"
# 3/8: "Debaryomyces_delbrueckii" and "Kluyveromyces_dobzhanskii"
cooccuring_pairs <- subset(results, p_lt >= 0.05 & p_gt < 0.05) 
# 3/17: "Debaryomyces_delbrueckii" and "Saccharomyces_paradoxus"
```


### Figure 5.A. [Map of sampling locations with] pie charts showing the proportion of ASVs found per yeast genus. 

```{r}
colnames(genus_tab_melted)[1:2] = c("Genus", "Sample")

# Remove pools
genus_tab_melted$group = metadata_clean$loc[match(genus_tab_melted$Sample,
                                                  metadata_clean$User.ID)]

# Save all the plots separately
loc <- read.xlsx("/your-path/bark_2024_yeasts.xlsx", sheetName = "temperature")
loc <- loc[,2:5]
loc$field_2[loc$field_1=="STH"] <- "Strömsholden"
loc$field_2[loc$field_1=="STR"] <- "Strömsrum"

genus_tab_melted <- merge(genus_tab_melted, loc, by.x="group", by.y="field_1")
samples <- unique(genus_tab_melted$field_2)

df = data.frame()
groups <- unique(genus_tab_melted$group)

for (grup in groups) {
  sub_tab <- subset(genus_tab_melted, group == grup)
  sub_tab$percent <- sub_tab$value/sum(sub_tab$value)
  sub_tab$percent = sub_tab$percent*100
  df <- rbind(df, sub_tab)
}

# mean of percentage, not cumulated nb 
# (to account for the fact that different samples have different number of total ASV)
df2 <- df %>%
  dplyr::group_by(group,Genus,field_2) %>%
  dplyr::summarize(percent2 = mean(percent), val2 = sum(value))

# Loop through each sample, plot, and save as SVG
for (sample in samples) {

   sample_data <- subset(df2, field_2 == sample)
    
   p2 <- ggplot(sample_data, aes(x="", y=percent2, fill=Genus)) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0)+
      theme_void()+
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),  
            axis.title.x = element_blank(),  axis.title.y = element_blank(), legend.position = "none") +
      ggtitle(unique(sample_data$field_2))
    
    ggsave(filename = paste0("Pie_", sample, ".svg"), plot = p2, device = "svg", width = 2, height = 2)
}
```


### Figure 5.B. Longitude correlated with the abundance of the three most common genera isolated across all locations (Kluyveromyces, Saccharomyces, and Pichia).

```{r}
head(metadata_clean)
head(TAXtab)
head(ASVtab)

genus_tab = aggregate(t(ASVtab) ~ TAXtab$Genus, FUN = "sum") # Sum all the ASV together in the same sample
maj_genus_tab <- subset(genus_tab, 
                       `TAXtab$Genus` %in% c("Kluyveromyces","Saccharomyces","Pichia"))

library(reshape2)
library(dplyr)

maj_genus_tab_long <- melt(maj_genus_tab, id.vars = 'TAXtab$Genus', 
                           variable.name = 'User.ID', value.name = 'value')
colnames(maj_genus_tab_long) <- c('Genus', 'User.ID', 'value')

merged_data <- merge(maj_genus_tab_long, subset(metadata_clean, select = -Genus), by = 'User.ID')
merged_data <- subset(merged_data, value!=0)

genus_sum_by_loc <- merged_data %>%
  dplyr::group_by(loc, Genus, Long) %>%
  dplyr::summarise(total_abundance = sum(value, na.rm = TRUE))

ggplot(genus_sum_by_loc, aes(x = Long, y = total_abundance, 
                             color = Genus, shape = Genus)) +
  geom_point(alpha = 0.7, size=2) +
  geom_smooth(method = "lm", se = FALSE) +  
  scale_color_manual(values = c("#00bdc2ff", "#c77affff", "#ff61ccff")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_y_continuous(labels = scales::scientific) +  
  labs(x = "Longitude", 
       y = "Total number of ASVs per genus") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"), 
    panel.grid = element_blank()               
  )

# Or formula = y ~ poly(x, 3),

# Linear model
library(dplyr)
genus_results <- data.frame(Genus = character(), R_squared = numeric(), P_value = numeric(), stringsAsFactors = FALSE)

for (genus in unique(genus_sum_by_loc$Genus)) {
  subset_data <- genus_sum_by_loc %>% filter(Genus == genus)
  model <- lm(total_abundance ~ Long, data = subset_data)
  
  model_summary <- summary(model)
  r_squared <- model_summary$r.squared
  p_value <- model_summary$coefficients[2, 4]  
  
  genus_results <- rbind(genus_results, data.frame(Genus = genus, R_squared = r_squared, P_value = p_value))
}
print(genus_results)

simple_lm <- lm(total_abundance ~  Genus * Long , data = genus_sum_by_loc)
summary(simple_lm) # p-value: 0.05198

slope_comparisons <- emtrends(simple_lm, pairwise ~ Genus, var = "Long")
summary(slope_comparisons$contrasts) # Kluyveromyces - Pichia: p-value = 0.0095

# Third order polynomial model
library(dplyr)
genus_results <- data.frame(Genus = character(), R_squared = numeric(), P_value = numeric(), stringsAsFactors = FALSE)

for (genus in unique(genus_sum_by_loc$Genus)) {
  subset_data <- genus_sum_by_loc %>% filter(Genus == genus)
  model <- lm(total_abundance ~ poly(Long, 3), data = subset_data)
  
  model_summary <- summary(model)
  r_squared <- model_summary$r.squared
  p_value <- model_summary$coefficients[2, 4] 
  
  genus_results <- rbind(genus_results, data.frame(Genus = genus, R_squared = r_squared, P_value = p_value))
}
print(genus_results)

third_order_model <- lm(total_abundance ~ poly(Long, 3) * Genus, 
                      data = genus_sum_by_loc) 
summary(third_order_model) # p-value: 0.04379

library(emmeans)
slope_comparisons <- emtrends(third_order_model, pairwise ~ Genus, var = "Long")
summary(slope_comparisons$contrasts) # Kluyveromyces - Pichia: p-value = 0.0095
```


### Figure 6. Correlation patterns between environmental variables, tree metrics and yeast diversity indices.
Significance levels of Pearson correlations are indicated by asterisks (p < 0.05*, 0.01**, and 0.001***).

```{r}
# Create a binomial score for island and tree species
non_pool_data$island_num[non_pool_data$island=="yes"] <- 1
non_pool_data$island_num[non_pool_data$island=="no"] <- 0

non_pool_data$tree_sp_num[non_pool_data$tree_sp=="QuRo"] <- 0
non_pool_data$tree_sp_num[non_pool_data$tree_sp=="QuPe"] <- 1

library(Hmisc)
library(ggplot2)
library(reshape2)

cor_results <- rcorr(as.matrix(non_pool_data[,c("richness","shannon","evenness","simpson","sp_nb",
                                                "Long","Lat","island_num","tree_sp_num",
                                                "Tree_diameter_cm","Tree_trunk_radius_mm",
                                                "Barkdepth_mm","Tree_height_m",
                                                "Tree_growth_mm_year","Tree_age")]), type = "pearson")

cor_matrix <- cor_results$r   
p_matrix <- cor_results$P
names <- c("Richness", "Shannon", "Evenness", "simpson","n species",
           "Longitude","Latitude","Insularity","Tree species",
           "Tree diameter","Tree trunk radius","Bark depth","Tree height",
           "Tree growth","Tree age")

rownames(cor_matrix) <- names
colnames(cor_matrix) <- names

rownames(p_matrix) <- names
colnames(p_matrix) <- names

cor_long <- melt(cor_matrix)
p_long <- melt(p_matrix)

cor_data <- merge(cor_long, p_long, by = c("Var1", "Var2"))
colnames(cor_data) <- c("Var1", "Var2", "correlation", "p_value")

cor_data$significance <- cut(cor_data$p_value,
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", ""))
cor_data$significance[is.na(cor_data$p_value)] <- ""

cor_data$label <- paste0(round(cor_data$correlation, 2), cor_data$significance)

library(reshape2)
ggheatmap <- ggplot(cor_data, aes(Var2, Var1, fill = correlation))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#101282", high = "#9c1e28", mid = "white", 
   midpoint = 0, limit = c(-1,1), name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
 theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
 coord_fixed()

ggheatmap <- ggheatmap + 
geom_text(aes(Var2, Var1, label = round(correlation, 2)), 
          color = "black", size=3) +
theme(
  axis.title.x = element_blank(), axis.title.y = element_blank(),
  panel.grid.major = element_blank(), panel.border = element_blank(),
  panel.background = element_blank(), axis.ticks = element_blank(),
  legend.justification = c(1, 0), legend.position = c(0.6, 0.7),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5))

ggheatmap + 
  geom_text(aes(Var2, Var1, label = significance), 
            color = "black", size = 3, vjust = 2)
```


### Figure 7.A. Relative abundance of yeast genera depending on host tree species.

```{r}
data_tree <- subset(non_pool_data, non_pool_data$tree_sp!="NA")

genus_tab_melted$Sample <- as.character(genus_tab_melted$Sample)
genus_tab_melted <- merge(genus_tab_melted, data_tree[,c("User.ID","tree_sp")], by.x="Sample", by.y="User.ID")

ggplot(genus_tab_melted,aes(x = tree_sp, y = value, fill = Genus))+
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank()) +
  ylab("Relative abundance")+
  scale_y_continuous(breaks = c(0, 0.5, 1))
```


### Figure 7.B. Relative abundance of yeast genera depending on insularity.

```{r}
genus_tab_melted <- merge(genus_tab_melted, non_pool_data[,c("User.ID","island")], by.x="Sample", by.y="User.ID")
ggplot(genus_tab_melted,aes(x = island, y = value, fill = Genus))+
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank()) + # to remove the 'ticks' on x-axis
  ylab("Relative abundance")+
  scale_y_continuous(breaks = c(0, 0.5, 1))
```


## 6 - Plotting the data - supplementary

### Supplementary Figure 2. Percentage of ASVs mapping to each family across all samples.

```{r}
# Do a plot to see which families are represented
family_tab = aggregate(t(ASVtab) ~ TAXtab$Family, FUN = "sum") # Sum all the ASV together in the same sample
family_tab_melted = melt(family_tab, value = family_tab$`TAXtab$Family`) 
rownames(family_tab) = family_tab$`TAXtab$Family`
family_tab = family_tab[,2:82]
total_ab_family = rowSums(family_tab)/sum(rowSums(family_tab))
total_ab_family = total_ab_family*100
names(total_ab_family) = rownames(family_tab)
total_ab_family <- as.data.frame(total_ab_family)
total_ab_family$names <- row.names(total_ab_family)

num_colors <- 7
library(RColorBrewer)
custom_palette <- colorRampPalette(brewer.pal(9, "Spectral"))(num_colors)

lbls <- paste(total_ab_family$names, " (",
              round(total_ab_family$total_ab_family, 2), "%)", sep = "")

ggplot(total_ab_family, aes(x="", y=total_ab_family, fill=names)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)+
  scale_fill_manual(values = scales::hue_pal()(length(unique(total_ab_family$names))), labels = lbls) +
  theme_void()+
  theme(legend.text = element_text()) +
  labs(fill = "Family")
```


### Supplementary Figure 3. Number of species in the three clusters that group the dominant genera.  
Boxplot extremities represent minimum and maximum values, whereas the box itself is composed of the first quartile, median (thick line), and third quartile. A jitter effect was added for better visibility of data points.

```{r}
head(metadata_clean)
maj_genus_metadata <- subset(metadata_clean, Genus %in% c("Kluyveromyces","Saccharomyces","Pichia"))

ggplot(data = maj_genus_metadata, aes(x = Genus, y = sp_nb, color = Genus)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  scale_color_manual(values = c("#00bdc2ff", "#c77affff","#ff61ccff"))+
  geom_jitter() +
  ylab("Number of yeast species") + 
  theme_bw()    

median(maj_genus_metadata$sp_nb[maj_genus_metadata$Genus=="Saccharomyces"])
median(maj_genus_metadata$sp_nb[maj_genus_metadata$Genus=="Kluyveromyces"])
median(maj_genus_metadata$sp_nb[maj_genus_metadata$Genus=="Pichia"])

kruskal_test <- kruskal.test(sp_nb ~ Genus, data = maj_genus_metadata)
print(kruskal_test) # Kruskal-Wallis chi-squared = 1.6876, df = 2, p-value = 0.4301
```


### Supplementary Figure 4. Co-occurrence heatmap between yeast species detected across all trees and sites.

```{r}
# Create the heatmap using ggplot
ggplot(data = co_occurrence_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("lightgrey", "#FFEAEB", "#AA2F33"), 
                       values = scales::rescale(c(0, 1, max(co_occurrence_long$value))),
                       breaks = c(0, 1, max(co_occurrence_long$value)),
                       labels = c("0", "1", "34")) + # max(co_occurrence_long$value)
  geom_tile(data = co_occurrence_long[co_occurrence_long$Var1 == co_occurrence_long$Var2, ],
            aes(x = Var1, y = Var2), fill = "grey") +
  geom_text(aes(label = value), color = "black", size = 3) + # Adjust size if needed
  theme_minimal() +
  labs(title = "Co-occurrence Heatmap", x = "Species", y = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        axis.text.y = element_text(face = "italic"))
```


### Supplementary Figure 5. Mean temperature in August along a south-Sweden longitudinal gradient. 
Non-linear regression using a third order polynomial model (R2 =0.159, p-value=xx). Climatic data from Meineri & Hylander (2017).

```{r}
library(xlsx)
library(tidyr)
library(ggplot2)

climatic_data <- read.xlsx("/your-path/bark_2024_yeasts.xlsx",sheetName = "temperature")

metadata_clean <- merge(metadata_clean, climatic_data, by.x="loc", by.y="field_1", all.x=T)

ggplot(data = metadata_clean, aes(x = Long,
                                  y = August.mean)) +
  geom_point(alpha = 0.5, color="darkgrey") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, color="grey") +
  ylim(15,18)+
  ylab("Mean temperature in August (°C)") +
  xlab("Longitude") +
  theme_bw()

third_order_model <- lm(August.mean ~ poly(Long, 2), data = metadata_clean) 
summary(third_order_model) # p-value: 0.001138
```

