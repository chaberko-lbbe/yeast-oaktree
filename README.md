# yeast-oaktree

Bioinformatic Pipeline associated to the paper "Fermentative yeast diversity at the northern range limit of their oak tree hosts" (Pinto, Haberkorn et al., 2025)
Available on BioRxiv: https://www.biorxiv.org/content/10.1101/2025.01.30.635702v2

[Update 16/04/2025: the associated files ASV_table.tsv, ASV_tax.unite-fungi.tsv, DADA2_table.rds and metadata_clean.xlsx will be shared here upon publication of the paper]

Javier Pinto<sup>†1</sup>, Chloé Haberkorn<sup>†1</sup>, Markus Franzén<sup>2</sup>, Ayco J. M. Tack<sup>3</sup>, Rike Stelkens<sup>1*</sup>

<sup>†</sup>J.P. and C.H. contributed equally to this work ;
<sup>1</sup>Department of Zoology, Stockholm University, Svante Arrheniusväg 18 B, 106 91 Stockholm, Sweden ;
<sup>2</sup>Center for Ecology and Evolution in Microbial Model Systems, Department of Biology and Environmental Science, Linnaeus University, Kalmar, Sweden ;
<sup>3</sup>Department of Ecology, Environment and Plant Sciences, Stockholm University, SE-106 91 Stockholm, Sweden ;
<sup>*</sup>Corresponding author: Rike Stelkens, Department of Zoology, Stockholm University, Sweden; email: rike.stelkens@zoologi.su.se


## 1 - Download the data

Raw reads were deposited on NCBI's Sequence Read Archive (SRA) database under BioProject PRJNA1114957:


## 2 - Download the software(s)

You will need: nf-core, R (Nextflow)


## 3 - Run nf-core/ampliseq

The data needs to be detailed in a samplesheet.

Our samples are named with numbers from 101 to 194
So, we are creating a file with those numbers:
```{bash}
seq -f "%02g" 01 94 > List.txt # 01 to 94
```

Then, we create the file Samplesheet.tsv:
```{python}
input_file = 'List.txt'
output_file = 'Samplesheet.tsv'

with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
    for line in in_file:
        line = line.strip()
        output_line = "S1{}\t/your-path/MetaOak/DataDelivery_2024-02-13_15-39-23_ngisthlm00569/files/P30009/P30009_1{}/02-FASTQ/240209_VH00203_376_AAFFLMLM5/P30009_1{}_S{}_L001_R1_001.fastq.gz\t/your-path/MetaOak/DataDelivery_2024-02-13_15-39-23_ngisthlm00569/files/P30009/P30009_1{}/02-FASTQ/240209_VH00203_376_AAFFLMLM5/P30009_1{}_S{}_L001_R2_001.fastq.gz".format(line, line, line, line, line, line,line,line)
        out_file.write(output_line + '\n')
```

Giving something like that:
> /your-path/MetaOak/DataDelivery_2024-02-13_15-39-23_ngisthlm00569/files/P30009/P30009_101/02-FASTQ/240209_VH00203_376_AAFFLMLM5/P30009_101_S01_L001_R1_001.fastq.gz

We have to mannually correct samples 1 to 9!
> /your-path/MetaOak/DataDelivery_2024-02-13_15-39-23_ngisthlm00569/files/P30009/P30009_101/02-FASTQ/240209_VH00203_376_AAFFLMLM5/P30009_101_S1_L001_R1_001.fastq.gz

Add this line on top (using nano):
sampleID	forwardReads	reverseReads

Now we can create and run the script with the true dataset:
```{bash}
nextflow run nf-core/ampliseq\
 -r dev -profile your-profile --project your-project-id\
 --input "/your-path/MetaOak/Samplesheet.tsv"\
 --illumina_pe_its\ 
 --FW_primer TCCGTAGGTGAACCTGCGG\
 --RV_primer GCTGCGTTCTTCATCGATGC\
 --dada_ref_taxonomy unite-fungi\
 --outdir "/your-path/MetaOak/results"
```

Default setting for taxonomic classification is DADA2 with the SILVA reference taxonomy database. Instead, we use 'unite-fungi' (eukaryotic nuclear ribosomal ITS region)

Running with nf-core/ampliseq v2.9.0dev-g2a23b1b

## 4 - Analysis of the output files produced:

First, load R packages:

```{r}
library(circlize)
library(cooccur)
library(corrplot)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Hmisc)
library(igraph)
library(network)
library(openxlsx)
library(RColorBrewer)
library(reshape2)
library(sna)
library(stringr)
library(tidyr)
library(vegan)
```

```{r}
ASV_tax.unite.fungi <- read.delim("/your-path/ASV_tax.unite-fungi.tsv") 
ASV_table <- read.delim("/your-path/ASV_table.tsv")
samples_infos <- read.delim("/your-path/R.Stelkens_23_01_sample_info.txt")

# Trimming with 97% confidence threshold
ASV_tax.unite.fungi_trim <- subset(ASV_tax.unite.fungi, confidence>=0.97)
# from 353 to 149
rm(ASV_tax.unite.fungi)

ASV_table <- ASV_table[ASV_table$ASV_ID %in% ASV_tax.unite.fungi_trim$ASV_ID,] 

samples_infos <- separate(samples_infos, col=NGI.ID, into=c('lib', 'ID'), sep='_')
samples_infos$ID <- paste0("S", samples_infos$ID)
samples_infos <- samples_infos[,c(2,3)]
```

We can first check where are the hits located.

The resulting abundance and taxonomy tables can be imported into R studio for further analysis.
Also see this tutorial: https://rpubs.com/mrgambero/taxa_alpha_beta

```{r}
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

# The metrics below will be all provided in a clean table, but just to show how it was computed
metadata_clean$richness = vegan::specnumber(ASVtab2) # calculate ASV richness
metadata_clean$shannon = vegan::diversity(ASVtab2, index = "shannon")
metadata_clean$evenness = metadata_clean$shannon / log(metadata_clean$richness)
metadata_clean$simpson <- vegan::diversity(ASVtab2, "simpson")

# Also add the true number of species per sample
ASVtab_t <- as.data.frame(t(ASVtab))
ASVtab_t$ASV_ID <- rownames(ASVtab_t)
ASVtab_tax <- merge(ASVtab_t, TAXtab[, c("ASV_ID", "Species")], by = "ASV_ID")

species_tab = aggregate(t(ASVtab) ~ TAXtab$Species, FUN = "sum") # Sum all the ASV together in the same sample
species_tab_melted <- melt(species_tab, id.vars = "Species", 
                           variable.name = "Sample", value.name = "Abundance")
species_tab_melted <- species_tab_melted %>%
  filter(Species != "" & !is.na(Species))

species_tab_clean <- species_tab_melted %>%
  filter(Abundance > 0) %>%
  group_by(Sample) %>%
  summarise(unique_species_count = n())

# Final output: each sample with its count of unique species
head(species_tab_clean)
```

### Load samples information
Containing samples ID, location data (longitude/latitude), diversity indexes, information on host trees, on insularity, on temperature and precipitation data, and dominant genus

```{r}
metadata_clean <- read.xlsx("/your-path/metadata_clean_160425.xlsx")
```


## 5 - Plotting the data - main text

### Figure 2.A. Percentage of ASVs mapping to each genus across all samples.

```{r}
ASVtab_t <- t(ASVtab)
ASV_with_genus <- cbind(Genus = TAXtab$Genus[match(rownames(ASVtab_t), TAXtab$ASV_ID)], ASVtab_t)
ASV_with_genus <- as.data.frame(ASV_with_genus)

ASV_with_genus_clean <- ASV_with_genus %>%
  mutate(across(.cols = -Genus, .fns = ~ as.numeric(.))) %>%
  mutate(Genus = as.character(Genus))

genus_tab <- aggregate(. ~ Genus, data = ASV_with_genus_clean, FUN = sum)

rownames(genus_tab) <- genus_tab$Genus
genus_tab <- genus_tab[, -1]  

total_ab_genus <- rowSums(genus_tab) / sum(rowSums(genus_tab)) * 100
total_ab_genus <- as.data.frame(total_ab_genus)
total_ab_genus$names <- rownames(total_ab_genus)

total_ab_genus <- subset(total_ab_genus, names != "NA")
genus_tab <- genus_tab[rownames(genus_tab) %in% total_ab_genus$names, ]

genus_tab$Genus <- rownames(genus_tab)
genus_tab_melted <- melt(genus_tab, id.vars = "Genus", variable.name = "Sample", value.name = "Abundance")

genus_identified <- genus_tab_melted %>%
  filter(Abundance > 0) %>%
  group_by(Genus) %>%
  summarise(places_identified = n_distinct(Sample))

total_ab_genus$formatted_percentage <- ifelse(total_ab_genus$total_ab_genus < 0.1, "<0.01",
                                               round(total_ab_genus$total_ab_genus, 2))

total_ab_genus <- subset(total_ab_genus, total_ab_genus>0)

lbls <- sapply(1:nrow(total_ab_genus), function(i) {
  bquote(italic(.(total_ab_genus$names[i])) ~ "(" ~ .(total_ab_genus$formatted_percentage[i]) ~ "%)")
})

custom_palette <- c(
  "Citeromyces" = "darkgreen",
  "Hanseniaspora" = "#6ab100",
  "Hyphopichia" = "#01bc66",
  "Kluyveromyces" = "#00bdc2ff",  #to keep (1)
  "Kregervanrija" = "#1283a3",  
  "Lachancea" = "#00a8ff",  
  "Monosporozyma" = "dodgerblue4",  #to keep (2)
  "Ogataea" = "darkorchid4",  
  "Oligophagozyma" = "purple",
  "Pichia" = "#d570ff",  #to keep (1)
  "Saccharomyces" = "#ff4cd1", #to keep (1)
  "Saccharomycodes" = "#ff6b66", 
  "Torulaspora" = "#d59300"   #to keep (1)
)

ggplot(total_ab_genus, aes(x = "", y = total_ab_genus, fill = names)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(
    values = custom_palette,
    labels = lbls
  ) +
  theme_void() +
  theme(legend.text = element_text(size = 10)) +
  labs(fill = "Genus")
```


### Figure 3.B. Projection of samples using non-metric multidimensional scaling (NMDS).

Ordination is based on ASVs richness with colours indicating the dominant genus. Outlier samples LANY-1, TAN-1 and VASV-1 were removed for display purposes. 

```{r}
### Figure 2.B. Projection of samples using non-metric multidimensional scaling (NMDS).

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

filtered_data <- subset(metadata_clean, cum_val > 0)
filtered_data <- subset(filtered_data, cum_val!='NA')


custom_palette <- c("darkgreen","#6ab100","#01bc66","#00bdc2ff","#1283a3","dodgerblue4","darkorchid4",
                    "#d570ff","#ff4cd1", "#ff6b66", "#d59300")

filtered_data$Genus <- as.factor(filtered_data$Genus) # 10 levels

genus_colors <- setNames(custom_palette, levels(filtered_data$Genus))
filtered_data <- filtered_data %>%
  mutate(Genus_Color = genus_colors[Genus])

colnames(filtered_data)[1] <- "Sample"
filtered_data <- subset(filtered_data,!(Sample %in% c("LANY_1","TAN_1","VASV_1")))

nmds <- metaMDS(asv_matrix, distance = "bray", k = 2, trymax = 1000)
# Run 357 stress 0.1707845 
# ... Procrustes: rmse 0.001342651  max resid 0.006971928 

nmds_scores <- scores(nmds)
nmds_scores <- as.data.frame(nmds_scores$sites)
nmds_scores$Sample <- rownames(nmds_scores)

plot_data <- left_join(nmds_scores, filtered_data, by = "Sample")

p <- ggplot(plot_data, aes(x = NMDS1, y = NMDS2, color = Genus_Color, size = cum_val)) +
  geom_point(alpha = 0.7) +
  scale_color_identity() +
  theme_minimal() +
  labs(title = '', x = 'NMDS1', y = 'NMDS2', size = "Number of\ndifferent species")
p
```


### Figure 3. Chord diagram depicting yeast species co-occurrence, based on shared presence across sampled trees. 
The connections (arcs) between species indicate how frequently they co-occur, with the thickness of the arcs representing the strength of co-occurrence ( thicker lines corresponding to more frequent co-occurrence). 

```{r}
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
  "Citeromyces matritensis" = "darkgreen",
  "Hanseniaspora osmophila" = "#6ab100",
  "Hyphopichia burtonii" = "#01bc66", #ok
  "Monosporozyma servazzii" = "dodgerblue4", #ok
  "Kluyveromyces dobzhanskii" = "#00bdc2ff", #ok
  "Kluyveromyces lactis" = "#7FDEE0",
  "Kregervanrija fluxuum" = "#1283a3", #ok
  "Lachancea kluyveri" = "#00a8ff",
  "Ogataea dorogensis" = "darkorchid4", #ok
  "Ogataea nonfermentans" = "#6C7AD8",
  "Oligophagozyma castellii" = "purple",
  "Pichia mandshurica" = "#9545B7",
  "Pichia kudriavzevii" = "#ab4dd2",
  "Pichia membranifaciens" = "#d570ff", #ok
  "Pichia sp" = "#AF61D0",
  "Saccharomyces paradoxus" = "#ff4cd1",
  "Saccharomycodes ludwigii" = "#ff6b66",
  "Torulaspora delbrueckii" = "#d59300" #ok
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
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1] + 0.1,
      labels = parse(text = paste0("italic('", sector.name, "')")),  # Italicize name
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.5
    )
  }
)

# Run statistics on co-occurrence 
cooccurrence_results <- cooccur(mat = presence_absence_clean, type = "spp_site")
summary(cooccurrence_results)
results <- cooccurrence_results$results

avoiding_pairs <- subset(results, p_gt >= 0.05 & p_lt < 0.05) 
# 4/18: "Kluyveromyces_dobzhanskii" and "Torulaspora_delbrueckii"
# 15/16: "Pichia_sp" and "Saccharomyces_paradoxus"
cooccuring_pairs <- subset(results, p_lt >= 0.05 & p_gt < 0.05) 
# 16/18: "Saccharomyces_paradoxus" and "Torulaspora_delbrueckii"
```


### Figure 4.A. [Map of sampling locations with] pie charts showing the proportion of ASVs found per yeast genus. 

```{r}
custom_palette <- c(
  "Citeromyces" = "darkgreen",
  "Hanseniaspora" = "#6ab100",
  "Hyphopichia" = "#01bc66",
  "Kluyveromyces" = "#00bdc2ff",  #to keep (1)
  "Kregervanrija" = "#1283a3",  
  "Lachancea" = "#00a8ff",  
  "Monosporozyma" = "dodgerblue4",  #to keep (2)
  "Ogataea" = "darkorchid4",  
  "Oligophagozyma" = "purple",
  "Pichia" = "#d570ff",  #to keep (1)
  "Saccharomyces" = "#ff4cd1", #to keep (1)
  "Saccharomycodes" = "#ff6b66", 
  "Torulaspora" = "#d59300"   #to keep (1)
)
colnames(genus_tab_melted)[1:2] = c("Genus", "Sample")

# Remove pools
genus_tab_melted$group = metadata_clean$loc[match(genus_tab_melted$Sample,
                                                  metadata_clean$User.ID)]

# Save all the plots separately
samples <- unique(metadata_clean$Place_name)

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
  dplyr::group_by(group,Genus,Place_name) %>%
  dplyr::summarize(percent2 = mean(percent), val2 = sum(value))

# Loop through each sample, plot, and save as SVG
for (sample in samples) {

   sample_data <- subset(df2, Place_name == sample)
    
   p2 <- ggplot(sample_data, aes(x="", y=percent2, fill=Genus)) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0)+
      scale_fill_manual(values = custom_palette) + 
      theme_void()+
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),  
            axis.title.x = element_blank(),  axis.title.y = element_blank(),
            legend.position = "none") +
      ggtitle(unique(sample_data$Place_name))
    
    ggsave(filename = paste0("Pie_", sample, ".svg"), plot = p2, device = "svg", width = 2, height = 2)
}
```


### Figure 4.B. Longitude correlated with the abundance of the three most common genera isolated across all locations (Kluyveromyces, Saccharomyces, and Pichia).

```{r}
head(metadata_clean)
head(TAXtab)
head(ASVtab)

ASVtab_df <- as.data.frame(t(ASVtab))
ASVtab_df$Species <- TAXtab$Species[match(rownames(ASVtab_df), TAXtab$ASV_ID)]

species_tab <- aggregate(. ~ Species, data = ASVtab_df, FUN = sum)
species_tab <- species_tab[!is.na(species_tab$Species), ]

maj_species_tab <- subset(species_tab, Species %in% c("Kluyveromyces_dobzhanskii","Saccharomyces_paradoxus","Pichia_membranifaciens"))

library(reshape2)
library(dplyr)

maj_species_tab_long <- melt(maj_species_tab, id.vars = 'Species', 
                           variable.name = 'User.ID', value.name = 'value')
colnames(maj_species_tab_long) <- c('Species', 'User.ID', 'value')

merged_data <- merge(maj_species_tab_long, metadata_clean, by = 'User.ID')
merged_data <- subset(merged_data, value!=0)

# Full model
genus_sum_by_loc <- merged_data %>%
  dplyr::group_by(loc, Genus, Long, Lat) %>%
  dplyr::summarise(total_abundance = sum(value, na.rm = TRUE))

full_lm <- lm(total_abundance ~  Species * Long * Lat - Long:Lat - Species:Long:Lat, 
              data = species_sum_by_loc)
# SpeciesPichia_membranifaciens       -5127594    2428978  -2.111   0.0386 *
# Long                                  -62323      25740  -2.421   0.0182 *
# SpeciesPichia_membranifaciens:Long     72061      40077   1.798   0.0767 .
# SpeciesSaccharomyces_paradoxus:Long    93595      36775   2.545   0.0133 *
# Multiple R-squared:  0.2315,	Adjusted R-squared:  0.1383 
# F-statistic: 2.485 on 8 and 66 DF,  p-value: 0.0202

# With longitude only - do the same for latitude only
species_sum_by_loc <- merged_data %>%
  dplyr::group_by(loc, Species, Long) %>%
  dplyr::summarise(total_abundance = sum(value, na.rm = TRUE))

ggplot(species_sum_by_loc, aes(x = Long, y = total_abundance, 
                             color = Species, shape = Species)) +
  geom_point(alpha = 0.7, size=2) +
  geom_smooth(method = "lm", se = FALSE) +  
  scale_color_manual(values = c("#00bdc2ff", "#c77affff", "#ff61ccff")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_y_continuous(labels = scales::scientific) +  
  labs(x = "Longitude", 
       y = "Total number of ASVs per species") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"), 
    panel.grid = element_blank()               
  )

# Linear model
library(dplyr)
library(emmeans)

species_results <- data.frame(Species = character(), R_squared = numeric(), P_value = numeric(), stringsAsFactors = FALSE)

for (species in unique(species_sum_by_loc$Species)) {
  subset_data <- species_sum_by_loc %>% filter(Species == species)
  model <- lm(total_abundance ~ Long, data = subset_data)
  
  model_summary <- summary(model)
  r_squared <- model_summary$r.squared
  p_value <- model_summary$coefficients[2, 4]  
  
  species_results <- rbind(species_results, data.frame(Species = species, R_squared = r_squared, P_value = p_value))
}

print(species_results) # For longitude
# Species                   R_squared  P_value
# Kluyveromyces_dobzhanskii 0.20425939 0.02045415
# Pichia_membranifaciens    0.11157560 0.11929367
# Saccharomyces_paradoxus   0.03429596 0.36508193

print(species_results) # For latitude
#                  Species  R_squared  P_value
# Kluyveromyces_dobzhanskii 0.05529702 0.2475235
#    Pichia_membranifaciens 0.26543745 0.0118734
#   Saccharomyces_paradoxus 0.01754168 0.5189422
```


### Figure 5. Correlation patterns between environmental variables, tree metrics and yeast diversity indices.
Significance levels of Pearson correlations are indicated by asterisks (p < 0.05*, 0.01**, and 0.001***).

```{r}
# Create a binomial score for island and tree species
metadata_clean$island_num[metadata_clean$island=="yes"] <- 1
metadata_clean$island_num[metadata_clean$island=="no"] <- 0

metadata_clean$tree_sp_num[metadata_clean$tree_sp=="QuRo"] <- 0
metadata_clean$tree_sp_num[metadata_clean$tree_sp=="QuPe"] <- 1
metadata_clean$evenness <- as.numeric(as.character(metadata_clean$evenness))

cor_results <- rcorr(as.matrix(metadata_clean[,c("tree_sp_num",
                                                "Tree_diameter_cm",
                                                "Barkdepth_mm","Tree_height_m",
                                                "Tree_growth_mm_year","Tree_age", # tree
                                                "cum_val","evenness","shannon","richness", # div idx
                                                "max_precipitation","min_precipitation",
                                                "maximum_temp","minimum_temp","Lat","Long" # envir
                                                )]), type = "pearson")
cor_matrix <- cor_results$r   
p_matrix <- cor_results$P
names <- c("Tree species","Tree diameter","Bark depth","Tree height",
           "Tree growth","Tree age",
           "n species","Evenness","Shannon","ASV richness",
           "Max Precipitation","Min Precipitation",
           "Max Temperature","Min Temperature",
           "Latitude","Longitude")

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


### Figure 6. Annual A) maximum temperature, B) minimum temperature, C) maximum precipitation and D) minimum precipitation correlated with the ASV abundance (total number of ASV reads) of the three most common genera isolated across all locations (Kluyveromyces, Saccharomyces, and Pichia)
R-squared and p-values are indicated for each linear model per genus if  significant (or approaching significance).

```{r}
species_sum_by_loc <- merged_data %>%
  dplyr::group_by(loc, Species, maximum_temp) %>%
  dplyr::summarise(total_abundance = sum(value, na.rm = TRUE))
species_sum_by_loc <- subset(species_sum_by_loc, total_abundance>0)

plotA <- ggplot(species_sum_by_loc, aes(x = maximum_temp, y = total_abundance, 
                             color = Species, shape = Species)) +
  geom_point(alpha = 0.7, size=2) +
  geom_smooth(method = "lm", se = FALSE) +  
  scale_color_manual(values = c("#00bdc2ff", "#c77affff", "#ff61ccff")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_y_continuous(labels = scales::scientific) +  
  labs(x = "Annual Maximum Temperature", 
       y = "ASV Abundance per Species") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"), 
    panel.grid = element_blank()               
  )

species_sum_by_loc <- merged_data %>%
  dplyr::group_by(loc, Species, minimum_temp) %>%
  dplyr::summarise(total_abundance = sum(value, na.rm = TRUE))
species_sum_by_loc <- subset(species_sum_by_loc, total_abundance>0)

plotB <- ggplot(species_sum_by_loc, aes(x = minimum_temp, y = total_abundance, 
                             color = Species, shape = Species)) +
  geom_point(alpha = 0.7, size=2) +
  geom_smooth(method = "lm", se = FALSE) +  
  scale_color_manual(values = c("#00bdc2ff", "#c77affff", "#ff61ccff")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_y_continuous(labels = scales::scientific) +  
  labs(x = "Annual Minimum Temperature", 
       y = "ASV Abundance per Species") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"), 
    panel.grid = element_blank()               
  )

species_sum_by_loc <- merged_data %>%
  dplyr::group_by(loc, Species, max_precipitation) %>%
  dplyr::summarise(total_abundance = sum(value, na.rm = TRUE))
species_sum_by_loc <- subset(species_sum_by_loc, total_abundance>0)

plotC <- ggplot(species_sum_by_loc, aes(x = max_precipitation, y = total_abundance, 
                             color = Species, shape = Species)) +
  geom_point(alpha = 0.7, size=2) +
  geom_smooth(method = "lm", se = FALSE) +  
  scale_color_manual(values = c("#00bdc2ff", "#c77affff", "#ff61ccff")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_y_continuous(labels = scales::scientific) +  
  labs(x = "Annual Maximum Precipitation", 
       y = "ASV Abundance per Species") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"), 
    panel.grid = element_blank()               
  )

species_sum_by_loc <- merged_data %>%
  dplyr::group_by(loc, Species, min_precipitation) %>%
  dplyr::summarise(total_abundance = sum(value, na.rm = TRUE))
species_sum_by_loc <- subset(species_sum_by_loc, total_abundance>0)

plotD <- ggplot(species_sum_by_loc, aes(x = min_precipitation, y = total_abundance, 
                             color = Species, shape = Species)) +
  geom_point(alpha = 0.7, size=2) +
  geom_smooth(method = "lm", se = FALSE) +  
  scale_color_manual(values = c("#00bdc2ff", "#c77affff", "#ff61ccff")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_y_continuous(labels = scales::scientific) +  
  labs(x = "Annual Minimum Precipitation", 
       y = "ASV Abundance per Species") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"), 
    panel.grid = element_blank()               
  )

# plot all together: 
library(ggpubr)
ggarrange(plotA, plotB, plotC, plotD, ncol = 2, nrow = 2)

# Linear model
library(dplyr)
library(emmeans)

species_sum_by_loc <- merged_data %>%
  dplyr::group_by(loc, Species, max_precipitation) %>%
  dplyr::summarise(total_abundance = sum(value, na.rm = TRUE))
species_sum_by_loc <- subset(species_sum_by_loc, total_abundance>0)

species_results <- data.frame(Species = character(), R_squared = numeric(), P_value = numeric(), stringsAsFactors = FALSE)

for (species in unique(species_sum_by_loc$Species)) {
  subset_data <- species_sum_by_loc %>% filter(Species == species)
  model <- lm(total_abundance ~ max_precipitation, data = subset_data)
  
  model_summary <- summary(model)
  r_squared <- model_summary$r.squared
  p_value <- model_summary$coefficients[2, 4]  
  
  species_results <- rbind(species_results, data.frame(Species = species, R_squared = r_squared, P_value = p_value))
}

species_results
# Example for Maximum Temperature:
#      Species  R_squared    P_value
# Kluyveromyces 0.15678420 0.05006936 .
#        Pichia 0.05319734 0.25696539
# Saccharomyces 0.00417995 0.75881740
```


### Figure 7. Relative abundance of yeast genera depending on A) host tree species and B) insularity. 
Significant differences in community composition depending on x-axis are indicated by asterisks (p < 0.05*).

```{r}
# With tree species: 
data_tree <- subset(metadata_clean, metadata_clean$tree_sp!="NA")
species_tab_melted <- merge(species_tab_melted, data_tree[,c("User.ID","tree_sp")], by.x="Sample", by.y="User.ID")

# or for island, re-load 'species_tab_melted' and do:
species_tab_melted <- merge(species_tab_melted, metadata_clean[,c("User.ID","island")], 
                          by.x="Sample", by.y="User.ID")

species_summary <- species_tab_melted %>%
  filter(Abundance > 0) %>%
  group_by(Species, Species_pretty, tree_sp) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance))

plot7A <- ggplot(species_summary, aes(x = tree_sp, y = TotalAbundance, fill = Species_pretty)) +
  geom_bar(stat = "identity", position = "fill", width = 0.9) +
  scale_fill_manual(values = species_colors) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_x_discrete(labels = c(
    "QuRo" = expression(italic("Quercus robur")),
    "QuPe" = expression(italic("Quercus petraea"))
  )) +
  ylab("Relative abundance") +
  xlab("Tree species") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3),
    legend.text = element_text(face = "italic")
  )

# Prepare the genera abundance table (wide format) and metadata
species_abundance <- species_tab_melted %>%
  dcast(Sample ~ Genus, value.var = "value", fun.aggregate = sum)

# Set 'Sample' as rownames and remove it from the table
rownames(species_abundance) <- species_abundance$Sample
species_abundance <- species_abundance[, -1]

# Prepare the metadata (ensure unique rows for each Sample)
metadata <- species_tab_melted %>%
  select(Sample, island) %>%
  distinct()

# Calculate Bray-Curtis distance matrix
dist_matrix <- vegdist(species_abundance, method = "bray")

# Run adonis2
adonis_result <- adonis2(dist_matrix ~ island, data = metadata, permutations = 999)

# View the results
print(adonis_result)
```

## 6 - Plotting the data - supplementary

### Supplementary Figure 1. Percentage of ASVs mapping to each family across all samples.

```{r}
ASVtab_df <- as.data.frame(t(ASVtab))  # Transpose ASVtab so ASVs are columns
ASVtab_df$SampleID <- rownames(ASVtab_df)

# Transpose ASVtab back to long format with ASV IDs as a column
ASVtab_long <- as.data.frame(t(ASVtab))
ASVtab_long$ASV_ID <- rownames(ASVtab_long)

# Merge ASV counts with taxonomy (match by ASV_ID)
ASV_with_tax <- merge(ASVtab_long, TAXtab[, c("ASV_ID", "Family")], by = "ASV_ID")

# Clean data: convert counts to numeric
ASV_with_tax[, 2:(ncol(ASV_with_tax)-1)] <- lapply(ASV_with_tax[, 2:(ncol(ASV_with_tax)-1)], as.numeric)

# Aggregate ASV counts by Family
family_tab <- ASV_with_tax %>%
  group_by(Family) %>%
  summarise(across(where(is.numeric), sum))

# Remove rows with NA Family
family_tab <- filter(family_tab, !is.na(Family))

# Prepare data for plotting
family_matrix <- as.data.frame(family_tab)
rownames(family_matrix) <- family_matrix$Family
family_matrix <- family_matrix[, -1]  # Drop Family column

# Compute relative abundance
total_ab_family <- rowSums(family_matrix) / sum(rowSums(family_matrix)) * 100
total_ab_family <- as.data.frame(total_ab_family)
total_ab_family$Family <- rownames(total_ab_family)
total_ab_family <- subset(total_ab_family, total_ab_family>0)

# Create labels
lbls <- paste0(total_ab_family$Family, " (", round(total_ab_family$total_ab_family, 2), "%)")

# Plot
ggplot(total_ab_family, aes(x = "", y = total_ab_family, fill = Family)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = scales::hue_pal()(length(lbls)), labels = lbls) +
  theme_void() +
  theme(legend.text = element_text(size = 9)) +
  labs(fill = "Family")
```


### Supplementary Figure 2. Number of species in the three clusters that group the dominant genera.  
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

kruskal_test <- kruskal.test(sp_nb ~ Genus, data = maj_genus_metadata)
print(kruskal_test)
```


### Supplementary Figure 3. Co-occurrence heatmap between yeast species detected across all trees and sites.

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


### Supplementary Figure 4. ASV richness of species in the three dominant genera, A) correlated with ASV abundance and B) represented as boxplots.

```{r}
head(metadata_clean)
maj_genus_metadata <- subset(metadata_clean, Genus %in% c("Kluyveromyces","Saccharomyces","Pichia"))

ggplot(data = maj_genus_metadata, aes(x = Genus, y = richness, color = Genus)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  scale_color_manual(values = c("#00bdc2ff", "#c77affff","#ff61ccff"))+
  geom_jitter() +
  ylab("ASV richness") + 
  theme_bw()    

kruskal_test <- kruskal.test(richness ~ Genus, data = maj_genus_metadata)
print(kruskal_test) 

pairwise.wilcox.test(maj_genus_metadata$richness, maj_genus_metadata$Genus, 
                     p.adjust.method = "bonferroni")
```
