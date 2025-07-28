#load packages

library(tidyverse)
library(vegan)
library(reshape2)
library(iNEXT.beta3D)
library(iNEXT.3D)
library(ape)
library(ggpubr)
library(RColorBrewer)
library(wesanderson)
library(ggrepel)
library(codyn)
library(mFD)
library(patchwork)



###############################################################################
######################## Upload data ##########################################
###############################################################################

# Upload species abundance data for each depth, site, and period combination
abun_combined_sites <- read.csv("abundance_sites_alpha.csv")
str(abun_combined_sites)

# Remove the first column 
abun_combined_sites <- abun_combined_sites %>% select(-1)


###############################################################################
#################### Format data for use in iNEXT #############################
###############################################################################

# Get column names of abundance data
col_names <- names(abun_combined_sites)

# Function to extract site, depth, and year from column name
parse_name <- function(name) {
  parts <- unlist(strsplit(name, "_"))
  list(site = parts[1], depth = parts[2], year = paste(parts[3], parts[4], sep = "_"))
}

# make a list of unique site and year combinations

data_5_15m <- list()
data_25_40m <- list()

sites_years <- unique(sapply(col_names, function(name) {
  parsed <- parse_name(name)
  paste(parsed$site, parsed$year, sep = "_")
}))

# Loop through each site_year combo to  depths
for (sy in sites_years) {
  site <- strsplit(sy, "_")[[1]][1]
  year <- paste(strsplit(sy, "_")[[1]][2:3], collapse = "_")
  
  name_5  <- paste(site, "5", year, sep = "_")
  name_15 <- paste(site, "15", year, sep = "_")
  name_25 <- paste(site, "25", year, sep = "_")
  name_40 <- paste(site, "40", year, sep = "_")
  
  # Combine 5m and 15m depths
  if (name_5 %in% col_names && name_15 %in% col_names) {
    new_name <- paste(site, "5_15m", year, sep = "_")
    data_5_15m[[new_name]] <- abun_combined_sites[[name_5]] + abun_combined_sites[[name_15]]
  }
  
  # Combine 25m and 40m depths
  if (name_25 %in% col_names && name_40 %in% col_names) {
    new_name <- paste(site, "25_40m", year, sep = "_")
    data_25_40m[[new_name]] <- abun_combined_sites[[name_25]] + abun_combined_sites[[name_40]]
  }
}

# Convert lists to data frames
df_5_15m <- as.data.frame(data_5_15m)
df_25_40m <- as.data.frame(data_25_40m)

# Combine depth categories into a single data frame
combined_depths <- bind_cols(df_5_15m, df_25_40m)

# Check how it looks. looks in correct format.
str(combined_depths)







###############################################################################
############### Estimate Taxonomic Alpha Diversity (Hill q = 0 and 2) #############
###############################################################################

alpha_diversity <- estimate3D(
  combined_depths,
  diversity = 'TD',
  q = c(0, 2),
  datatype = "abundance",
  base = "coverage"
)

# NOTE: Output was manually cleaned (and site, depth, period separated) using Excel.
# no idea how to do this in r and spent too long trying to work it out.
# write.csv(alpha_diversity, "alpha_diversity_output_sites_depthscombined.csv")
# level/SC of 0.9715554 was taken from the outputs for q = 0 in the excel document


###############################################################################
########################### Load  output ######################################
###############################################################################

#read in data
alpha_data <- read.csv("alpha_diversity_output_sites_depthscombined.csv")

# Convert columns to factors for plotting
alpha_data$Period <- factor(alpha_data$Period)
alpha_data$Site <- factor(alpha_data$Site)
alpha_data$Depth <- factor(alpha_data$Depth, levels = c("5_15", "25_40"))
alpha_data$Order.q <- factor(alpha_data$Order.q)

# Rename periods for clarity and for plots
alpha_data$Period <- recode(alpha_data$Period,
                            "2014_2015" = "2014/15",
                            "2022_2023" = "2022/23")

# Split data by Hill number (q = 0 and q = 2)
data_q0 <- alpha_data %>% filter(Order.q == "0") %>% as_tibble()
data_q2 <- alpha_data %>% filter(Order.q == "2") %>% as_tibble()

# Calculate confidence interval widths
data_q0 <- data_q0 %>%
  mutate(qTD.LCL.adj = qTD - qTD.LCL,
         qTD.UCL.adj = qTD.UCL - qTD)

data_q2 <- data_q2 %>%
  mutate(qTD.LCL.adj = qTD - qTD.LCL,
         qTD.UCL.adj = qTD.UCL - qTD)




###############################################################################
############################## Plot tax alpha diversity #######################
###############################################################################

# Define color palette
colors <- wes_palette("Cavalcanti1", n = 5, type = "discrete")


### Alpha diversity plots (q = 0) 

alpha_5_15_0 <- data_q0 %>%
  filter(Depth == "5_15") %>%
  ggplot(aes(x = Period, y = qTD, color = Site, group = Site)) +
  annotate("rect", xmin = 1.5, xmax = Inf, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "grey50") +
  geom_line(position = position_dodge(0.45), linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 2, position = position_dodge(0.45)) +
  geom_errorbar(aes(ymin = qTD - qTD.LCL.adj, ymax = qTD + qTD.UCL.adj),
                width = 0.4, position = position_dodge(0.45), lwd = 0.4) +
  scale_color_manual(values = colors) +
  labs(x = "", y = "Taxonomic diversity") +
  ylim(0, 64) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "right",
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

alpha_25_40_0 <- data_q0 %>%
  filter(Depth == "25_40") %>%
  ggplot(aes(x = Period, y = qTD, color = Site, group = Site)) +
  annotate("rect", xmin = 1.5, xmax = Inf, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "grey50") +
  geom_line(position = position_dodge(0.45), linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 2, position = position_dodge(0.45)) +
  geom_errorbar(aes(ymin = qTD - qTD.LCL.adj, ymax = qTD + qTD.UCL.adj),
                width = 0.4, position = position_dodge(0.45), lwd = 0.4) +
  scale_color_manual(values = colors) +
  labs(x = "", y = "Taxonomic diversity") +
  ylim(0, 64) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "right",
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


#### Alpha diversity plots (q = 2)

alpha_5_15_2 <- data_q2 %>%
  filter(Depth == "5_15") %>%
  ggplot(aes(x = Period, y = qTD, color = Site, group = Site)) +
  annotate("rect", xmin = 1.5, xmax = Inf, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "grey50") +
  geom_line(position = position_dodge(0.45), linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 2, position = position_dodge(0.45)) +
  geom_errorbar(aes(ymin = qTD - qTD.LCL.adj, ymax = qTD + qTD.UCL.adj),
                width = 0.4, position = position_dodge(0.45), lwd = 0.4) +
  scale_color_manual(values = colors) +
  labs(x = "", y = "Taxonomic diversity") +
  ylim(0, 64) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "right",
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

alpha_25_40_2 <- data_q2 %>%
  filter(Depth == "25_40") %>%
  ggplot(aes(x = Period, y = qTD, color = Site, group = Site)) +
  annotate("rect", xmin = 1.5, xmax = Inf, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "grey50") +
  geom_line(position = position_dodge(0.45), linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 2, position = position_dodge(0.45)) +
  geom_errorbar(aes(ymin = qTD - qTD.LCL.adj, ymax = qTD + qTD.UCL.adj),
                width = 0.4, position = position_dodge(0.45), lwd = 0.4) +
  scale_color_manual(values = colors) +
  labs(x = "", y = "Taxonomic diversity") +
  ylim(0, 64) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "right",
    text = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


#### combine the above plots (half of first results figure)

# Remove  y-axis labels and adjust legends
alpha_5_15_0  <- alpha_5_15_0  + theme(axis.title.y = element_blank()) +
  guides(color = guide_legend(title = NULL, nrow = 2))
alpha_5_15_2  <- alpha_5_15_2  + theme(axis.title.y = element_blank()) +
  guides(color = guide_legend(title = NULL, nrow = 2))
alpha_25_40_0 <- alpha_25_40_0 + theme(axis.title.y = element_blank()) +
  guides(color = guide_legend(title = NULL, nrow = 2))
alpha_25_40_2 <- alpha_25_40_2 + theme(axis.title.y = element_blank()) +
  guides(color = guide_legend(title = NULL, nrow = 2))

# Arrange four plots into a single panel
arranged_plot <- ggarrange(
  alpha_5_15_0, alpha_5_15_2,
  alpha_25_40_0, alpha_25_40_2,
  ncol = 4, nrow = 1,
  common.legend = TRUE,
  legend = "none",
  labels = c("5–15 m (q = 0)", "5–15 m (q = 2)", "25–40 m (q = 0)", "25–40 m (q = 2)"),
  vjust = 2.5,
  hjust = -0.4,
  font.label = list(size = 11),
  align = "v"
)

# Add overall y-axis label
final_plot <- annotate_figure(
  arranged_plot,
  left = text_grob(expression(Taxonomic~alpha-diversity), rot = 90, size = 14, hjust = 0.3)
)

# Half of first results figure done
print(final_plot)










###############################################################################
################### Functional alpha  diversity ########################
###############################################################################

# Load Fish Trait Data
fish_traits <- read.csv("species_traits.csv", header = TRUE, row.names = 1)

# Correct column types
fish_traits$diet <- factor(fish_traits$diet)
fish_traits$depth.range <- as.numeric(fish_traits$depth.range)
fish_traits$body.size <- as.numeric(fish_traits$body.size)
fish_traits$water.column.position <- factor(fish_traits$water.column.position)
fish_traits$Gregariousness <- factor(fish_traits$Gregariousness)

glimpse(fish_traits)

# Explore Trait Correlations
cor.test(fish_traits$body.size, fish_traits$depth.range, method = "pearson")
cor.test(fish_traits$body.size, fish_traits$trophic.level, method = "pearson")
cor.test(fish_traits$depth.range, fish_traits$trophic.level, method = "pearson")

# Trophic level correlates with others and is probably redundant so i'll remove it
fish_traits <- fish_traits %>% dplyr::select(-trophic.level)
glimpse(fish_traits)

#### Create Trait Category Data Frame
# Define the type of each trait: 
# "N" = nominal (categorical), "Q" = quantitative (numerical)
fish_traits_cat <- data.frame(names(fish_traits), c("N", "Q", "N", "Q", "N"))
colnames(fish_traits_cat) <- c("trait_name", "trait_type")
fish_traits_cat

#  Calculate Functional Distance Matrix
dist_fISH <- mFD::funct.dist(
  sp_tr        = fish_traits,
  tr_cat       = fish_traits_cat,
  metric       = "gower",          # Gower distance for mixed trait types
  scale_euclid = "scale_center",   # Scale quantitative traits
  ordinal_var  = "classic",        # Treat ordinal variables with classic method
  weight_type  = "equal",          # Equal weighting of traits
  stop_if_NA   = TRUE              # Stop if NA values present
)

# Summary of distance matrix
summary(as.matrix(dist_fISH))

# Define  Species List
# These are the species to extract for a sub-matrix
ex_species <- c(
  "Abudefduf saxatilis", "Acanthostracion polygonius", "Acanthostracion quadricornis", 
  "Acanthurus bahianus", "Acanthurus chirurgus", "Acanthurus coeruleus",
  "Aluterus scriptus", "Anisotremus surinamensis", "Anisotremus virginicus",
  "Apsilus dentatus", "Aulostomus maculatus", "Balistes vetula", "Bodianus rufus",
  "Calamus bajonado", "Calamus calamus", "Cantherhines pullus", "Canthidermis sufflamen",
  "Canthigaster rostrata", "Caranx bartholomaei", "Caranx crysos", "Caranx ruber",
  "Cephalopholis cruentata", "Cephalopholis fulva", "Chaetodon capistratus",
  "Chaetodon ocellatus", "Chaetodon sedentarius", "Chaetodon striatus",
  "Chromis cyanea", "Chromis insolata", "Chromis multilineata", "Chromis scotti",
  "Clepticus parrae", "Diodon hystrix", "Echeneis neucratoides",
  "Epinephelus adscensionis", "Epinephelus guttatus", "Epinephelus striatus",
  "Equetus lanceolatus", "Equetus punctatus", "Gramma loreto", "Gramma melacara",
  "Gymnothorax funebris", "Haemulon album", "Haemulon aurolineatum",
  "Haemulon bonariense", "Haemulon carbonarium", "Haemulon flavolineatum",
  "Haemulon macrostomum", "Haemulon parra", "Haemulon plumierii",
  "Haemulon sciurus", "Haemulon vittatum", "Halichoeres bivittatus",
  "Halichoeres cyanocephalus", "Halichoeres garnoti", "Halichoeres maculipinna",
  "Halichoeres pictus", "Halichoeres radiatus", "Holacanthus ciliaris",
  "Holacanthus tricolor", "Holocentrus adscensionis", "Holocentrus rufus",
  "Hypoplectrus aberrans", "Hypoplectrus gummigutta", "Hypoplectrus indigo",
  "Hypoplectrus nigricans", "Hypoplectrus puella", "Kyphosus sectatrix",
  "Kyphosus vaigiensis", "Lachnolaimus maximus", "Lactophrys bicaudalis",
  "Lactophrys trigonus", "Lactophrys triqueter", "Lutjanus analis",
  "Lutjanus apodus", "Lutjanus griseus", "Lutjanus jocu", "Lutjanus mahogoni",
  "Lutjanus synagris", "Malacanthus plumieri", "Melichthys niger",
  "Microspathodon chrysurus", "Mulloidichthys martinicus", "Mycteroperca bonaci",
  "Mycteroperca interstitialis", "Mycteroperca tigris", "Mycteroperca venenosa",
  "Neoniphon marianus", "Ocyurus chrysurus", "Pomacanthus arcuatus",
  "Pomacanthus paru", "Prognathodes aculeatus", "Pseudupeneus maculatus",
  "Pterois volitans", "Sargocentron bullisi", "Scarus coeruleus", "Scarus iseri",
  "Scarus taeniopterus", "Scarus vetula", "Serranus tabacarius", "Serranus tigrinus",
  "Serranus tortugarum", "Sparisoma atomarium", "Sparisoma aurofrenatum",
  "Sparisoma chrysopterum", "Sparisoma rubripinne", "Sparisoma viride",
  "Sphyraena barracuda", "Stegastes adustus", "Stegastes diencaeus",
  "Stegastes leucostictus", "Stegastes partitus", "Stegastes planifrons",
  "Stegastes variabilis", "Synodus intermedius", "Synodus saurus",
  "Synodus synodus", "Thalassoma bifasciatum", "Trachinotus falcatus"
)

#  Extract Subset of Distance Matrix for Selected Species
distance_matrix <- round(as.matrix(dist_fISH)[ex_species, ex_species], 2)

########### format abundance data to use in inext. functional slightly different to tax #################

# Load site abundance data (rows = species, columns = site-depth-period)
abun_combined_sites <- read.csv("abundance_sites_alpha.csv", row.names = 1)
str(abun_combined_sites)

# take all column names from abundance data
col_names <- names(abun_combined_sites)

# split column name into site, depth, and year components
parse_name <- function(name) {
  parts <- unlist(strsplit(name, "_"))
  list(site = parts[1], depth = parts[2], year = paste(parts[3], parts[4], sep = "_"))
}

# make empty list to store new depth-combined abundance columns
data_5_15m <- list()
data_25_40m <- list()

# take all unique site-year combinations
sites_years <- unique(sapply(col_names, function(name) {
  parsed <- parse_name(name)
  paste(parsed$site, parsed$year, sep = "_")
}))

# For each site-year combo, sum depths 5+15m and 25+40m
for (sy in sites_years) {
  site <- strsplit(sy, "_")[[1]][1]
  year <- paste(strsplit(sy, "_")[[1]][2:3], collapse = "_")
  
  name_5 <- paste(site, "5", year, sep = "_")
  name_15 <- paste(site, "15", year, sep = "_")
  name_25 <- paste(site, "25", year, sep = "_")
  name_40 <- paste(site, "40", year, sep = "_")
  
  # Combine 5m + 15m
  if (name_5 %in% col_names && name_15 %in% col_names) {
    new_name <- paste(site, "5_15m", year, sep = "_")
    data_5_15m[[new_name]] <- abun_combined_sites[[name_5]] + abun_combined_sites[[name_15]]
  }
  
  # Combine 25m + 40m
  if (name_25 %in% col_names && name_40 %in% col_names) {
    new_name <- paste(site, "25_40m", year, sep = "_")
    data_25_40m[[new_name]] <- abun_combined_sites[[name_25]] + abun_combined_sites[[name_40]]
  }
}

# change combined depth lists to data frames, keeping species row names
df_5_15m  <- as.data.frame(data_5_15m,  row.names = rownames(abun_combined_sites))
df_25_40m <- as.data.frame(data_25_40m, row.names = rownames(abun_combined_sites))

# Combine both depth layers into a single abundance matrix
abun_combined_sites_functional <- bind_cols(df_5_15m, df_25_40m)
rownames(abun_combined_sites_functional) <- rownames(abun_combined_sites)

# Check structure of final functional abundance matrix. all looks good!
str(abun_combined_sites_functional)


############### Estimate functional diversity (FD) at q = 0 and q = 2 ##########

data <- abun_combined_sites_functional
distM <- distance_matrix

functional_diversity <- estimate3D(data, diversity = 'FD', datatype = "abundance", 
                                   base = "coverage", q = 0, level = 0.9715554,
                                   FDdistM = distM, FDtype = 'AUC')

functional_diversity_1 <- estimate3D(data, diversity = 'FD', datatype = "abundance", 
                                     base = "coverage", q = 2, level = 1,
                                     FDdistM = distM, FDtype = 'AUC')

write.csv(functional_diversity, "functional_alpha_data_q0_revised.csv")
write.csv(functional_diversity_1, "functional_alpha_data_q2_revised.csv")

# datasets were combined in excel and columns added for sites, depths and periods
# still can't workout how best to do it in r. output from inext hard to subset

#### Load  combine both datasets and format for plots 

functional_alpha_data <- read.csv("functional_alpha_data_q0_revised.csv")
functional_alpha_data$Period <- factor(functional_alpha_data$Period)
functional_alpha_data$Site <- factor(functional_alpha_data$Site)
functional_alpha_data$Depth <- factor(functional_alpha_data$Depth, levels = c("5_15", "25_40"))
functional_alpha_data$Order.q <- factor(functional_alpha_data$Order.q)

# Rename periods for plots
functional_alpha_data$Period <- recode(functional_alpha_data$Period,
                                       "2014_2015" = "2014/15",
                                       "2022_2023" = "2022/23")

# Split by q order
functional_data_q0 <- functional_alpha_data %>% filter(Order.q == "0") %>% as_tibble()
functional_data_q2 <- functional_alpha_data %>% filter(Order.q == "2") %>% as_tibble()

# Calculate error bars
functional_data_q0 <- functional_data_q0 %>%
  mutate(qFD.LCL.adj = qFD - qFD.LCL,
         qFD.UCL.adj = qFD.UCL - qFD)

functional_data_q2 <- functional_data_q2 %>%
  mutate(qFD.LCL.adj = qFD - qFD.LCL,
         qFD.UCL.adj = qFD.UCL - qFD)



# Define the color palette
colors <- wes_palette("Cavalcanti1", n = 5, type = "discrete")

##### functional alpha diversity plots ######################

functional_alpha_5_15_0 <-functional_data_q0 %>% filter(Depth == "5_15")%>%
  ggplot(aes(x = as.factor(Period), y = qFD, color = Site, group = Site)) +
  annotate("rect", xmin=1.5, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.15, fill="grey50")+
  # Add lines connecting the same sites
  geom_line(position = position_dodge(width = 0.45), linewidth = 0.8 , alpha = 0.7) +
  # Plot means with standard error bars
  geom_point(aes(color = Site), size = 2, position = position_dodge(width = 0.45)) +
  geom_errorbar(aes(ymin = qFD - qFD.LCL.adj, ymax = qFD + qFD.UCL.adj, color = Site),
                position = position_dodge(width = 0.45), width = 0.4, lwd = 0.4) +
  scale_color_manual(values = colors) +
  labs(x = "", y = "Functional diversity") +
  ylim(0, 18) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right",
        legend.text=element_text(size=14),
        text = element_text(size = 12, family = "sans", color = "black"),
        axis.title.y = element_text(size = 12, family = "sans", color = "black"),
        axis.text = element_text(size = 12, family = "sans", color = "black"),
        strip.text = element_text(size = 12, family = "sans", face = "bold", color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

functional_alpha_25_40_0 <-functional_data_q0 %>% filter(Depth == "25_40")%>%
  ggplot(aes(x = as.factor(Period), y = qFD, color = Site, group = Site)) +
  annotate("rect", xmin=1.5, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.15, fill="grey50")+
  # Add lines connecting the same sites
  geom_line(position = position_dodge(width = 0.45), linewidth = 0.8 , alpha = 0.7) +
  # Plot means with standard error bars
  geom_point(aes(color = Site), size = 2, position = position_dodge(width = 0.45)) +
  geom_errorbar(aes(ymin = qFD - qFD.LCL.adj, ymax = qFD + qFD.UCL.adj, color = Site),
                position = position_dodge(width = 0.45), width = 0.4, lwd = 0.4) +
  scale_color_manual(values = colors) +
  labs(x = "", y = "Functional diversity") +
  ylim(0, 18) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right",
        legend.text=element_text(size=14),
        text = element_text(size = 12, family = "sans", color = "black"),
        axis.title.y = element_text(size = 12, family = "sans", color = "black"),
        axis.text = element_text(size = 12, family = "sans", color = "black"),
        strip.text = element_text(size = 12, family = "sans", face = "bold", color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

functional_alpha_5_15_2 <-functional_data_q2 %>% filter(Depth == "5_15")%>%
  ggplot(aes(x = as.factor(Period), y = qFD, color = Site, group = Site)) +
  annotate("rect", xmin=1.5, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.15, fill="grey50")+
  # Add lines connecting the same sites
  geom_line(position = position_dodge(width = 0.45), linewidth = 0.8 , alpha = 0.7) +
  # Plot means with standard error bars
  geom_point(aes(color = Site), size = 2, position = position_dodge(width = 0.45)) +
  geom_errorbar(aes(ymin = qFD - qFD.LCL.adj, ymax = qFD + qFD.UCL.adj, color = Site),
                position = position_dodge(width = 0.45), width = 0.4, lwd = 0.4) +
  scale_color_manual(values = colors) +
  labs(x = "", y = "Functional diversity") +
  ylim(0, 18) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right",
        legend.text=element_text(size=14),
        text = element_text(size = 12, family = "sans", color = "black"),
        axis.title.y = element_text(size = 12, family = "sans", color = "black"),
        axis.text = element_text(size = 12, family = "sans", color = "black"),
        strip.text = element_text(size = 12, family = "sans", face = "bold", color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

functional_alpha_25_40_2 <-functional_data_q2 %>% filter(Depth == "25_40")%>%
  ggplot(aes(x = as.factor(Period), y = qFD, color = Site, group = Site)) +
  annotate("rect", xmin=1.5, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.15, fill="grey50")+
  # Add lines connecting the same sites
  geom_line(position = position_dodge(width = 0.45), linewidth = 0.8 , alpha = 0.7) +
  # Plot means with standard error bars
  geom_point(aes(color = Site), size = 2, position = position_dodge(width = 0.45)) +
  geom_errorbar(aes(ymin = qFD - qFD.LCL.adj, ymax = qFD + qFD.UCL.adj, color = Site),
                position = position_dodge(width = 0.45), width = 0.4, lwd = 0.4) +
  scale_color_manual(values = colors) +
  labs(x = "", y = "Functional diversity") +
  ylim(0, 18) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right",
        legend.text=element_text(size=14),
        text = element_text(size = 12, family = "sans", color = "black"),
        axis.title.y = element_text(size = 12, family = "sans", color = "black"),
        axis.text = element_text(size = 12, family = "sans", color = "black"),
        strip.text = element_text(size = 12, family = "sans", face = "bold", color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# combine these functional plots into one panel 

functional_alpha_5_15_0 <- functional_alpha_5_15_0 +
  theme(axis.title.y = element_blank()) +
  guides(color = guide_legend(title = NULL, nrow = 2))
functional_alpha_5_15_2 <- functional_alpha_5_15_2 +
  theme(axis.title.y = element_blank()) +
  guides(color = guide_legend(title = NULL,nrow = 2))

functional_alpha_25_40_0 <- functional_alpha_25_40_0 +
  theme(axis.title.y = element_blank()) +
  guides(color = guide_legend(title = NULL,nrow = 2))

functional_alpha_25_40_2 <- functional_alpha_25_40_2 +
  theme(axis.title.y = element_blank()) +
  guides(color = guide_legend(title = NULL,nrow = 2))

# Remove y-axis title from each plot
functional_alpha_5_15_0 <- functional_alpha_5_15_0 + theme(axis.title.y = element_blank())
functional_alpha_5_15_2 <- functional_alpha_5_15_2 + theme(axis.title.y = element_blank())
functional_alpha_25_40_0 <- functional_alpha_25_40_0 + theme(axis.title.y = element_blank())
functional_alpha_25_40_2 <- functional_alpha_25_40_2 + theme(axis.title.y = element_blank())

# Arrange the plots
functional_arranged_plot <- ggarrange(functional_alpha_5_15_0, functional_alpha_5_15_2, functional_alpha_25_40_0, functional_alpha_25_40_2,
                                      ncol = 4, nrow = 1,
                                      common.legend = TRUE,
                                      legend = "bottom",
                                      labels = c("5-15 m (q = 0)", "5-15 m (q = 2)", "25-40 m (q = 0)", "25-40 m (q = 2)"),
                                      vjust = 2.5,
                                      hjust = -0.4,
                                      font.label = list(size = 11),
                                      align = c("v"))

# Annotate the arranged plot with the new titles
functional_final_plot <- annotate_figure(functional_arranged_plot,
                                         left = text_grob(expression(Functional~alpha-diversity), rot = 90, size = 14, hjust = 0.3),
                                         bottom = text_grob("Time period", size = 14, vjust = -6, hjust = 0.3,))


# Print the final plot
print(functional_final_plot)





################################################################################
########### combine tax and functional plots to make alpha diversity figure ####
################################################################################

# Combine the two plots with a common legend below
combined_plot <- ggarrange(final_plot, functional_final_plot, 
                           ncol = 1, nrow = 2,
                           heights = c(0.9, 1.15),
                           labels = c("(a)", "(b)"),
                           common.legend = TRUE, 
                           legend = "bottom",
                           align = c("hv"))

# need to add top margin to final_plot
final_plot <- final_plot + 
  theme(plot.margin = margin(t = 30, r = 5, b = 5, l = 5))  # t = top margin in points

# Combine the two plots
combined_plot <- ggarrange(final_plot, functional_final_plot, 
                           ncol = 1, nrow = 2,
                           heights = c(0.9, 1.15),
                           common.legend = TRUE, 
                           legend = "bottom",
                           align = "hv")

print(combined_plot) # this is the alpha diversity plot in results












################################################################################
###########  Beta Diversity Analysis - Same Site & Depth, Different Years #########
################################################################################


##### Taxonomic Beta Diversity

#  compare beta diversity across time periods for the same site and depth combination.
# First, extract columns corresponding to each site-depth pair across both time periods.

# 5–15 m depth
CV_5_15m <- combined_depths[, c("CV_5_15m_2014_2015", "CV_5_15m_2022_2023")]
LB_5_15m <- combined_depths[, c("LB_5_15m_2014_2015", "LB_5_15m_2022_2023")]
RC_5_15m <- combined_depths[, c("RC_5_15m_2014_2015", "RC_5_15m_2022_2023")]
RP_5_15m <- combined_depths[, c("RP_5_15m_2014_2015", "RP_5_15m_2022_2023")]
TM_5_15m <- combined_depths[, c("TM_5_15m_2014_2015", "TM_5_15m_2022_2023")]

# 25–40 m depth
CV_25_40m <- combined_depths[, c("CV_25_40m_2014_2015", "CV_25_40m_2022_2023")]
LB_25_40m <- combined_depths[, c("LB_25_40m_2014_2015", "LB_25_40m_2022_2023")]
RC_25_40m <- combined_depths[, c("RC_25_40m_2014_2015", "RC_25_40m_2022_2023")]
RP_25_40m <- combined_depths[, c("RP_25_40m_2014_2015", "RP_25_40m_2022_2023")]
TM_25_40m <- combined_depths[, c("TM_25_40m_2014_2015", "TM_25_40m_2022_2023")]


# Create Named Lists for iNEXTbeta3D Input

CV_5_15m_list <- list(CV_5_15m = CV_5_15m)
LB_5_15m_list <- list(LB_5_15m = LB_5_15m)
RC_5_15m_list <- list(RC_5_15m = RC_5_15m)
RP_5_15m_list <- list(RP_5_15m = RP_5_15m)
TM_5_15m_list <- list(TM_5_15m = TM_5_15m)

CV_25_40m_list <- list(CV_25_40m = CV_25_40m)
LB_25_40m_list <- list(LB_25_40m = LB_25_40m)
RC_25_40m_list <- list(RC_25_40m = RC_25_40m)
RP_25_40m_list <- list(RP_25_40m = RP_25_40m)
TM_25_40m_list <- list(TM_25_40m = TM_25_40m)


# calculate beta diversity (TD, q = 0 and 2) between time periods for each site-depth combo

# This must be done one by one. cannot be automated in this case.

beta_CV_5_15m <- iNEXTbeta3D(data = CV_5_15m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95)
beta_LB_5_15m <- iNEXTbeta3D(data = LB_5_15m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95)
beta_RC_5_15m <- iNEXTbeta3D(data = RC_5_15m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95)
beta_RP_5_15m <- iNEXTbeta3D(data = RP_5_15m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95)
beta_TM_5_15m <- iNEXTbeta3D(data = TM_5_15m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95)

beta_CV_25_40m <- iNEXTbeta3D(data = CV_25_40m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95)
beta_LB_25_40m <- iNEXTbeta3D(data = LB_25_40m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95)
beta_RC_25_40m <- iNEXTbeta3D(data = RC_25_40m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95)
beta_RP_25_40m <- iNEXTbeta3D(data = RP_25_40m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95)
beta_TM_25_40m <- iNEXTbeta3D(data = TM_25_40m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95)

beta_CV_5_15m0 <- iNEXTbeta3D(data = CV_5_15m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95)
beta_LB_5_15m0 <- iNEXTbeta3D(data = LB_5_15m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95)
beta_RC_5_15m0 <- iNEXTbeta3D(data = RC_5_15m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95)
beta_RP_5_15m0 <- iNEXTbeta3D(data = RP_5_15m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95)
beta_TM_5_15m0 <- iNEXTbeta3D(data = TM_5_15m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95)

beta_CV_25_40m0 <- iNEXTbeta3D(data = CV_25_40m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95)
beta_LB_25_40m0 <- iNEXTbeta3D(data = LB_25_40m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95)
beta_RC_25_40m0 <- iNEXTbeta3D(data = RC_25_40m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95)
beta_RP_25_40m0 <- iNEXTbeta3D(data = RP_25_40m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95)
beta_TM_25_40m0 <- iNEXTbeta3D(data = TM_25_40m_list, diversity = "TD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95)


# again output put into excel sheet because i can't workout how to do it automatically. Ask r group for help with this maybe
# Create Plots of Taxonomic Beta Diversity

tax_beta <- read.csv("taxonomic_beta_groupeddepth_byyear.csv")
tax_beta$Site <- factor(tax_beta$Site)

# reanme depths and order depth
tax_beta$Depth <- recode(tax_beta$Depth,
                         "5_15" = "5-15",
                         "25_40" = "25-40")
tax_beta$Depth <- factor(tax_beta$Depth, levels = c("5-15", "25-40"))

# Plot for q = 0
tax_beta_0_plot <- ggplot(tax_beta, aes(x = as.factor(Site), y = Beta_q0, color = Depth, group = interaction(Site, Depth))) +
  geom_point(size = 2, position = position_dodge(width = -0.45)) +
  geom_errorbar(aes(ymin = Beta_LCL_q0, ymax = Beta_UCL_q0),
                position = position_dodge(width = -0.45),
                width = 0.3, lwd = 0.5) +
  labs(x = NULL, y = expression(Taxonomic~beta-diversity)) +
  scale_color_manual(values = c("5-15" = "red", "25-40" = "black")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    text = element_text(size = 12, family = "sans", color = "black"),
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_flip() +
  scale_y_continuous(limits = c(0.8, 1.9), breaks = c(1.0, 1.4, 1.8))

# Plot for q = 2
tax_beta_2_plot <- ggplot(tax_beta, aes(x = as.factor(Site), y = Beta_q2, color = Depth, group = interaction(Site, Depth))) +
  geom_point(size = 2, position = position_dodge(width = -0.45)) +
  geom_errorbar(aes(ymin = Beta_LCL_q2, ymax = Beta_UCL_q2),
                position = position_dodge(width = -0.45),
                width = 0.3, lwd = 0.5) +
  labs(x = NULL, y = expression(Taxonomic~beta-diversity)) +
  scale_color_manual(name = "Depth (m)", values = c("5-15" = "red", "25-40" = "black")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    text = element_text(size = 12, family = "sans", color = "black"),
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_flip() +
  scale_y_continuous(limits = c(0.8, 1.9), breaks = c(1.0, 1.4, 1.8, 1.8))


# Combine Plots

temporal_site_tax_plot <- tax_beta_0_plot + tax_beta_2_plot +
  plot_layout(ncol = 2, nrow = 1, axis_titles = "collect")










################################################################################################
############################### FUNCTIONAL DISSIMILARITY ANALYSIS ###############################
################################################################################################

# Functional beta-diversity comparison: same site-depths, different time periods
# Subset abundance data for each site-depth combination (2014–2015 vs 2022–2023)

CV_5_15m_functional <- abun_combined_sites_functional %>% select(CV_5_15m_2014_2015, CV_5_15m_2022_2023)
LB_5_15m_functional <- abun_combined_sites_functional %>% select(LB_5_15m_2014_2015, LB_5_15m_2022_2023)
RC_5_15m_functional <- abun_combined_sites_functional %>% select(RC_5_15m_2014_2015, RC_5_15m_2022_2023)
RP_5_15m_functional <- abun_combined_sites_functional %>% select(RP_5_15m_2014_2015, RP_5_15m_2022_2023)
TM_5_15m_functional <- abun_combined_sites_functional %>% select(TM_5_15m_2014_2015, TM_5_15m_2022_2023)

CV_25_40m_functional <- abun_combined_sites_functional %>% select(CV_25_40m_2014_2015, CV_25_40m_2022_2023)
LB_25_40m_functional <- abun_combined_sites_functional %>% select(LB_25_40m_2014_2015, LB_25_40m_2022_2023)
RC_25_40m_functional <- abun_combined_sites_functional %>% select(RC_25_40m_2014_2015, RC_25_40m_2022_2023)
RP_25_40m_functional <- abun_combined_sites_functional %>% select(RP_25_40m_2014_2015, RP_25_40m_2022_2023)
TM_25_40m_functional <- abun_combined_sites_functional %>% select(TM_25_40m_2014_2015, TM_25_40m_2022_2023)

# Calculate functional beta diversity for each site-depth using iNEXTbeta3D
# First with q = 2 and level = 1

functional_beta_CV_5_15m_2 <- iNEXTbeta3D(data = CV_5_15m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_LB_5_15m_2 <- iNEXTbeta3D(data = LB_5_15m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_RC_5_15m_2 <- iNEXTbeta3D(data = RC_5_15m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_RP_5_15m_2 <- iNEXTbeta3D(data = RP_5_15m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_TM_5_15m_2 <- iNEXTbeta3D(data = TM_5_15m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)

functional_beta_CV_25_40m_2 <- iNEXTbeta3D(data = CV_25_40m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_LB_25_40m_2 <- iNEXTbeta3D(data = LB_25_40m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_RC_25_40m_2 <- iNEXTbeta3D(data = RC_25_40m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_RP_25_40m_2 <- iNEXTbeta3D(data = RP_25_40m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_TM_25_40m_2 <- iNEXTbeta3D(data = TM_25_40m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)

# Do it again for q = 0 and level = 0.9715554

functional_beta_CV_5_15m_0 <- iNEXTbeta3D(data = CV_5_15m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_LB_5_15m_0 <- iNEXTbeta3D(data = LB_5_15m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_RC_5_15m_0 <- iNEXTbeta3D(data = RC_5_15m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_RP_5_15m_0 <- iNEXTbeta3D(data = RP_5_15m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_TM_5_15m_0 <- iNEXTbeta3D(data = TM_5_15m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)

functional_beta_CV_25_40m_0 <- iNEXTbeta3D(data = CV_25_40m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_LB_25_40m_0 <- iNEXTbeta3D(data = LB_25_40m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_RC_25_40m_0 <- iNEXTbeta3D(data = RC_25_40m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_RP_25_40m_0 <- iNEXTbeta3D(data = RP_25_40m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)
functional_beta_TM_25_40m_0 <- iNEXTbeta3D(data = TM_25_40m_functional, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = 'tau_value', FDcut_number = 30)

# again had problems with extracting data from the output. so did it manually in excel. ask r group!

# Read in summary CSV, adjust factors

functional_tax_beta <- read.csv("functional_beta_groupeddepth_byyear.csv")
functional_tax_beta$Site <- factor(functional_tax_beta$Site)
functional_tax_beta$Depth <- recode(functional_tax_beta$Depth, "5_15" = "5-15", "25_40" = "25-40")
functional_tax_beta$Depth <- factor(functional_tax_beta$Depth, levels = c("5-15", "25-40"))

# Plot for q = 0
functional_tax_beta_0_plot <- ggplot(functional_tax_beta, aes(x = as.factor(Site), y = Beta_q0, color = Depth, group = interaction(Site, Depth))) +
  geom_point(size = 2, position = position_dodge(width = -0.45)) +
  geom_errorbar(aes(ymin = Beta_LCL_q0, ymax = Beta_UCL_q0), position = position_dodge(width = -0.45), width = 0.3, lwd = 0.5) +
  labs(x = NULL, y = NULL) +
  scale_color_manual(values = c("5-15" = "red", "25-40" = "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right",
        text = element_text(size = 12, family = "sans", color = "black"),
        axis.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_flip() +
  scale_y_continuous(limits = c(0.8, 1.9), breaks = c(1.0, 1.4, 1.8))

# Plot for q = 2
functional_tax_beta_2_plot <- ggplot(functional_tax_beta, aes(x = as.factor(Site), y = Beta_q2, color = Depth, group = interaction(Site, Depth))) +
  geom_point(size = 2, position = position_dodge(width = -0.45)) +
  geom_errorbar(aes(ymin = Beta_LCL_q2, ymax = Beta_UCL_q2), position = position_dodge(width = -0.45), width = 0.3, lwd = 0.5) +
  labs(x = NULL, y = NULL) +
  scale_color_manual(name = "Depth (m)", values = c("5-15" = "red", "25-40" = "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right",
        text = element_text(size = 12, family = "sans", color = "black"),
        axis.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_flip() +
  scale_y_continuous(limits = c(0.8, 1.9), breaks = c(1.0, 1.4, 1.8))

# Side-by-side summary plot
temporal_site_functional_plot <- functional_tax_beta_0_plot + functional_tax_beta_2_plot +
  plot_layout(ncol = 2, nrow = 1, axis_titles = "collect")










################################################################################
##################### Dissimilarity Matrices (functional) ######################
################################################################################

############### Coral View Functional dissim matrices ##########################


# Select paired datasets


CV_5_15m_2014_2015_vs_5_15m_2022_2023 <- abun_combined_sites_functional %>% select(CV_5_15m_2014_2015, CV_5_15m_2022_2023)
CV_5_15m_2014_2015_vs_25_40m_2014_2015 <- abun_combined_sites_functional %>% select(CV_5_15m_2014_2015, CV_25_40m_2014_2015)
CV_5_15m_2014_2015_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>% select(CV_5_15m_2014_2015, CV_25_40m_2022_2023)
CV_5_15m_2022_2023_vs_25_40m_2014_2015 <- abun_combined_sites_functional %>% select(CV_5_15m_2022_2023, CV_25_40m_2014_2015)
CV_5_15m_2022_2023_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>% select(CV_5_15m_2022_2023, CV_25_40m_2022_2023)
CV_25_40m_2014_2015_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>% select(CV_25_40m_2014_2015, CV_25_40m_2022_2023)
  


##### q = 2 


functional_beta_CV_5_15m_2014_2015_vs_5_15m_2022_2023_2 <- iNEXTbeta3D(
  data = CV_5_15m_2014_2015_vs_5_15m_2022_2023, diversity = "FD", datatype = "abundance",
  base = "coverage", q = 2, level = 1, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_CV_5_15m_2014_2015_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data = CV_5_15m_2014_2015_vs_25_40m_2014_2015, diversity = "FD", datatype = "abundance",
  base = "coverage", q = 2, level = 1, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_CV_5_15m_2014_2015_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data = CV_5_15m_2014_2015_vs_25_40m_2022_2023, diversity = "FD", datatype = "abundance",
  base = "coverage", q = 2, level = 1, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_CV_5_15m_2022_2023_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data = CV_5_15m_2022_2023_vs_25_40m_2014_2015, diversity = "FD", datatype = "abundance",
  base = "coverage", q = 2, level = 2, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_CV_5_15m_2022_2023_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data = CV_5_15m_2022_2023_vs_25_40m_2022_2023, diversity = "FD", datatype = "abundance",
  base = "coverage", q = 2, level = 1, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_CV_25_40m_2014_2015_vs_2022_2023_2 <- iNEXTbeta3D(
  data = CV_25_40m_2014_2015_vs_25_40m_2022_2023, diversity = "FD", datatype = "abundance",
  base = "coverage", q = 2, level = 1, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)


##### q = 0 

functional_beta_CV_5_15m_2014_2015_vs_5_15m_2022_2023_0 <- iNEXTbeta3D(
  data = CV_5_15m_2014_2015_vs_5_15m_2022_2023, diversity = "FD", datatype = "abundance",
  base = "coverage", q = 0, level = 0.9715554, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_CV_5_15m_2014_2015_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data = CV_5_15m_2014_2015_vs_25_40m_2014_2015, diversity = "FD", datatype = "abundance",
  base = "coverage", q = 0, level = 0.9715554, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_CV_5_15m_2014_2015_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data = CV_5_15m_2014_2015_vs_25_40m_2022_2023, diversity = "FD", datatype = "abundance",
  base = "coverage", q = 0, level = 0.9715554, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_CV_5_15m_2022_2023_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data = CV_5_15m_2022_2023_vs_25_40m_2014_2015, diversity = "FD", datatype = "abundance",
  base = "coverage", q = 0, level = 0.9715554, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_CV_5_15m_2022_2023_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data = CV_5_15m_2022_2023_vs_25_40m_2022_2023, diversity = "FD", datatype = "abundance",
  base = "coverage", q = 0, level = 0.9715554, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_CV_25_40m_2014_2015_vs_2022_2023_0 <- iNEXTbeta3D(
  data = CV_25_40m_2014_2015_vs_25_40m_2022_2023, diversity = "FD", datatype = "abundance",
  base = "coverage", q = 0, level = 0.9715554, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)


# Dissimilarity Plotfor q = 2


cv_q2_dis <- read.csv("cv_q2_dissimilarity.csv", row.names = 1)

row_new_names <- c("5-15m (2014/15)", "5-15m (2022/23)", "25-40m (2014/15)")
column_new_names <- c("5-15m (2022/23)", "25-40m (2014/15)", "25-40m (2022/23)")

rownames(cv_q2_dis) <- row_new_names
colnames(cv_q2_dis) <- column_new_names

library(reshape2)
melted_cv_q2_dis <- melt(as.matrix(cv_q2_dis), na.rm = TRUE)

melted_cv_q2_dis$Var1 <- factor(melted_cv_q2_dis$Var1, levels = row_new_names)
melted_cv_q2_dis$Var2 <- factor(melted_cv_q2_dis$Var2, levels = rev(column_new_names))

cv_q2_func_dis_plot <- ggplot(data = melted_cv_q2_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "lightblue1", high = "black", limits = c(1, 2),
    name = expression(atop(Functional ~ beta-diversity, (q == 2)))
  ) +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(0, 1),
    legend.position = c(0.9, 0.7),
    legend.direction = "vertical",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7, title.position = "top", title.hjust = 0.5))


# Dissimilarity Plot — q = 0


cv_q0_dis <- read.csv("cv_q0_dissimilarity.csv", row.names = 1)

rownames(cv_q0_dis) <- row_new_names
colnames(cv_q0_dis) <- column_new_names

melted_cv_q0_dis <- melt(as.matrix(cv_q0_dis), na.rm = TRUE)

melted_cv_q0_dis$Var1 <- factor(melted_cv_q0_dis$Var1, levels = row_new_names)
melted_cv_q0_dis$Var2 <- factor(melted_cv_q0_dis$Var2, levels = rev(column_new_names))

cv_q0_func_dis_plot <- ggplot(data = melted_cv_q0_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "lightblue1", high = "black", limits = c(1, 2),
    name = expression(atop(Functional ~ beta-diversity, (q == 0)))
  ) +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(0, 1),
    legend.position = c(0.9, 0.7),
    legend.direction = "vertical",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7, title.position = "top", title.hjust = 0.5))


############### Little bight Functional dissim matrices ###########################################

# Data prep
LB_5_15m_2014_2015_vs_5_15m_2022_2023 <- abun_combined_sites_functional %>% select(LB_5_15m_2014_2015, LB_5_15m_2022_2023)
LB_5_15m_2014_2015_vs_25_40m_2014_2015 <- abun_combined_sites_functional %>% select(LB_5_15m_2014_2015, LB_25_40m_2014_2015)
LB_5_15m_2014_2015_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>% select(LB_5_15m_2014_2015, LB_25_40m_2022_2023)
LB_5_15m_2022_2023_vs_25_40m_2014_2015 <- abun_combined_sites_functional %>% select(LB_5_15m_2022_2023, LB_25_40m_2014_2015)
LB_5_15m_2022_2023_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>% select(LB_5_15m_2022_2023, LB_25_40m_2022_2023)
LB_25_40m_2014_2015_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>% select(LB_25_40m_2014_2015, LB_25_40m_2022_2023)

# q = 2
functional_beta_LB_5_15m_2014_2015_vs_5_15m_2022_2023_2 <- iNEXTbeta3D(data = LB_5_15m_2014_2015_vs_5_15m_2022_2023, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = "tau_value", FDcut_number = 30)
functional_beta_LB_5_15m_2014_2015_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(data = LB_5_15m_2014_2015_vs_25_40m_2014_2015, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = "tau_value", FDcut_number = 30)
functional_beta_LB_5_15m_2014_2015_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(data = LB_5_15m_2014_2015_vs_25_40m_2022_2023, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = "tau_value", FDcut_number = 30)
functional_beta_LB_5_15m_2022_2023_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(data = LB_5_15m_2022_2023_vs_25_40m_2014_2015, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = "tau_value", FDcut_number = 30)
functional_beta_LB_5_15m_2022_2023_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(data = LB_5_15m_2022_2023_vs_25_40m_2022_2023, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = "tau_value", FDcut_number = 30)
functional_beta_LB_25_40m_2014_2015_vs_2022_2023_2 <- iNEXTbeta3D(data = LB_25_40m_2014_2015_vs_25_40m_2022_2023, diversity = "FD", datatype = "abundance", base = "coverage", q = 2, level = 1, conf = 0.95, FDdistM = distM, FDtype = "tau_value", FDcut_number = 30)

# q = 0
functional_beta_LB_5_15m_2014_2015_vs_5_15m_2022_2023_0 <- iNEXTbeta3D(data = LB_5_15m_2014_2015_vs_5_15m_2022_2023, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = "tau_value", FDcut_number = 30)
functional_beta_LB_5_15m_2014_2015_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(data = LB_5_15m_2014_2015_vs_25_40m_2014_2015, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = "tau_value", FDcut_number = 30)
functional_beta_LB_5_15m_2014_2015_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(data = LB_5_15m_2014_2015_vs_25_40m_2022_2023, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = "tau_value", FDcut_number = 30)
functional_beta_LB_5_15m_2022_2023_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(data = LB_5_15m_2022_2023_vs_25_40m_2014_2015, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = "tau_value", FDcut_number = 30)
functional_beta_LB_5_15m_2022_2023_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(data = LB_5_15m_2022_2023_vs_25_40m_2022_2023, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = "tau_value", FDcut_number = 30)
functional_beta_LB_25_40m_2014_2015_vs_2022_2023_0 <- iNEXTbeta3D(data = LB_25_40m_2014_2015_vs_25_40m_2022_2023, diversity = "FD", datatype = "abundance", base = "coverage", q = 0, level = 0.9715554, conf = 0.95, FDdistM = distM, FDtype = "tau_value", FDcut_number = 30)

# q = 2 plot
lb_q2_dis <- read.csv("lb_q2_dissimilarity.csv", row.names = 1)
row_new_names <- c("5-15m (2014/15)", "5-15m (2022/23)", "25-40m (2014/15)")
column_new_names <- c("5-15m (2022/23)", "25-40m (2014/15)", "25-40m (2022/23)")
rownames(lb_q2_dis) <- row_new_names
colnames(lb_q2_dis) <- column_new_names
melted_lb_q2_dis <- melt(as.matrix(lb_q2_dis), na.rm = TRUE)
melted_lb_q2_dis$Var1 <- factor(melted_lb_q2_dis$Var1, levels = row_new_names)
melted_lb_q2_dis$Var2 <- factor(melted_lb_q2_dis$Var2, levels = rev(column_new_names))
lb_q2_func_dis_plot <- ggplot(data = melted_lb_q2_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue1", high = "black", limits = c(1, 2), name = "q = 2 Beta") +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.ticks = element_blank(), legend.justification = c(1, 0), legend.position = c(0.9, 0.3), legend.direction = "horizontal", legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))

# q = 0 plot
lb_q0_dis <- read.csv("lb_q0_dissimilarity.csv", row.names = 1)
rownames(lb_q0_dis) <- row_new_names
colnames(lb_q0_dis) <- column_new_names
melted_lb_q0_dis <- melt(as.matrix(lb_q0_dis), na.rm = TRUE)
melted_lb_q0_dis$Var1 <- factor(melted_lb_q0_dis$Var1, levels = row_new_names)
melted_lb_q0_dis$Var2 <- factor(melted_lb_q0_dis$Var2, levels = rev(column_new_names))
lb_q0_func_dis_plot <- ggplot(data = melted_lb_q0_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue1", high = "black", limits = c(1, 2), name = "q = 0 Beta") +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.ticks = element_blank(), legend.justification = c(1, 0), legend.position = c(0.9, 0.3), legend.direction = "horizontal", legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))



############### Raggedy Cay Functional dissim matrices ###########################################

########## Same for Rggedy Cay ########################

# Select paired site data
RC_5_15m_2014_2015_vs_5_15m_2022_2023 <- abun_combined_sites_functional %>%
  select(RC_5_15m_2014_2015, RC_5_15m_2022_2023)

RC_5_15m_2014_2015_vs_25_40m_2014_2015 <- abun_combined_sites_functional %>%
  select(RC_5_15m_2014_2015, RC_25_40m_2014_2015)

RC_5_15m_2014_2015_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>%
  select(RC_5_15m_2014_2015, RC_25_40m_2022_2023)

RC_5_15m_2022_2023_vs_25_40m_2014_2015 <- abun_combined_sites_functional %>%
  select(RC_5_15m_2022_2023, RC_25_40m_2014_2015)

RC_5_15m_2022_2023_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>%
  select(RC_5_15m_2022_2023, RC_25_40m_2022_2023)

RC_25_40m_2014_2015_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>%
  select(RC_25_40m_2014_2015, RC_25_40m_2022_2023)


# Functional beta diversity for q = 2
functional_beta_RC_5_15m_2014_2015_vs_5_15m_2022_2023_2 <- iNEXTbeta3D(
  data = RC_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "FD", datatype = "abundance", base = "coverage",
  q = 2, level = 1, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_RC_5_15m_2014_2015_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data = RC_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "FD", datatype = "abundance", base = "coverage",
  q = 2, level = 1, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_RC_5_15m_2014_2015_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data = RC_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "FD", datatype = "abundance", base = "coverage",
  q = 2, level = 1, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_RC_5_15m_2022_2023_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data = RC_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "FD", datatype = "abundance", base = "coverage",
  q = 2, level = 1, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_RC_5_15m_2022_2023_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data = RC_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "FD", datatype = "abundance", base = "coverage",
  q = 2, level = 1, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_RC_25_40m_2014_2015_vs_2022_2023_2 <- iNEXTbeta3D(
  data = RC_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "FD", datatype = "abundance", base = "coverage",
  q = 2, level = 1, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)



# Functional beta diversity for q = 0

functional_beta_RC_5_15m_2014_2015_vs_5_15m_2022_2023_0 <- iNEXTbeta3D(
  data = RC_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "FD", datatype = "abundance", base = "coverage",
  q = 0, level = 0.9715554, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_RC_5_15m_2014_2015_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data = RC_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "FD", datatype = "abundance", base = "coverage",
  q = 0, level = 0.9715554, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_RC_5_15m_2014_2015_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data = RC_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "FD", datatype = "abundance", base = "coverage",
  q = 0, level = 0.9715554, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_RC_5_15m_2022_2023_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data = RC_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "FD", datatype = "abundance", base = "coverage",
  q = 0, level = 0.9715554, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_RC_5_15m_2022_2023_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data = RC_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "FD", datatype = "abundance", base = "coverage",
  q = 0, level = 0.9715554, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

functional_beta_RC_25_40m_2014_2015_vs_2022_2023_0 <- iNEXTbeta3D(
  data = RC_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "FD", datatype = "abundance", base = "coverage",
  q = 0, level = 0.9715554, conf = 0.95,
  FDdistM = distM, FDtype = "tau_value", FDcut_number = 30
)

# Load and prepare q = 2 dissimilarity matrix for heatmap
rc_q2_dis <- read.csv("rc_q2_dissimilarity.csv", row.names = 1)

row_new_names <- c("5-15m (2014/15)", "5-15m (2022/23)", "25-40m (2014/15)")
column_new_names <- c("5-15m (2022/23)", "25-40m (2014/15)", "25-40m (2022/23)")

rownames(rc_q2_dis) <- row_new_names
colnames(rc_q2_dis) <- column_new_names

library(reshape2)
melted_rc_q2_dis <- melt(as.matrix(rc_q2_dis), na.rm = TRUE)

# Set order
melted_rc_q2_dis$Var1 <- factor(melted_rc_q2_dis$Var1, levels = new_names)
melted_rc_q2_dis$Var2 <- factor(melted_rc_q2_dis$Var2, levels = rev(new_names))

# Plot q = 2
rc_q2_func_dis_plot <- ggplot(data = melted_rc_q2_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue1", high = "black", limits = c(1, 2), name = "q = 2 Beta") +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.3),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))

# q = 0 plot
rc_q0_dis <- read.csv("rc_q0_dissimilarity.csv", row.names = 1)

rownames(rc_q0_dis) <- row_new_names
colnames(rc_q0_dis) <- column_new_names

melted_rc_q0_dis <- melt(as.matrix(rc_q0_dis), na.rm = TRUE)
melted_rc_q0_dis$Var1 <- factor(melted_rc_q0_dis$Var1, levels = new_names)
melted_rc_q0_dis$Var2 <- factor(melted_rc_q0_dis$Var2, levels = rev(new_names))

rc_q0_func_dis_plot <- ggplot(data = melted_rc_q0_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue1", high = "black", limits = c(1, 2), name = "q = 0 Beta") +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.3),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))



############### Functional dissim matrices for rocky point #########################################

RP_5_15m_2014_2015_vs_5_15m_2022_2023 <- abun_combined_sites_functional %>%
  select(RP_5_15m_2014_2015, RP_5_15m_2022_2023)

RP_5_15m_2014_2015_vs_25_40m_2014_2015 <- abun_combined_sites_functional %>%
  select(RP_5_15m_2014_2015, RP_25_40m_2014_2015)

RP_5_15m_2014_2015_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>%
  select(RP_5_15m_2014_2015, RP_25_40m_2022_2023)

RP_5_15m_2022_2023_vs_25_40m_2014_2015 <- abun_combined_sites_functional %>%
  select(RP_5_15m_2022_2023, RP_25_40m_2014_2015)

RP_5_15m_2022_2023_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>%
  select(RP_5_15m_2022_2023, RP_25_40m_2022_2023)

RP_25_40m_2014_2015_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>%
  select(RP_25_40m_2014_2015, RP_25_40m_2022_2023)

functional_beta_RP_5_15m_2014_2015_vs_5_15m_2022_2023_2 <- iNEXTbeta3D(
  data = RP_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_RP_5_15m_2014_2015_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data = RP_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_RP_5_15m_2014_2015_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data = RP_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_RP_5_15m_2022_2023_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data = RP_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_RP_5_15m_2022_2023_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data = RP_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_RP_25_40m_2014_2015_vs_2022_2023_2 <- iNEXTbeta3D(
  data = RP_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)


# q=0 

functional_beta_RP_5_15m_2014_2015_vs_5_15m_2022_2023_0 <- iNEXTbeta3D(
  data = RP_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_RP_5_15m_2014_2015_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data = RP_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_RP_5_15m_2014_2015_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data = RP_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_RP_5_15m_2022_2023_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data = RP_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_RP_5_15m_2022_2023_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data = RP_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_RP_25_40m_2014_2015_vs_2022_2023_0 <- iNEXTbeta3D(
  data = RP_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)



rp_q2_dis <-read.csv("rp_q2_dissimilarity.csv", row.names = 1) 

row_new_names <- c("5-15m (2014/15)", "5-15m (2022/23)", "25-40m (2014/15)")  # Make sure this matches the dimensions
column_new_names <- c("5-15m (2022/23)", "25-40m (2014/15)", "25-40m (2022/23)")  # Make sure this matches the dimensions


rownames(rp_q2_dis) <- row_new_names
colnames(rp_q2_dis) <- column_new_names

library(reshape2)

melted_rp_q2_dis <- melt(as.matrix(rp_q2_dis), na.rm = TRUE)

# Set the order explicitly
melted_rp_q2_dis$Var1 <- factor(melted_rp_q2_dis$Var1, levels = new_names)
melted_rp_q2_dis$Var2 <- factor(melted_rp_q2_dis$Var2, levels = rev(new_names))

rp_q2_func_dis_plot <-ggplot(data = melted_rp_q2_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "lightblue1",
    high = "black",
    limits = c(1, 2),
    name = "q = 2 Beta"
  ) +
  theme_minimal(base_size = 12) +  # Controls base font size for all text
  coord_fixed() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.3),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(
    barwidth = 7, barheight = 1,
    title.position = "top", title.hjust = 0.5
  ))

# q = 0 plot 

rp_q0_dis <-read.csv("rp_q0_dissimilarity.csv", row.names = 1) 

row_new_names <- c("5-15m (2014/15)", "5-15m (2022/23)", "25-40m (2014/15)")  # Make sure this matches the dimensions
column_new_names <- c("5-15m (2022/23)", "25-40m (2014/15)", "25-40m (2022/23)")  # Make sure this matches the dimensions


rownames(rp_q0_dis) <- row_new_names
colnames(rp_q0_dis) <- column_new_names


library(reshape2)

melted_rp_q0_dis <- melt(as.matrix(rp_q0_dis), na.rm = TRUE)

# Set the order explicitly
melted_rp_q0_dis$Var1 <- factor(melted_rp_q0_dis$Var1, levels = new_names)
melted_rp_q0_dis$Var2 <- factor(melted_rp_q0_dis$Var2, levels = rev(new_names))

rp_q0_func_dis_plot <-ggplot(data = melted_rp_q0_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "lightblue1",
    high = "black",
    limits = c(1, 2),
    name = "q = 0 Beta"
  ) +
  theme_minimal(base_size = 12) +  # Controls base font size for all text
  coord_fixed() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.3),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(
    barwidth = 7, barheight = 1,
    title.position = "top", title.hjust = 0.5
  ))

####### Functional dissim matrices for the maze ################################

TM_5_15m_2014_2015_vs_5_15m_2022_2023 <- abun_combined_sites_functional %>%
  select(TM_5_15m_2014_2015, TM_5_15m_2022_2023)

TM_5_15m_2014_2015_vs_25_40m_2014_2015 <- abun_combined_sites_functional %>%
  select(TM_5_15m_2014_2015, TM_25_40m_2014_2015)

TM_5_15m_2014_2015_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>%
  select(TM_5_15m_2014_2015, TM_25_40m_2022_2023)

TM_5_15m_2022_2023_vs_25_40m_2014_2015 <- abun_combined_sites_functional %>%
  select(TM_5_15m_2022_2023, TM_25_40m_2014_2015)

TM_5_15m_2022_2023_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>%
  select(TM_5_15m_2022_2023, TM_25_40m_2022_2023)

TM_25_40m_2014_2015_vs_25_40m_2022_2023 <- abun_combined_sites_functional %>%
  select(TM_25_40m_2014_2015, TM_25_40m_2022_2023)

functional_beta_TM_5_15m_2014_2015_vs_5_15m_2022_2023_2 <- iNEXTbeta3D(
  data = TM_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_TM_5_15m_2014_2015_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data = TM_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_TM_5_15m_2014_2015_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data = TM_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_TM_5_15m_2022_2023_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data = TM_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_TM_5_15m_2022_2023_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data = TM_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_TM_25_40m_2014_2015_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data = TM_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

# q = 0


functional_beta_TM_5_15m_2014_2015_vs_5_15m_2022_2023_0 <- iNEXTbeta3D(
  data = TM_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_TM_5_15m_2014_2015_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data = TM_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_TM_5_15m_2014_2015_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data = TM_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_TM_5_15m_2022_2023_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data = TM_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_TM_5_15m_2022_2023_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data = TM_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

functional_beta_TM_25_40m_2014_2015_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data = TM_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "FD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95,
  FDdistM = distM,
  FDtype = "tau_value",
  FDcut_number = 30
)

tm_q2_dis <-read.csv("tm_q2_dissimilarity.csv", row.names = 1) 

row_new_names <- c("5-15m (2014/15)", "5-15m (2022/23)", "25-40m (2014/15)")  # Make sure this matches the dimensions
column_new_names <- c("5-15m (2022/23)", "25-40m (2014/15)", "25-40m (2022/23)")  # Make sure this matches the dimensions


rownames(tm_q2_dis) <- row_new_names
colnames(tm_q2_dis) <- column_new_names

library(reshape2)

melted_tm_q2_dis <- melt(as.matrix(tm_q2_dis), na.rm = TRUE)

# Set the order explicitly
melted_tm_q2_dis$Var1 <- factor(melted_tm_q2_dis$Var1, levels = new_names)
melted_tm_q2_dis$Var2 <- factor(melted_tm_q2_dis$Var2, levels = rev(new_names))

tm_q2_func_dis_plot <-ggplot(data = melted_tm_q2_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "lightblue1",
    high = "black",
    limits = c(1, 2),
    name = "q = 2 Beta"
  ) +
  theme_minimal(base_size = 12) +  # Controls base font size for all text
  coord_fixed() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.3),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(
    barwidth = 7, barheight = 1,
    title.position = "top", title.hjust = 0.5
  ))

# q = 0 

tm_q0_dis <-read.csv("tm_q0_dissimilarity.csv", row.names = 1) 

row_new_names <- c("5-15m (2014/15)", "5-15m (2022/23)", "25-40m (2014/15)")  # Make sure this matches the dimensions
column_new_names <- c("5-15m (2022/23)", "25-40m (2014/15)", "25-40m (2022/23)")  # Make sure this matches the dimensions


rownames(tm_q0_dis) <- row_new_names
colnames(tm_q0_dis) <- column_new_names


library(reshape2)

melted_tm_q0_dis <- melt(as.matrix(tm_q0_dis), na.rm = TRUE)

# Set the order explicitly
melted_tm_q0_dis$Var1 <- factor(melted_tm_q0_dis$Var1, levels = new_names)
melted_tm_q0_dis$Var2 <- factor(melted_tm_q0_dis$Var2, levels = rev(new_names))

tm_q0_func_dis_plot <-ggplot(data = melted_tm_q0_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "lightblue1",
    high = "black",
    limits = c(1, 2),
    name = "q = 0 Beta"
  ) +
  theme_minimal(base_size = 12) +  # Controls base font size for all text
  coord_fixed() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.3),
    legend.direction = "horizontal",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(
    barwidth = 7, barheight = 1,
    title.position = "top", title.hjust = 0.5
  ))


############# combine plots from the functional dissim matrices  ################
############### Load Required Libraries ###############
library(patchwork)
library(cowplot)

############### q = 2 Functional Dissimilarity Plots ###############

# Add titles and bold styling
cv_q2_func_dis_plot <- cv_q2_func_dis_plot +
  ggtitle("Coral View") +
  theme(plot.title = element_text(face = "bold"))

lb_q2_func_dis_plot <- lb_q2_func_dis_plot +
  ggtitle("Little Bight") +
  theme(plot.title = element_text(face = "bold"))

rc_q2_func_dis_plot <- rc_q2_func_dis_plot +
  ggtitle("Raggedy Cay") +
  theme(plot.title = element_text(face = "bold"))

rp_q2_func_dis_plot <- rp_q2_func_dis_plot +
  ggtitle("Rocky Point") +
  theme(plot.title = element_text(face = "bold"))

tm_q2_func_dis_plot <- tm_q2_func_dis_plot +
  ggtitle("The Maze") +
  theme(plot.title = element_text(face = "bold"))

# Remove legends for combining
cv_plot_noleg <- cv_q2_func_dis_plot + theme(legend.position = "none")
lb_plot_noleg <- lb_q2_func_dis_plot + theme(legend.position = "none")
rc_plot_noleg <- rc_q2_func_dis_plot + theme(legend.position = "none")
rp_plot_noleg <- rp_q2_func_dis_plot + theme(legend.position = "none")
tm_plot_noleg <- tm_q2_func_dis_plot + theme(legend.position = "none")

# Extract legend from one plot
legend_plot <- cv_q2_func_dis_plot + theme(legend.position = "right")
legend <- get_legend(legend_plot)
legend_patch <- ggdraw(legend)

# Combine plots and legend
combined_plot_f2 <- (
  cv_plot_noleg + lb_plot_noleg + rc_plot_noleg +
    rp_plot_noleg + tm_plot_noleg + legend_patch
) + plot_layout(ncol = 3, nrow = 2)



############### q = 0 Functional Dissimilarity Plots ###############

# Add titles and bold styling
cv_q0_func_dis_plot <- cv_q0_func_dis_plot +
  ggtitle("Coral View") +
  theme(plot.title = element_text(face = "bold"))

lb_q0_func_dis_plot <- lb_q0_func_dis_plot +
  ggtitle("Little Bight") +
  theme(plot.title = element_text(face = "bold"))

rc_q0_func_dis_plot <- rc_q0_func_dis_plot +
  ggtitle("Raggedy Cay") +
  theme(plot.title = element_text(face = "bold"))

rp_q0_func_dis_plot <- rp_q0_func_dis_plot +
  ggtitle("Rocky Point") +
  theme(plot.title = element_text(face = "bold"))

tm_q0_func_dis_plot <- tm_q0_func_dis_plot +
  ggtitle("The Maze") +
  theme(plot.title = element_text(face = "bold"))

# Remove legends for combining
cv_plot_noleg <- cv_q0_func_dis_plot + theme(legend.position = "none")
lb_plot_noleg <- lb_q0_func_dis_plot + theme(legend.position = "none")
rc_plot_noleg <- rc_q0_func_dis_plot + theme(legend.position = "none")
rp_plot_noleg <- rp_q0_func_dis_plot + theme(legend.position = "none")
tm_plot_noleg <- tm_q0_func_dis_plot + theme(legend.position = "none")

# Extract legend from one plot
legend_plot <- cv_q0_func_dis_plot + theme(legend.position = "right")
legend <- get_legend(legend_plot)
legend_patch <- ggdraw(legend)

# Combine plots and legend
combined_plot_f0 <- (
  cv_plot_noleg + lb_plot_noleg + rc_plot_noleg +
    rp_plot_noleg + tm_plot_noleg + legend_patch
) + plot_layout(ncol = 3, nrow = 2)



#####################################################
# COMBINE ALL FUNCTIONAL BETA PLOTS
#####################################################

# Adjust margins for top and bottom plots
functional_tax_beta_0_plot <- functional_tax_beta_0_plot & 
  theme(plot.margin = unit(c(0, 0, 20, 0), "pt"))  # Top, right, bottom, left

functional_tax_beta_2_plot <- functional_tax_beta_2_plot & 
  theme(plot.margin = unit(c(20, 0, 0, 0), "pt"))  # Top, right, bottom, left

# Combine taxonomic beta diversity plots (q = 0 on top of q = 2)
combined_func_diverse <- (
  functional_tax_beta_0_plot / functional_tax_beta_2_plot
)

# Adjust margins for functional dissimilarity plots
combined_plot_f0 <- combined_plot_f0 & 
  theme(plot.margin = unit(c(0, 0, 20, 0), "pt"))  # Top, right, bottom, left

combined_plot_f2 <- combined_plot_f2 & 
  theme(plot.margin = unit(c(20, 0, 0, 0), "pt"))  # Top, right, bottom, left

# Combine functional dissimilarity plots (q = 0 on top of q = 2)
combined_funcbeta <- combined_plot_f0 / combined_plot_f2

# Final combined plot: functional beta temporal (left) and functional beta matrices (right)
combined_func_all <- (
  combined_func_diverse |
    combined_funcbeta
) +
  plot_layout(
    ncol = 2,
    widths = c(1, 5),
    heights = c(1, 1)
  )

# Final dimensions: 11.5 x 11.5 (portrait)








################################################################################
##################### Dissimilarity Matrices (taxonomic) ######################################
################################################################################

############### Coral View taxonomic dissim matrices ###########################################

# Create depth/time subset datasets


CV_5_15m_2014_2015_vs_5_15m_2022_2023 <- combined_depths %>%
  select(CV_5_15m_2014_2015, CV_5_15m_2022_2023)

CV_5_15m_2014_2015_vs_25_40m_2014_2015 <- combined_depths %>%
  select(CV_5_15m_2014_2015, CV_25_40m_2014_2015)

CV_5_15m_2014_2015_vs_25_40m_2022_2023 <- combined_depths %>%
  select(CV_5_15m_2014_2015, CV_25_40m_2022_2023)

CV_5_15m_2022_2023_vs_25_40m_2014_2015 <- combined_depths %>%
  select(CV_5_15m_2022_2023, CV_25_40m_2014_2015)

CV_5_15m_2022_2023_vs_25_40m_2022_2023 <- combined_depths %>%
  select(CV_5_15m_2022_2023, CV_25_40m_2022_2023)

CV_25_40m_2014_2015_vs_25_40m_2022_2023 <- combined_depths %>%
  select(CV_25_40m_2014_2015, CV_25_40m_2022_2023)



# Taxonomic Dissimilarity (q = 2) 


tax_beta_CV_5_15m_2014_2015_vs_5_15m_2022_2023_2 <- iNEXTbeta3D(
  data      = CV_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_CV_5_15m_2014_2015_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data      = CV_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_CV_5_15m_2014_2015_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data      = CV_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_CV_5_15m_2022_2023_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data      = CV_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_CV_5_15m_2022_2023_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data      = CV_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_CV_25_40m_2014_2015_vs_2022_2023_2 <- iNEXTbeta3D(
  data      = CV_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)



# Taxonomic Dissimilarity (q = 0) - Species Richness Based


tax_beta_CV_5_15m_2014_2015_vs_5_15m_2022_2023_0 <- iNEXTbeta3D(
  data      = CV_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_CV_5_15m_2014_2015_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data      = CV_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_CV_5_15m_2014_2015_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data      = CV_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_CV_5_15m_2022_2023_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data      = CV_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_CV_5_15m_2022_2023_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data      = CV_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_CV_25_40m_2014_2015_vs_2022_2023_0 <- iNEXTbeta3D(
  data      = CV_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)


############### Little bight taxonomic dissim matrices ###########################################

# Data Selection


LB_5_15m_2014_2015_vs_5_15m_2022_2023 <- combined_depths %>%
  select(LB_5_15m_2014_2015, LB_5_15m_2022_2023)

LB_5_15m_2014_2015_vs_25_40m_2014_2015 <- combined_depths %>%
  select(LB_5_15m_2014_2015, LB_25_40m_2014_2015)

LB_5_15m_2014_2015_vs_25_40m_2022_2023 <- combined_depths %>%
  select(LB_5_15m_2014_2015, LB_25_40m_2022_2023)

LB_5_15m_2022_2023_vs_25_40m_2014_2015 <- combined_depths %>%
  select(LB_5_15m_2022_2023, LB_25_40m_2014_2015)

LB_5_15m_2022_2023_vs_25_40m_2022_2023 <- combined_depths %>%
  select(LB_5_15m_2022_2023, LB_25_40m_2022_2023)

LB_25_40m_2014_2015_vs_25_40m_2022_2023 <- combined_depths %>%
  select(LB_25_40m_2014_2015, LB_25_40m_2022_2023)


# Taxonomic Beta Diversity (q = 2)


tax_beta_LB_5_15m_2014_2015_vs_5_15m_2022_2023_2 <- iNEXTbeta3D(
  data      = LB_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_LB_5_15m_2014_2015_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data      = LB_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_LB_5_15m_2014_2015_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data      = LB_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_LB_5_15m_2022_2023_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data      = LB_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_LB_5_15m_2022_2023_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data      = LB_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_LB_25_40m_2014_2015_vs_2022_2023_2 <- iNEXTbeta3D(
  data      = LB_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)



# Taxonomic Beta Diversity (q = 0)


tax_beta_LB_5_15m_2014_2015_vs_5_15m_2022_2023_0 <- iNEXTbeta3D(
  data      = LB_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_LB_5_15m_2014_2015_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data      = LB_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_LB_5_15m_2014_2015_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data      = LB_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_LB_5_15m_2022_2023_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data      = LB_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_LB_5_15m_2022_2023_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data      = LB_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_LB_25_40m_2014_2015_vs_2022_2023_0 <- iNEXTbeta3D(
  data      = LB_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)



############### Raggedy cay taxonomic dissim matrices ###########################################

#  Data Selection 
RC_5_15m_2014_2015_vs_5_15m_2022_2023 <- combined_depths %>%
  select(RC_5_15m_2014_2015, RC_5_15m_2022_2023)

RC_5_15m_2014_2015_vs_25_40m_2014_2015 <- combined_depths %>%
  select(RC_5_15m_2014_2015, RC_25_40m_2014_2015)

RC_5_15m_2014_2015_vs_25_40m_2022_2023 <- combined_depths %>%
  select(RC_5_15m_2014_2015, RC_25_40m_2022_2023)

RC_5_15m_2022_2023_vs_25_40m_2014_2015 <- combined_depths %>%
  select(RC_5_15m_2022_2023, RC_25_40m_2014_2015)

RC_5_15m_2022_2023_vs_25_40m_2022_2023 <- combined_depths %>%
  select(RC_5_15m_2022_2023, RC_25_40m_2022_2023)

RC_25_40m_2014_2015_vs_25_40m_2022_2023 <- combined_depths %>%
  select(RC_25_40m_2014_2015, RC_25_40m_2022_2023)

#Taxonomic Beta Diversity (q = 2) 

tax_beta_RC_5_15m_2014_2015_vs_5_15m_2022_2023_2 <- iNEXTbeta3D(
  data      = RC_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_RC_5_15m_2014_2015_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data      = RC_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_RC_5_15m_2014_2015_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data      = RC_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_RC_5_15m_2022_2023_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data      = RC_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_RC_5_15m_2022_2023_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data      = RC_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

tax_beta_RC_25_40m_2014_2015_vs_2022_2023_2 <- iNEXTbeta3D(
  data      = RC_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 2,
  level     = 1,
  conf      = 0.95
)

# Taxonomic Beta Diversity (q = 0) 

tax_beta_RC_5_15m_2014_2015_vs_5_15m_2022_2023_0 <- iNEXTbeta3D(
  data      = RC_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_RC_5_15m_2014_2015_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data      = RC_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_RC_5_15m_2014_2015_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data      = RC_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_RC_5_15m_2022_2023_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data      = RC_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_RC_5_15m_2022_2023_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data      = RC_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)

tax_beta_RC_25_40m_2014_2015_vs_2022_2023_0 <- iNEXTbeta3D(
  data      = RC_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype  = "abundance",
  base      = "coverage",
  q         = 0,
  level     = 0.9715554,
  conf      = 0.95
)


############## Tax dissim matrix data for rocky point ############################################

RP_5_15m_2014_2015_vs_5_15m_2022_2023 <- combined_depths %>%
  select(RP_5_15m_2014_2015, RP_5_15m_2022_2023)

RP_5_15m_2014_2015_vs_25_40m_2014_2015 <- combined_depths %>%
  select(RP_5_15m_2014_2015, RP_25_40m_2014_2015)

RP_5_15m_2014_2015_vs_25_40m_2022_2023 <- combined_depths %>%
  select(RP_5_15m_2014_2015, RP_25_40m_2022_2023)

RP_5_15m_2022_2023_vs_25_40m_2014_2015 <- combined_depths %>%
  select(RP_5_15m_2022_2023, RP_25_40m_2014_2015)

RP_5_15m_2022_2023_vs_25_40m_2022_2023 <- combined_depths %>%
  select(RP_5_15m_2022_2023, RP_25_40m_2022_2023)

RP_25_40m_2014_2015_vs_25_40m_2022_2023 <- combined_depths %>%
  select(RP_25_40m_2014_2015, RP_25_40m_2022_2023)

tax_beta_RP_5_15m_2014_2015_vs_5_15m_2022_2023_2 <- iNEXTbeta3D(
  data = RP_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95
)

tax_beta_RP_5_15m_2014_2015_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data = RP_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95
)

tax_beta_RP_5_15m_2014_2015_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data = RP_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95
)

tax_beta_RP_5_15m_2022_2023_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data = RP_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95
)

tax_beta_RP_5_15m_2022_2023_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data = RP_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95
)

tax_beta_RP_25_40m_2014_2015_vs_2022_2023_2 <- iNEXTbeta3D(
  data = RP_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95
)

# q = 0

tax_beta_RP_5_15m_2014_2015_vs_5_15m_2022_2023_0 <- iNEXTbeta3D(
  data = RP_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95
)

tax_beta_RP_5_15m_2014_2015_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data = RP_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95
)

tax_beta_RP_5_15m_2014_2015_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data = RP_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95
)

tax_beta_RP_5_15m_2022_2023_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data = RP_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95
)

tax_beta_RP_5_15m_2022_2023_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data = RP_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95
)

tax_beta_RP_25_40m_2014_2015_vs_2022_2023_0 <- iNEXTbeta3D(
  data = RP_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95
)


############## tax dissim matrix data for the maze ####################################

TM_5_15m_2014_2015_vs_5_15m_2022_2023 <- combined_depths %>%
  select(TM_5_15m_2014_2015, TM_5_15m_2022_2023)

TM_5_15m_2014_2015_vs_25_40m_2014_2015 <- combined_depths %>%
  select(TM_5_15m_2014_2015, TM_25_40m_2014_2015)

TM_5_15m_2014_2015_vs_25_40m_2022_2023 <- combined_depths %>%
  select(TM_5_15m_2014_2015, TM_25_40m_2022_2023)

TM_5_15m_2022_2023_vs_25_40m_2014_2015 <- combined_depths %>%
  select(TM_5_15m_2022_2023, TM_25_40m_2014_2015)

TM_5_15m_2022_2023_vs_25_40m_2022_2023 <- combined_depths %>%
  select(TM_5_15m_2022_2023, TM_25_40m_2022_2023)

TM_25_40m_2014_2015_vs_25_40m_2022_2023 <- combined_depths %>%
  select(TM_25_40m_2014_2015, TM_25_40m_2022_2023)

tax_beta_TM_5_15m_2014_2015_vs_5_15m_2022_2023_2 <- iNEXTbeta3D(
  data = TM_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95
)

tax_beta_TM_5_15m_2014_2015_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data = TM_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95
)

tax_beta_TM_5_15m_2014_2015_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data = TM_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95
)

tax_beta_TM_5_15m_2022_2023_vs_25_40m_2014_2015_2 <- iNEXTbeta3D(
  data = TM_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95
)

tax_beta_TM_5_15m_2022_2023_vs_25_40m_2022_2023_2 <- iNEXTbeta3D(
  data = TM_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95
)

tax_beta_TM_25_40m_2014_2015_vs_2022_2023_2 <- iNEXTbeta3D(
  data = TM_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 2,
  level = 1,
  conf = 0.95
)

# q = 0

tax_beta_TM_5_15m_2014_2015_vs_5_15m_2022_2023_0 <- iNEXTbeta3D(
  data = TM_5_15m_2014_2015_vs_5_15m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95
)

tax_beta_TM_5_15m_2014_2015_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data = TM_5_15m_2014_2015_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95
)

tax_beta_TM_5_15m_2014_2015_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data = TM_5_15m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95
)

tax_beta_TM_5_15m_2022_2023_vs_25_40m_2014_2015_0 <- iNEXTbeta3D(
  data = TM_5_15m_2022_2023_vs_25_40m_2014_2015,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95
)

tax_beta_TM_5_15m_2022_2023_vs_25_40m_2022_2023_0 <- iNEXTbeta3D(
  data = TM_5_15m_2022_2023_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95
)

tax_beta_TM_25_40m_2014_2015_vs_2022_2023_0 <- iNEXTbeta3D(
  data = TM_25_40m_2014_2015_vs_25_40m_2022_2023,
  diversity = "TD",
  datatype = "abundance",
  base = "coverage",
  q = 0,
  level = 0.9715554,
  conf = 0.95
)

###############################################################################

# Again, outputs from inextbeta3d were tricky to automatically put into a data frame (i couldn't work it out anyway)
# so i added them into individual excel spreadsheets 

###############################################################################


############## Make heatmaps for the Maze q = 2 #####################

#  Load Dissimilarity Matrix 
tax_tm_q2_dis <- read.csv("tax_tm_q2_dissimilarity.csv", row.names = 1) 

# Rename rows and columns (ensure dimensions match)
row_new_names    <- c("5-15m (2014/15)", "5-15m (2022/23)", "25-40m (2014/15)")  
column_new_names <- c("5-15m (2022/23)", "25-40m (2014/15)", "25-40m (2022/23)")  

rownames(tax_tm_q2_dis) <- row_new_names
colnames(tax_tm_q2_dis) <- column_new_names

#  Melt for ggplot 
library(reshape2)
melted_tax_tm_q2_dis <- melt(as.matrix(tax_tm_q2_dis), na.rm = TRUE)

# Explicit factor ordering for plot axes
melted_tax_tm_q2_dis$Var1 <- factor(melted_tax_tm_q2_dis$Var1, levels = new_names)
melted_tax_tm_q2_dis$Var2 <- factor(melted_tax_tm_q2_dis$Var2, levels = rev(new_names))

#  make Heatmap 
tm_q2_tax_dis_plot <- ggplot(data = melted_tax_tm_q2_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low    = "lightblue1",
    high   = "black",
    limits = c(1, 2),
    name   = "q = 2 Beta"
  ) +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y      = element_blank(),
    axis.title.x     = element_blank(),
    axis.title.y     = element_blank(),
    panel.grid.major = element_blank(),
    panel.border     = element_blank(),
    panel.background = element_blank(),
    axis.ticks       = element_blank(),
    legend.justification = c(1, 0),
    legend.position      = c(0.9, 0.3),
    legend.direction     = "horizontal",
    legend.title         = element_text(size = 12),
    legend.text          = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 7,
    barheight      = 1,
    title.position = "top",
    title.hjust    = 0.5
  ))


############## Make heatmaps for coral view q = 2 #####################

tax_cv_q2_dis <- read.csv("tax_cv_q2_dissimilarity.csv", row.names = 1) 

# Rename rows and columns
new_names <- c("5-15m (2014/15)", "5-15m (2022/23)", 
               "25-40m (2014/15)", "25-40m (2022/23)")  
rownames(tax_cv_q2_dis) <- new_names
colnames(tax_cv_q2_dis) <- new_names

# Melt for ggplot
library(reshape2)
melted_tax_cv_q2_dis <- melt(as.matrix(tax_cv_q2_dis), na.rm = TRUE)

# Set explicit order
melted_tax_cv_q2_dis$Var1 <- factor(melted_tax_cv_q2_dis$Var1, levels = new_names)
melted_tax_cv_q2_dis$Var2 <- factor(melted_tax_cv_q2_dis$Var2, levels = rev(new_names))

# Plot
cv_q2_tax_dis_plot <- ggplot(data = melted_tax_cv_q2_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low    = "lightblue1",
    high   = "black",
    limits = c(1, 2),
    name   = expression(atop(Taxonomic~beta-diversity, (q==2)))
  ) +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x       = element_blank(),
    axis.text.y       = element_text(size = 10, colour = "black"),
    axis.title.x      = element_blank(),
    axis.title.y      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.border      = element_blank(),
    panel.background  = element_blank(),
    axis.ticks        = element_blank(),
    legend.justification = c(0, 1),
    legend.position      = c(0.9, 0.7),
    legend.direction     = "vertical",
    legend.title         = element_text(size = 12),
    legend.text          = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 1,
    barheight      = 7,
    title.position = "top",
    title.hjust    = 0.8,
    title.vjust    = 5
  ))


############## Make heatmaps for the little bight q = 2 #####################

tax_lb_q2_dis <- read.csv("tax_lb_q2_dissimilarity.csv", row.names = 1) 

rownames(tax_lb_q2_dis) <- new_names
colnames(tax_lb_q2_dis) <- new_names

melted_tax_lb_q2_dis <- melt(as.matrix(tax_lb_q2_dis), na.rm = TRUE)
melted_tax_lb_q2_dis$Var1 <- factor(melted_tax_lb_q2_dis$Var1, levels = new_names)
melted_tax_lb_q2_dis$Var2 <- factor(melted_tax_lb_q2_dis$Var2, levels = rev(new_names))

lb_q2_tax_dis_plot <- ggplot(data = melted_tax_lb_q2_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low    = "lightblue1",
    high   = "black",
    limits = c(1, 2),
    name   = "q = 2 Beta"
  ) +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x       = element_blank(),
    axis.text.y       = element_blank(),
    axis.title.x      = element_blank(),
    axis.title.y      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.border      = element_blank(),
    panel.background  = element_blank(),
    axis.ticks        = element_blank(),
    legend.justification = c(1, 0),
    legend.position      = c(0.9, 0.3),
    legend.direction     = "horizontal",
    legend.title         = element_text(size = 12),
    legend.text          = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 7,
    barheight      = 1,
    title.position = "top",
    title.hjust    = 0.5
  ))


# ---- RC ----
tax_rc_q2_dis <- read.csv("tax_rc_q2_dissimilarity.csv", row.names = 1) 

rownames(tax_rc_q2_dis) <- new_names
colnames(tax_rc_q2_dis) <- new_names

melted_tax_rc_q2_dis <- melt(as.matrix(tax_rc_q2_dis), na.rm = TRUE)
melted_tax_rc_q2_dis$Var1 <- factor(melted_tax_rc_q2_dis$Var1, levels = new_names)
melted_tax_rc_q2_dis$Var2 <- factor(melted_tax_rc_q2_dis$Var2, levels = rev(new_names))

rc_q2_tax_dis_plot <- ggplot(data = melted_tax_rc_q2_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low    = "lightblue1",
    high   = "black",
    limits = c(1, 2),
    name   = "q = 2 Beta"
  ) +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y       = element_blank(),
    axis.title.x      = element_blank(),
    axis.title.y      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.border      = element_blank(),
    panel.background  = element_blank(),
    axis.ticks        = element_blank(),
    legend.justification = c(1, 0),
    legend.position      = c(0.9, 0.3),
    legend.direction     = "horizontal",
    legend.title         = element_text(size = 12),
    legend.text          = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 7,
    barheight      = 1,
    title.position = "top",
    title.hjust    = 0.5
  ))


############## Make heatmaps for rocky point q = 2 #####################
tax_rp_q2_dis <- read.csv("tax_rp_q2_dissimilarity.csv", row.names = 1) 

rownames(tax_rp_q2_dis) <- new_names
colnames(tax_rp_q2_dis) <- new_names

melted_tax_rp_q2_dis <- melt(as.matrix(tax_rp_q2_dis), na.rm = TRUE)
melted_tax_rp_q2_dis$Var1 <- factor(melted_tax_rp_q2_dis$Var1, levels = new_names)
melted_tax_rp_q2_dis$Var2 <- factor(melted_tax_rp_q2_dis$Var2, levels = rev(new_names))

rp_q2_tax_dis_plot <- ggplot(data = melted_tax_rp_q2_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low    = "lightblue1",
    high   = "black",
    limits = c(1, 2),
    name   = "q = 2 Beta"
  ) +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y       = element_text(size = 10, colour = "black"),
    axis.title.x      = element_blank(),
    axis.title.y      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.border      = element_blank(),
    panel.background  = element_blank(),
    axis.ticks        = element_blank(),
    legend.justification = c(1, 0),
    legend.position      = c(0.9, 0.3),
    legend.direction     = "horizontal",
    legend.title         = element_text(size = 12),
    legend.text          = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(
    barwidth       = 7,
    barheight      = 1,
    title.position = "top",
    title.hjust    = 0.5
  ))


##############################################################################
####### now make the heatmaps for q = 0. exact same process ##################
##############################################################################

# TM q = 0
tax_tm_q0_dis <- read.csv("tax_tm_q0_dissimilarity.csv", row.names = 1)

new_names <- c("5-15m (2014/15)", "5-15m (2022/23)",
               "25-40m (2014/15)", "25-40m (2022/23)")

rownames(tax_tm_q0_dis) <- new_names
colnames(tax_tm_q0_dis) <- new_names

melted_tax_tm_q0_dis <- melt(as.matrix(tax_tm_q0_dis), na.rm = TRUE)
melted_tax_tm_q0_dis$Var1 <- factor(melted_tax_tm_q0_dis$Var1, levels = new_names)
melted_tax_tm_q0_dis$Var2 <- factor(melted_tax_tm_q0_dis$Var2, levels = rev(new_names))

tm_q0_tax_dis_plot <- ggplot(melted_tax_tm_q0_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue1", high = "black",
                      limits = c(1, 2), name = "q = 0 Beta") +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major   = element_blank(),
    panel.border       = element_blank(),
    panel.background   = element_blank(),
    axis.ticks         = element_blank(),
    legend.justification = c(1, 0),
    legend.position      = c(0.9, 0.3),
    legend.direction     = "horizontal",
    legend.title         = element_text(size = 12),
    legend.text          = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

# CV q = 0
tax_cv_q0_dis <- read.csv("tax_cv_q0_dissimilarity.csv", row.names = 1)

rownames(tax_cv_q0_dis) <- new_names
colnames(tax_cv_q0_dis) <- new_names

melted_tax_cv_q0_dis <- melt(as.matrix(tax_cv_q0_dis), na.rm = TRUE)
melted_tax_cv_q0_dis$Var1 <- factor(melted_tax_cv_q0_dis$Var1, levels = new_names)
melted_tax_cv_q0_dis$Var2 <- factor(melted_tax_cv_q0_dis$Var2, levels = rev(new_names))

cv_q0_tax_dis_plot <- ggplot(melted_tax_cv_q0_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue1", high = "black",
                      limits = c(1, 2),
                      name = expression(atop(Taxonomic~beta-diversity, (q==0)))) +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(size = 10, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major   = element_blank(),
    panel.border       = element_blank(),
    panel.background   = element_blank(),
    axis.ticks         = element_blank(),
    legend.justification = c(0, 1),
    legend.position      = c(0.9, 0.7),
    legend.direction     = "vertical",
    legend.title         = element_text(size = 12),
    legend.text          = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,
                               title.position = "top", title.hjust = 0.8, title.vjust = 5))


# LB q = 0
tax_lb_q0_dis <- read.csv("tax_lb_q0_dissimilarity.csv", row.names = 1)

rownames(tax_lb_q0_dis) <- new_names
colnames(tax_lb_q0_dis) <- new_names

melted_tax_lb_q0_dis <- melt(as.matrix(tax_lb_q0_dis), na.rm = TRUE)
melted_tax_lb_q0_dis$Var1 <- factor(melted_tax_lb_q0_dis$Var1, levels = new_names)
melted_tax_lb_q0_dis$Var2 <- factor(melted_tax_lb_q0_dis$Var2, levels = rev(new_names))

lb_q0_tax_dis_plot <- ggplot(melted_tax_lb_q0_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue1", high = "black",
                      limits = c(1, 2), name = "q = 0 Beta") +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major   = element_blank(),
    panel.border       = element_blank(),
    panel.background   = element_blank(),
    axis.ticks         = element_blank(),
    legend.justification = c(1, 0),
    legend.position      = c(0.9, 0.3),
    legend.direction     = "horizontal",
    legend.title         = element_text(size = 12),
    legend.text          = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


# RC q = 0
tax_rc_q0_dis <- read.csv("tax_rc_q0_dissimilarity.csv", row.names = 1)

rownames(tax_rc_q0_dis) <- new_names
colnames(tax_rc_q0_dis) <- new_names

melted_tax_rc_q0_dis <- melt(as.matrix(tax_rc_q0_dis), na.rm = TRUE)
melted_tax_rc_q0_dis$Var1 <- factor(melted_tax_rc_q0_dis$Var1, levels = new_names)
melted_tax_rc_q0_dis$Var2 <- factor(melted_tax_rc_q0_dis$Var2, levels = rev(new_names))

rc_q0_tax_dis_plot <- ggplot(melted_tax_rc_q0_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue1", high = "black",
                      limits = c(1, 2), name = "q = 0 Beta") +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major   = element_blank(),
    panel.border       = element_blank(),
    panel.background   = element_blank(),
    axis.ticks         = element_blank(),
    legend.justification = c(1, 0),
    legend.position      = c(0.9, 0.3),
    legend.direction     = "horizontal",
    legend.title         = element_text(size = 12),
    legend.text          = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


# RP q = 0
tax_rp_q0_dis <- read.csv("tax_rp_q0_dissimilarity.csv", row.names = 1)

rownames(tax_rp_q0_dis) <- new_names
colnames(tax_rp_q0_dis) <- new_names

melted_tax_rp_q0_dis <- melt(as.matrix(tax_rp_q0_dis), na.rm = TRUE)
melted_tax_rp_q0_dis$Var1 <- factor(melted_tax_rp_q0_dis$Var1, levels = new_names)
melted_tax_rp_q0_dis$Var2 <- factor(melted_tax_rp_q0_dis$Var2, levels = rev(new_names))

rp_q0_tax_dis_plot <- ggplot(melted_tax_rp_q0_dis, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue1", high = "black",
                      limits = c(1, 2), name = "q = 2 Beta") +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y  = element_text(size = 10, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major   = element_blank(),
    panel.border       = element_blank(),
    panel.background   = element_blank(),
    axis.ticks         = element_blank(),
    legend.justification = c(1, 0),
    legend.position      = c(0.9, 0.3),
    legend.direction     = "horizontal",
    legend.title         = element_text(size = 12),
    legend.text          = element_text(size = 12)
  ) +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))



############### Combine Taxonomic Dissimilarity Plots #######################

library(patchwork)
library(cowplot)

### q = 0 Plots 

# Add titles to q = 0 plots
cv_q0_tax_dis_plot <- cv_q0_tax_dis_plot + ggtitle("Coral View") + theme(plot.title = element_text(face = "bold"))
lb_q0_tax_dis_plot <- lb_q0_tax_dis_plot + ggtitle("Little Bight") + theme(plot.title = element_text(face = "bold"))
rc_q0_tax_dis_plot <- rc_q0_tax_dis_plot + ggtitle("Raggedy Cay") + theme(plot.title = element_text(face = "bold"))
rp_q0_tax_dis_plot <- rp_q0_tax_dis_plot + ggtitle("Rocky Point") + theme(plot.title = element_text(face = "bold"))
tm_q0_tax_dis_plot <- tm_q0_tax_dis_plot + ggtitle("The Maze") + theme(plot.title = element_text(face = "bold"))

# Remove legends for combined layout
cv_plot_noleg <- cv_q0_tax_dis_plot + theme(legend.position = "none")
lb_plot_noleg <- lb_q0_tax_dis_plot + theme(legend.position = "none")
rc_plot_noleg <- rc_q0_tax_dis_plot + theme(legend.position = "none")
rp_plot_noleg <- rp_q0_tax_dis_plot + theme(legend.position = "none")
tm_plot_noleg <- tm_q0_tax_dis_plot + theme(legend.position = "none")

# take legend from one plot
legend_plot   <- cv_q0_tax_dis_plot + theme(legend.position = "right")
legend        <- get_legend(legend_plot)
legend_patch  <- ggdraw(legend)

# Combine q = 0 plots with legend
combined_plot_t0 <- (
  cv_plot_noleg + lb_plot_noleg + rc_plot_noleg +
    rp_plot_noleg + tm_plot_noleg + legend_patch
) +
  plot_layout(ncol = 3, nrow = 2)


###  q = 2 Plots 

# Add titles to q = 2 plots
cv_q2_tax_dis_plot <- cv_q2_tax_dis_plot + ggtitle("Coral View") + theme(plot.title = element_text(face = "bold"))
lb_q2_tax_dis_plot <- lb_q2_tax_dis_plot + ggtitle("Little Bight") + theme(plot.title = element_text(face = "bold"))
rc_q2_tax_dis_plot <- rc_q2_tax_dis_plot + ggtitle("Raggedy Cay") + theme(plot.title = element_text(face = "bold"))
rp_q2_tax_dis_plot <- rp_q2_tax_dis_plot + ggtitle("Rocky Point") + theme(plot.title = element_text(face = "bold"))
tm_q2_tax_dis_plot <- tm_q2_tax_dis_plot + ggtitle("The Maze") + theme(plot.title = element_text(face = "bold"))

# Remove legends for combined layout
cv_plot_noleg <- cv_q2_tax_dis_plot + theme(legend.position = "none")
lb_plot_noleg <- lb_q2_tax_dis_plot + theme(legend.position = "none")
rc_plot_noleg <- rc_q2_tax_dis_plot + theme(legend.position = "none")
rp_plot_noleg <- rp_q2_tax_dis_plot + theme(legend.position = "none")
tm_plot_noleg <- tm_q2_tax_dis_plot + theme(legend.position = "none")

# take legend from one plot
legend_plot   <- cv_q2_tax_dis_plot + theme(legend.position = "right")
legend        <- get_legend(legend_plot)
legend_patch  <- ggdraw(legend)

# Combine q = 2 plots with legend
combined_plot_t2 <- (
  cv_plot_noleg + lb_plot_noleg + rc_plot_noleg +
    rp_plot_noleg + tm_plot_noleg + legend_patch
) +
  plot_layout(ncol = 3, nrow = 2)


### Final Combined Layout!

# function to remove margins from plots (apply this to other plots to save time!)
clean_plot <- function(p) {
  p + theme(plot.margin = margin(0, 0, 0, 0))
}

# Clean margins
tax_beta_0_plot <- clean_plot(tax_beta_0_plot)
tax_beta_2_plot <- clean_plot(tax_beta_2_plot)
combined_plot_t0 <- clean_plot(combined_plot_t0)
combined_plot_t2 <- clean_plot(combined_plot_t2)

# Final combined plot: Left column = taxonomic beta plots, Right column = combined site dissimilarity plots
final_tax_dissim_combined_plot <- (
  (tax_beta_0_plot / tax_beta_2_plot) |
    (combined_plot_t0 / combined_plot_t2)
) +
  plot_layout(widths = c(0.5, 4)) # Adjust left-right width ratio

# Display
print(final_tax_dissim_combined_plot)