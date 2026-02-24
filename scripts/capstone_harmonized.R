# Setting up renv (project-local environment) ----
getRversion() # renv tested with R v4.5.1
if (require("renv") == FALSE) install.packages("renv")
renv::restore()

# Loading packages ----
library(ccrepe)
library(curatedMetagenomicData)
library(ggbeeswarm)
library(ggpubr)
library(janitor)
library(Maaslin2)
library(MMUPHin)
library(rstatix)
library(scales)
library(tidyverse)
library(vegan)

# Loading data ----
load("data/DIABIMMUNE_Karelia_metadata_sub.RData", verbose = TRUE)

###############
#Q1 ----
##Q1A Doing some initial eda on meta data ----
dim(metadata_all)               # 1946 rows 18 col
view(metadata_all)
colnames(metadata_all)          # 18 colnames ....

##Q1B Filtering the samples by shotgun (non-NA gid_wgs, by age and first ever collected sample) ----
metadata_filt <- metadata_all %>%
  filter(!is.na(gid_wgs), age_at_collection <= 365 ) %>%
  group_by(subjectID) %>%
  slice_min(age_at_collection , n = 1) %>%
  column_to_rownames("gid_wgs") # 162 individuals , first collected sample of each only

dim(metadata_filt) # 162  17


##Q1C Identifying species level taxonomy profiles corresponding to the filtered samples ----
### Loading taxa data from CMData from Bioconductor Experiment Hub  -----
cdata <- curatedMetagenomicData::curatedMetagenomicData(
  "VatanenT_2016.relative_abundance",
  dryrun = FALSE ,
  counts = FALSE )              #Large list (7.3 MB)

cdata <- cdata[[1]]
taxa_data_all <- assay(cdata)
rm(cdata)

dim(taxa_data_all)                  # 618 rows 785 cols
taxa_data_all[1:5 , 1:5]            # rownames = taxanomy  , colnames = gid_wgs ID

### Changing the taxa name to species level
ranks <- c("Kingdom", "Phylum", "Class", "Order",
           "Family", "Genus", "Species", "Strain")

# keep only species level data for ALL samples
species_data_all <- taxa_data_all %>%
  as.data.frame() %>%
  rownames_to_column("Org") %>%
  separate(
    Org,
    into = ranks,
    sep = "\\|",
    fill = "right",
    extra = "drop"
  ) %>%
  filter(!is.na(Species), is.na(Strain)) %>%
  select(Species, everything(), -Kingdom, -Phylum,
         -Class, -Order, -Family, -Genus, -Strain) %>%
  mutate(Species = sub("^s_+", "", Species)) %>%
  mutate(Species = make.unique(Species)) %>%
  column_to_rownames("Species") %>%
  as.matrix()

dim(species_data_all)
dim(taxa_data_all)
rm(taxa_data_all)

### Checking for overlap between species_data_all & metadata_filt ----
ggvenn::ggvenn(
  list(
    metadata_filt = rownames(metadata_filt),
    species_data_all   = colnames(species_data_all)
  )
)
ggsave("results/1c_samples_venn_diagram.png", width = 6, height = 4, dpi = 600)

#--#--#---#--#--#--#--#--#--#--#--#--#
keep_ids <- intersect(colnames(species_data_all) , rownames(metadata_filt))
length(keep_ids)                # 785 samples ->  147 samples (gid_wgs)

metadata_filt <- metadata_filt[keep_ids , ]
taxa_data <- species_data_all[ , keep_ids] # get first year subset of taxa data, retain full version for Q6
# if u want to check go and check the venn diagram

### Checking if both taxa_data and metadata_filt are identical ----
identical(rownames(metadata_filt) , colnames(taxa_data))    # TRUE
rm(keep_ids)

## Q1D Part 1 ----
###Filter low abundance taxa (threshold = 0.01% ----

taxa_data[ taxa_data < 0.01 ] <- 0  #logical indexing # [] used for subsetting

## Q1D Part 2 ----
### How many species remain? (335)----

keep_species <- which( rowSums(taxa_data) > 0 )
taxa_data <- taxa_data[ keep_species, ]    # all rows with net 0 removed

dim(taxa_data)                          # (335rows 147 cols) , 618 rows dropped to 335 (species)
rm(keep_species)

##Q1E ----
### Apply total sum scaling (TSS) normalization ----
col_totals <- colSums(taxa_data)
taxa_tss_data <- sweep(
  taxa_data[, col_totals > 0, drop = FALSE],
  2,
  col_totals[col_totals > 0], # avoid div0 error for samples with no species detected
  "/"
)

###############
# Qn 2  ----

## Q2A Milk and egg allergies across countries ----

ggplot(metadata_filt , aes(x = country , fill = allergy_milk)) +
  geom_bar(position = "fill") +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = scales::percent)
ggsave("results/2a_milk_allergy_by_country_v1.png", width = 6, height = 6, dpi = 600)

##better looking version of 2A. ----
ggplot(metadata_filt %>% tibble::as_tibble(),
       aes(x = country, fill = factor(allergy_milk))) +
  geom_bar(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_brewer(palette = "Set1", name = "Milk allergy") +
  labs(x = "Country", y = "Number of children") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("results/2a_milk_allergy_by_country_v2.png", width = 6, height = 6, dpi = 600)


##Q2B Covariates across countries  ----

metadata_covariates <- metadata_filt %>%
  select(country , delivery , Exclusive_breast_feeding , age_at_collection)

#age distribution
ggplot(metadata_covariates, aes(x = country, y = age_at_collection, fill = country)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
  labs(x = "country", y = "age at collection (days)") +
  theme(legend.position = "none")
ggsave("results/2b_age_at_collection_by_country.png", width = 6, height = 6, dpi = 600)

#delivery distribution
ggplot(metadata_covariates , aes(x = country  , fill = delivery)) +
  geom_bar( ) + labs( x = "country" , y = " type of delivery ")
ggsave("results/2b_delivery_method_by_country_v1.png", width = 6, height = 6, dpi = 600)

##better version for delivery distribution ----
ggplot(metadata_covariates %>% tibble::as_tibble(),
       aes(x = country, fill = delivery)) +
  geom_bar(position = position_dodge(width = 0.8),
           width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    position = position_dodge(width = 0.8),
    vjust = -0.35,
    size = 4
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Delivery") +
  labs(title = "Delivery mode by country",
       x = "Country", y = "Number of children") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
ggsave("results/2b_delivery_method_by_country_v2.png", width = 6, height = 6, dpi = 600)


#breast feeding distribution 
ggplot(metadata_covariates , aes(x = country  , fill = Exclusive_breast_feeding)) + 
  geom_bar( ) + labs( x = "country" , y = "Breastfed children") # this is the simpler version of the chart 
ggsave("results/2b_exclusive_breast_feeding_by_country_v1.png", width = 6, height = 6, dpi = 600)

##better version for breast feeding distribution ----
df_bf <- metadata_covariates %>%
  tibble::as_tibble() %>%
  mutate(Exclusive_breast_feeding = as.logical(Exclusive_breast_feeding)) %>%
  dplyr::count(country, Exclusive_breast_feeding, name = "n")
print(df_bf)

ggplot(df_bf, aes(x = country, y = n, color = Exclusive_breast_feeding)) +
  geom_linerange(aes(ymin = 0, ymax = n),
                 position = position_dodge(width = 0.5),
                 linewidth = 1) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.5),
            vjust = -0.6, size = 4) +
  labs(x = "Country", y = "Number of children", color = "Exclusive breastfeeding") +
  theme_minimal(base_size = 14)
ggsave("results/2b_exclusive_breast_feeding_by_country_v2.png", width = 6, height = 6, dpi = 600)


###############
# Qn 3  ----

## Q3a alpha diversity across countries ----
alpha_diversity_plot <- function(metadata,
                                 species,
                                 group_var,
                                 facet_ncol = 3,
                                 add_stats = TRUE) {
  # ---- calculate alpha diversity ----
  metadata_out <- metadata %>%
    mutate(
      richness = specnumber(species, MARGIN = 2),
      shannon_diversity = diversity(species, index = "shannon", MARGIN = 2),
      simpson_diversity = diversity(species, index = "simpson", MARGIN = 2)
    )

  # ---- reshape for plotting ----
  plot_df <- metadata_out %>%
    pivot_longer(
      cols = c(richness, shannon_diversity, simpson_diversity),
      names_to = "metric",
      values_to = "value"
    )

  # ---- base plot ----
  p <- ggplot(
    plot_df,
    aes(x = .data[[group_var]], y = value, fill = .data[[group_var]])
  ) +
    geom_boxplot() +
    facet_wrap(~ metric, scales = "free", ncol = facet_ncol) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_blank()
    )

  # ---- optional stats ----
  if (add_stats) {
    p <- p +
      stat_compare_means(
        aes(group = .data[[group_var]]),
        label = "p.format",
        hjust = -0.25
      )
  }

  # ---- return both data + plot ----
  return(
    list(
      metadata = metadata_out,
      plot = p
    )
  )
}


res <- alpha_diversity_plot(
  metadata = metadata_filt,
  species  = taxa_tss_data,
  group_var = "country"
)

# updated metadata with diversity columns
metadata_filt <- res$metadata

# plot
print(res$plot)
ggsave("results/3a_alpha_diversity_metrics_by_country.png", res$plot, width = 6, height = 6, dpi = 600)


## Q3c beta diversity across countries ----
beta_diversity_pcoa <- function(metadata,
                                species,
                                group_var,
                                method = c("bray", "jaccard"),
                                binary = TRUE,
                                k = 2,
                                ellipse = TRUE) {

  method <- match.arg(method)

  # ---- confirm sample name matches ----

  stopifnot(
    identical(rownames(metadata),
              colnames(species))
  )

  # ---- distance matrix ----
  dist_mat <- vegdist(
    t(species),
    method = method,
    binary = ifelse(method == "jaccard", binary, FALSE)
  )


  # ---- PCoA ----
  pcoa <- cmdscale(dist_mat, k = k, eig = TRUE)

  # ---- variance explained ----
  ve <- round(100 * pcoa$eig / sum(pcoa$eig), 1)
  axis_labels <- paste0(
    "PCoA", seq_len(k), " (", ve[seq_len(k)], "%)"
  )

  # ---- plotting dataframe ----
  plot_df <- data.frame(
    PCoA1 = pcoa$points[, 1],
    PCoA2 = pcoa$points[, 2],
    country = metadata[[group_var]],
    age = metadata$age_at_collection
  )

  # ---- PCoA plot ----
  p_country <- ggplot(plot_df,
                      aes(x = PCoA1,
                          y = PCoA2,
                          colour = country)) +
    geom_point(size = 3, alpha = 0.75) +
    labs(
      x = axis_labels[1],
      y = axis_labels[2],
      colour = group_var
    ) +
    theme_bw()

  # Ellipses should ONLY apply to country, ellipses on numeric age are meaningless statistically.

  if (ellipse) {
    p_country <- p_country + stat_ellipse(show.legend = FALSE)
  }

  p_age <- ggplot(plot_df,
                  aes(x = PCoA1,
                      y = PCoA2,
                      colour = age)) +
    geom_point(size = 3, alpha = 0.75) +
    scale_colour_gradient(name = "Age") +
    labs(
      x = axis_labels[1],
      y = axis_labels[2]
    ) +
    theme_bw()

  # ---- return everything ----
  return(
    list(
      distance = dist_mat,
      metadata_for_permanova = metadata,
      pcoa = pcoa,
      plot_data = plot_df,
      plot_country = p_country,
      plot_age = p_age
    )
  )
}


res_bray <- beta_diversity_pcoa(
  metadata = metadata_filt,
  species  = taxa_tss_data,
  group_var = "country",
  method = "bray"
)

print(res_bray$plot_country)
ggsave("results/3c_bray_pcoa_by_country.png", res_bray$plot_country, width = 6, height = 6, dpi = 600)

print(res_bray$plot_age)
ggsave("results/3c_bray_pcoa_by_age.png", res_bray$plot_age, width = 6, height = 6, dpi = 600)



res_jaccard <- beta_diversity_pcoa(
  metadata = metadata_filt,
  species  = taxa_tss_data,
  group_var = "country",
  method = "jaccard"
)

print(res_jaccard$plot_country)
ggsave("results/3c_jaccard_pcoa_by_country.png", res_jaccard$plot_country, width = 6, height = 6, dpi = 600)

print(res_jaccard$plot_age)
ggsave("results/3c_jaccard_pcoa_by_age.png", res_jaccard$plot_age, width = 6, height = 6, dpi = 600)


## Q3d PERMANOVA ----

dist_mat <- res_bray$distance
meta     <- res_bray$metadata_for_permanova

set.seed(123)

permanova <- adonis2(
  dist_mat ~
    age_at_collection +
    gender +
    delivery +
    Exclusive_breast_feeding +
    abx_first_year +
    country,
  data = meta,
  permutations = 999,
  by = "terms",
  na.action = na.omit
)
print(permanova)

sink("results/3d_bray_permanova_adonis2.txt") # save results
print(permanova)
sink(NULL)

rm(dist_mat, meta, res, res_bray, res_jaccard, permanova)


###############
# Qn 4  ----

## Q4a Top 15 species for each country ----

top_species_barplot <- function(metadata,
                                taxa_data,
                                group_var,
                                top_n = 15){

  # ---- sanity check sample matching ----
  stopifnot(
    identical(rownames(metadata),
              colnames(taxa_data))
  )

  groups <- levels(as.factor(metadata[[group_var]]))

  # ---- Compute top species per group ----
  top <- lapply(groups, function(g){

    idx <- which(metadata[[group_var]] == g)

    taxa_data[, idx] %>%
      rowMeans(na.rm = TRUE) %>%
      sort(decreasing = TRUE) %>%
      head(top_n) %>%
      enframe(name = "species",
              value = "mean_abundance") %>%
      mutate(group = g)

  }) %>%
    bind_rows() %>%
    group_by(group) %>%
    mutate(
      species = fct_reorder(species,
                            mean_abundance,
                            .desc = TRUE),
      rank_id = row_number()
    ) %>%
    ungroup()

  # ---- Plot ----
  p <- ggplot(top,
              aes(x = group,
                  y = mean_abundance,
                  fill = species,
                  label = rank_id)) +
    geom_bar(stat = "identity",
             colour = "black",
             width = 0.6) +
    geom_text(position = position_stack(vjust = 0.5),
              size = 3,
              colour = "black") +
    labs(
      x = NULL,
      y = "Relative abundance (mean)",
      fill = "Species"
    ) +
    theme_bw()

  return(
    list(
      data = top,
      plot = p
    )
  )
}

res_top <- top_species_barplot(
  metadata = metadata_filt,
  taxa_data = taxa_tss_data,
  group_var = "country",
  top_n = 15
)

print(res_top$plot)
ggsave("results/4a_top_15_species_by_country.png", res_top$plot, width = 8, height = 4, dpi = 600)

rm(res_top)

## Q4b Masslin2 ----

run_maaslin_country <- function(metadata,
                                taxa_data,
                                include_hla = FALSE,
                                ref_country = "RUS",
                                prevalence = 0.10,
                                abundance = 0.001,
                                q_cutoff = 0.01,
                                plot_scatter = TRUE,
                                res_folder = "results/maaslin_res"){
  
  # ---- sanity check ----
  stopifnot(
    identical(rownames(metadata),
              colnames(taxa_data))
  )

  # ---- set reference level ----
  metadata$country <- relevel(
    as.factor(metadata$country),
    ref = ref_country
  )

  # ---- fixed effects ----
  base_effects <- c("country",
                    "age_at_collection",
                    "gender",
                    "delivery",
                    "Exclusive_breast_feeding",
                    "abx_first_year")

  if(include_hla){
    fixed_effects <- c(base_effects, "hla_risk_class")
  } else {
    fixed_effects <- base_effects
  }

  # ---- create output folder ----
  if(!dir.exists(res_folder)){
    dir.create(res_folder, recursive = TRUE)
  }

  # ---- run Maaslin2 ----
  fit <- Maaslin2(
    input_data = taxa_data,
    input_metadata = metadata,
    output = res_folder,
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    min_prevalence = prevalence,
    min_abundance = abundance,
    fixed_effects = fixed_effects,
    correction = "BH",
    standardize = FALSE,
    reference = c(paste0("country,", ref_country)),
    max_significance = q_cutoff,
    plot_scatter = plot_scatter
  )

  # ---- extract country effects only ----
  results <- fit$results

  diff_taxa <- results %>%
    filter(metadata == "country") %>%
    filter(qval < q_cutoff)

  return(
    list(
      fit = fit,
      results = diff_taxa
    )
  )
}


res_maas <- run_maaslin_country(
  metadata = metadata_filt,
  taxa_data = taxa_tss_data,
  res_folder = "results/4b_maaslin_res",
  plot_scatter = FALSE # deprecated since ggplot2 v3.4.0
)

## Q4d visualization ----

plot_maaslin_diff_abundance <- function(maaslin_results,
                                        taxa_data,
                                        metadata,
                                        group_var = "country",
                                        q_cutoff = 0.01,
                                        pseudocount = 1e-6){
  
  # ---- sanity check ----
  stopifnot(
    identical(rownames(metadata),
              colnames(taxa_data))
  )

  # ---- filter sig taxa for group ----
  diff_taxa_group <- maaslin_results %>%
    filter(metadata == group_var) %>%
    filter(qval < q_cutoff)

  if(nrow(diff_taxa_group) == 0){
    stop("No significant taxa found for this group.")
  }

  sig_taxa <- unique(diff_taxa_group$feature)

  # ---- reshape taxa ----
  taxa_long <- taxa_data %>%
    as.data.frame() %>%
    rownames_to_column("taxa") %>%
    filter(taxa %in% sig_taxa) %>%
    pivot_longer(cols = -taxa,
                 names_to = "sample_id",
                 values_to = "abundance")

  # ---- reshape metadata ----
  metadata_needed <- metadata %>%
    rownames_to_column("sample_id") %>%
    select(sample_id, !!sym(group_var))

  # ---- merge ----
  taxa_long <- taxa_long %>%
    left_join(metadata_needed, by = "sample_id")

  # ---- raw abundance plot ----
  p_raw <- ggplot(taxa_long,
                  aes(x = get(group_var),
                             y = get("abundance"),
                             fill = get(group_var))) +
    geom_boxplot() +
    geom_jitter(width = 0.2,
                alpha = 1,
                size = 1) +
    facet_wrap(~taxa,
               scales = "free_y") +
    stat_compare_means(method = "anova",
                       label = "p.format",
                       vjust = 1) +
    theme(legend.position = "none") +
    labs(x = group_var,
         y = "Relative abundance")

  # ---- log plot ----
  p_log <- ggplot(taxa_long, aes(x = get(group_var), y = get("abundance"), fill = get(group_var))) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 1, size = 1) +
    facet_wrap(~taxa, scales = "free_y") +
    stat_compare_means(method = "anova", label = "p.format", vjust = 1) +
    scale_y_continuous(trans = "log10") +   # <-- log scale applied to axis
    theme(legend.position = "none") +
    labs(x = group_var, y = "Relative abundance (log10)")

  return(
    list(
      sig_table = diff_taxa_group,
      raw_plot = p_raw,
      log_plot = p_log
    )
  )
}


vis_maas <- plot_maaslin_diff_abundance(
  maaslin_results = res_maas$results,
  taxa_data = taxa_tss_data,
  metadata = metadata_filt,
  group_var = "country"
)

print(vis_maas$raw_plot)
ggsave("results/4d_maaslin_diff_abundance_raw.png", vis_maas$raw_plot, width = 6, height = 6, dpi = 600)

print(vis_maas$log_plot)
ggsave("results/4d_maaslin_diff_abundance_log.png", vis_maas$log_plot, width = 6, height = 6, dpi = 600)

print(vis_maas$sig_table)       # vector of significant taxa
write.csv(vis_maas$sig_table, "results/4d_maaslin_diff_abundance_significant_species.csv", row.names = FALSE, quote = FALSE)

rm(res_maas, vis_maas)

###############
# 5a ----------------------------------------------------------------------
#### 5. Functional Pathway Analysis
### a. Subset the pathway relative abundance profiles to correspond with your filtered samples (similar to what you did for species profiles)

# Load VatanenT_2016 metabolic pathway abundance
capstone_pathway_data <- curatedMetagenomicData(
  "VatanenT_2016.pathway_abundance",
  dryrun = FALSE,
  counts = FALSE
)

# Returns a list with one TSE object
capstone_pathway_tse <- capstone_pathway_data[[1]]

## Extract Pathway Matrix
# Get the assay (pathway matrix)
pathway_matrix <- assay(capstone_pathway_tse)
# Dimensions: features × samples
dim(pathway_matrix) # 24605 x 785
# Convert to data frame
pathway_matrix_df <- pathway_matrix %>% data.frame()

pathway <- pathway_matrix_df %>% rownames_to_column("feature") %>% filter(
  !grepl("\\.g__|\\.s__|\\||\\.unclassified", feature, ignore.case= TRUE) &
    !grepl("^UNINTEGRATED", feature, ignore.case= TRUE) &
    !grepl("^UNCLASSIFIED", feature, ignore.case= TRUE) &
    !grepl("^UNMAPPED", feature, ignore.case= TRUE)
) %>% column_to_rownames("feature")

pathway_filtered <- pathway
pathway_filtered[pathway_filtered < 0.00001] <- 0 # Removes pathways less than 0.001% abundant
pathway_filtered <- pathway_filtered[rowSums(pathway_filtered > 0) > 0, ] # Removes pathways where total row sum adds up to 0

nrow(pathway_filtered) # 452 rows

# Apply TSS
pathway_tss <- sweep(
  pathway_filtered,
  2,
  colSums(pathway_filtered),
  FUN = "/"
) * 100

summary(colSums(pathway_tss)) # should all be 100

metadata_filt_gidwgs <- metadata_filt %>%
  rownames_to_column("gid_wgs")

common_ids_pathway <- intersect(metadata_filt_gidwgs$gid_wgs, colnames(pathway_tss))
length(common_ids_pathway) # 147

# Taking common pathways
pathway_first_year <- pathway_tss %>%
  select(all_of(common_ids_pathway))

length(pathway_first_year) # should be 147

# Simplifying pathway names (e.g. KDO-NAGLIPASYN-PWY: superpathway of (Kdo)2-lipid A biosynthesis -> KDO-NAGLIPASYN-PWY)
# Also keep an original copy for future lookup table
ori_pathway_filt <- pathway_first_year
rownames(pathway_first_year) <- sub(":.*", "", rownames(pathway_first_year))

metadata_pathway <- metadata_filt_gidwgs[match(common_ids_pathway,
                                               metadata_filt_gidwgs$gid_wgs), ]

rm(capstone_pathway_data, capstone_pathway_tse, pathway, pathway_matrix, pathway_tss, pathway_filtered, common_ids_pathway, metadata_filt_gidwgs)

# 5b Task 3 ---------------------------------------------------------------
### b. Repeat the analyses from Tasks 3 and 4, addressing the same questions but from a functional pathway perspective
#### 3. Community Diversity and Structure (use 147 post-filtering)
### a. Calculate alpha diversity and visualize its distribution by country.
### Interpret any between-country differences.
### b. Bonus: can you statistically test differences in alpha diversity across countries?

## Richness
metadata_pathway$richnesspathway <- specnumber(pathway_first_year, MARGIN = 2)

ggplot(metadata_pathway, aes(x = country, y = richnesspathway, fill = country)) +
  geom_boxplot() +
  stat_compare_means(label = "p.format", hjust=0) +
  geom_jitter(width = 0.2, alpha = 1, size = 2) +
  labs(x = "Country", y = "Richness") +
  theme_bw() + theme(legend.position = "none")
ggsave("results/5b_3a_pathway_richness_by_country.png", width = 6, height = 6, dpi = 600)

## Simpson
metadata_pathway$simpson_diversity_pathway <- diversity(pathway_first_year, index = "simpson", MARGIN = 2)

ggplot(metadata_pathway, aes(x = country, y = simpson_diversity_pathway, fill = country)) +
  geom_boxplot() +
  stat_compare_means(label = "p.format", hjust=0) +
  geom_jitter(width = 0.2, alpha = 1, size = 2) +
  labs(x = "Country", y = "Simpson") +
  theme_bw() + theme(legend.position = "none")
ggsave("results/5b_3a_pathway_simpson_diversity_by_country.png", width = 6, height = 6, dpi = 600)

## Shannon
metadata_pathway$shannon_diversity_pathway <- diversity(pathway_first_year, index = "shannon", MARGIN = 2)

ggplot(metadata_pathway, aes(x = country, y = shannon_diversity_pathway, fill = country)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label = "p.format", hjust=0) +
  geom_jitter(width = 0.2, alpha = 1, size = 2) +
  labs(x = "Country", y = "Shannon") +
  theme_bw() + theme(legend.position = "none")
ggsave("results/5b_3a_pathway_shannon_diversity_by_country.png", width = 6, height = 6, dpi = 600)

## Multiple alpha-diversity
tmp_pathway <- metadata_pathway %>%
  pivot_longer(cols = c(richnesspathway, simpson_diversity_pathway, shannon_diversity_pathway)) %>%
  select(country, name, value)

head(tmp_pathway)

ggplot(tmp_pathway, aes(x = country, y = value, fill = country)) +
  facet_wrap( ~ name, scales = "free", ncol = 3) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label = "p.format", hjust = -0.5) +
  geom_jitter(width = 0.2, alpha = 1, size = 1) +
  theme_bw() + theme(legend.position = "none") + labs(x = NULL, y = NULL)
ggsave("results/5b_3a_pathway_ALL_diversity_by_country.png", width = 6, height = 4, dpi = 600)

rm(tmp_pathway)

### c. Compute beta-diversity using Bray-Curtis dissimilarity
### and perform Principal Coordinate Analysis.
### Visualize PCoA of samples colored by country and age.
### Hint: consider the approach in Figure 2A of the original manuscript.

## Bray-Curtis
bray_dist_pathway <- vegdist(t(pathway_first_year), method = "bray")
bray_pcoa_pathway <- cmdscale(bray_dist_pathway, k = 2, eig = TRUE)

mdf_pathway <- data.frame(
  PCoA1 = bray_pcoa_pathway$points[ , 1],
  PCoA2 = bray_pcoa_pathway$points[ , 2],
  country = metadata_pathway$country,
  age = metadata_pathway$age_at_collection) %>%
  mutate(age_bin = cut(age,
                       breaks = c(0, 90, 180, 270, 365),
                       labels = c("0–3 mo", "3–6 mo", "6–9 mo", "9-12 mo")))

VE_pathway <- round(100*bray_pcoa_pathway$eig / sum(bray_pcoa_pathway$eig), 1)
VE_pathway <- paste0("PCoA", 1:length(VE_pathway), " (", VE_pathway, "%)")

head(mdf_pathway, 3)

# Country
ggplot(mdf_pathway, aes(x = PCoA1, y = PCoA2, colour = country)) +
  geom_point(size = 3, alpha = 0.75) +
  stat_ellipse(show.legend = FALSE) +
  labs(x = VE_pathway[1], y = VE_pathway[2]) + theme_bw()
ggsave("results/5b_3c_pathway_bray_pcoa_by_country.png", width = 6, height = 6, dpi = 600)

# Age
ggplot(mdf_pathway, aes(x = PCoA1, y = PCoA2, colour = age)) +
  geom_point(size = 3, alpha = 0.75) +
  scale_colour_viridis_c(option = "plasma") +
  labs(x = VE_pathway[1], y = VE_pathway[2], colour = "Age (days)") +
  theme_bw()
ggsave("results/5b_3c_pathway_bray_pcoa_by_age.png", width = 6, height = 6, dpi = 600)

# If binning age into groups
ggplot(mdf_pathway, aes(PCoA1, PCoA2, colour = age_bin)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse() +
  labs(x = VE_pathway[1], y = VE_pathway[2], colour = "Age group") +
  theme_bw()
ggsave("results/5b_3c_pathway_bray_pcoa_by_age_groups.png", width = 6, height = 6, dpi = 600)

rm(mdf_pathway, VE_pathway, bray_dist_pathway, bray_pcoa_pathway)

### d. Statistically test between-country differences while
### adjusting for the following covariates:
### age at collection, gender, delivery mode,
### breastfeeding status, antibiotic usage in the first year.
### Interpret your results. Hint: PERMANOVA

## PERMANOVA
set.seed(123)

bray_dist_pathway <- vegdist(t(pathway_first_year), method = "bray")

fit_adonis_pathway <- adonis2(bray_dist_pathway ~ country + age_at_collection + gender + delivery + Exclusive_breast_feeding + abx_first_year,
                              data = metadata_pathway, by="terms", na.action = na.omit)
broom::tidy(fit_adonis_pathway) %>%
  knitr::kable()

sink("results/5b_3d_pathway_bray_permanova_adonis2.txt") # save results
print(fit_adonis_pathway)
sink(NULL)

rm(bray_dist_pathway, fit_adonis_pathway)

### Repeat of Task 4
### a. Visualize the abundance of the top 15 species (ranked by mean abundance)
### for each country using stacked bar plots. What visual differences do you observe between countries?

# Top 15 pathways by mean abundance
path_means <- rowMeans(pathway_first_year)

top15_path <- names(sort(path_means, decreasing = TRUE))[1:15]

plot_path <- pathway_first_year[top15_path, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("pathway") %>%
  pivot_longer(
    cols = -pathway,
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  left_join(metadata_pathway, by = c("sample_id" = "gid_wgs"))

# Plot
plot_path %>%
  group_by(country, pathway) %>%
  summarise(mean_abundance = mean(abundance), .groups = "drop") %>%
  ggplot(aes(x = country, y = mean_abundance, fill = pathway)) +
  geom_bar(stat = "identity", position = "stack", colour = "darkslategrey") +
  theme_bw() +
  ylab("Mean pathway relative abundance")
ggsave("results/5b_4a_top_15_pathways_by_country.png", width = 6, height = 6, dpi = 600)


rm(plot_path,top15_path, path_means)

# Comparing the barplot by country, Russia have less ileusyn PWY (L-isoleucine biosynthesis),
# but more PWY7222 (guanosine deoxyribonucleotides de novo biosynthesis II),
# PWY7220 (adenosine deoxyribonucleotides de novo biosynthesis II)
# and valsyn PWY (L-valine biosynthesis).
# These pathways contribited to amino acid metabolism and DNA synthesis.

### b. Identify species with statistically significant differential abundance between countries
## [Use Russia as the reference group; significance threshold = 0.01;
## focus on species with abundance >0% in at least 10% of samples]

# Differentially abundant pathways (MaAsLin2)
# Extract clean pathway IDs
path_ids <- sub("^([^:|]+).*", "\\1", rownames(ori_pathway_filt))

# Create MaAsLin-style sanitized feature names
mas_features <- path_ids %>%
  gsub("[^A-Za-z0-9]", ".", .) %>%   # replace non-alphanumeric with .
  gsub("\\.+", ".", .)               # collapse multiple dots

# Build lookup table
pathway_map <- data.frame(
  pathway_id   = path_ids,           # PWY-7111
  mas_feature  = mas_features,        # PWY.7111
  original_name = rownames(ori_pathway_filt),
  stringsAsFactors = FALSE
)

# Sanity check
identical(rownames(pathway_first_year), pathway_map$pathway_id)

# Create modeling metadata (subset only needed variables)
sub_metadata <- metadata_pathway %>%
  dplyr::select(
    gid_wgs,
    country,
    age_at_collection,
    gender,
    delivery,
    Exclusive_breast_feeding,
    abx_first_year
  )

# Set rownames using gid_wgs for MaAsLin2
rownames(sub_metadata) <- sub_metadata$gid_wgs
sub_metadata$gid_wgs <- NULL

# Align ordering with pathway matrix
sub_metadata <- sub_metadata[colnames(pathway_first_year), ]
all(colnames(pathway_first_year) == rownames(sub_metadata))

# Run Maaslin2
diff_path <- Maaslin2(
  input_data = pathway_first_year,
  input_metadata = sub_metadata,
  output = "results/5b_4b_maaslin_pathways",
  fixed_effects = c(
    "country",
    "age_at_collection",
    "gender",
    "delivery",
    "Exclusive_breast_feeding",
    "abx_first_year"
  ),
  reference = "country,RUS",
  normalization = "TSS",
  transform = "LOG",
  correction = "BH",
  standardize = FALSE,
  min_prevalence = 0.1, min_abundance = 0,
  max_significance = 0.01,
  plot_scatter = FALSE
)

rm(mas_features, path_ids)

### c. Compare with original publication-> moved after d.


### d. Visualize differentially abundant pathways across countries
maaslin_results <- read.delim(
  "results/5b_4b_maaslin_pathways/all_results.tsv"
) %>%
  filter(metadata == "country", qval < 0.01) %>%
  left_join(pathway_map, by = c("feature" = "mas_feature")) %>%
  filter(!is.na(pathway_id)) %>%
  distinct(pathway_id, .keep_all = TRUE) %>%   # one row per pathway
  arrange(qval)

# Pick top 8 paths
top8_paths <- maaslin_results %>%
  slice_head(n = 8) %>%
  pull(pathway_id) %>%
  intersect(rownames(pathway_first_year))

# Subset filter_path using pathway_id
plot_sig_path <- pathway_first_year[top8_paths, , drop = FALSE] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(
    cols = -sample_id,
    names_to = "pathway_id",
    values_to = "abundance"
  ) %>%
  left_join(metadata_pathway, by = c("sample_id" = "gid_wgs")) %>%
  left_join(pathway_map, by = "pathway_id") %>%
  mutate(
    pathway_label = ifelse(
      is.na(original_name),
      pathway_id,
      original_name
    )
  ) %>%
  filter(!is.na(country))

# Plot
ggplot(plot_sig_path,
       aes(x = country, y = abundance, fill = country)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.1, height = 0.0004, size = 0.1, alpha = 0.15) +
  facet_wrap(~ pathway_label, scales = "free_y", nrow = 2) +
  scale_y_continuous(labels = scales::label_number(suffix = "%")) +
  theme_bw() +
  theme(legend.position = "none", strip.text = element_text(hjust = 0)) +
  ylab("Pathway abundance (%)")
ggsave("results/5b_4d_maaslin_diff_abundance_raw.png", width = 8, height = 4, dpi = 600)

# Print q val table
top8_table <- maaslin_results %>%
  select(original_name, qval, coef) %>%
  slice_head(n = 8)
print(top8_table)

rm(sub_metadata, diff_path, plot_sig_path, top8_table, top8_paths)

### c. Compare with original publication

# Top 8 differential pathways by country (RUS as ref) suggestes:
# Lower folate biosynthesis (PWY-7539, PWY-6147)
# Lower preQ0 biosynthesis (PWY-6703) -> lower tRNA modification
# Lower flavin biosynthesis (PWY-6168) -> lower B2 production
# Lower CDP-diacylglycerol biosynthesis (PWY-5667, PWY0-1319) -> lower membrane lipid precursor

# In which, membrane lipid precursor may be relevant with the LPS and lipid A biothesis mentioned in Fig3A.

# Also higher sucrose degradation IV (PWY-5384) -> higher carbohydrate breakdown
# Higher UDP-N-acetyl-D-glucosamine biosynthesis (UDPNAGSYN-PWY) -> higher peptidoglycan (lipid A) precursor

# In which, glycolysis process and lipid A biosynthetic process are also mentioned in the paper Fig3A.

### Bonus: In the original publication, the authors analyzed metagenomic data-derived GO terms
### and noted that "lipid A biosynthetic process (GO: 0009245) showed a striking difference in abundance
### between countries." Investigate whether you can detect similar differences using metaCyc pathway data

# Select lipid A–related pathways
lipidA_paths <- pathway_map %>%
  filter(str_detect(original_name, regex("lipid|Lipid", ignore_case = TRUE))) %>%
  pull(pathway_id) %>%
  intersect(rownames(pathway_first_year))

# Intersect with maaslin_results
lipidA_paths <- intersect(maaslin_results$pathway_id, lipidA_paths)

# Make plotting data frame.
lipidA_df <- pathway_first_year[lipidA_paths, , drop = FALSE] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(
    cols = -sample_id,
    names_to = "pathway_id",
    values_to = "abundance"
  ) %>%
  left_join(metadata_pathway, by = c("sample_id" = "gid_wgs")) %>%
  filter(!is.na(country))

# Plot with stats
ggplot(lipidA_df,
       aes(x = country, y = abundance * 100, fill = country)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_quasirandom(size = 0.5, alpha = 0.35) +
  facet_wrap(~ pathway_id, scales = "free_y") +
  scale_y_continuous(labels = scales::label_number(suffix = "%")) +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  ylab("Pathway abundance (%)")
ggsave("results/5c_lipid_pathway_abundance_by_country.png", width = 6, height = 6, dpi = 600)

rm(lipidA_df)

# NAGLIPASYN-PWY is significantly different between countries (lower in Russia)
# This pathway is relevent with cell wall peptidoglycan synthesis, a precursor
# For Lipid A biosynthetic process. That is similar to the authors' conclusion.

### d. Bonus: If yes for 5c, please identify the top species or
### genera contributing to the lipid A biosynthesis pathway.

# Focus on the NAGLIPASYN-PWY
lipidA_taxa_df <- pathway_matrix_df %>%
  rownames_to_column("pathway_taxa") %>%
  filter(grepl(paste0("^", lipidA_paths, ":"), pathway_taxa) |
           grepl(paste0("^", lipidA_paths, "\\|"), pathway_taxa))

lipidA_taxa_df <- lipidA_taxa_df %>%
  filter(grepl("\\|", pathway_taxa))

lipidA_long <- lipidA_taxa_df %>%
  pivot_longer(
    cols = -pathway_taxa,
    names_to = "sample_id",
    values_to = "abundance"
  )

# Extract genus and species from the row names
lipidA_long <- lipidA_long %>%
  mutate(
    pathway_genus = ifelse(grepl("\\|", pathway_taxa),
                           gsub("^.*\\|g__([^.]+)\\..*$", "\\1", pathway_taxa),
                           NA),
    pathway_species = ifelse(grepl("\\|", pathway_taxa),
                             gsub("^.*\\.s__([^.]+).*$", "\\1", pathway_taxa),
                             NA)
  )

# 4. Add metadata (country)
lipidA_long <- lipidA_long %>%
  left_join(metadata_pathway %>% select(gid_wgs, country),
            by = c("sample_id" = "gid_wgs")) %>%
  filter(!is.na(country))

# Aggregate mean abundance per genus/species per country
taxa_summary <- lipidA_long %>%
  # remove rows where genus or species is NA or unclassified
  filter(
    !is.na(pathway_genus), !is.na(pathway_species),
    !grepl("unclassified", pathway_genus, ignore.case = TRUE),
    !grepl("unclassified", pathway_species, ignore.case = TRUE)
  ) %>%
  group_by(country, pathway_species) %>%
  summarise(total_abundance = sum(abundance), .groups = "drop")

# Order by abundance
taxa_summary$pathway_species <- factor(
  taxa_summary$pathway_species,
  levels = taxa_summary %>%
    group_by(pathway_species) %>%
    summarise(total = sum(total_abundance)) %>%
    arrange(desc(total)) %>%
    pull(pathway_species)
)

# Plot stacked bar by country
ggplot(taxa_summary, aes(x = country, y = total_abundance, fill = pathway_species)) +
  geom_bar(stat = "identity", position = "stack", colour = "darkslategrey") +
  labs(
    x = "Country",
    y = "",
    title = "Total contribution to NAGLIPASYN-PWY",
    fill = "Species"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("results/5d_NAGLIPASYN_species_by_country.png", width = 12, height = 6, dpi = 600)

rm(lipidA_taxa_df, lipidA_long, taxa_summary,maaslin_results, lipidA_paths)


#### Q6 Longitudinal data ----
## Q6a. Return to the original metadata and taxonomy files. Some subjects were sampled at multiple timepoints. Summarize the number of samples collected from each subject. *Hint: group_by subjectID* ----
# keep only samples with corresponding row in metadata and column in species table
valid_samples <- intersect(colnames(species_data_all), metadata_all$gid_wgs)
longitudinal_species <- species_data_all %>% as.data.frame() %>% select(all_of(valid_samples))
dim(longitudinal_species)
longitudinal_meta <- metadata_all %>% filter(!is.na(gid_wgs), gid_wgs %in% valid_samples)
dim(longitudinal_meta)
rm(valid_samples)

# clean all species-level longitudinal abundance data
longitudinal_species <- longitudinal_species %>%
  mutate(across(everything(), ~ ifelse(. < 0.01, 0, .))) %>% # remove low abundance "noise"
  filter(rowSums(.) != 0) # remove species with abundance = 0 for all samples
longitudinal_species <- sweep(longitudinal_species, 2, colSums(longitudinal_species), FUN = "/") # TSS

# plot number of samples per subject
ggplot(
  longitudinal_meta %>% group_by(subjectID) %>%
    summarise(num_collect = n()),
  aes(x = num_collect)
) +
  geom_histogram(binwidth = 1, boundary = 0) +
  labs(title = "Number of samples collected per subject", x = "", y = "") +
  theme_bw() + theme(legend.position = "none")
ggsave("results/6a_num_samples_per_subject.png", width = 6, height = 4, dpi = 600)

# plot number of samples per subject, for each country
ggplot(longitudinal_meta %>% group_by(country, subjectID) %>% summarise(num_collect = n(), .groups = "keep"), aes(x = num_collect, fill = country)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0) +
  labs(x = "# samples collected per subject", y = "") +
  facet_wrap(~ country, ncol = 1) +
  theme_bw() + theme(legend.position = "none")
ggsave("results/6a_num_samples_per_subject_by_country.png", width = 6, height = 4, dpi = 600)


## Q6b. Select the subject with the most samples and visualize community composition across all their samples, ordered by age at collection ----
# count samples per subject, excluding invalid sample IDs with no taxa data
sample_count <- longitudinal_meta %>%
  group_by(subjectID) %>% summarise(num_samples = n()) %>%
  arrange(desc(num_samples))
print(head(sample_count))

# most samples in relative abundance table = P014839 (13)
subj_max_samples <- (sample_count %>% slice_max(num_samples))$subjectID
max_samples_meta <- longitudinal_meta %>% filter(subjectID == subj_max_samples, !is.na(gid_wgs)) %>% arrange(collection_month)
print(max_samples_meta)

max_samples_species <- longitudinal_species %>% select(all_of(max_samples_meta$gid_wgs))

# plot top 5 species by relative abundance at each time
max_samples_plot_df <- merge(
  max_samples_meta %>%
    mutate(collection_month = factor(collection_month)) %>%
    select(gid_wgs, collection_month),
  t(max_samples_species),
  by.x = "gid_wgs", by.y = 0
) %>%
  pivot_longer(cols = -(1:2), names_to = "species", values_to = "rel_ab")
ggplot(
  max_samples_plot_df %>% group_by(collection_month) %>% slice_max(rel_ab, n = 5),
  aes(x = collection_month, y = rel_ab, fill = species)
) +
  geom_bar(colour = "darkslategrey", stat = "identity", position = "stack") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  labs(y = "Relative abundance", x = "Month of collection", fill = "Species") +
  theme(legend.position = "right")
ggsave(paste0("results/6b_longitudinal_top_5_species_", subj_max_samples, ".png"), width = 8, height = 6, dpi = 600)


# plot top 5 genus by relative abundance at each time
print(
  ggplot(
    max_samples_plot_df %>% mutate(
      genus = sapply(max_samples_plot_df$species, function(s) str_split_1(s, "_")[[1]])
    ) %>% group_by(gid_wgs, collection_month, genus) %>% summarise(total_relab = sum(rel_ab), .groups = "keep") %>%
      group_by(collection_month) %>% slice_max(total_relab, n = 5),
    aes(x = collection_month, y = total_relab, fill = genus)
  ) +
    geom_bar(colour = "darkslategrey", stat = "identity", position = "stack") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    labs(y = "Relative abundance", x = "Month of collection", fill = "Genus") +
    theme(legend.position = "right")
)
ggsave(paste0("results/6b_longitudinal_top_5_genus_", subj_max_samples, ".png"), width = 8, height = 6, dpi = 600)


## Q6c. Using all collected samples, identify differentially abundant species between countries while accounting for repeated sampling. *Hint: consider subject ID as a random effect* ----
longitudinal_species_maaslin_fit <- Maaslin2(
  input_data = longitudinal_species,
  input_metadata = longitudinal_meta %>% column_to_rownames("gid_wgs"),
  output = "results/6c_longitudinal_maaslin_res",
  normalization ="TSS", transform = "log",
  analysis_method = "LM",
  min_prevalence = 0.1,
  min_abundance = 0.001,
  fixed_effects = c("country", "age_at_collection", "gender", "delivery", "Exclusive_breast_feeding", "abx_first_year"),
  random_effects = c("subjectID"),
  reference = "country,RUS",
  correction = "BH",
  standardize = FALSE,
  max_significance = 0.01,
  plot_scatter = FALSE # deprecated
)

# subset only country-related differences at 0.01 significance threshold
diff_taxa_longit <- longitudinal_species_maaslin_fit$results %>% filter(qval < 0.01, metadata == "country")
print(diff_taxa_longit %>% arrange(qval))

# plot only 8 most significant features (ranked based on min. qval)
plot_species <- diff_taxa_longit %>% group_by(feature) %>%
  summarise(min_qval = min(qval)) %>%
  slice_min(min_qval, n = 8)
plot_df <- longitudinal_species %>%
  rownames_to_column("feature") %>%
  filter(feature %in% plot_species$feature) %>%
  pivot_longer(cols = -feature, names_to = "gid_wgs", values_to = "relab") %>%
  left_join(
    longitudinal_meta %>% select(subjectID, gid_wgs, country),
    by = "gid_wgs"
  )

ggplot(plot_df, aes(
  y = relab,
  x = country,
  fill = country
)) +
  geom_violin(scale = "count", quantile.linetype = 1) +
  facet_wrap(~ feature, scales="free_y", nrow = 2) +
  labs(x= "", y= "Relative abundance") +
  theme_bw() + theme(legend.position = "none")
ggsave("results/6c_longitudinal_maaslin_diff_abundance_raw.png", width = 12, height = 6, dpi = 600)

ggplot(plot_df, aes(
  y=(relab + 1e-50), # add pseudo-count for plotting 0 values on log scale
  x=country,
  fill=country
)) +
  geom_violin(scale = "count", quantile.linetype = 1) +
  facet_wrap(~ feature, scales="free_y", nrow = 2) +
  scale_y_log10() +
  labs(x= "", y= "Relative abundance (log10)") +
  theme_bw() + theme(legend.position = "none")
ggsave("results/6c_longitudinal_maaslin_diff_abundance_log.png", width = 12, height = 6, dpi = 600)
