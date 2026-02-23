# Setting up and loading data ----
setwd("C:/previous desktop folders/Desktop/CAPSTONE PROJECT")
pacman::p_load(tidyverse, janitor, ggpubr, vegan, Maaslin2)
rm(list = ls())

load("Data/DIABIMMUNE_Karelia_metadata_sub.RData", verbose = TRUE)
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
taxa_data <- assay(cdata)
rm(cdata)

dim(taxa_data)                  # 618 rows 785 cols
taxa_data[1:5 , 1:5]            # rownames = taxanomy  , colnames = gid_wgs ID

### Checking for overlap between taxa_data & metadata_filt ----
ggvenn::ggvenn(
  list(
    metadata_filt = rownames(metadata_filt),
    taxa_data   = colnames(taxa_data)
  )
)
#--#--#---#--#--#--#--#--#--#--#--#--#
keep_ids <- intersect(colnames(taxa_data) , rownames(metadata_filt))
length(keep_ids)                # 785 samples ->  147 samples (gid_wgs) 

metadata_filt <- metadata_filt[keep_ids , ]
taxa_data <- taxa_data[ , keep_ids]
# if u want to check go and check the venn diagram 

### Changing the taxa name to species level
ranks <- c("Kingdom", "Phylum", "Class", "Order",
           "Family", "Genus", "Species", "Strain")

taxa_data_byspecies <- taxa_data %>%
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

dim(taxa_data_byspecies)
dim(taxa_data)
rm(taxa_data)

### Checking if both taxa_data and metadata_filt are identical ----
identical(rownames(metadata_filt) , colnames(taxa_data_byspecies))    # TRUE 
rm(keep_ids)

## Q1D Part 1 ----
###Filter low abundance taxa (threshold = 0.01% ----

taxa_data_byspecies[ taxa_data_byspecies < 0.01 ] <- 0  #logical indexing # [] usede for subsetting 

## Q1D Part 2 ----
### How many species remain? (335)----

keep_species <- which( rowSums(taxa_data_byspecies) > 0 )
taxa_data_byspecies <- taxa_data_byspecies[ keep_species, ]    # all rows with net 0 removed 

dim(taxa_data_byspecies)                      # (335rows 147 cols) , 618 rows dropped to 335 (species)  
rm(keep_species) 

taxa_data <- taxa_data_byspecies              #long name shortened back because dont really need to specify after this qn 
rm(taxa_data_byspecies)

##Q1E ----
### Apply total sum scaling (TSS) normalization ----
col_totals <- colSums(taxa_data)
taxa_tss_data <- sweep(
  taxa_data[, col_totals > 0, drop = FALSE],
  2,
  col_totals[col_totals > 0],
  "/"
)

# Qn 2  ----

## Q2A Milk and egg allergies across countries ----

ggplot(metadata_filt , aes(x = country , fill = allergy_milk)) +
  geom_bar(position = "fill") +
  labs(x = NULL, y = NULL) + 
  scale_y_continuous(labels = scales::percent)

##better looking version of 2A. ----
library(tidyverse)
library(scales)

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


##Q2B Covariates across countries  ----

metadata_covariates <- metadata_filt %>% 
  select(country , delivery , Exclusive_breast_feeding , age_at_collection)

#age distribution 
ggplot(metadata_covariates, aes(x = country, y = age_at_collection, fill = country)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
  labs(x = "country", y = "age at collection (days)") +
  theme(legend.position = "none")

#delivery distribution 
ggplot(metadata_covariates , aes(x = country  , fill = delivery)) + 
  geom_bar( ) + labs( x = "country" , y = " type of delivery ")
##better version for delivery distribution ----
library(tidyverse)

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


#breast feeding distribution 
ggplot(metadata_covariates , aes(x = country  , fill = Exclusive_breast_feeding)) + 
  geom_bar( ) + labs( x = "country" , y = "Breastfed children") # this is the simpler version of the chart 

##better version for breast feeding distribution ----
library(tidyverse)
df_bf <- metadata_covariates %>%
  tibble::as_tibble() %>%
  mutate(Exclusive_breast_feeding = as.logical(Exclusive_breast_feeding)) %>%
  dplyr::count(country, Exclusive_breast_feeding, name = "n")

df_bf

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




# Qn 3  ----

## Q3a alpha diversity across countries ----
alpha_diversity_plot <- function(metadata,
                                 species,
                                 group_var,
                                 facet_ncol = 3,
                                 add_stats = TRUE) {
  # packages
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(vegan)
  require(rstatix)
  
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
res$plot


## Q3c beta diversity across countries ----
beta_diversity_pcoa <- function(metadata,
                                species,
                                group_var,
                                method = c("bray", "jaccard"),
                                binary = TRUE,
                                k = 2,
                                ellipse = TRUE) {
  # packages
  require(vegan)
  require(ggplot2)
  require(dplyr)
  
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

res_bray$plot_country
res_bray$plot_age



res_jaccard <- beta_diversity_pcoa(
  metadata = metadata_filt,
  species  = taxa_tss_data,
  group_var = "country",
  method = "jaccard"
)

res_jaccard$plot_country
res_jaccard$plot_age


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

rm(dist_mat, meta, res, res_bray, res_jaccard, permanova)

# Qn 4  ----


## Q4a Top 15 species for each country ----

top_species_barplot <- function(metadata,
                                taxa_data,
                                group_var,
                                top_n = 15){
  
  require(dplyr)
  require(ggplot2)
  require(tibble)
  require(forcats)
  
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

res_top$plot

rm(res_top)

## Q4b Masslin2 ----

run_maaslin_country <- function(metadata,
                                taxa_data,
                                include_hla = FALSE,
                                ref_country = "RUS",
                                prevalence = 0.10,
                                abundance = 0.001,
                                q_cutoff = 0.01){
  
  require(Maaslin2)
  require(dplyr)
  
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
  if(!dir.exists("results")){
    dir.create("results")
  }
  
  if(!dir.exists("results/maaslin_res")){
    dir.create("results/maaslin_res")
  }
  
  # ---- run Maaslin2 ----
  fit <- Maaslin2(
    input_data = taxa_data,
    input_metadata = metadata,
    output = "results/maaslin_res",
    normalization = "TSS",
    transform = "log",
    analysis_method = "LM",
    min_prevalence = prevalence,
    min_abundance = abundance,
    fixed_effects = fixed_effects,
    correction = "BH",
    standardize = FALSE,
    reference = c(paste0("country,", ref_country)),
    max_significance = q_cutoff
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
  taxa_data = taxa_tss_data
)

## Q4d visualization ----

plot_maaslin_diff_abundance <- function(maaslin_results,
                                        taxa_data,
                                        metadata,
                                        group_var = "country",
                                        q_cutoff = 0.01,
                                        pseudocount = 1e-6){
  
  require(dplyr)
  require(tidyr)
  require(tibble)
  require(ggplot2)
  require(ggpubr)
  
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
                  aes_string(x = group_var,
                             y = "abundance",
                             fill = group_var)) +
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
  p_log <- ggplot(taxa_long, aes_string(x = group_var, y = "abundance", fill = group_var)) +
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

vis_maas$raw_plot
vis_maas$log_plot
vis_maas$sig_table       # vector of significant taxa

rm(res_maas, vis_maas)
