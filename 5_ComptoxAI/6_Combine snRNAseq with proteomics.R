# Sandbox analysis
source(here::here("!libraries.r"))
source(here::here("!directories.r"))
source(here::here("!load_clean_data.R"))
library(igraph)

# 1. Read in data ----


## a. snRNAseq data -------
# Read and combine all CSV files into a single data frame, adding a column for file names
# Read and combine all CSV files into a single data frame, adding a column for file names
scRNAseq_deg <- list.files(path = fs::path(dir_data, "scRNAseq DEGs"),
                           pattern = "*.csv", 
                           full.names = TRUE) %>%
  set_names(.) %>%
  map_df(~ read_csv(.x) %>% 
           mutate(file_name = tools::file_path_sans_ext(basename(.x))))

# Rename first col
scRNAseq_deg <- scRNAseq_deg |> 
  dplyr::rename("gene" = 1) |> 
  dplyr::mutate(
    file_name = str_remove(file_name, ".de.markers"),
    cell_type_time = str_remove(file_name, "PFNAvsCTR_"),
    week = if_else(str_detect(file_name, "week1"), "week 1", "week 2"),
    
    cell_type = str_remove(cell_type_time, "week1_") |> 
      str_remove("week2_"))

length(unique(scRNAseq_deg$gene))

table(scRNAseq_deg$pct.2>scRNAseq_deg$pct.1)
# # Identify top DEGs
scRNAseq_deg <- scRNAseq_deg |>
  dplyr::mutate(
    top_scrnaseq = ifelse(abs(avg_log2FC) > quantile(abs(avg_log2FC), 0.5),
                          1, 0))

scRNAseq_deg_sig <- scRNAseq_deg |>
  tidylog::filter(p_val<0.05)
  # tidylog::filter(p_val_adj<0.05)

length(unique(scRNAseq_deg_sig$gene))

# 2. ComptoxAI graph generated in python ----
# 35 protein graph
g <- read_graph(fs::path(dir_results,
                         "ComptoxAI",
                         "PFAS_prot_in_vitro_sig_091224.graphml"),
                format = "graphml")

# g <- read_graph(fs::path(dir_results,
#                          "ComptoxAI",
#                          "PFAS_prot_dkd_all_prot_032724.graphml"),
#                 format = "graphml")

graph_with_clusters <- read_csv(fs::path(dir_results, "network_graph_analysis.csv"))


# Calculate graph data
node_metadata <- tibble(
  geneSymbol = V(g)$id,
  eigen_centrality = igraph::eigen_centrality(g)$vector,
  authority = igraph::authority_score(g)$vector,
  harmonic_centrality = harmonic_centrality(g, mode = "out"),
  hub_score = hub_score(g)$vector,
  page_rank = page_rank(g)$vector,
  strength = strength(g),
  # power_centrality = power_centrality(g),
  # subgraph_centrality = subgraph_centrality(g),
  # alpha_centrality = alpha_centrality(g),
  betweenness = igraph::betweenness(g),
  degree_all   = igraph::degree(g, mode = "all"),
  degree_in    = igraph::degree(g, mode = "in"),
  degree_out   = igraph::degree(g, mode = "out"),
  degree_total = igraph::degree(g, mode = "total"),
  eccentricity = eccentricity(g),
  coreness = coreness(g),
  closeness = closeness(g)
)

mean(node_metadata$degree_in)

# ## a. Calculate top comptoxai hubs -------
node_metadata <- node_metadata |>
  mutate(top_comptox = ifelse(degree_out > quantile(degree_out, 0.5), 1, 0))

# 3. Proteomics meet in middle results ----
# Meet in middle results:
meet_in_middle_effects <- read.csv(fs::path(dir_results, "med_res_df_020924.csv"))

# Results of all mediation analyses
pfas_proteomics_res <- read.csv(
  fs::path(dir_results,
           "pairwise_mediation_result.csv"))

# Combine protein metadata and proteomics results
pfas_proteomics_all <- prot_metadata |> 
  # tidylog::select(-c(category:ratio_salivary_gland)) |>
  tidylog::select(-c(SeqId:EntrezGeneID),
                  -c(Organism:maxrownum)) |>
  tidylog::left_join(pfas_proteomics_res) 

# For now, since we didnt include an interaction, we can exclude the .T and 
# .C vars from the mediation results
pfas_proteomics_all <- pfas_proteomics_all |>
  tidylog::select(-contains(".T"), 
                  -contains(".C"), 
                  # -contains(".pval"),
                  -contains(".hi"), 
                  -contains(".lo")) |>
  rename_all(~str_remove(., ".avg"))

## a. Filter to select only 35 sig from meet in middle -------
pfas_proteomics_35 <- pfas_proteomics_all |>
  tidylog::filter(AptName %in% meet_in_middle_effects$feature_name)



# Clean proteomics gene names to merge better
pfas_proteomics_35 <- pfas_proteomics_35 %>%
  dplyr::mutate(
    EntrezGeneSymbol = case_when(EntrezGeneSymbol == "FGF8_3" ~ "FGF8", 
                                 EntrezGeneSymbol == "HSPA1A_2" ~ "HSPA1A",
                                 EntrezGeneSymbol == "CLEC4G_2" ~ "CLEC4G",
                                 EntrezGeneSymbol == "HEY1_2" ~ "HEY1",
                                 EntrezGeneSymbol == "GPNMB_4" ~ "GPNMB",
                                 EntrezGeneSymbol == "RSPO2_1" ~ "RSPO2",
                                 EntrezGeneSymbol == "IL27|EBI3" ~ "IL27",
                                 EntrezGeneSymbol == "GDF11|MSTN" ~ "GDF11", 
                                 TRUE ~ EntrezGeneSymbol)) |>
  # Remove the "_#" pattern from EntrezGeneSymbol
  mutate(EntrezGeneSymbol = gsub("_\\d+", "", EntrezGeneSymbol)) %>%
  group_by(EntrezGeneSymbol) |> 
  tidylog::filter(abs(PMed) == max(abs(PMed))) |> 
  slice_head() |> 
  ungroup() 

# ## b. Calculate top 10 -----
pfas_proteomics_35 <- pfas_proteomics_35 |>
  mutate(top_prop_med = ifelse(PMed > quantile(PMed, 0.5), 1, 0))

# 4. Combine data from proteomics, comptoxai, and scRNAseq ----
# proteomics and comptoxai
ppw <-  tidylog::right_join(
  node_metadata,
  pfas_proteomics_35,
  by = c("geneSymbol" = "EntrezGeneSymbol"))

# Combine with scRNAseq
ppw2 <- tidylog::left_join(ppw,
  scRNAseq_deg_sig,
  # by = c("EntrezGeneSymbol" = "gene"))
  by = c("geneSymbol" = "gene"))

length(unique(ppw2$geneSymbol))


# Examine the results
# temp <- ppw2 |> 
#   dplyr::select(geneSymbol, cell_type, week, avg_log2FC, 
#                 pct.1, pct.2, PMed, PMed.pval, ACME, ACME_rank) |>
#   dplyr::mutate(same_dir = (avg_log2FC/ACME)/abs(avg_log2FC/ACME))

# 5. Analyze combined data -----------
length(unique(ppw2$geneSymbol))
table(ppw2$week, ppw2$cell_type)


# collapse data across scRNseq cell types and time by selecting the most 
# significant p_value from the scRNAseq data
ppw_summarized <- ppw2 |> 
  group_by(geneSymbol) |> 
  tidylog::filter(p_val == min(p_val)) |> 
  ungroup()
  # dplyr::summarise(
  #   across(where(is.numeric), mean, na.rm = TRUE),
  #   across(where(is.character), function(x){str_c(unique(x), collapse = "; ")}))
# rm(ppw2)

# Determine top 10% of mediation, comptoxai, and scRNAseq effects:
ppw_summarized <- ppw_summarized |>
  mutate(
    top_comptox = ifelse(degree_out > quantile(degree_out, 0.66), 1, 0),
    top_prop_med = ifelse(PMed > quantile(PMed, 0.66), 1, 0),
    top_scrnaseq = ifelse(abs(avg_log2FC) > quantile(abs(avg_log2FC), 0.66), 1, 0))

# Filter top features across all approaches
ppw_top_pct <- ppw_summarized |> 
  tidylog::filter(top_comptox == 1, 
                  top_prop_med == 1, 
                  top_scrnaseq == 1)


# Remove RSPO2, CHRNA5, FSCN1 because effect estimates are not possible?
# pfas_prot_network <- pfas_prot_network |>
#   tidylog::filter(!(geneSymbol %in% c ("CHRNA5", "FSCN1")))

# pivot data wider on mediation effect estimates
# ppw <- ppw |> 
#   pivot_wider(names_from = Effect, 
#               values_from = pe, 
#               id_cols = c(geneSymbol:closeness)) 


## a. Filter to only positive percent mediated ---------------------------------
ppw_pos <- ppw_summarized |> 
  tidylog::filter(PMed > 0) |>
  mutate(ACME_abs = abs(ACME))

# Rank the data values
ppw_rank <- ppw_pos |> 
  tidylog::select(-geneSymbol) |> #, -c(AptName:ensembl_gene_id)
  tidylog::mutate(
    across(c(eigen_centrality:closeness, TE:PMed.pval), 
      ~rank(.))) %>%
  bind_cols(ppw_pos[,1], .) 


## b. correlation plots ----
# # (cor_matrix <- cor(janitor::remove_constant(ppw_rank[,-1]), 
# #                    method = "spearman"))
# # corrplot::corrplot(cor_matrix, method = 'ellipse', type = 'upper')
# # PerformanceAnalytics::chart.Correlation(as.matrix(ppw_rank[,-1] |> 
# #                                                     dplyr::select(-nobs)),
# #                                         histogram=TRUE, 
# #                                         cex.labels=50,
# #                                         method = "spearman")
# 
# # Get name of all mediation variables of interest
# mediation_res_colnames <- ppw_rank |> 
#   dplyr::select(-nobs, -geneSymbol, -AptName, #-Target,
#                 -all_of(colnames(node_metadata[,-1])), 
#                 -contains(".pval")) |>
#   colnames()
# 
# # Test all associations 
# resout <- epiomics::owas(
#   df = ppw_rank, 
#   var = colnames(node_metadata[,-1]), 
#   omics = mediation_res_colnames, 
#   var_exposure_or_outcome = "exposure")
# 
# cor(ppw_rank$eigen_centrality, ppw_rank$ACME, method = "spearman")
# # pm_quartile = gtools::quantcut(ppw_rank$pm, 4) |> as.numeric()

# get colors
ppw_rank1 <- ppw_rank %>% 
  mutate(Group = case_when(top_comptox + top_scrnaseq + top_prop_med == 3 ~ "Identified across all methods",
                           top_comptox == 1 & top_scrnaseq == 1 ~ "ComptoxAI and in-vitro", 
                           top_prop_med == 1 & top_scrnaseq == 1 ~ "TODAY and in-vitro", 
                           top_prop_med == 1 & top_comptox == 1 ~ "TODAY and ComptoxAI", 
                           top_comptox  == 1 ~ "ComptoxAI",
                           top_prop_med == 1 ~ "TODAY",
                           top_scrnaseq == 1 ~ "In-vitro",
                           TRUE ~ "None"))
table(ppw_rank1$top_comptox, ppw_rank1$top_prop_med, ppw_rank1$top_scrnaseq)
table(ppw_rank1$Group)
# Color = case_when(top_comptox + top_scrnaseq + top_prop_med == 3 ~ "#f54278",
#                   top_comptox == 1 & top_scrnaseq == 1 ~ ""
#                   top_30_hc == 1 ~ "#fc8c03",
#                   top_30_prop_med == 1 ~ "#4287f5",
#                   top_30_hc != 1&top_30_prop_med != 1 ~ "white"),





## c. Plot relationships of graph characteristics with mediation results ----
# (cor_plot <- ggplot(ppw_rank, aes(x = hub_score, y = Rpnie, label = geneSymbol)) +
(cor_plot <- ggplot(ppw_rank1, aes(x = harmonic_centrality, 
                                   y = PMed, 
                                   color = Group, 
                                   label = geneSymbol)) +
   ggrepel::geom_label_repel(force = 1, box.padding = .5, min.segment.length = 0) +
   # annotate("rect", xmin = -Inf, xmax = 15, ymin = -Inf, ymax = 15, fill = "lightgreen", alpha = 0.3) +
   # annotate("rect", xmin = 15, xmax = Inf, ymin = -Inf, ymax = 15, fill = "lightblue", alpha = 0.3) +
   # annotate("rect", xmin = -Inf, xmax = 15, ymin = 15, ymax = Inf, fill = "lightyellow", alpha = 0.3) +
   # annotate("rect", xmin = 15, xmax = Inf, ymin = 15, ymax = Inf, fill = "lightpink", alpha = 0.3) +
   geom_point(aes(size = abs(avg_log2FC))) + 
   # geom_point(aes(size = pm), shape = 21, color = "black", fill = "grey50") +
   # geom_point(aes(size = pm_quartile, fill = pm_quartile), shape = 21, color = "black") +
   scale_size(name = "Average\nLog2 FC\n") +
   xlab("Harmonic Centrality from ComptoxAI (Rank)") +
   ylab("Proportion Mediated\nfrom TODAY Study (Rank)") + 
   xlim(c(2, 32)) 
) 

quantile(ppw_rank$harmonic_centrality, 2/3)


ggsave(cor_plot, 
       filename = fs::path(dir_figure, "Fig3b_Centrality vs PM.jpg"), 
       height = 6, width = 6)





# 5. Plot Graph ------


V(g)$node_color <- case_when(
  V(g)$type == "Disease"                ~ "purple",  
  V(g)$type == "mediating-protein"      ~ "Light Blue",
  V(g)$type == "non-identified protein" ~ "grey50", 
  V(g)$type == "PFAS"                   ~ "red")


str(g)

plot.igraph(g, layout = layout_with_fr(g), 
            vertex.color=V(g)$node_color, 
            vertex.frame.color="#555555", 
            vertex.label=V(g)$id, vertex.label.color="black")


# Save 
?read.graph
read.graph(fs::path(dir_results, 
                    "ComptoxAI",
                    "PFAS_prot_dkd_expanded_020924_with_nodes.graphml"),
           format = "graphml")


# With NetworkD3


g_d3 <- igraph_to_networkD3(g, group = members)

# Create force directed network plot
forceNetwork(Links = karate_d3$links, Nodes = karate_d3$nodes, 
             Source = 'source', Target = 'target', 
             NodeID = 'name', Group = 'group')


# With GGraph
library(ggraph)
