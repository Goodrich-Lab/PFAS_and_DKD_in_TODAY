# Combine data across comptoxai, scRNAseq, and TODAY
source(here::here("!libraries.r"))
source(here::here("!directories.r"))
source(here::here("!load_clean_data.R"))
library(igraph)

# 1. Read data from python ----

# ComptoxAI graph generated in python
g <- read_graph(fs::path(dir_results,
                         "ComptoxAI",
                         "PFAS_prot_in_vitro_sig_fdr_trimmed_073024.graphml"),
                format = "graphml")

# Read in network analysis from python
graph_dat <- read_csv(fs::path(dir_results, 'network_graph_analysis.csv')) |> 
  tidylog::select(-EntrezGeneSymbol) |> 
  rename(geneSymbol = gene)

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


# 2. (new) Combine with network analysis data from python ------------
# 13 wont match because they were trimmed from the network
full_data <- node_metadata |> 
  tidylog::inner_join(graph_dat, by = "geneSymbol")


graph_dat_l_on_centrality <- pivot_longer(
  full_data, 
  cols = c(eigen_centrality:closeness), 
  values_to = "metric_value", names_to = "metric")


graph_dat_l_on_neo <- pivot_longer(
  full_data, 
  cols = c(eigen_centrality:closeness), 
  values_to = "metric_value", names_to = "metric")


## A. Plot ----
graph_dat |> 
  tidylog::filter(sig_overall != "scRNAseq only" ) |>
  ggplot(aes(x = tsne_dim1, y = tsne_dim2, 
             color = as.factor(sig_overall))) + 
  geom_point(size = 2, alpha = .5)

# By specific significance
graph_dat |> 
  tidylog::filter(!is.na(sig_overall)) |>
  ggplot(aes(x = sig_overall,
             y = tsne_dim2, 
             fill = sig_overall_simplified)) +
  geom_jitter(width = .2, alpha = .5) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_violin(alpha = .8) + 
  geom_boxplot(width = .1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

graph_dat$sig_overall

# 
graph_dat_l |> 
  filter(!is.na(sig_overall), 
         metric_value != 0) |> 
  ggplot(aes(x = sig_overall,
             y = log2(metric_value+.001), 
             fill = sig_overall_simplified)) +
  geom_jitter(width = .2, alpha = .5) + 
  # geom_hline(yintercept = 0, linetype = 2) + 
  geom_violin(alpha = .8) + 
  geom_boxplot(width = .1) + 
  facet_wrap(~metric, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


# 3. -----------------

# 2. (OLD) Proteomics meet in middle results ----
# Meet in middle results:
study_res <- read_csv(
  # fs::path(dir_results, "Combined_proteomics_scRNAseq_sig_results_073024.csv"))
  fs::path(dir_results, "Combined_proteomics_scRNAseq_sig_fdr_results_073024.csv"))


## a. Filter to select only 35 sig from meet in middle -------
# pfas_proteomics_35 <- pfas_proteomics_all |>
#   tidylog::filter(AptName %in% meet_in_middle_effects$feature_name)


## b. Combine data from proteomics and from comptoxai ----
ppw <- node_metadata |> 
  tidylog::filter(geneSymbol != "PFNA", 
                  geneSymbol != "C0011881") |> 
  tidylog::inner_join(study_res,
                      by = c("geneSymbol" = "EntrezGeneSymbol")) |>
  dplyr::select(geneSymbol, sig_overall, sig_overall_simplified, everything())


table(ppw$sig_overall_simplified)
table(study_res$sig_overall, study_res$PMed>0)

# # Determine top 10% of effects:
# ppw2 <- ppw |> 
#   mutate(#top_10_hub = ifelse(hub_score > quantile(PMed, 0.66), 1, 0), 
#          top_10_comptox = ifelse(degree_out > quantile(degree_out, 0.66), 1, 0), 
#          top_10_prop_med = ifelse(PMed > quantile(PMed, 0.66), 1, 0))
# 
# ppw_top_10_pct <- ppw2 |> 
#   tidylog::filter(top_10_comptox == 1, 
#                   top_10_prop_med == 1)
# ppw_top_10_pct
# hist(ppw$authority)
# Remove RSPO2, CHRNA5, FSCN1 because effect estimates are not possible?
# pfas_prot_network <- pfas_prot_network |>
#   tidylog::filter(!(geneSymbol %in% c ("CHRNA5", "FSCN1")))

# pivot data wider on mediation effect estimates
# ppw <- ppw |> 
#   pivot_wider(names_from = Effect, 
#               values_from = pe, 
#               id_cols = c(geneSymbol:closeness)) 


## c. Filter to only positive percent mediated ---------------------------------
ppw_pos <- ppw |> 
  tidylog::filter(PMed > 0) #|>
# mutate(ACME_abs = abs(ACME))

# Rank the data values
ppw_rank <- ppw_pos |> 
  tidylog::select(-geneSymbol) |> #, -c(AptName:ensembl_gene_id)
  tidylog::mutate_all(~rank(.)) %>%
  bind_cols(ppw_pos[,1], .) 

## d. correlation plots ----
(cor_matrix <- cor(janitor::remove_constant(ppw_rank[,-1]), 
                   method = "spearman"))
corrplot::corrplot(cor_matrix, method = 'ellipse', type = 'upper')
PerformanceAnalytics::chart.Correlation(as.matrix(ppw_rank[,-1] |> 
                                                    dplyr::select(-nobs)),
                                        histogram=TRUE, 
                                        cex.labels=50,
                                        method = "spearman")

# Get name of all mediation variables of interest
mediation_res_colnames <- ppw_rank |> 
  dplyr::select(-nobs, -geneSymbol, -AptName, -Target,
                -all_of(colnames(node_metadata[,-1])), 
                -contains(".pval")) |>
  colnames()

# Test all associations 
resout <- epiomics::owas(
  df = ppw_rank, 
  var = colnames(node_metadata[,-1]), 
  omics = mediation_res_colnames, 
  var_exposure_or_outcome = "exposure")

hist(ppw$page_rank)
cor(ppw_rank$eigen_centrality, ppw_rank$ACME, method = "spearman")
# pm_quartile = gtools::quantcut(ppw_rank$pm, 4) |> as.numeric()

# Plot relationships between graph characteristics and mediation results
# (cor_plot <- ggplot(ppw_rank, aes(x = hub_score, y = Rpnie, label = geneSymbol)) +
(cor_plot <- ggplot(ppw_rank, aes(x = harmonic_centrality, 
                                  y = PMed, 
                                  label = geneSymbol)) +
    ggrepel::geom_label_repel(force = 1, box.padding = .5, min.segment.length = 0) +
    # annotate("rect", xmin = -Inf, xmax = 15, ymin = -Inf, ymax = 15, fill = "lightgreen", alpha = 0.3) +
    # annotate("rect", xmin = 15, xmax = Inf, ymin = -Inf, ymax = 15, fill = "lightblue", alpha = 0.3) +
    # annotate("rect", xmin = -Inf, xmax = 15, ymin = 15, ymax = Inf, fill = "lightyellow", alpha = 0.3) +
    # annotate("rect", xmin = 15, xmax = Inf, ymin = 15, ymax = Inf, fill = "lightpink", alpha = 0.3) +
    geom_point(size = 2, shape = 21, color = "black", fill = "grey50") + 
    # geom_point(aes(size = pm), shape = 21, color = "black", fill = "grey50") +
    # geom_point(aes(size = pm_quartile, fill = pm_quartile), shape = 21, color = "black") +
    scale_size(name = "Percent\nMediated\n(Rank)") +
    xlab("Harmonic Centrality from ComptoxAI (Rank)") +
    ylab("Proportion Mediated\nfrom TODAY Study (Rank)") + 
    xlim(c(2, 32)) 
) 

quantile(ppw_rank$harmonic_centrality, 2/3)


ggsave(cor_plot, 
       filename = fs::path(dir_figure, "Fig3b_Centrality vs PM.jpg"), 
       height = 6, width = 6)





# 3. Plot Graph ------


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
