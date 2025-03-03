# Pathway analysis of graph data

source(here::here("!libraries.r"))
source(here::here("!directories.r"))
# library(topGO)
# library(org.Hs.eg.db)

graph_with_clusters <- read_csv(fs::path(dir_results, "network_graph_analysis.csv"))
graph_with_clusters$order = c(1:nrow(graph_with_clusters))

# 1. IPA analysis ----------
## 1.a. Format data for IPA ---------
ipa_data <- graph_with_clusters |>
  tidylog::filter(gene != "PFNA", 
                  gene != "C0011881") |>
  tidylog::select(EntrezGeneSymbol, 
                  leidenalg_community_8) |> 
  mutate(leidenalg_community_8 = leidenalg_community_8+1, 
         value = 0) |>
  arrange(leidenalg_community_8) |>
  # Pivot data wider so that each community is represented as a column
  pivot_wider(names_from = leidenalg_community_8,
              values_from = value,
              values_fill = list(value = 1),
              names_prefix = "community_") #%>%

# Here, the values represent p-values, with 0 if included and 1 if not included
# Write data for IPA
write_csv(ipa_data, fs::path(dir_results, "ipa_data.csv"))


# 1.b. Read back in IPA results ---------
# Get a list of all .tsv files in the folder
tsv_files <- list.files(path = fs::path(dir_results, "comptoxai", "pathway analysis"), 
                        pattern = "\\.txt", full.names = TRUE)

# Read all TSV files into a list
tsv_list <- map(tsv_files, ~read_tsv(., skip = 1) |> 
                  janitor::clean_names() %>%
                  dplyr::rename(., "name" = 1))

# Print the names of the files for confirmation
names(tsv_list) <- basename(tsv_files) |> str_remove(".txt")

# Inspect the first few rows of each data frame (optional)
ipa_pw <- bind_rows(tsv_list, .id = "type") |> 
  rename_all(~str_replace(., "observation", "community"))


# Pivot longer
ipa_pw_l <- ipa_pw |> 
  pivot_longer(cols = community_1:community_8, names_to = "community") |> 
  mutate(community_numeric = str_remove(community, "community_") %>% as.numeric(), 
         community_numeric = community_numeric-1)


# 2. Read in graph to get names of communities -----------
library(igraph)

g <- read_graph(fs::path(dir_results,
                         "ComptoxAI",
                         "PFAS_prot_in_vitro_sig_100924.graphml"),
                format = "graphml")

## 2.a. Calculate graph data ------
node_metadata <- tibble(
  gene = V(g)$id,
  order = 1:length(V(g)$id),
  leidenalg_community = V(g)$leidenalg_community,
  leidenalg_community_8 = V(g)$leidenalg_community_8,
  eigen_centrality = igraph::eigen_centrality(g)$vector,
  authority = igraph::authority_score(g)$vector,
  harmonic_centrality = harmonic_centrality(g, mode = "out"),
  hub_score = hub_score(g)$vector,
  page_rank = page_rank(g)$vector,
  strength = strength(g),
  betweenness = igraph::betweenness(g),
  degree_all   = igraph::degree(g, mode = "all"),
  degree_in    = igraph::degree(g, mode = "in"),
  degree_out   = igraph::degree(g, mode = "out"),
  degree_total = igraph::degree(g, mode = "total"),
  eccentricity = eccentricity(g),
  coreness = coreness(g),
  closeness = closeness(g)
)

# # List of column names to be used in the regression
# metrics <- c("eigen_centrality", "authority", "harmonic_centrality", "hub_score", 
#              "page_rank", "strength", "betweenness", "degree_all", 
#              "degree_in", "degree_out", "degree_total", "eccentricity", 
#              "closeness")
# # Run the linear models and tidy the results using lapply
# res_lst <- lapply(metrics, function(metric) {
#   formula <- as.formula(paste0("scale(", metric, ") ~ -1 + leidenalg_community_8"))
#   tidy(lm(formula, data = node_metadata))
# })
# # Name the list elements with corresponding metric names
# names(res_lst) <- paste0(metrics, "_res")
# # Combine the results into a single data frame
# lm_res <- bind_rows(res_lst, .id = "name")
# 
# ggplot(lm_res, aes(x = term, y = estimate, color = name, group = name)) + 
#   geom_point() + geom_line() + 
#   facet_wrap(~name) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


## 2.b. Merge with node_metadata ------
# Sanity check:
table(V(g)$id == graph_with_clusters$gene) # The order is the same
table(V(g)$leidenalg_community_8 == graph_with_clusters$leidenalg_community_8) 


# Combine data from graph analysis with data from python on clusters
node_data <- tidylog::full_join(node_metadata, graph_with_clusters) 


## 2.c. Calculate top comptoxai hubs -------
select_top_hubs <- node_data |> 
  group_by(leidenalg_community_8) |>
  arrange(-hub_score, -harmonic_centrality) |>
  tidylog::slice_head(n = 3) |>
  summarise(genes = str_c(gene, collapse = "; ")) |>
  ungroup()


## 2.d. Summarise cell types in top  nodes ----
sig_cells <- node_data %>%
  pivot_longer(
    cols = c("week1_cycling":"week2_mesangial_cells_2"), 
    names_to = "week_cellname_number", 
    values_to = "p_value"
  ) %>%
  # Filter out rows with NA p-values
  tidylog::filter(!is.na(p_value)) |>
  mutate(cellname = str_remove(week_cellname_number, "cells_") |> 
           str_replace_all("_", " ") |>
           str_remove("week1 ") |> str_remove("week2 ") |>
            tools::toTitleCase())

cell_type_summary <- sig_cells |> 
  group_by(leidenalg_community_8, cellname) |> 
  summarise(n_sig = length(cellname)) |>
  ungroup()


num_sig_genescommunities <- sig_cells |>
  group_by(leidenalg_community_8) |> 
  summarise(total_sig = length(cellname)) |>
  ungroup()

cell_type_sig <- tidylog::full_join(num_sig_genescommunities, 
                                    cell_type_summary)

cell_type_sig <- cell_type_sig |> 
  mutate(pct_sig = n_sig/total_sig)


ggplot(cell_type_sig, 
       aes(x = leidenalg_community_8, y = pct_sig)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~cellname, scales = "free")

# No real difference in the cell types between the communities



# 3. Plot Pathway Results ----------

## 3.a Prep Data ----------
# rename clusters with decriptions of top clusters
ipa_pw_l <- ipa_pw_l |> 
  mutate(community_name = case_when(
    community_numeric == 2 ~ "Grey Community (PFNA)", 
    community_numeric == 1 ~ "Blue Community (DKD)", 
    community_numeric == 0 ~ "Green Community (THOC1)", 
    community_numeric == 3 ~ "Orange Community (SPRED1/MFN1)", 
    community_numeric == 4 ~ "Pink Community", 
    community_numeric == 5 ~ "Turqiose Community", 
    community_numeric == 6 ~ "Lilac Community", 
    community_numeric == 7 ~ "Purple Community")) |> 
  mutate(community_name = factor(
    community_name,
    levels = c("Grey Community (PFNA)", 
               "Blue Community (DKD)", 
               "Green Community (THOC1)", 
               "Orange Community (SPRED1/MFN1)", 
               "Pink Community", 
               "Turqiose Community", 
               "Lilac Community", 
               "Purple Community")))


# Get dataframe for cluster colors for ggplot figure
ipa_pw_l <- ipa_pw_l |> 
  mutate(community_color = case_when(
    community_name == "Green Community (THOC1)" ~ "#73C000",
    community_name == "Blue Community (DKD)" ~ "#00C4FF",
    community_name == "Grey Community (PFNA)" ~ "#4C463E",
    community_name == "Orange Community (SPRED1/MFN1)" ~ "#FF8805",
    community_name == "Pink Community" ~ "#FF5584",
    community_name == "Turqiose Community" ~ "#00BD94",
    community_name == "Lilac Community" ~ "#D3B3B0",
    community_name == "Purple Community" ~ "#DF89FF") )


# Select top 10 and order by factor levels
top_5 <- ipa_pw_l |> 
  tidylog::group_by(community_name) |> 
  tidylog::filter(type == "canonical pathways") |>
  tidylog::top_n(3, value) |>
  ungroup() |>
  mutate(name = fct_reorder(name, -value))




## 3.b. Plot ----
# Bargraph of top results using ggplot, facet wrap for communities
(plotout <- ggplot(top_5, aes(y = value, 
                  x = fct_rev(name), 
                  fill = community_color)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, color = "grey50") +
  coord_flip() +
  facet_wrap(~community_name, scales = "free_y", ncol = 1) +
  labs(x = NULL, y = "-log10 p-value") + 
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        strip.background = element_blank())) 


# scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
# theme(legend.position = "none")

ggsave(plotout, filename = fs::path(dir_figure, "Pathway Results Long.svg"), 
       width = 7.5, height = 9)


#
#
#
#
#
#



















# 2. (Not Run): topGO analysis --------
# Run analysis across all communities --------
# List of community columns to iterate over
community_columns <- paste0("community_", 1:length(unique(graph_with_clusters$leidenalg_community_8)))

# Initialize an empty list to store results for each community
all_results <- list()

# Loop through each community
for (community in community_columns) {
  
  cat("Running pathway analysis for:", community, "\n")
  
  # Creating a topGO compatible list (factor)
  gene_list <- factor(as.integer(ipa_data[[community]] == 0))
  names(gene_list) <- ipa_data$EntrezGeneSymbol
  
  # Define the topGO Data Object
  GOdata <- new("topGOdata",
                ontology = "BP", # BP for Biological Process, or "MF" or "CC"
                allGenes = gene_list,
                geneSelectionFun = function(x) x == 1,
                annot = annFUN.org,
                mapping = "org.Hs.eg.db", # For human genes
                ID = "symbol") # Depending on your gene identifier type
  
  # Run the Enrichment Analysis
  resultFisher_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  
  # Extract and interpret the results
  resultTable <- GenTable(GOdata,
                          classicFisher = resultFisher_classic,
                          weight01Fisher = resultFisher_weight01,
                          orderBy = "classicFisher",
                          topNodes = 100)
  
  # Calculate p-value
  resultTable <- resultTable %>% 
    mutate(p_fisher = str_remove(classicFisher, "< ") %>% as.numeric(),
           community = community) # Add community information
  
  # # Visualize Results 
  # par(cex = .3)
  # showSigOfNodes(GOdata, score(resultFisher), 
  #                firstSigNodes = 10, useInfo = "all")
  # dev.off()
  
  # Store the result in the list
  all_results[[community]] <- resultTable
}



# Combine all results into a single data frame
final_results <- bind_rows(all_results) 

# FDR
final_results <- final_results |> 
  mutate(FDR = p.adjust(p_fisher, method = "fdr"), 
         enrichment = Significant/Annotated)

# Significant Results 
sig_results <- final_results |> 
  tidylog::filter(FDR < 0.01) |>
  arrange(FDR, -1*enrichment)


# Compare unique pathways from GO analysis -----------
library(UpSetR)
sig_w <- sig_results |> 
  tidylog::select(community, GO.ID, FDR) |> 
  tidylog::pivot_wider(names_from = community, 
                       values_from = FDR, 
                       values_fill = NA) |> 
  column_to_rownames("GO.ID") |> 
  mutate(across(starts_with("community_"), ~ ifelse(is.na(.), 0, 1)))

table(sig_results$community)


UpSetR::upset(sig_w, 
              sets = community_columns, 
              order.by = "freq")
empty.intersections = "on")

# Plot top results -----------
# Bargraph of top results using ggplot, facet wrap for communities
sig_results |> 
  group_by(community) |> 
  top_n(10, FDR) |> 
  ggplot(aes(x = reorder(Term, enrichment), y = enrichment, )) +
  geom_col() +
  coord_flip() +
  facet_wrap(~community, scales = "free_y", ncol = 2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Top 10 Pathways by Community",
       x = "GO.ID",
       y = "FDR") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(legend.position = "none")
