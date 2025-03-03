# Clean data

# 1. Clean putcomes, proteins, and metabolites ------
prot_metadata <- read_csv(fs::path(dir_data, 
                                   "clean_data",
                                   "Proteome metadata.csv"))

outcomes <- read_csv(fs::path(dir_data, 
                              "raw_data",
                              "TODAY data for ViCTER.csv")) |>
  janitor::remove_empty(which = c("cols"))

# Proteomics
proteome <- outcomes |> 
  select(SampleId, all_of(prot_metadata$AptName)) |>
  rename(sample_id = SampleId)

# Metabolomics
plasma_metabolome <- outcomes |>
  tidylog::select(SampleId, contains(".uM"), -contains("mM.Creatinine"))  |>
  rename_all(~str_remove(., ".in.uM")) |>
  rename_all(~paste0("met_", .)) |>
  janitor::clean_names() |>
  rename(sample_id = met_sample_id)

# Outcomes
outcomes_only <- outcomes |> 
  tidylog::select(-contains(".uM"),
                  -contains("mM.Creatinine"),
                  -starts_with("seq.")) |>
  janitor::remove_constant() |>
  tidylog::select(-releaseid, -index, -c(SampleDescription:Subarray),
                  -c(SampleMatrix:bsi_id.x), 
                  -c(current_label.y:Date.Drawn.y)) |>
  select(SampleId, everything()) |>
  janitor::clean_names()

# Clean outcomes 
outcomes_only <- outcomes_only |> 
  mutate(tx = as.factor(tx), 
         sex = if_else(sex == 1, "female", "male"), 
         monthstoglyc = daystoglyc/12) |>
  group_by(case1_control0) |>
  mutate(time_to_glyc_scld = scale(monthstoglyc) |> as.numeric()) |>
  ungroup()

# 1. Clean PFAS data ---- 
pfas_data <- read_csv(fs::path(dir_data, "raw_data",
                               "TODAY PFAS 320-95494-1_TalStandard.csv"), 
                      na = c("ND", "")) %>% 
  janitor::clean_names() |>
  tidylog::filter(!is.na(client_sample_id)) |>
  # janitor::remove_empty(which = c("rows", "cols")) |>
  janitor::remove_constant() 

# Two samples were entered incorrectly: 
pfas_data <- pfas_data |>
  tidylog::mutate(client_sample_id = case_when(
    client_sample_id == "353016922-6218313" ~ "3533016922-6218313",
    client_sample_id == "353301612-6196721" ~ "3533016612-6196721", 
    TRUE ~ client_sample_id))


## Remove Internal Standards (indicated with "13C", "18O2", d3, d5) ----
pfas_data_reduced <- pfas_data |>
  tidylog::filter(!is.na(high_limit))


## Impute 1/sqrt(lod) for those below the reporting limit 
pfas_data_reduced <- pfas_data_reduced |>
  mutate(pfas_imputed = if_else(is.na(result), 
                                high_limit/sqrt(2), 
                                result))

length(unique(pfas_data_reduced$analyte))
table(pfas_data_reduced$analyte)

## Get PFAS Metadata -----
pfas_metadata <- pfas_data_reduced |>
  tidylog::select(lab_sample_id, analyte, cas, result)  |>
  group_by(analyte) |>
  summarise(cas = cas[1], 
            n_na = sum(is.na(result)),
            pct_na = n_na/length(result), 
            n_unique = length(unique(result)), 
            pct_unique = n_unique/length(result), 
            geomean_ci = jag2::fungm(result, 
                                      n.digits = 2, 
                                      na.rm = TRUE)) |>
  # tidylog::filter(pct_na < .8) |>
  ungroup()

# Get common name 
pfas_metadata <- pfas_metadata |> 
  mutate(pfas_name = case_when(
    analyte == "Br-Perfluorooctanesulfonic acid" ~ "br_pfos",
    analyte == "L-Perfluorooctanesulfonic acid" ~ "l_pfos",
    analyte == "L-Perfluorooctanoic acid" ~ "l_pfoa",
    analyte == "M2-4:2 FTS" ~ "4_2_fts",
    analyte == "M2-8:2 FTS" ~ "8_2_fts",
    analyte == "NMeFOSAA" ~ "nmefosaa",
    analyte == "Perfluorodecanoic acid (PFDA)" ~ "pfda",
    analyte == "Perfluoroheptanesulfonic acid (PFHpS)" ~ "pfhps",
    analyte == "Perfluoroheptanoic acid (PFHpA)" ~ "pfhpa",
    analyte == "Perfluorohexanesulfonic acid (PFHxS)" ~ "pfhxs",
    analyte == "Perfluorononanoic acid (PFNA)" ~ "pfna",
    analyte == "Perfluoroundecanoic acid (PFUnA)" ~ "pfuna",
    analyte == "Total PFOA" ~ "pfoa",
    analyte == "Total PFOS" ~ "pfos"))


# PFAS metadata 
write_csv(pfas_metadata, fs::path(dir_data, "clean_data", "PFAS metadata.csv"))


# Join long data with metadata, and remove those which are not in metadata
pfas_data_reduced_2 <- pfas_data_reduced |> 
  tidylog::left_join(pfas_metadata |> 
                       tidylog::select(analyte, pfas_name)) |>
  tidylog::filter(!is.na(pfas_name))

## Pivot PFAS data wider ----
pfas <- pfas_data_reduced_2 |>  
  tidylog::select(client_sample_id, pfas_name, pfas_imputed) |>
  tidylog::pivot_wider(id_cols = c(client_sample_id), 
                       names_from = pfas_name, 
                       values_from = pfas_imputed) |>
  rename_all(~paste0("pfas_",. )) |> 
  rename(sample_id = pfas_client_sample_id)


# Save PFAS data 
write_csv(pfas, fs::path(dir_data, "clean_data", "TODAY_PFAS.csv"))

# Combine all data ------
full_data <- outcomes_only |>
  tidylog::full_join(pfas) |>
  tidylog::full_join(plasma_metabolome) |>
  tidylog::full_join(proteome)


write_rds(full_data, 
          fs::path(dir_data, 
                   "clean_data", 
                   "TODAY_complete_data_0417.RDS"))

