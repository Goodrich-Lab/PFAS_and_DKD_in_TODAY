# ! Load clean data 

# 1. Load exposure and outcome data ---------
## A. Read in data --------
original_data <- read_rds(fs::path(dir_data, "clean_data", "TODAY_complete_data_0417.RDS"))


# load protein meta data
meta_pro <- read_csv(fs::path(dir_data, "clean_data", "Proteome metadata.csv"))%>%
  filter(Organism == "Human")


# Set var names
pfas_names <- original_data |> tidylog::select(contains("pfas_")) |> colnames()
met_names <- original_data |> tidylog::select(contains("met_")) |> colnames()
prot_names <- original_data |> tidylog::select(all_of(meta_pro$AptName)) |> colnames()




# Get informative metabolites/protein names
omic_names <- c(prot_names, met_names)
var_uni_value <- original_data %>% 
  pivot_longer(cols = all_of(omic_names)) %>%
  group_by(name) %>%
  summarise(n = length(unique(value))) %>%
  filter(n == 1)

# Get all omic names
omic_names <- original_data %>% 
  dplyr::select(all_of(omic_names), 
                -all_of(var_uni_value$name)) |> 
  colnames()

# Log transform and scale Metabolites and Proteins
original_data <- original_data |>
  mutate_at(.vars = vars(contains("seq.")), .funs = ~log(.) %>% scale(.) %>% as.numeric(.)) %>% 
  mutate_at(.vars = vars(contains("met_")), .funs = ~log(.) %>% scale(.) %>% as.numeric(.)) 

## B. Set exposure and outcome names ----
outcome_glu <- c("codi", "si_1_ins0", "hb_a1c")

outcome_biomaker <- c("seq.2836.68",
                      "seq.9021.1",
                      "seq.5661.15",
                      "seq.11516.7",
                      "seq.17138.8",
                      "seq.12446.49",
                      "seq.15509.2")

pfas_names_all <- c("pfas_pfna", 
                    "pfas_pfda", 
                    "pfas_pfhpa", 
                    "pfas_pfoa",
                    "pfas_pfos", 
                    "pfas_pfhps",
                    "pfas_pfhxs",
                    "pfas_nmefosaa", 
                    "pfas_pfuna")

covars <- c("sex_male", "agebase", "eGFR", "dxtime")

levels <- c("PFUnA", "PFHpS","PFHpA","NMeFOSAA","PFDA","PFNA",
            "PFHxS", "PFOA", "PFOS", 
            "PFAS Burden Score",
            "PFSA Burden Score",
            "PFCA Burden Score")


## C. Calculate eGFR

# Function to calculate eGFR, based on CKD-EPI 2021 formula
calculate_eGFR <- function(Scr, Age, Gender) {
  # Constants for males and females
  if (Gender == 'female') {
    K <- 0.7
    alpha <- -0.241
    gender_factor <- 1.012
  } else if (Gender == 'male') {
    K <- 0.9
    alpha <- -0.302
    gender_factor <- 1
  } else {
    stop("Gender must be 'female' or 'male'")
  }
  
  # eGFR calculation
  eGFR <- 142 * (min(Scr/K, 1)^alpha) * (max(Scr/K, 1)^-1.200) * (0.9938^Age) * gender_factor
  
  return(eGFR)
}


original_data <- original_data |>
  rowwise() |>
  mutate(eGFR = calculate_eGFR(serum_creat, Age = agebase, Gender = sex)) |> 
  ungroup()

# 2. Calculate PFAS burden scores --------------
## A. Create categorical PFAS which are needed for the score ------
data_analysis <- original_data %>%
  mutate_at(.vars = vars(all_of(pfas_names_all[c(1,4,5,7,8)])), 
            .funs = list(
              quartile = ~as.integer(cut(., quantile(., probs = seq(0, 1, 0.25)), 
                                         include.lowest = TRUE
              )))) %>%
  mutate(pfas_pfda_detected = ifelse(pfas_pfda<0.05|pfas_pfda==0.2/sqrt(2), 1, 2)) %>%
  mutate(pfas_pfuna_detected = ifelse(pfas_pfuna<0.05|pfas_pfuna==0.1/sqrt(2), 1, 2)) %>%
  # mutate_at(.vars = vars(c("pfas_pfhpa","pfas_pfuna")),
  #           .funs = list(detected = ~ifelse(.<0.05|.==0.1/sqrt(2), 1, 2))) %>%
  tidylog::mutate_at(.vars = vars(c("pfas_pfhps", "pfas_pfhpa")), #, 
                     .funs = list(tert = ~cut(., 
                                                     quantile(., probs = seq(0, 1, 1/3)), 
                                                     include.lowest = TRUE) %>% 
                                    as.integer()))

# Select just the categorical PFAS
burden_score_pfas <- data_analysis %>% 
  dplyr::select(contains("detected"), 
                contains("quartile"), 
                contains("tert")) %>% 
  as.data.frame()

# Select PFCAs
burden_score_pfcas <- burden_score_pfas |> 
  dplyr::select(pfas_pfoa_quartile, 
                pfas_pfna_quartile, 
                pfas_pfda_detected, 
                pfas_pfuna_detected, 
                pfas_pfhpa_tert)

# Select PFSAs
burden_score_pfsas <- burden_score_pfas |> 
  dplyr::select(pfas_pfos_quartile, 
                pfas_pfhxs_quartile, 
                pfas_pfhps_tert)

## B. Calculate Scores ------
eap.pfas <- ltm::factor.scores(grm(burden_score_pfas), 
                               method="EAP", 
                               resp.patterns = burden_score_pfas)$score.dat$z1

eap.pfcas <- ltm::factor.scores(grm(burden_score_pfcas), 
                                method="EAP", 
                                resp.patterns = burden_score_pfcas)$score.dat$z1
eap.pfsas <- ltm::factor.scores(grm(burden_score_pfsas), 
                                method="EAP", 
                                resp.patterns = burden_score_pfsas)$score.dat$z1

# Plot scores
# plot(grm(burden_score_pfcas), type = "IIC",legend = TRUE,
#      xlab = "PFAS Burden", main = "", cx = "topright", cex = 0.6)

# Add burden scores back into data
data_raw <- data_analysis %>% 
  bind_cols(score = eap.pfas) %>%  
  bind_cols(score_pfcas = eap.pfcas) %>% 
  bind_cols(score_pfsas = eap.pfsas) 

rm(eap.pfsas, eap.pfcas, eap.pfas, 
   burden_score_pfas, burden_score_pfcas, burden_score_pfsas)

## C. Create Categorical Scores  ------
data2 <- data_raw %>% 
  mutate_at(.vars = vars(c(all_of(pfas_names_all), score)),
            .funs = list(median = ~ifelse(.<median(.),"0", "1"))) %>%
  mutate(score_tertile  = cut(score, quantile(score, probs = seq(0, 1, 1/3)), include.lowest = TRUE) %>% as.integer(), 
         score_quartile = cut(score, quantile(score, probs = seq(0, 1, 1/4)), include.lowest = TRUE) %>% as.integer(), 
         score_quintile = cut(score, quantile(score, probs = seq(0, 1, 1/5)), include.lowest = TRUE) %>% as.integer() )

## D. Create dummy vars, scale outcomes -----
data_log2 <- data2 |>
  fastDummies::dummy_cols(select_columns = c('sex'),
                          remove_selected_columns = FALSE, 
                          remove_most_frequent_dummy = TRUE) |>
  mutate_at(.vars = vars(c(outcome_glu)), 
            .funs = ~scale(.) %>% as.numeric(.)) %>%
  mutate_at(.vars = vars(all_of(pfas_names_all)), 
            .funs = ~ log2(.) %>% as.numeric(.)) # log transform here 

data_scaled <- data_log2 |> 
  tidylog::mutate_at(.vars = vars(all_of(pfas_names_all)), 
          .funs = ~ scale(.) %>% as.numeric(.)) # log transform here 

# Clean up data environment
rm(data_analysis, data2)

# 3. Proteomics metadata -------------------------------
# proteomics metadata specific to this project
prot_metadata1 <- read_csv(fs::path(dir_data, "raw_data", "protein_analytes.csv"), 
                           show_col_types = FALSE) |>
  dplyr::select(-c(Dilution:Dilution2)) 

prot_names_reduced <- prot_metadata1 |> 
  tidylog::filter(Organism == "Human")
prot_names_human <- prot_names_reduced$AptName

# read in gene expression by tissue information from Fagerberg et al
gene_exp_by_tissue <- read_csv(
  fs::path(dir_data,
           "protein_expression_by_organ",
           "Protein Atlas Fagerberg et al",
           "Fagerberg et al MCP 2014_cleaned.csv"), 
  show_col_types = FALSE)


# 4. Determine which tissue has highest expression of the genes --------
# Define a function to extract the highest and second highest values and their column names
extract_values <- function(row) {
  sorted_indices <- order(row, decreasing = TRUE)
  highest_value <- row[sorted_indices[1]]
  second_highest_value <- row[sorted_indices[2]]
  highest_column <- names(gene_exp_by_tissue)[c(4:30)][sorted_indices[1]]
  second_highest_column <- names(gene_exp_by_tissue)[c(4:30)][sorted_indices[2]]
  return(c(highest_column, highest_value, second_highest_column, second_highest_value))
}

# Apply the function to each row
extracted_values <- t(apply(gene_exp_by_tissue[,c(4:30)], 1, extract_values))

# Add the results to the dataframe
gene_exp_by_tissue$highest_expression <- extracted_values[, 1]
gene_exp_by_tissue$highest_expression_value <- as.numeric(extracted_values[, 2])
gene_exp_by_tissue$second_highest_expression <- extracted_values[, 3]
gene_exp_by_tissue$second_highest_expression_value <- as.numeric(extracted_values[, 4])
gene_exp_by_tissue <- gene_exp_by_tissue |> 
  mutate(high_second_high_ratio = highest_expression_value/second_highest_expression_value)



# Merge metadata and gene expression by tissue information
prot_metadata2 <- prot_metadata1 |> 
  filter(!is.na(EntrezGeneID)) |>
  tidylog::left_join(gene_exp_by_tissue, by = c("EntrezGeneID", "EntrezGeneSymbol")) #|>
  # dplyr::select(-c(colon:salivary_gland))

# Fix issues with missing Entrez gene name
prot_metadata2 <- prot_metadata2 |> 
  group_by(EntrezGeneSymbol) |> 
  mutate(rownum = row_number(), 
         maxrownum = max(rownum)) |> 
  ungroup() |>
  mutate(EntrezGeneSymbol = if_else(maxrownum > 1,
                                     str_c(EntrezGeneSymbol, "_", rownum), 
                                     EntrezGeneSymbol), 
    EntrezGeneSymbol = case_when(
    EntrezGeneID == "648791" ~ "PPP1R3G", # Missing Entrez gene name
    EntrezGeneID == "647087"  ~ "STMP1",
    TargetFullName == "Isthmin-1" ~ "ISM1",
    EntrezGeneID == "100134938" ~ "UPK3BL1",
    TRUE ~ EntrezGeneSymbol))

# Calculate ratio of expression in kidney versus other tissues. 
# Select key variables
prot_tissue_expression <- prot_metadata2 %>%
  dplyr::select(AptName, colon:salivary_gland) %>% 
  group_by(AptName) %>% 
  slice_head() %>%
  ungroup()

# Calculate sum of expression across tissues
prot_tissue_expression$sum_expression <-  prot_tissue_expression |> 
  dplyr::select(colon:salivary_gland) %>% 
  rowSums(na.rm = TRUE)

# Calculate the ratios
prot_tissue_expression_ratio <- prot_tissue_expression %>%
  mutate(across(where(is.numeric), ~(.)/sum_expression, 
                .names = "ratio_{.col}")) |> 
  dplyr::select(-c(colon:salivary_gland), -contains("sum_expression"))

# Add the ratios to the metadata
prot_metadata <- prot_metadata2 |> 
  tidylog::full_join(prot_tissue_expression_ratio) |>
  tidylog::filter(Organism == "Human") |> 
  dplyr::select(-c(colon:salivary_gland))


# Clean up data environment
# rm(prot_tissue_expression, prot_metadata1, prot_metadata2, prot_names_reduced,
#    prot_tissue_expression_ratio, gene_exp_by_tissue, 
#    extracted_values, extract_values,
#    var_uni_value)
