rename_pfas <- function(df){
  df <- df %>%
    mutate(pfas = case_when(grepl("pfos", pfas) ~ "PFOS",
                            grepl("pfoa", pfas) ~ "PFOA",
                            grepl("pfda", pfas) ~ "PFDA",
                            grepl("pfhpa", pfas) ~ "PFHpA",
                            grepl("pfhps", pfas) ~ "PFHpS",
                            grepl("pfhxs", pfas) ~ "PFHxS",
                            grepl("pfna", pfas) ~ "PFNA",
                            grepl("pfuna", pfas) ~ "PFUnA",
                            grepl("nmefosaa", pfas) ~ "NMeFOSAA",
                            grepl("score_pfsas", pfas) ~ "PFSA Burden Score",
                            grepl("score_pfcas", pfas) ~ "PFCA Burden Score",
                            grepl("score", pfas) ~ "PFAS Burden Score",
                            TRUE ~ pfas))
  return(df)
}

# Rename PFAS with Carbon Chain Length
rename_pfas_with_chainlength <- function(pfas_names_cleaned){
  pfas <- pfas_names_cleaned
  pfas = case_when(grepl("PFUnA", pfas) ~  "PFUnA (11)",
                   grepl("PFDA", pfas) ~  "PFDA (10)",
                   grepl("PFNA", pfas) ~  "PFNA (9)",
                   grepl("PFOA", pfas) ~  "PFOA (8)",
                   grepl("PFHpA", pfas) ~  "PFHpA (7)",
                   grepl("PFOS", pfas) ~  "PFOS (8)",
                   grepl("PFHpS", pfas) ~  "PFHpS (7)",
                   grepl("PFHxS", pfas) ~  "PFHxS (6)",
                   TRUE ~ pfas) 

}


# Get chain length
order_pfas_by_chain_length <- function(pfas_names_cleaned_chainlength){
  pfas <- pfas_names_cleaned_chainlength
  fct_relevel(pfas, 
              "PFUnA (11)",   # C 11 HF21O2
              "PFDA (10)",    # C 10 HF19O2
              "PFNA (9)",     # C 9 HF17O2
              "PFOA (8)",     # C 8 HF15O2
              "PFHpA (7)",    # C 7 HF13O2
              "PFOS (8)",    # C8HF17O3S
              "PFHpS (7)",   # C7F15SO3H
              "PFHxS (6)",   # C6F13SO3H	
              "NMeFOSAA", #C11H6F17NO4S
              "PFAS Burden Score",
              "PFCA Burden Score", 
              "PFSA Burden Score")
}


rename_outcomes <- function(x){
  nms <- case_when(
    x == "time_to_glyc_scld" ~ "Time to glycemic failure", 
    x == "codi" ~ "β-cell function (C-Peptide oDI)",
    x == "si_1_ins0" ~  "Ins. Sensitivity (1/fast. ins)",
    x == "hb_a1c" ~ "Hyperglycemia (HbA1c)", 
    x == "daystomic" ~ "Time to microalbuminuria", 
    x == "daystorapid" ~ "Time to rapid",
    x == "daystohyp" ~ "Time to hyperfiltration",
    x == "est_creat_clear" ~ "Estimated Creatinine Clearance",
    x == "daystomac" ~ "Time to macroalbuminuria",
    x == "uacid" ~ "Uric Acid",
    x == "u_alb_creat" ~ "u_alb_creat",
    x == "serum_creat" ~ "Serum Creatanine",
    x == "eGFR" ~ "Estimated Glomerular Filtration Rate",
    x == "seq.2836.68" ~"NGAL", 
    x == "seq.9021.1" ~ "KIM-1", 
    x == "seq.5661.15" ~ "IL18",
    x == "seq.11516.7" ~ "FABPL",
    x == "seq.17138.8" ~ "a-GST", 
    x == "seq.12446.49" ~ "a-GST1",
    x == "seq.15509.2" ~ "NAG")
}




# Rename PFAS -----------------------
rename_pfas_new <- function(pfas_names, include_asterisk = FALSE, arrange_by_class = FALSE){
  x <- tibble(pfas = pfas_names)
  pfas2 <-  x %>%
    mutate(pfas2 = case_when(
      grepl("pfos", pfas) ~ "PFOS",
      grepl("pfoa", pfas) ~ "PFOA",
      grepl("pfda", pfas) ~ "PFDA",
      grepl("pfhpa", pfas) ~ "PFHpA",
      grepl("pfhps", pfas) ~ "PFHpS",
      grepl("pfhxs", pfas) ~ "PFHxS",
      grepl("pfna", pfas) ~ "PFNA",
      grepl("pfuna", pfas) ~ "PFUnA",
      grepl("nmefosaa", pfas) ~ "NMeFOSAA",
      TRUE ~ toupper(pfas)) %>% 
        as.factor() %>% 
        fct_relevel(., 
                    "PFUnA", 
                    "PFDA", 
                    "PFNA", 
                    "PFOA", 
                    "PFHpA",
                    "PFOS", 
                    "PFHpS",
                    "PFHxS",
                    "NMeFOSAA")) 
  
  if(include_asterisk == TRUE){ 
    pfas2 <-  x %>%
      mutate(pfas2 = case_when(
        grepl("pfos", pfas) ~ "PFOS",
        grepl("pfoa", pfas) ~ "PFOA",
        grepl("pfda", pfas) ~ "PFDA",
        grepl("pfhpa", pfas) ~ "PFHpA",
        grepl("pfhps", pfas) ~ "PFHpS",
        grepl("pfhxs", pfas) ~ "PFHxS",
        grepl("pfna", pfas) ~ "PFNA",
        grepl("pfuna", pfas) ~ "PFUnA",
        grepl("nmefosaa", pfas) ~ "NMeFOSAA",
        TRUE ~ toupper(pfas)) %>% 
          as.factor() %>% 
          fct_relevel(., 
                      "PFUnA", 
                      "PFDA", 
                      "PFNA", 
                      "PFOA", 
                      "PFHpA",
                      "PFOS", 
                      "PFHpS",
                      "PFHxS",
                      "NMeFOSAA")) 
  }
  
  if(arrange_by_class == TRUE){ 
    pfas2 <-  pfas2 %>% 
      left_join(lod, by = "pfas") %>% 
      mutate(pfas2 = fct_reorder(pfas2, order_by_class))
  }
  
  return(pfas2$pfas2)
}

