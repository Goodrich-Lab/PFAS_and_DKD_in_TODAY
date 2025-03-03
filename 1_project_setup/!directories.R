# directory paths for file architecture

# home directory for project
dir_home <- here::here() |>
  dirname() |>
  dirname()


# data folder 
dir_data <- fs::path(dir_home,"0_data")

# report folder
dir_results <- fs::path(dir_home,"2_results")

# figure folder
dir_figure <- fs::path(dir_home, "3_figures")

