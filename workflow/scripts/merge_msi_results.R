log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

rlang::global_entrace()

library(tidyverse)

combined_files_with_header <- read_tsv(
    file = snakemake@input[["msi_results"]],
    id = "group"
  ) |>
  rename(
    n_all_sites = Total_Number_of_Sites,
    n_unstable_sites = Number_of_Unstable_Sites,
    msi_score =	`%`
  ) |>
  mutate(
    group = str_match(
        group,
        "results/[^/]+/(?<match>[^/]+)/.+"
    )[,"match"]
  )

write_tsv(
    x = combined_files_with_header,
    file = snakemake@output[["tsv"]]
)