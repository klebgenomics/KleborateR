require(tidyverse)

#' quant_kaptive_genes
#' @description
#' A general function to take a gene column from a Kaptive output and count
#' the number of each gene
#' It can take the `Genes in locus` column but also the `Missing genes` column
#' if the delimiter is set to a column.
#' `quant_kaptive_genes(df$K_locus_missing_genes, delim = ",")`
#' @param genes_column A list from a Kaptive output table
#' @param delim A string delimiter between each gene
#' @export
quant_kaptive_genes <- function(
    genes_column, 
    delim = ";",
    sep_cols = c("locus", "gene_number", "gene_name")
    ) {
  return(
    stringr::str_split(genes_column, delim, simplify = TRUE) |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(cols=everything()) |>
      dplyr::pull(value) %>%
      stringr::str_split_fixed("_", n=3) |>
      tibble::as_tibble() %>%
      rlang::set_names(sep_cols) |>
      dplyr::filter(!is.na(Locus), Locus!="") %>%
      dplyr::count(Locus, Gene_number, Gene_name) |>
      dplyr::mutate(perc = n/length(genes_column)*100) |>
      dplyr::arrange(-n)
  )
}

#' quant_missing_genes
#' @description
#' A more specific function to quantify the number of missing genes in a 
#' Pathogenwatch Kleborate output, and gives extra counts for groups passed
#' with the `group_vars` argument.
#' @param kaptive_data A tibble from a Kleborate/Kaptive run
#' @param antigen A string indicating either the K or O antigen
#' @param gene_col_delim A string delimiter between each gene
#' @param group_vars A character vector of column names to count groups
#' @export
quant_missing_genes <- function(
    kaptive_data,
    antigen = "K",
    gene_col_delim = ",",
    group_vars = c("Country", "ST")
    ) {
  if(!is_tibble(kaptive_data)){stop("kaptive_data is not a tibble")}
  
  locus_var = switch(
    match.arg(as.character(antigen), c("K", "O")), K = 'K_locus', O = 'O_locus'
  )
  missing_var = paste0(locus_var, "_missing_genes")
  
  for(i in c(locus_var, missing_var, group_vars)){
    if(!i %in% names(kaptive_data)){stop(paste(i, "not in kaptive_data"))}
    }
  
  count_vars = c(group_vars, locus_var)
  all_vars = c(count_vars, missing_var)
  count_var_name = paste(count_vars, collapse = "_per_")
  
  return(
    kaptive_data |>
      dplyr::mutate(index = 1:n()) |>
      dplyr::select(index, dplyr::all_of(all_vars)) |>
      dplyr::add_count(across(dplyr::all_of(count_vars)), name = count_var_name) |>
      tidyr::pivot_longer(dplyr::all_of(count_vars)) |>
      dplyr::add_count(value, name = 'total') |>
      tidyr::pivot_wider(names_from = name, values_from = c(value, total)) |>
      dplyr::rename_with(
        .fn = ~stringr::str_remove(.x, "^value_"), 
        .cols = tidyselect::starts_with("value_")
      ) |>
      dplyr::filter(if_any(dplyr::all_of(missing_var), ~!is.na(.x))) |>
      tidyr::separate_rows(dplyr::all_of(missing_var), sep = gene_col_delim) |>
      tidyr::pivot_longer(dplyr::all_of(count_vars)) |>
      dplyr::group_by(.data[[missing_var]]) |>
      dplyr::mutate(total_missing_gene = n()) |>
      dplyr::add_count(value, name = 'missing_in') |>
      tidyr::pivot_wider(names_from = name, values_from = c(value, missing_in)) |>
      dplyr::rename_with(
        .fn = ~stringr::str_remove(.x, "^value_"), 
        .cols = tidyselect::starts_with("value_")
      ) |> 
      dplyr::select(!index) |>
      dplyr::relocate(dplyr::all_of(c(missing_var, locus_var, group_vars, count_var_name))) |>
      dplyr::distinct() |>
      dplyr::arrange(-total_missing_gene)
    )
  }

