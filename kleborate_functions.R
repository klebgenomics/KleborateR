require(tidyverse)
require(fs)
require(glue)

#' Calculate the raw and outbreak adjusted counts of a Kleborate column
#' @description
#' `get_counts()` is used to calculate the raw and outbreak adjusted 
#' proportions of a Kleborate column
#' WARNING: due to the nature of the `all_of()` function, var1 must be removed 
#' from the outbreak adjustment vars, which may lead to inconsistencies in the 
#' outbreak adjustment method.
#' @param kleborate_data A tibble
#' @param var1 A string referring to a column to calculate counts for
#' @param large_cluster_size An integer for the minimum number of genomes in a 
#' large cluster
#' @param adj_vars A string vector of columns used to perform outbreak adjustment
#' @param epi_vars A string vector of columns to count distinct occurrences of
#' per var1 group
#' @export
get_counts <- function(
    kleborate_data, 
    var1, 
    large_cluster_size = 3,
    adj_vars=c('Cluster', 'Site', 'Country'),
    epi_vars = c("ST", "K_locus", "K_type", "O_locus", "O_type", "Cluster")
){
  
  # Check args
  if(!is_tibble(kleborate_data)){stop(paste(base::quote(kleborate_data), "must be a tibble"))}
  if(!is.character(var1)){stop("var1 must be a character string")}
  for(i in unique(c(var1, adj_vars, epi_vars))){
    if(!i %in% names(kleborate_data)){
      stop(paste(i, "not in", base::quote(kleborate_data)))
    }
  }
  if(large_cluster_size != round(large_cluster_size)){
    stop("large_cluster_size must be an integer")}
  
  # We NEED to remove var1 from *_vars for the all_of() functions
  if(var1 %in% epi_vars){
    message(paste('var1:', var1, 'in epi_vars, removing from epi_vars'))
    epi_vars = epi_vars[!epi_vars == var1]
  }
  if(var1 %in% adj_vars){
     message(paste('var1:', var1, 'in adj_vars, removing from adj_vars'))
    adj_vars = adj_vars[!adj_vars == var1]
  }
  # Add Large_cluster if Cluster is in adj_vars
  # Addresses: change if statement to determine whether large cluster calculation is needed #2
  # Only really applicable for plotting, consider depreciating
  if('Cluster' %in% adj_vars){  # Identify large clusters
    kleborate_data <- kleborate_data |>
      dplyr::add_count(Cluster, name = "cluster_size") |> 
      dplyr::mutate(
        Large_cluster = dplyr::if_else(cluster_size >= large_cluster_size, Cluster, NA_character_)
        ) |> 
      dplyr::select(!cluster_size)
    epi_vars = c(epi_vars, 'Large_cluster')
  }
  message(paste("Grouping var:", var1))
  message(paste("Epi vars:", paste(epi_vars, collapse = ", ")))
  message(paste("Adj vars:", paste(adj_vars, collapse = ", ")))

  #if ("K_locus" %in% var1){
  #	kleborate_data <- kleborate_data %>%
  #		filter(K_locus != "unknown (KL107)") %>%
  #		filter(!(K_locus =="KL107" & K_locus_confidence=="None")) %>%
  #		mutate(K_locus = str_replace(K_locus, "unknown \\(", "")) %>%
  #		mutate(K_locus = str_replace(K_locus, "\\)", ""))
 # }
    
  return(
    kleborate_data |>
      dplyr::reframe(  # Perform a raw and adjusted count
        .by = all_of(var1),  # Per group (mitigates a group_by function call)
        # n = nrow(kleborate_data),
        raw_count = dplyr::n(),                                                      # Raw
        adj_count = dplyr::n_distinct(dplyr::across(tidyselect::all_of(adj_vars))),  # Adjusted
        dplyr::across(
          tidyselect::all_of(epi_vars), n_distinct, .names = "n_{.col}"
        )
      ) |>
      dplyr::mutate(  # Perfrom proportion calculation
        raw_prop = raw_count/sum(raw_count),  # Raw
        adj_prop = adj_count/sum(adj_count)   # Adjusted
      ) |>
      dplyr::select(  # Re-order the cols
        tidyselect::all_of(var1), raw_count, adj_count, raw_prop, adj_prop, dplyr::everything()
      ) |>
      dplyr::arrange(-adj_prop) # Sort
  )}

#' Calculate the raw and outbreak adjusted counts and proportions for variable 
#' combinations of interest.
#' 
#' @description
#' `raw_adj_prop()` is used to calculate raw and outbreak-adjusted counts, stratified by 
#' combinations of variables. The output will include counts and proportions for all 
#' observed combinations of the grouping_vars, and adjusted proportions will use as the
#' denominator all unique combinations of the grouping_vars AND adj_vars (which would 
#' usually be clusters).
#' For example, to calculate K-locus counts and frequencies per country, adjusted for
#' outbreak clusters, we set grouping_vars=c("K_locus", "Country"), summarise_by="Country"
#' and leave the default setting of adj_vars="Cluster" in order to produce outbreak-adjusted
#' counts and proportions. Or, to see the frequency of STs per K-locus, we set 
#' grouping_vars=c("K_locus", "ST"), summarise_by="K_locus".
#' Note that all unique values for the specified variables are counted separately, and 
#' contribute to unique variable combinations that will count as separate clusters.
#' This means that two strains in the same cluster and country, one with K_locus='KL112' 
#' and the other with K_locus='unknown (KL112)' will be counted as 2 distinct clusters.
#' 
#' WARNING: due to the nature of the `all_of()` function, var1 must be removed 
#' from the outbreak adjustment vars, which may lead to inconsistencies in the 
#' outbreak adjustment method.
#' 
#' @param kleborate_data A tibble
#' @param grouping_vars A string vector of columns to calculate proportions for
#' @param summarise_by A string referring to the column used to summarise proportion
#' denominators by
#' @param adj_vars A string vector of columns used to perform outbreak adjustment
#' @export
raw_adj_prop <- function(
    kleborate_data,
    grouping_vars = c("K_locus", "Country"),
    summarise_by = "Country",
    adj_vars = c("Cluster")
    ){
  
  # Check args
  if(!is_tibble(kleborate_data)){stop(paste(base::quote(kleborate_data), "must be a tibble"))}
  if(!is.character(summarise_by)){stop("summarise_by must be a string (variable name)")}
  if(summarise_by=='default'){summarise_by=grouping_vars[1]}
  for(i in c(adj_vars, grouping_vars)){
    if(!i %in% names(kleborate_data)){
      stop(paste(i, "not in", base::quote(kleborate_data)))
    }
  }

  # Check grouping vars
  for (i in grouping_vars) { 
    if(i %in% adj_vars){
      message(paste('grouping_var:', i, 'in adj_vars, removing from adj_vars'))
      adj_vars = adj_vars[!adj_vars == i]
      }
  }
  message(paste("Grouping vars:", paste(grouping_vars, collapse = ", ")))
  message(paste("Summarising by:", summarise_by))
  message(paste("Adj vars:", paste(adj_vars, collapse = ", ")))
  
  #if ("K_locus" %in% grouping_vars){
  	kleborate_data <- kleborate_data %>%
  		filter(K_locus != "unknown (KL107)") %>%
  		filter(!(K_locus =="KL107" & K_locus_confidence=="None")) %>%
  		mutate(K_locus = str_replace(K_locus, "unknown \\(", "")) %>%
  		mutate(K_locus = str_replace(K_locus, "\\)", ""))
  #}
  
  return(
    kleborate_data |>
      dplyr::reframe(  # Perform a raw and adjusted count
        .by = tidyselect::all_of(grouping_vars),  # Per group (mitigates a group_by function call)
        raw_count = dplyr::n(),                                   # Raw
        adj_count = dplyr::n_distinct(across(all_of(adj_vars))),  # Adjusted
        ) |>
      dplyr::mutate(  # Perfrom proportion calculation
        .by = tidyselect::all_of(summarise_by),  # Group by the summary variable (mitigates a group_by function call)
        raw_prop = raw_count/sum(raw_count),  # Raw
        adj_prop = adj_count/sum(adj_count),   # Adjusted
        raw_sum = sum(raw_count),
        adj_sum = sum(adj_count)
      ) |>
      dplyr::distinct() |>
     dplyr:: arrange(-adj_count)
  )}  


#' Filter genomes in a Kleborate output tibble
#' @description
#' `genome_filter()` removes undesirable genomes from Kleborate results
#' using pre-defined parameters that can be tweaked.
#' Defaults for contig count and genome size are those set by KlebNET GSP, see
#' https://bigsdb.pasteur.fr/klebsiella/genome-quality-check/
#' @param kleborate_data A tibble
#' @param species A string vector of species to keep
#' @param k_typable Logical to drop K-locus confidence calls of "Low" and "None"
#' @param o_typable Logical to drop O-locus confidence calls of "Low" and "None"
#' @param max_contigs Integer for the max number of contigs in each genome
#' @param max_size Integer for the maximum size of each genome
#' @param min_size Integer for the minimum size of each genome
#' @export
genome_filter <- function(
    kleborate_data, 
    species = c("Klebsiella pneumoniae"), 
    k_typable=FALSE, 
    o_typable=TRUE,
    max_contigs=500, 
    max_size=6200000,
    min_size=5000000
    ){
  # Check args
  if(!is_tibble(kleborate_data)){stop("kleborate_data must be a tibble")}
  if(!is.character(species) | !is.vector(species)){
    stop("species must be a character vector")}
  for(i in c(k_typable, k_typable)){
    if(!is.logical(i)){
      stop(paste(deparse(substitute(i)), "must be logical"))}
  }
  for(i in c(max_contigs, max_size, min_size)){
    if(i != round(i)){
      stop(paste(deparse(substitute(i)), "must be an integer"))}
  }

  kleborate_data %>% 
    filter(species %in% species) %>% 
    filter(total_size <= max_size) %>% 
    filter(total_size >= min_size) %>% 
    filter(contig_count <= max_contigs) %>%
    {if(k_typable) filter(., !K_locus_confidence %in% c('Low', 'None')) else .} %>%
    {if(o_typable) filter(., !O_locus_confidence %in% c('Low', 'None')) else .}
}

#' Clean and homogenise Kleborate data, adapted from Kaptive-Web
#' @description
#' `clean_data()` takes a Kleboate result tibble and cleans the columns
#' to be compatible with the Klebsiella Sero-epi Shiny app.
#' Some country names will be changed to be used with the `join_world_data()`
#' function, which is called within this function, but really this should be 
#' done at the data pre-processing step.
#' As with the Kaptive-Web function, simplified ESBL and Carbapenemase columns
#' added.
#' @param kleborate_data A tibble
#' @export
clean_data <- function(kleborate_data) {
  return(
    kleborate_data %>%
      # Tidy up country data for plotting map
     dplyr::mutate(Country = case_when(
        Country == 'UK' ~ 'United Kingdom',
        Country == 'Laos' ~ 'Lao PDR',
        Country == 'Republic of Ireland' ~ 'Ireland',
        Country == 'USA' ~ 'United States',
        # For Eva's Caribbean study in collaboration with (CARPHA), 
        # its headquartered in Trinidad and Tobago so setting this as the country
        Country == 'Caribbean' ~ 'Trinidad and Tobago',
        TRUE ~ Country)) %>% 
      join_world_data() %>%
      
      # convert scientific name to common name
     dplyr::mutate(Source =dplyr::if_else(str_detect(Source, 'homo'), 'human', Source)) %>% 
      
      # unify age groups
     dplyr::mutate(
        Age_group = case_when(
         stringr::str_detect(Age_group, 'neo') ~ 'neonatal',
         stringr::str_detect(Age_group, 'adul') ~ 'adult',
          as.numeric(Age_group) >= 16 ~ 'adult',
          as.numeric(Age_group) < 16 ~ 'child',
          as.numeric(Age_group) == 0 ~ 'neonatal',
          is.na(Age_group) ~ 'unknown',
          TRUE ~ tolower(Age_group))) %>%
      
     dplyr::mutate(Sample = str_replace(Sample,'[()]|^other_', "")) %>% 
     dplyr::mutate(Sample = str_replace(Sample,'_', " ")) %>%
      
      # Clean the _locus/_type columns from [unknown (best match = )]
     dplyr::mutate(
       dplyr::across(
          matches("_locus$|_type$"),
          ~if_else(str_detect(.x, 'unknown'), str_extract(.x, LOCUS_TYPE_REGEX), .x)
          )
        ) %>%
      # simplify omp
     dplyr::mutate(Omp_mutations_simplified = str_replace_all(Omp_mutations, "-[0-9]+%", "-trunc"), 
             Omp_simple =dplyr::if_else(Omp_mutations == "-", "wt", "mut")) %>%
      
      # simplify carbapenemases and combine with omp
     dplyr::mutate(Bla_Carb_simplified = case_when(
       stringr::str_detect(Bla_Carb_acquired, "IMP") ~ "IMP", 
       stringr::str_detect(Bla_Carb_acquired, "KPC") ~ "KPC",
       stringr::str_detect(Bla_Carb_acquired, "OXA") ~ "OXA", 
       stringr::str_detect(Bla_Carb_acquired, "NDM") ~ "NDM",
       stringr::str_detect(Bla_Carb_acquired, "VIM") ~ "VIM", 
       stringr::str_detect(Bla_Carb_acquired, ";") ~ "multiple",
       stringr::str_detect(Bla_Carb_acquired, "[A-Z]+") ~ "other",
        TRUE ~ "-")) %>%
     dplyr::mutate(carbapenemase_omp_combination = paste(Bla_Carb_simplified, Omp_simple, sep = " ")) %>%
      
      # simplify ESBLs and combine with omp
     dplyr::mutate(Bla_ESBL_simplified = case_when(
       stringr::str_detect(Bla_ESBL_acquired, "CTX-M") ~ "CTX-M-other", 
        Bla_ESBL_acquired == "CTX-M-14" ~ "CTX-M-14",
        Bla_ESBL_acquired == "CTX-M-15" ~ "CTX-M-15",
        Bla_ESBL_acquired == "CTX-M-65" ~ "CTX-M-65",
       stringr::str_detect(Bla_ESBL_acquired, "SHV") ~ "SHV",
       stringr::str_detect(Bla_ESBL_acquired, "TEM") ~ "TEM",
       stringr::str_detect(Bla_ESBL_acquired, ";") ~ "multiple",
       stringr::str_detect(Bla_ESBL_acquired, "[A-Z]+") ~ "other",
        TRUE ~ "-")) %>%
     dplyr::mutate(ESBL_omp_combination = paste(Bla_ESBL_simplified, Omp_simple, sep = " ")) %>%
      
      # simplify bla acquired and combine with omp
     dplyr::mutate(Bla_acq_simplified = case_when(
       stringr::str_detect(Bla_acquired, "TEM") ~ "TEM",
       stringr::str_detect(Bla_acquired, "OXA") ~ "OXA",
       stringr::str_detect(Bla_acquired, "LAP") ~ "LAP",
       stringr::str_detect(Bla_acquired, "DHA") ~ "DHA",
       stringr::str_detect(Bla_acquired, ";") ~ "multiple", 
       stringr::str_detect(Bla_acquired, "[A-Z]+") ~ "other",
        TRUE ~ "-")) %>%
     dplyr::mutate(Bla_acquired_omp_combination = paste(Bla_acq_simplified, Omp_simple, sep = " ")) %>%
      
      #  # rmpADC lineage simplification
      # dplyr::mutate(rmpADC_simplified = case_when(
      #   stringr::str_detect(RmpADC, "rmp") ~ str_extract(RmpADC, "rmp [0-9]+"),
      #   stringr::str_detect(RmpADC, "rmp unknown") ~ "rmp unknown",
      #   stringr::str_detect(RmpADC, "rmp 2A") ~ "rmp 2A",
      #   stringr::str_detect(RmpADC, ",") &stringr::str_detect(RmpADC, "rmp") ~ "multiple rmp",
      #    TRUE ~ "-")) %>%
      # 
      #  # rmpADC truncations
      # dplyr::mutate(rmpADC_trunc = case_when(
    #   stringr::str_detect(RmpADC, "rmp") ~ "intact", 
    #   stringr::str_detect(RmpADC, "incomplete") ~ "truncated",
    #    TRUE ~ "-")) %>%
    # 
    #  # rmpA2 truncations
    # dplyr::mutate(rmpA2_trunc = case_when(
    #   stringr::str_detect(rmpA2, "rmp") ~ "intact",
    #   stringr::str_detect(rmpA2, "%") ~ "truncated",
    #    TRUE ~ "-")) %>%
    
    # convert sample and source to lowercase
   dplyr::mutate(
      Sample = str_to_lower(Sample), 
      Source = str_to_lower(Source), 
      Age_group = str_to_lower(Age_group)
      )
  )
}

#' Adds `rnaturalearth` geographic information for plotting with `sf`
#' @description
#' `join_world_data()` takes a Kleboate result tibble and a specified column of 
#' geographic information to join on additional geographical information
#' specified by `info_cols` using the `ne_col` column.
#' It will show a message if a variable in the specified geographic column is
#' not found in the specified column in `rnaturalearth::ne_countries`
#' @param kleborate_data A tibble
#' @param geo_col A string of the geographic col to join `rnaturalearth::ne_countries`
#' @param ne_col A string of the corresponding `rnaturalearth::ne_countries` col
#' @param info_cols A string vector of `rnaturalearth::ne_countries` cols to add
#' @export
join_world_data <- function(
    kleborate_data,
    geo_col = "Country",
    ne_col='name_long',
    info_cols = c("continent", "income_grp", "region_un", "subregion", "region_wb")
    ){
  if(!geo_col %in% names(kleborate_data)){
    stop(paste(geo_col, "not in", base::quote(kleborate_data)))
  }
  ne_data = sf::st_drop_geometry(
    rnaturalearth::ne_countries(
      returnclass = "sf", type = "countries", scale = "large"
      )
    )
  ne_cols = unique(c(ne_col, info_cols))
  for(col in ne_cols){
    if(!col %in% names(ne_data)){
      stop(glue::glue("{col} not in ne_data"))
    }
  }
  missing_countries = setdiff(kleborate_data[[geo_col]], ne_data[[ne_col]])
  if (length(missing_countries) > 0) {
    message(glue::glue('{missing_countries} were not found in {geo_col}'))
    }
  
 return(
   kleborate_data |> 
     left_join(
      dplyr::select(ne_data, all_of(ne_cols)),
       by = setNames(nm=geo_col, ne_col)
       ) |>
     rename_with(stringr::str_to_title, .cols = any_of(ne_cols))
 )
}

trait_stats <- function(kleborate_data, trait1, trait2, stat_test, 
                        conf_level = 0.95,  adjusted = FALSE) {
  
  stat_test = match.arg(
    stat_test, c('chisq.test', 'prop.test', 'cor.test', 'cor'))
  
  if(stat_test %in% c('xsq', 'prop')){
    count_type =dplyr::if_else(adjusted == TRUE, 'adj', 'raw')
    count_cols = paste(count_type, c('count', trait1), sep = "_")
    names(count_cols) <- c('successes', 'trials')
    props = raw_adj_prop(kleborate_data, c(trait1, trait2))
  }

  if(stat_test == 'prop.test'){
    props |>
     dplyr::select(all_of(c(trait1, trait2, count_cols))) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        tst = list(
          broom::tidy(prop.test(successes, trials, conf.level = conf_level))
        )
      ) |>
      tidyr::unnest(tst)
  } else if(stat_test == 'chisq.test'){
    props |>
      dplyr::select(all_of(c(trait1, trait2, count_cols[1]))) |>
      tidyr::pivot_wider(names_from = 2, values_from = 3, values_fill = 0) |>
      tibble::column_to_rownames(trait1) |>
      as.matrix() |>
      chisq.test()
  } else if(stat_test == 'cor.test'){
    cor.test(kleborate_data[[trait1]], kleborate_data[[trait2]])
  } else if(stat_test == 'cor'){
 dplyr::select(kleborate_data, where(is.numeric)) |>
    as.matrix() |>
    cor()
  }
}

kleborate_column_spec <- vroom::cols(
  `Genome ID` = vroom::col_character(),
  `Genome Name` = vroom::col_character(),
  Version = vroom::col_character(),
  `Kleborate version` = vroom::col_character(),
  strain = vroom::col_character(),
  species = vroom::col_character(),
  species_match = vroom::col_character(),
  contig_count = vroom::col_double(),
  N50 = vroom::col_double(),
  largest_contig = vroom::col_double(),
  total_size = vroom::col_double(),
  ambiguous_bases = vroom::col_character(),
  QC_warnings = vroom::col_character(),
  ST = vroom::col_character(),
  virulence_score = vroom::col_double(),
  resistance_score = vroom::col_double(),
  num_resistance_classes = vroom::col_double(),
  num_resistance_genes = vroom::col_double(),
  Yersiniabactin = vroom::col_character(),
  YbST = vroom::col_character(),
  Colibactin = vroom::col_character(),
  CbST = vroom::col_character(),
  Aerobactin = vroom::col_character(),
  AbST = vroom::col_character(),
  Salmochelin = vroom::col_character(),
  SmST = vroom::col_character(),
  RmpADC = vroom::col_character(),
  RmST = vroom::col_character(),
  rmpA2 = vroom::col_character(),
  wzi = vroom::col_character(),
  K_locus = vroom::col_character(),
  K_type = vroom::col_character(),
  K_locus_problems = vroom::col_character(),
  K_locus_confidence = vroom::col_character(),
  K_locus_identity = vroom::col_character(),
  K_locus_missing_genes = vroom::col_character(),
  O_locus = vroom::col_character(),
  O_type = vroom::col_character(),
  O_locus_problems = vroom::col_character(),
  O_locus_confidence = vroom::col_character(),
  O_locus_identity = vroom::col_character(),
  O_locus_missing_genes = vroom::col_character(),
  AGly_acquired = vroom::col_character(),
  Col_acquired = vroom::col_character(),
  Fcyn_acquired = vroom::col_character(),
  Flq_acquired = vroom::col_character(),
  Gly_acquired = vroom::col_character(),
  MLS_acquired = vroom::col_character(),
  Phe_acquired = vroom::col_character(),
  Rif_acquired = vroom::col_character(),
  Sul_acquired = vroom::col_character(),
  Tet_acquired = vroom::col_character(),
  Tgc_acquired = vroom::col_character(),
  Tmt_acquired = vroom::col_character(),
  Bla_acquired = vroom::col_character(),
  Bla_inhR_acquired = vroom::col_character(),
  Bla_ESBL_acquired = vroom::col_character(),
  Bla_ESBL_inhR_acquired = vroom::col_character(),
  Bla_Carb_acquired = vroom::col_character(),
  Bla_chr = vroom::col_character(),
  SHV_mutations = vroom::col_character(),
  Omp_mutations = vroom::col_character(),
  Col_mutations = vroom::col_character(),
  Flq_mutations = vroom::col_character(),
  truncated_resistance_hits = vroom::col_character(),
  spurious_resistance_hits = vroom::col_character(),
  Chr_ST = vroom::col_character(),
  gapA = vroom::col_double(),
  infB = vroom::col_double(),
  mdh = vroom::col_double(),
  pgi = vroom::col_double(),
  phoE = vroom::col_double(),
  rpoB = vroom::col_double(),
  tonB = vroom::col_double(),
  ybtS = vroom::col_character(),
  ybtX = vroom::col_character(),
  ybtQ = vroom::col_character(),
  ybtP = vroom::col_character(),
  ybtA = vroom::col_character(),
  irp2 = vroom::col_character(),
  irp1 = vroom::col_character(),
  ybtU = vroom::col_character(),
  ybtT = vroom::col_character(),
  ybtE = vroom::col_character(),
  fyuA = vroom::col_character(),
  clbA = vroom::col_character(),
  clbB = vroom::col_character(),
  clbC = vroom::col_character(),
  clbD = vroom::col_character(),
  clbE = vroom::col_character(),
  clbF = vroom::col_character(),
  clbG = vroom::col_character(),
  clbH = vroom::col_character(),
  clbI = vroom::col_character(),
  clbL = vroom::col_character(),
  clbM = vroom::col_character(),
  clbN = vroom::col_character(),
  clbO = vroom::col_character(),
  clbP = vroom::col_character(),
  clbQ = vroom::col_character(),
  iucA = vroom::col_character(),
  iucB = vroom::col_character(),
  iucC = vroom::col_character(),
  iucD = vroom::col_character(),
  iutA = vroom::col_character(),
  iroB = vroom::col_character(),
  iroC = vroom::col_character(),
  iroD = vroom::col_character(),
  iroN = vroom::col_character(),
  rmpA = vroom::col_character(),
  rmpD = vroom::col_character(),
  rmpC = vroom::col_character(),
  spurious_virulence_hits = vroom::col_character()
  )

read_kleborate_file <- function(path) {
  if(fs::is_file(path) && !fs::is_file_empty(path)){
    return(readr::read_csv(path, col_types = kleborate_column_spec, 
                           show_col_types = FALSE))
  } else {
    warning(paste(path, "is not a valid file"))
    return(NULL)
  }
}

read_kleborate_files <- function(paths, show_progress=FALSE) {
  purrr::list_rbind(purrr::map(paths, read_kleborate_file, 
                               .progress = show_progress))
}
