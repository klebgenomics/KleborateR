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
    adj_vars=c('Cluster', 'Site', 'Study', 'ST'),
    epi_vars = c("ST", "K_locus", "K_type", "O_locus", "O_type", "Country", 
                  "Study", "Site", "Cluster")
){
  
  # Check args
  if(!is_tibble(kleborate_data)){stop("kleborate_data must be a tibble")}
  if(!is.character(var1)){stop("var1 must be a string")}
  for(i in unique(c(var1, adj_vars, epi_vars))){
    if(!i %in% names(kleborate_data)){
      stop(paste(i, "not in", base::quote(kleborate_data)))
    }
  }
  if(large_cluster_size != round(large_cluster_size)){
    stop("large_cluster_size must be an integer")}
  
  # We NEED to remove var1 from *_vars for the all_of() functions
  if(var1 %in% epi_vars){
    message(paste('Removing', var1, 'from epi_vars'))
    epi_vars = epi_vars[!epi_vars == var1]
  }
  if(var1 %in% adj_vars){
    message(paste('Removing', var1, 'from adj_vars'))
    adj_vars = adj_vars[!adj_vars == var1]
  }
  # Add Large_cluster if Cluster column present
  if('Cluster' %in% names(kleborate_data)){  # Identify large clusters
    kleborate_data <- kleborate_data |>
      dplyr::add_count(Cluster) |> 
      dplyr::mutate(
        Large_cluster =dplyr::if_else(n >= large_cluster_size, Cluster, NA_character_)
        ) |> 
      dplyr::select(!n)
    epi_vars = c(epi_vars, 'Large_cluster')
  }
  message(paste("Grouping var:", var1))
  message(paste("Epi vars:", paste(epi_vars, collapse = ", ")))
  message(paste("Adj vars:", paste(adj_vars, collapse = ", ")))
  
  return(
    kleborate_data |>
      dplyr::reframe(
        .by = all_of(var1),
        # n = nrow(kleborate_data),
        raw_count = dplyr::n(),
        adj_count = dplyr::n_distinct(dplyr::across(tidyselect::all_of(adj_vars))),
        dplyr::across(
          tidyselect::all_of(epi_vars), n_distinct, .names = "n_{.col}"
        )
      ) |>
      dplyr::mutate(
        raw_prop = raw_count/sum(raw_count),
        adj_prop = adj_count/sum(adj_count)
      ) |>
      dplyr::select(
        tidyselect::all_of(var1), raw_count, adj_count, raw_prop, adj_prop, dplyr::everything()
      ) |>
      dplyr::arrange(-adj_prop)
  )}

#' Calculate the raw and outbreak adjusted proportions of Kleborate columns
#' @description
#' `raw_adj_prop()` is used to calculate the raw and outbreak adjusted 
#' proportions of multiple Kleborate columns relative to one another. For example,
#' if grouping_vars = c('K_locus', 'ST'), the proportions are the percent of each
#' K-locus that has a particular ST, with total/adjusted K-locus counts being
#' the denominator ('K_locus' by default).
#' If we wanted to reverse this, and see what proportion of each ST has a 
#' particular K-locus, we can change the denominator to 'ST', or simply
#' change the variable order c('ST', 'K_locus').
#' WARNING: due to the nature of the `all_of()` function, the grouping vars
#' must be removed from the outbreak adjustment vars, which may lead to
#' inconsistencies in the outbreak adjustment method.
#' @param kleborate_data A tibble
#' @param grouping_vars A string vector of columns to calculate proportions for
#' @param denominator A string referring to the column used as the proportion
#' denominator
#' @param adj_vars A string vector of columns used to perform outbreak adjustment
#' @export
raw_adj_prop <- function(
    kleborate_data,
    grouping_vars,
    denominator = 'default',
    adj_vars = c('Cluster', 'Site', 'Study', 'ST')
    ){
  
  # Check args
  if(!is_tibble(kleborate_data)){stop("kleborate_data must be a tibble")}
  if(!is.character(denominator)){stop("denominator must be a string")}
  if(denominator=='default'){denominator=grouping_vars[1]}
  for(i in c(adj_vars, grouping_vars)){
    if(!i %in% names(kleborate_data)){
      stop(paste(i, "not in"  base::quote(kleborate_data)))
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
  message(paste("Denominator:", denominator))
  message(paste("Adj vars:", paste(adj_vars, collapse = ", ")))
  
  return(
    kleborate_data |>
      dplyr::reframe(
        .by = tidyselect::all_of(grouping_vars),
        raw_count = dplyr::n(),
        adj_count = dplyr::n_distinct(across(all_of(adj_vars))),
        ) |>
      dplyr::mutate(
        .by = tidyselect::all_of(denominator),
        raw_prop = raw_count/sum(raw_count),
        adj_prop = adj_count/sum(adj_count)
      ) |>
      dplyr::distinct() |>
     dplyr:: arrange(-adj_count)
  )} 

#' Filter genomes in a Kleborate output tibble
#' @description
#' `genome_filter()` removes undesirable genomes from Kleborate results
#' using pre-defined parameters that can be tweaked.
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
    k_typable=TRUE, 
    o_typable=TRUE,
    max_contigs=275, 
    max_size=6500000,
    min_size=2500000
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
