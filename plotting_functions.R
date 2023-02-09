require(tidyverse)
require(glue)
require(gggenes)
require(ggrepel)
require(sf)
require(ggpol)
require(rnaturalearth)
require(rnaturalearthhires)

# TODO Replace all `aes_string()` usages ####
# TODO Make all functions 'safe' by adding column name checks ####
# TODO Make sure all functions refer to columns by strings that can be checked or passed as arguments ####

plot_raw_adj_pyramid <- function(
    counts, low_col = "#41b6c4", high_col = "#081d58", y_order = c('default'),
    fill_var = 'n_Large_cluster'
    ){
  
  # Check required columns --------------------------------------------------
  if ('default' %in% y_order){y_order = counts[[1]]}
  for(name in c(fill_var, 'Total')){
    if(!name %in% names(counts)){
      stop(glue::glue('name: {name} not in counts data'))
    }
  }
  # Prepare counts --------------------------------------------------
  counts %>%
    # This filter makes sure there are no empty columns after scaling by order
    filter(.[[1]] %in% y_order) %>% 
    mutate(raw_prop = raw_prop * -1) %>% # Reverse raw_props to make the pyramid
    pivot_longer(c(raw_prop, adj_prop), names_to = 'prop', values_to = 'percent') %>% 
    mutate(prop=as_factor(str_to_sentence(str_replace(prop, "_prop", " %")))) %>%
    
    # Plot data --------------------------------------------------
    ggplot(aes(x=.data[[names(.)[1]]], fill=.data[[fill_var]], label=Total, y=percent)) +
    geom_bar(stat="identity") +
    ggpol::facet_share(~prop, dir="h", scales="free", reverse_num=T) +
    coord_flip() +
    scale_x_discrete(limits = rev(y_order)) +
    scale_fill_gradient(
      str_replace_all(fill_var, "_", " "), low = low_col, high = high_col
      ) +
    labs(
      y=paste("Raw and outbreak-adjusted", 
              str_replace_all(names(counts)[1], "_", " "), "proportions"),
      x=NULL) +
    theme(
      legend.position = "top", legend.background = element_rect(color = NA),
      panel.border = element_blank(), panel.background = element_blank(),
      plot.background = element_blank()
      )
  }

plot_var_heatmap <- function(
    raw_adj_props, fill_var='adj_prop', low_col='yellow', high_col='red', 
    y_order=c('default'), max_x=20
    ){
  
  # Check required columns --------------------------------------------------
  for(name in c(fill_var, 'raw_count')){
    if(!name %in% names(raw_adj_props)){
      stop(glue::glue('name: {name} not in raw_adj_props data'))
    }
  }
  
  # X axis is ALWAYS the first column
  y = names(raw_adj_props)[1]
  
  # Default order is the first column
  if ('default' %in% y_order){
    y_order = raw_adj_props |>  pull(1) |>  unique()
    }
  
  # X axis is ALWAYS the second column
  x = names(raw_adj_props)[2]
  
  x_rank = raw_adj_props |> 
    group_by(.data[[x]]) |> 
    summarise(n=sum(adj_count)) |> 
    arrange(-n) |> 
    slice_head(n=max_x) |> # Slice the top adj values for ordered x-axis
    pull(1)
  
  raw_adj_props %>% 
    ggplot(
      aes(x=.data[[x]], y=.data[[y]], fill=.data[[fill_var]], label=raw_count)
           ) +
    geom_tile(colour='black', linewidth=0.4) + 
    geom_text(colour='black', fontface="bold") +
    theme_minimal() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_discrete(limits = rev(y_order)) +
    scale_x_discrete(limits = x_rank) +
    scale_fill_gradient("Adj %", low = low_col, high = high_col) +
    labs(
      caption='Numbers represent raw count\nColours represent adjusted percentage',
      x=str_replace_all(x, "_", " "), y=str_replace_all(y, "_", " ")
      )
}

summary_data <- function(kleborate_data) {
  return(
    tribble(
      ~name, ~value,
      'Samples:', n_distinct(kleborate_data$`Genome Name`),
      'Datasets:', n_distinct(kleborate_data$Study),
      'Countries:',n_distinct(kleborate_data$Country),
      'Sites:', n_distinct(kleborate_data$Site),
      'Clusters:', n_distinct(kleborate_data$Cluster),
      'Years:', n_distinct(kleborate_data$Year),
      'STs:', n_distinct(kleborate_data$ST),
      'KLs:', n_distinct(kleborate_data$K_locus),
      'OLs:', n_distinct(kleborate_data$O_locus))
  )}


plot_distinct_bars <- function(counts, var1, y_order=c('default')){
  if(!is.vector(y_order)){
    stop("y_order is not a vector but a ", typeof(y_order))}
  if ('default' %in% y_order){y_order = counts[[1]]}
  g <- counts %>% 
    filter(.[[1]] %in% y_order) %>%
    ggplot(aes(x = .data[[names(.)[[1]]]], y = .data[[var1]])) + 
    geom_bar(stat="identity", fill = "#000000") +
    theme_minimal() +
    theme(legend.position = "top") +
    scale_x_discrete(limits=rev(y_order)) +
    coord_flip() +
    labs(y = sprintf("`%s`", var1))
  return(g)}


plot_country_coverage <- function(raw_adj_props){
  sf::sf_use_s2(FALSE)  # Prevents MULTIPOINT error in R Shiny
  world <- rnaturalearth::ne_countries(returnclass = "sf", scale = "large")
  countries <- world %>% 
    right_join(raw_adj_props, by=c("name_long"="Country")) %>% 
    rename("Country"="name_long")
  world_cropped <- st_crop(world, st_bbox(countries)) %>% 
    # Prevents this bug https://github.com/plotly/plotly.R/issues/1785#issuecomment-643563527
    st_cast("MULTIPOLYGON")
  g <- ggplot(world_cropped) + 
    geom_sf(fill = "white", size=.1) +
    geom_sf(
      data=countries, size=0.4, col='black',
      aes(
        group = Country, label = raw_prop, fill = adj_prop,
        text = paste0("Accounting for ", raw_count, "/", raw_Country, " infections"))
      ) +
    scale_fill_viridis_c("Adj %", na.value='grey70', alpha = .8) +
    labs(
      x="Longitude",
      y="Latitude", 
      title="Adjusted coverage of infections by a single locus") +
    theme(legend.position = "right",
          panel.background = element_rect(fill = "aliceblue"))
  return(g)}

#' Plot the cumulative adjusted proportion of the 1st column, for each of the 
#' variables in the 2nd columns of the `tibble` output of the `raw_adj_prop()` function.
#' The x-axis and y-axis will always be the 1st and 2nd columns respectively. 
#' This means that the `denominator` argument should be set to the value of the 
#' `second` variable like so:
#' `plot_cumulative_coverage(raw_adj_prop(data, c("x", "y"), denominator="y"))`
#' @param raw_adj_props A tibble output of `raw_adj_prop()`
#' @param x_order A string vector
#' @param valency_intercept An integer vector
plot_cumulative_coverage <-  function(
    raw_adj_props, 
    x_order=c('default'),
    adjusted_prop=TRUE,
    valency_intercept=c(10, 20, 30)
    ){
  if(!is.vector(x_order)){
    stop("x_order is not a vector but a ", typeof(x_order))
  }
  prop_col = if_else(adjusted_prop == TRUE, 'adj_prop', 'raw_prop')
  if(!prop_col %in% names(raw_adj_props)){
    stop(paste(prop_col, "not found in raw_adj_props"))
    }
  
  valency_intercept = tibble(valen_x=valency_intercept) |> mutate(valen_y=100)
  x = names(raw_adj_props)[1]
  y = names(raw_adj_props)[2]
  prop_col_idx = which(names(raw_adj_props) == prop_col)
  if ('default' %in% x_order){
    x_order = raw_adj_props |> pull(1) |> unique()
  }
  raw_adj_props |> 
    select(all_of(c(1, 2, prop_col_idx))) |> 
    # Add 0% back in with two pivots and replace NA with 0
    pivot_wider(names_from = 2, values_from = 3) |> 
    mutate(across(!1, ~replace_na(., 0))) |> 
    pivot_longer(cols=!x, names_to = y, values_to = prop_col) |>
    right_join(   # Order with a join and add rank for x-axis
      tibble(x_order, seq(1, length(x_order)), .name_repair = ~ c(x, 'rank')),
      by=x
      ) |> 
    arrange(rank) |>
    group_by(across(2)) |> 
    mutate(!!prop_col := unlist(cumsum(across(all_of(prop_col))))) |>
    drop_na() %>%
    ggplot(aes(x=rank, y=.data[[prop_col]], colour=.data[[y]], label=.data[[x]])) +
    geom_line(alpha=0.4, linewidth=1.5) +
    theme_minimal() +
    labs(title=glue("Cumulative adjusted coverage per {y}"),
         subtitle=glue("Ordered by adjusted frequency of {x}"),
         y=glue("% coverage of {y}"), x=glue("{x}: most frequent > less frequent")) +
    geom_vline(
      xintercept=valency_intercept$valen_x, linetype="dashed", color = "grey", 
      linewidth=0.5
      )
}

#' Calculate the raw and outbreak adjusted proportions over time
#' @param kleborate_data A tibble of Kleborate output data
#' @param var1 A string corresponding to a column header to quantify over time
#' @param var2 A string corresponding to a column header to quantify over time
#' @export
plot_temporal_for_var <- function(kleborate_data, var1, var2, top_n=10) {
  
  # First, calculate the proportions of the variable for each country
  props <- kleborate_data %>% raw_adj_prop(var2, var1) %>%
    slice(1:top_n) # Limit the top n of var based on user input
  
  g <- kleborate_data %>% 
    # Next, calculate the proportions of the variable for each Year
    raw_adj_prop("Year", var1) %>% 
    right_join(props, by=names(props)) %>% # Right join the country props to order the data
    ggplot(aes_string(x="Year", y="adj_prop", group=var1, col=var1)) +
    geom_line(alpha=0.4, linewidth=1.5) +
    geom_point(aes(y=raw_prop), shape=3) +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    labs(y="Adjusted %", caption="+ = raw %") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(g)}

#' Calculate the raw and outbreak adjusted proportions over time
#' @param kleborate_data A tibble of Kleborate output data
#' @param var1 A string corresponding to a column header to quantify over time
#' @export
plot_temporal <- function(kleborate_data, var1, prop='adj_prop', top_n=20) {
  g <- raw_adj_prop(kleborate_data, 'Year', var1) %>%
    filter(.[[var1]] %in% get_counts(kleborate_data, var1)[[var1]][1:top_n]) %>% 
    ggplot(aes(x=Year, y=.data[[prop]], group=.data[[var1]], col=.data[[var1]])) +
    geom_line(alpha=.4, linewidth=1.5) +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    theme_minimal() +
    labs(y=paste(str_to_sentence(str_replace(prop, '_prop', ' proportion'))), "%") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(g)
  }

#' Plot a Kaptive locus with `gggenes`
#' @param LOCUS_GENES A `tibble` imported from `antigen_db.xlsx`
#' @export
plot_locus <- function(LOCUS_GENES){
  LOCUS_GENES %>%
    mutate(middle = START+(GENE_LENGTH/2)) %>% # Add midpoint to each gene
    ggplot(
      aes(xmin=START, xmax=END, y=LOCUS, x=middle, fill=GENE_TYPE,
          label=GENE_NAME, forward=STRAND)
      ) +
    scale_fill_manual(
      values =  c(
        'Core genes'='#0c7fb8',
        'Initiating glycosyltransferase'='#17a658', 
        'Other sugar synthesis and processing'='#a0cc91',
        'Flippase'='#dc4b43',
        'Rhamonse synthesis and processing' = '#C6B3D3',
        'Mannose synthesis and processing' = '#643F95',
        'Capsule repeat unit polymerase'='#e68747', 
        'Hypothetical/unknown protein'='#898a8b',
        'Transposase'='#ffffff',
        'Insertion sequence'='#EFC000FF',
        'O-antigen biosynthesis'="#B03060"
          )
      ) +
    gggenes::geom_gene_arrow(
      arrowhead_height = unit(10, "mm"), arrow_body_height=unit(10, "mm"), 
      arrowhead_width = unit(2, "mm"), show.legend = TRUE
      ) +
    geom_text(angle = 45, hjust = 0.2, nudge_y=0.2, size = 4, fontface = "bold") +
    facet_wrap(vars(LOCUS), ncol = 1, scales = "free_y") + 
    labs(y = NULL, x='Length (basepairs)', fill='Gene Type') +
    guides(fill = guide_legend(ncol = 4)) + gggenes::theme_genes() + 
    theme(legend.position='bottom')
}

#' Plot a SNP network graph calculated from from a distance matrix
#' @param snp_graph A ggnetwork object
#' @param kleborate_data A `tibble` of Kleborate output data
#' @param var1 A string corresponding to a column header to colour nodes
#' @export
plot_transmission_network <- function(snp_graph, kleborate_data, var1) {
  snp_graph %>% 
    left_join(
      kleborate_data %>% select(all_of(c('Genome Name', var1))),
      by=c('name'='Genome Name')) %>%
    ggplot(aes(x=x, y=y, xend=xend, yend=yend, label=name)) +
    geom_edges() +
    geom_nodes(
      aes(col=.data[[var1]]), 
      size=3.5, shape=16, alpha=.5
    ) +
    theme_blank() +
    guides(col=FALSE, fill=FALSE)
}

clear_plot_sides <- function(plot){
  plot + theme(
    axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
    axis.line.y = element_blank(), axis.text.y = element_blank())
}

plot_antigen_structure <- function(pdb_path) {
  r3dmol::r3dmol() |>
    r3dmol::m_add_model(
      data = readr::read_file(pdb_path),
      format = "pdb", keepH = TRUE
    ) |>
    r3dmol::m_set_style(style = r3dmol::m_style_stick()) |>
    r3dmol::m_add_surface(style = r3dmol::m_style_surface(opacity = 0.4)) |>
    r3dmol::m_add_outline() |>
    r3dmol::m_add_label(
      fs::path_ext_remove(fs::path_file(pdb_path)),
      style = r3dmol::m_style_label(backgroundOpacity = 0.4)
    ) |>
    r3dmol::m_zoom_to()
}

#' plot_antigen_image
#' @description A function to plot mutliple jpeg files as a patchwork plot
#' by reading the jpegs, converting them to rasters and plotting with ggplot2
#' @param molecule_file_df A tibble with columns `locus` and `path` which point 
#' to jpeg files
#'
#' @return A patchwork plot of the jpeg files
#' @export
plot_antigen_image <- function(molecule_file_df) {
  patchwork::wrap_plots(
    purrr::map2(
      molecule_file_df$locus, molecule_file_df$path, 
      ~ggplot2::ggplot() +
        ggplot2::annotation_custom(grid::rasterGrob(jpeg::readJPEG(.y, native = TRUE))) +
        ggplot2::ggtitle(.x) + 
        ggplot2::theme_void()
    ), ncol = 1)
}

