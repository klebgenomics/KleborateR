require(tidyverse)
require(igraph)

get_cluster_membership_from_graph <- function(kleborate_data, snp_graph) {
  if(!is_tibble(kleborate_data)){stop('kleborate_data is not a tibble')}
  if(!is_igraph(snp_graph)){stop('snp_graph is not an igraph')}
  igraph::components(snp_graph)$membership %>% 
    as_tibble(rownames = "Genome Name") %>% 
    mutate(Cluster = as.character(value)) %>% 
    select(!value) %>% 
    right_join(kleborate_data, by="Genome Name") %>% 
    mutate(Cluster=coalesce(as.character(.$Cluster), `Genome Name`))
}

pw_distmat_to_graph <- function(snp_data, snp_dist=10, directed=FALSE) {
  if(!is_tibble(snp_data)){stop('snp_data is not a tibble')}
  if(!'dist' %in% names(snp_data)){stop('dist column not in snp_data')}
  snp_data %>% 
    filter(dist <= snp_dist) %>% 
    select(-dist) %>% 
    as.matrix() %>% 
    igraph::graph_from_edgelist(., directed = directed)
}

get_cluster_membership_from_distmat <- function(
    kleborate_data, snp_data, snp_dist=10
){
  snp_data %>% 
    pw_distmat_to_graph(snp_dist) %>% 
    get_cluster_membership_from_graph(kleborate_data, .)
}

read_snp_diff <- function(distance_matrix_path){
  if(!fs::is_file(distance_matrix_path)){
    stop(paste(distance_matrix_path, 'is not a valid file'))
  }
  snp_data = read_csv(distance_matrix_path)
  if(nrow(snp_data) + 1 != ncol(snp_data)){
    stop('Number of rows and cols must be the same')
    }
  if(!'Name' %in% names(snp_data)){
    stop('Name column not in snp_data, is it from Pathogenwatch?')
    }
  return(
    snp_data %>% 
      pivot_longer(cols = !Name, values_to = 'dist', names_to = 'iso2') %>%
      rename(iso1=Name)
  )}
