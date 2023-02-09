require(tidyverse)
require(rentrez)

#' Search pubmed with a query and return a tibble of the results
#' @param search_term string of a term to search pubmed
#' @param type string argument of 'full', 'authors', 'history' or 'articleids'
search_pubmed <- function(search_term, type='full', .retmax=20, link_pmid=FALSE){
  type <- match.arg(type, choices = c('full', 'authors', 'history', 'articleids'))
  pubmed_search <- rentrez::entrez_search(
    db="pubmed", term=search_term, retmax=.retmax, use_history=TRUE
  )
  if(length(pubmed_search$ids) < 1){stop('No publications found')}
  
  summary <- rentrez::entrez_summary(
    db="pubmed", id=pubmed_search$ids, web_history = pubmed_search$WebEnv,
    always_return_list=TRUE
  )
  out_table = switch(
    type,
    full=list_rbind(map(summary, ~as_tibble(keep(.x, is.character)))) %>% rename('pmid'=uid),
    authors=list_rbind(map(summary, ~as_tibble(.x$authors)), names_to = 'pmid'),
    history=list_rbind(map(summary, ~as_tibble(.x$history)), names_to = 'pmid'),
    articleids= list_rbind(map(summary, ~as_tibble(.x$articleids)), names_to = 'pmid')
    )
  return(
    out_table %>% 
      select_if(~!(all(is.na(.)) | all(. == ""))) %>% 
      {if(link_pmid) mutate(., pmid=sprintf("<a href='https://pubmed.ncbi.nlm.nih.gov/%s/' target='_blank'>%s</a>", pmid, pmid)) else .}
    )
}

