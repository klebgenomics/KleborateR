## testing_ground.R
##
## This file is to test new functions or other code snippets which will find
## their way into 'kleborate_functions.R' once they are complete. 


weighting_country_and_study <- function(df, countries){
  
  # This function requires functions get_counts and raw_adj_prop
  
  new_counts <- get_counts(df, 'Country')

  # SECTION ONE
  # If there is more than one country then we will have to add a weighting 
  # for each country. This section of code does that. 
  
  if(length(countries) > 1){
    
    raw_adj_prop_output <- raw_adj_prop(df, c("Country", "Study", "K_locus"))
    
    # Calculate national weighting
    
    country_proportions <- count(raw_adj_prop_output, Country)
    country_proportions$country_prop <- country_proportions$n/sum(country_proportions$n)
    raw_adj_prop_output$country_weight <- NA
    
    raw_adj_prop_output <- inner_join(raw_adj_prop_output, country_proportions)
    
    
  }
  
  olddf <- tibble(K_locus =c(NA), raw_weight =c(NA), adj_weight =c(NA), country=c(NA))
  
  # SECTION TWO
  # Once that has been completed or skipped, this section of code deals with
  # applying weights to each study within the country or countries. 
  
  for (c in 1:length(countries)){
    
    country_filtered <- raw_adj_prop_output %>% filter(., Country == countries[c])
    
    # Section 2 Step 1: Assigns weights to studies per country, unless there is
    #                   only 1 study
    
    if(length(unique(country_filtered$Study)) > 1){
      
      study_proportions <- count(country_filtered, Study)
      study_proportions$study_prop <- study_proportions$n/sum(study_proportions$n)
      country_filtered <- inner_join(country_filtered, study_proportions, by = "Study")
      
      # Section 2 Step 2: Generate weighted counts
     
      country_filtered <- country_filtered %>% 
        mutate(.,raw_count_weight = raw_count * study_prop * country_prop) %>% 
        mutate(.,adj_count_weight = adj_count * study_prop * country_prop) 
      sum_raw <- sum(country_filtered$raw_count_weight)
      sum_adj <- sum(country_filtered$adj_count_weight)
      country_filtered <- country_filtered %>% 
        mutate(.,raw_weight = raw_count_weight/sum_raw *100) %>% 
        mutate(.,adj_weight = adj_count_weight/sum_adj *100) 
      
      newdf <- tibble(K_locus = unique(country_filtered$K_locus))
      
      
      for(i in 1:length(unique(country_filtered$K_locus))){
        Kloc <- as.character(newdf[i,"K_locus"])
        country_filtered_KL <- country_filtered %>% filter(., K_locus == paste(Kloc))
        newdf[i,"raw_weight"] <- sum(country_filtered_KL$raw_weight)
        newdf[i,"adj_weight"] <- sum(country_filtered_KL$adj_weight)
        newdf[i,"country"] <- countries[c]
      }
      
      
      olddf <- rbind(newdf, olddf)
      
    }
  }
  
  return(olddf)
}






  







#newdf <- newdf %>% arrange(., desc(adj_weight))


# Load data and define scope of countries to explore
source('kleborate_functions.R')
library(dplyr)

df <- readr::read_csv('/Users/LSHSK42/My Drive/SeroEpi/klebsiella-seroepi/data/combined_kleborate_100722.csv')
countries <-  unique(df$Country)
#countries <- "Australia"

# Quantify run time
# a is the start time
a <- Sys.time()

weighting_country_and_study(df, countries)


# Quantify run time
# b is the end time
# b-a is the total time elapsed
b <- Sys.time()
b-a
