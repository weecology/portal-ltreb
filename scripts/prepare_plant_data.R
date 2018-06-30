#' @title get winter or summer plant data: crosstab
#'
#' @param selected_plots plot numbers: 1-24
#' @param plant_community which plant community (e.g. annuals or perennials; input 'type' in plant_abundance function)
#' @param summer_winter 'summer' or 'winter' census; or "both" for perennials
#' @param threshold value used to remove "rare" or "transient" species from the table. For example, a threshold of .1 removes species that are present in less than 10% of sampling years
#'
#' @return table of plant counts by species
#'
make_plant_table = function(selected_plots,plant_community,summer_winter,threshold=0.1) {
  plant_data = portalr::plant_abundance('repo',level='Plot',type=plant_community,
                               correct_sp=T,unknowns=F,plots=selected_plots,
                               shape='flat')
  
  # remove data before 1983 to avoid having to adjust by quadrat area per plot
  season_data = dplyr::filter(plant_data,year>1982)
  # filter by season
  if (tolower(summer_winter) %in% c('summer','winter')) {
    season_data = dplyr::filter(plant_data,season==summer_winter)
  }
 
  
  # find species that occurred in >10% of years
  transients = find_transient_species(season_data,threshold)
  
  # filter based on selected_plots, remove transient species
  seasonplants = filter(season_data,plot %in% selected_plots, !(species %in% transients$species))
  
  # sum by year (effort doesn't vary by year)
  seasontotal = aggregate(seasonplants$abundance,by=list(year=seasonplants$year,season=seasonplants$season,species=seasonplants$species,plot=seasonplants$plot),FUN=sum)
  seasontable = portalr:::make_crosstab(seasontotal,variable_name='x')
  seasontable[is.na(seasontable)] <- 0
  

  return(seasontable)
}

#' @title find transient species
#'
#' @description identifies species found in less than a given percentage of sampling events
#'
#' @param data_table table of data including species and year
#' @param threshold value for determining transient. Default 10%
#' 
#' @return vector of species names that are rare
#'
find_transient_species = function(data_table,threshold=0.1) {
  ntimesteps = dplyr::select(data_table,year) %>% unique() %>% nrow()
  sp_time = dplyr::select(data_table,year,species) %>% unique()
  yrs_present = dplyr::count(sp_time,species)
  transient_names = dplyr::filter(yrs_present,n <= threshold*ntimesteps) %>% select(species)
  return(transient_names)
}

# -----------------------------------------------------------------------------------------------------------------

