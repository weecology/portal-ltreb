# NMDS for Plants
# EKB
# June 2018


# LIBRARIES #

library(dplyr)
library(reshape2)
library(vegan)
library(ggplot2)
#devtools::install_github("dgrtwo/gganimate")
library(gganimate)

#=============================================================================#
# DATA PREP #

plants <- read.csv("WinterAnnuals.csv", header = TRUE, stringsAsFactors = FALSE)
plants <- filter(plants, year > 1987)

# add a column for data before or after the switch

plants$switch_time <- plants$treatment
for (i in 1:length(plants$plot)){
  if (plants$year[i] > 1987 & plants$year[i] <= 2015){
    plants$switch_time[i] = 'Before'
  } else {
    plants$switch_time[i] = 'After' 
  }
}

# add columns for plot treatment and switch type

plants$treatment <- plants$plot

for (i in 1:length(plants$plot)){
  
  # treatment type (dependent on year)
  if (plants$year[i] <= 2015){
    
    if (plants$plot[i] %in% c(1, 2, 4, 8, 9, 11, 12, 14, 17, 22)){
      plants$treatment[i] = 'Control'
    } else if (plants$plot[i] %in% c(3, 6, 13, 15, 18, 19, 20, 21)){
      plants$treatment[i] = 'Dipo_Exclosure' 
    } else {
      plants$treatment[i] = 'Removal'
    }
    
  } else {
    
    if (plants$plot[i] %in% c(4, 5, 6, 7, 11, 13, 14, 17, 18, 24)){
      plants$treatment[i] = 'Control'
    } else if (plants$plot[i] %in% c(2, 3, 8, 15, 19, 20, 21, 22)){
      plants$treatment[i] = 'Dipo_Exclosure' 
    } else {
      plants$treatment[i] = 'Removal'
    }
    
  }
}  

plants$switch_type <- plants$plot

for (i in 1:length(plants$plot)){
  
  # switch type
  if(plants$plot[i] %in% c(4, 11, 14, 17)){
    plants$switch_type[i] = 'CC'
  } else if (plants$plot[i] %in% c(2, 8, 22)){
    plants$switch_type[i] = 'CE'
  } else if (plants$plot[i] %in% c(1, 9, 12)){
    plants$switch_type[i] = 'CR'
  } else if (plants$plot[i] %in% c(3, 15, 19, 20, 21)){
    plants$switch_type[i] = 'EE'
  } else if (plants$plot[i] %in% c(6, 18, 13)){
    plants$switch_type[i] = 'EC'
  } else if (plants$plot[i] %in% c(5, 7, 24)){
    plants$switch_type[i] = 'RC'
  } else {
    plants$switch_type[i] = 'RR'
  }
  
}


# make a matrix with relative abundances

## new data frame
species <- as.list(colnames(plants)[4:43])
rel_abund <- as.data.frame(matrix(0, ncol = length(species), nrow = length(plants$plot)))
colnames(rel_abund) <- c(species)

## get rel abundances
for (i in 1:length(plants$plot)){
  temp <- plants[i, 4:43]
  total_abund <- rowSums(temp)
  #plants$total_abund[i] <- total_abund
  temp2 <- apply(temp, 2, function(x) x/total_abund)
  rel_abund[i,] <- temp2
}

## bind columns for new dataframe
row_data <- plants[,c(1:3,44:46)]
rel_abund_matrix <- bind_cols(row_data, rel_abund)

#=============================================================================#
# RUN NMDS #

years <- unique(plants$year) # list of years

# dataframe for storing NMDS values
NMDS_df <- data.frame(group = character(0), 
                      MDS1 = numeric(0), 
                      MDS2 = numeric(0), 
                      year = integer(0))

# dataframe for storing means
NMDS_mean_df <- data.frame(group = character(0), 
                           MDS1 = numeric(0), 
                           MDS2 = numeric(0), 
                           year = integer(0))

# dataframe for storing ellipse elements
df_ell_all <- data.frame(group = character(0), 
                         MDS1 = numeric(0), 
                         MDS2 = numeric(0), 
                         year = integer(0))

### STOPPED HERE - 6/28/2018 ###
# To do from here:
#   - try to figure out how to run *same* NMDS through time 
#   - and probably amongst different groups
#===============================================================================#

for (i in 1:length(years)){                
  
  # create species x site communities, by year
  #   - run NMDS on each
  #   - get the mean NMDS oridation location by group type
  #   - put this data into a dataframe for plotting
  
  yr <- years[i]                            # select each year
  
  data <- rel_abund_matrix %>% filter(year == yr)         # filter rodent data by year
  matrix <- data[,-c(1,2,4,5)]                         # remove year column
  matrix <- tibble::column_to_rownames(as.data.frame(matrix), "plot")
  
  # run NMDS using vegan
  dist <- metaMDS(matrix, trymax = 100)
  
  # make dataframe with with MDS values and grouping variable
  NMDS <- data.frame(MDS1 = dist$points[,1], 
                     MDS2 = dist$points[,2], 
                     group = rel_abund_matrix$switch_type[yr], #changed this, should work
                     year = yr)
  
  # add to NMDS dataframe
  NMDS_df <- rbind(NMDS_df, NMDS)
  
  # get mean point for each group
  NMDS.mean <- aggregate(NMDS[,1:2], list(group = occupied$treatment), mean)
  NMDS.mean$year <- year # add year to dataframe
  
  # add to NMDS mean dataframe
  NMDS_mean_df <- rbind(NMDS_mean_df, NMDS.mean)
  
  # get ellipse elements
  df_ell <- data.frame()
  for(g in levels(NMDS$group)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g, ], 
                                                     vegan:::veganCovEllipse(cov.wt(cbind(MDS1, MDS2), wt=rep(1/length(MDS1),
                                                                                                              length(MDS1)))$cov, center=c(mean(MDS1), mean(MDS2))))), group = g))
  }
  df_ell$year <- year
  
  df_ell_all <- rbind(df_ell_all, df_ell)
  
}

#===================================================
# PLOT RESULTS

# use ggplot2 to plot the data
plot <- ggplot(data = NMDS_df, aes(MDS1, MDS2, frame = year)) + 
  geom_point(aes(color = group)) +
  geom_path(data = NMDS_mean_df, aes(MDS1, MDS2, color = group)) +
  #geom_path(data = df_ell_all, aes(x = MDS1, y = MDS2, colour = group), size = 1) +
  theme_bw()
gganimate(plot, "All_species.gif")