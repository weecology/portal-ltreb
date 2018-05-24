# NMDS for plot composition through time
# April 2017


# LIBRARIES AND DATA #
library(RCurl)
library(dplyr)
library(reshape2)
library(vegan)
library(ggplot2)
#devtools::install_github("dgrtwo/gganimate")
library(gganimate)

# get data
rodents <- getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv") 
rdat <- read.csv(text = rodents, header = TRUE, na.strings = c(""))

#=====================================================
# DATA PREP #

### GROUPING DATA ###

plot_type <- data.frame(plot = c(1,2,3,4,6,8,9,11,12,13,14,15,17,18,19,20,21,22), 
                        treatment = NA)

for (i in 1:18){
  if (plot_type$plot[i] %in% c(1, 2, 4, 8, 9, 11, 12, 14, 17, 22)){
    plot_type$treatment[i] = 'Control'
  } else {
    plot_type$treatment[i] = 'Dipo_Exclosure' 
  }
}

### RODENT DATA ###

rspec <- c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO')

rdat <- select(rdat, period, year, plot, species) %>%
  filter(period > 0,
         year > 1987 & year < 2015,
         plot %in% c(1,2,3,4,6,8,9,11,12,13,14,15,17,18,19,20,21,22),
         species %in% rspec) %>% 
  group_by(year, plot, species) %>% 
  summarise(count = n())

#=====================================================
# RUN NMDS #

years <- unique(rdat$year) # list of years

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

for (i in 1:length(years)){
  
  # create species x site communities, by year
  #   - run NMDS on each
  #   - get the mean NMDS oridation location by group type
  #   - put this data into a dataframe for plotting
  
  yr <- years[i]                            # select each year
  
  data <- rdat %>% filter(year == yr)         # filter rodent data by year
  matrix <- data[,-1]                         # remove year column
  matrix <- dcast(matrix,                     # make presence-absence matrix
                  formula = plot ~ species, 
                  fun.aggregate = length, 
                  value.var = "count")
  matrix <- tibble::column_to_rownames(as.data.frame(matrix), "plot")

  plot <- as.numeric(row.names(matrix))
  plot <- as.data.frame(plot)
  occupied <- semi_join(plot_type, plot)
  
  # run NMDS using vegan
  dist <- metaMDS(matrix, trymax = 100)
  
  # make dataframe with with MDS values and grouping variable
  NMDS <- data.frame(MDS1 = dist$points[,1], 
                     MDS2 = dist$points[,2], 
                     group = occupied$treatment)
  NMDS$year <- year
  
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
