library(dplyr)
library(vegan)
library(portalr)
library(ggplot2)
source('scripts/prepare_plant_data.R')

# data ---------------------------------------------------------------------------------------------------------
data=load_data(clean=FALSE)

longterm_plots = c(3,4,10,11,14,15,16,17,19,20,21,23)
flipped_plots = c(1,2,5,6,7,8,9,12,13,18,20,24)

plots = data$plots_table %>%
  select(year,plot,treatment) %>%
  group_by(year,plot) %>%
  summarise(treatment=first(treatment))

plants = make_plant_table(selected_plots=1:24,
                          plant_community='Annuals',
                          summer_winter='winter',
                          threshold = 0.1) %>%
  select(-season) %>%
  mutate_at(vars(-one_of("year", "plot")), sqrt) %>%
  left_join(plots) %>%
  mutate(treatment = as.character(treatment),plotflip = ifelse(year<2016," Pre Flip","Post Flip"),
         plottype = ifelse(plot %in% longterm_plots, "longterm", "flipped")) %>%
  filter(year>2011,year!=2013)

###########All plots##########################
pcca <- cca(select(plants, -c(year,plot,treatment,plotflip,plottype)),
            select(plants, c(treatment)),
            select(plants, c(year)))

# plot with year, plot, treatment
results=data.frame(scores(pcca)$sites,
                   year=plants$year,
                   plot=plants$plot,
                   treatment=plants$treatment,
                   plotflip = plants$plotflip,
                   plottype = plants$plottype)

ggplot(results,aes(CCA1, CCA2)) + 
  stat_ellipse(aes(color = treatment)) +
  geom_text(aes(label = plot, color = as.factor(year))) +
  facet_grid(plottype~plotflip) +
  coord_cartesian(xlim = c(-5, 4), ylim = c(-7, 7)) +
  theme(legend.title=element_blank()) +
  scale_color_discrete(direction=-1)

#######################
# Other figures, analyzing flipped and longterm plots separately
longterm = make_plant_table(selected_plots=longterm_plots,
                            plant_community='Annuals',
                            summer_winter='winter',
                            threshold = 0.1) %>%
  select(-season) %>%
  mutate_at(vars(-one_of("year", "plot")), sqrt) %>%
  left_join(plots) %>%
  mutate(treatment = as.character(treatment),
         plotflip = ifelse(year<2016,"Longterm 2012 - 2015","Longterm 2016 - 2018")) %>%
  filter(year>2011,year!=2013)


flipped = make_plant_table(selected_plots=flipped_plots,
                           plant_community='Annuals',
                           summer_winter='winter',
                           threshold = 0.1) %>%
  select(-season) %>%
  mutate_at(vars(-one_of("year", "plot")), sqrt) %>%
  left_join(plots) %>%
  mutate(treatment = as.character(treatment),plotflip = ifelse(year<2016," Pre Flip","Post Flip")) %>%
  filter(year>2011,year!=2013)

###########Flipped plots##########################
preflip = filter(flipped,year<2016)
preflip.pcca <- cca(select(preflip, -c(year,plot,treatment,plotflip)),
                     select(preflip, c(treatment)),
                     select(preflip, c(year,plot)))
postflip = filter(flipped,year>2015)
postflip.pcca <- cca(select(postflip, -c(year,plot,treatment,plotflip)),
                    select(postflip, c(treatment)),
                    select(postflip, c(year,plot)))


# plot with year, plot, treatment
flip_results=data.frame(rbind(scores(preflip.pcca)$sites,scores(postflip.pcca)$sites),
                           year=flipped$year,
                           plot=flipped$plot,
                           treatment=flipped$treatment,
                           plotflip = flipped$plotflip)

ggplot(flip_results,aes(CCA1, CCA2)) + 
  stat_ellipse(aes(color = treatment)) +
  geom_text(aes(label = plot, color = as.factor(year))) +
  facet_grid(.~plotflip) +
  theme(legend.title=element_blank()) +
  scale_color_discrete(direction=-1)


###########Longterm plots##########################
longterm1 = filter(longterm,year<2016)
longterm.pcca1 <- cca(select(longterm1, -c(year,plot,treatment,plotflip)),
                    select(longterm1, c(treatment)),
                    select(longterm1, c(year,plot)))

longterm2 = filter(longterm,year>2015)
longterm.pcca2 <- cca(select(longterm2, -c(year,plot,treatment,plotflip)),
                     select(longterm2, c(treatment)),
                     select(longterm2, c(year,plot)))

# plot with year, plot, treatment
longterm_results=data.frame(rbind(scores(longterm.pcca1)$sites,scores(longterm.pcca2)$sites),
                        year=longterm$year,
                        plot=longterm$plot,
                        treatment=longterm$treatment,
                        plotflip = longterm$plotflip)

ggplot(longterm_results,aes(CCA1, CCA2)) + 
  stat_ellipse(aes(color = treatment)) +
  geom_text(aes(label = plot, color = as.factor(year))) +
  facet_grid(.~plotflip) +
  theme(legend.title=element_blank()) +
  scale_color_discrete(direction=-1)




  