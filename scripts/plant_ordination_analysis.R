# code for doing ordination analysis on Portal plant communities (summer and winter annuals)
# I followed Sarah's code from Supp et al 2012 https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/12-0370.1
# Sarah did a pcca on the plant community (summer and winter annuals)

# this is also a useful example: https://rgriff23.github.io/2017/05/23/mosquito-community-ecology-in-vegan.html


library(dplyr)
source('scripts/prepare_plant_data.R')

# data ---------------------------------------------------------------------------------------------------------
winterannuals = make_plant_table(selected_plots=1:24,
                                 plant_community='Annuals',
                                 summer_winter='winter',
                                 threshold = 0.33)

#write.csv(winterannuals,'WinterAnnuals.csv',row.names=F)

#dat = read.csv('WinterAnnuals.csv')
dat = winterannuals
# filter: 3 years before switch to present
dat.filtered = dplyr::filter(dat,year >= 2012)
dat.filtered = dat

# add treatment
dat.filtered$treatment = rep('')
dat.filtered$treatment[dat.filtered$plot %in% c(4,11,14,17)] <- 'control'
dat.filtered$treatment[dat.filtered$plot %in% c(3,15,19,20,21)] <- 'exclosure'
dat.filtered$treatment[dat.filtered$plot %in% c(10,16,23)] <- 'removal'
dat.filtered$treatment[dat.filtered$plot %in% c(2,8,22)] <- 'control-exclosure'
dat.filtered$treatment[dat.filtered$plot %in% c(6,13,18)] <- 'exclosure-control'
dat.filtered$treatment[dat.filtered$plot %in% c(5,7,24)] <- 'removal-control'
dat.filtered$treatment[dat.filtered$plot %in% c(1,12,9)] <- 'control-removal'

dat.filtered$treat_before = rep('')
dat.filtered$treat_before[dat.filtered$plot %in% c(1,2,4,8,9,11,12,14,17,22)] <- 'control'
dat.filtered$treat_before[dat.filtered$plot %in% c(3,6,13,15,18,19,20,21)] <- 'exclosure'
dat.filtered$treat_before[dat.filtered$plot %in% c(5,7,10,16,23,24)] <- 'removal'


dat.filtered$treat_after = rep('')
dat.filtered$treat_after[dat.filtered$plot %in% c(4,5,6,7,11,13,14,17,18,24)] <- 'control'
dat.filtered$treat_after[dat.filtered$plot %in% c(2,3,8,15,19,20,21,22)] <- 'exclosure'
dat.filtered$treat_after[dat.filtered$plot %in% c(1,9,10,12,16,23)] <- 'removal'

dat.winter = dat.filtered
dat.winter = filter(dat.filtered,year>2010,year<2016)
dat.winter = filter(dat.filtered,year>2015)
# pcca model ------------------------------------------------------------------------------------------------------
# sqrt transofrm abundance data to account for huge differences in abundance year to year
dat.species = sqrt(as.matrix(dat.winter[,!names(dat.filtered) %in% c('year','season','plot','treatment','treat_before','treat_after')]))

# rest of variables
dat.year = as.factor(dat.winter$year)
dat.trt = as.factor(dat.winter$treat_before)
dat.plot = as.factor(dat.winter$plot)

colorpalette = RColorBrewer::brewer.pal(3,"Set1")

# pcca
win.pcca <- cca(dat.species ~ dat.trt + Condition(dat.year))
# check variance inflation factor
vif.cca(win.pcca) # all good
# plot sites
plot(win.pcca, display=c("sites", "bp"), type="n", main="Sites", scaling="sites")
text(win.pcca, display="sites", col=site.cols, scaling="sites")
#text(win.pcca, display="bp", col="red")
legend("topright", legend = levels(dat.trt), bty = "n",
       col = colorpalette, pch = 21, pt.bg = colorpalette)
# plot species
plot(win.pcca, display=c("species", "bp"), type="n", ylab="", main="Species", scaling="species")
text(win.pcca, display="species", col="black", scaling="species")
text(win.pcca, display="bp", col="red")
# proportion of variance explained
win.pcca$CCA$tot.chi/win.pcca$tot.chi
# significance test
anova(win.pcca)
permutest(win.pcca,permutations=500) # should be similar to anova on pcca
anova(win.pcca,strata=dat.year) # more conservative test

# =================================================================================================
# other analyses

# just a correspondence analysis
win.ca <- cca(dat.species)
barplot(win.ca$CA$eig/win.ca$tot.chi, names.arg = 1:win.ca$CA$rank, cex.names = 0.5, ylab="Proportion of variance explained", xlab="CA axis")
site.cols = colorpalette[as.numeric(dat.trt)]
plot(win.ca, display="sites", type="n", main="Sites", scaling="sites")
text(win.ca, display="sites", col=site.cols, scaling="sites")
legend("topright", legend = levels(dat.trt), bty = "n",
       col = colorpalette, pch = 21, pt.bg = colorpalette)


# cca
win.cca <- cca(dat.species ~ dat.trt)
# check variance inflation factor (>10 is bad)
vif.cca(win.cca) # all good
barplot(win.cca$CA$eig/win.cca$tot.chi, names.arg = 1:win.cca$CA$rank, cex.names = 0.5, ylab="Proportion of variance explained", xlab="CCA axis")
# plot sites
plot(win.cca, display=c("sites", "bp"), type="n", main="Sites", scaling="sites")
text(win.cca, display="sites", col=site.cols, scaling="sites")
text(win.cca, display="bp", col="red")
legend("topright", legend = levels(dat.trt), bty = "n",
       col = colorpalette, pch = 21, pt.bg = colorpalette)
# plot species
plot(win.cca, display=c("species", "bp"), type="n", ylab="", main="Species", scaling="species")
text(win.cca, display="species", col="black", scaling="species")
text(win.cca, display="bp", col="red")
# proportion of variance explained
win.cca$CCA$tot.chi/win.cca$tot.chi
# significance test for individual predictors (type 3 test) (Not Applicable here because I only have one predictor)
anova(win.cca, by="margin")
# significance test for entire model
anova(win.cca)


### ADONIS - another test to look for compositional differences, with similar results to above.
# winter
win.spp.canb = vegdist(dat.species, method = "canb")
win.canb = adonis(win.spp.canb ~ dat.trt, permutation=1000)
win.canb
