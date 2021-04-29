

# THIS IS THE LOOP FOR RUNNING PHYLO-CORRECTED LINEAR MODELS TO SEE IF POLLINATION SYNDROME IS A SIGNIFICANT PREDICTOR OF CV, AC1, AND DIP
# since there are many possible configurations of the phylogeny (3 species in lemon and limes, etc.) we run the models on all possible combinations (216).
# then we take the mean chi-sq and p-values and coefficients from those models. That is what we report in the text and what we plot for figures 3B and 4B. 


library(ape) 
library(Matrix) 
library(lme4) 
library(MASS) 
library(glmmTMB) 
library(coda) 
library(lattice) 
library(broom) 
library(dplyr)
library(geiger)
library(phytools)
library(tidyverse)
library(mgcv)
library(lme4)
library(car)
library(bbmle)
library(diptest)
library(plotrix)

# set your working directory
setwd("")
# you need to save these source files in the whatever working directory you have
source("lme4_phylo_setup.R")  # this is the one we're using
#source("glmmTMB_phylo_setup.R")

outdat = read.csv("outdat.csv") # this is generated in the "FAO_AB_timeseries-loop" R code

# full phylogenetic tree
new_tree = read.tree("FAO_AB_tree_update.tre")
plotTree(new_tree)
attributes(new_tree)

# make a copy of the data to fiddle with
outdat_fiddle<- outdat
outdat_fiddle$sp_latin

# make an object with the Latin names of the FAO "lemons and limes" category
lemons <- c("Citrus_aurantiifolia", "Citrus_medica","Citrus_limon")
# find where the citrus species named above are in outdat
outdat_fiddle$sp_latin %in% lemons
# replace the initial Latin name with the label of lemons so that we can cycle through all lemon/lime Latin names in the loop 
outdat_fiddle$sp_latin[outdat_fiddle$sp_latin %in% lemons] = "lemons"

# repeat for all the grouped categories
plums <- c("Prunus_spinosa", "Prunus_domestica")
outdat_fiddle$sp_latin %in% plums
outdat_fiddle$sp_latin[outdat_fiddle$sp_latin %in% plums] = "plums"

walnuts <- c("Juglans_cinerea" , "Juglans_regia", "Juglans_nigra" )
outdat_fiddle$sp_latin[outdat_fiddle$sp_latin %in% walnuts] = "walnuts"

mangos <- c( "Garcinia_x_mangostana", "Psidium_guajava", "Mangifera_indica")
outdat_fiddle$sp_latin[outdat_fiddle$sp_latin %in% mangos] = "mangos"

oranges <- c("Citrus_aurantium", "Citrus_xsinensis")
outdat_fiddle$sp_latin[outdat_fiddle$sp_latin %in% oranges] = "oranges"

persimmons <- c( "Diospyros_kaki" ,  "Diospyros_virginiana" )
outdat_fiddle$sp_latin[outdat_fiddle$sp_latin %in% persimmons] = "persimmons"

#make unique combinations of all Latin names in the grouped categories
crop_repeats <- expand.grid(lemons, plums, walnuts, mangos, oranges, persimmons)
head(crop_repeats)

# fixing the headers in the dataframe with all possible combos
crop_repeats = as.data.frame(crop_repeats, stringsAsFactors = F)
crop_repeats <- data.frame(lapply(crop_repeats, as.character), stringsAsFactors=FALSE)
names(crop_repeats) = c("lemons", "plums", "walnuts", "mangos", "oranges", "persimmons")
head(crop_repeats)

# making the results table. this is where we will save the model outputs for each of the 216 iterations of the loop
results_phylo = crop_repeats

results_phylo$ac1_chisq = rep(-99, nrow(results_phylo))
results_phylo$ac1_p = rep(-99, nrow(results_phylo))
results_phylo$ac_lik = rep(-99, nrow(results_phylo))
results_phylo$ac_AIC = rep(-99, nrow(results_phylo))

results_phylo$ac1_coef_wind= rep(-99, nrow(results_phylo))
results_phylo$ac1_025_wind= rep(-99, nrow(results_phylo))
results_phylo$ac1_975_wind= rep(-99, nrow(results_phylo))
results_phylo$ac1_coef_insect= rep(-99, nrow(results_phylo))
results_phylo$ac1_025_insect= rep(-99, nrow(results_phylo))
results_phylo$ac1_975_insect= rep(-99, nrow(results_phylo))

results_phylo$cv_chisq= rep(-99, nrow(results_phylo))
results_phylo$cv_p= rep(-99, nrow(results_phylo))
results_phylo$cv_lik = rep(-99, nrow(results_phylo))
results_phylo$cv_AIC = rep(-99, nrow(results_phylo))

results_phylo$cv_coef_insect= rep(-99, nrow(results_phylo))
results_phylo$cv_025_insect= rep(-99, nrow(results_phylo))
results_phylo$cv_975_insect= rep(-99, nrow(results_phylo))
results_phylo$cv_coef_wind= rep(-99, nrow(results_phylo))
results_phylo$cv_025_wind= rep(-99, nrow(results_phylo))
results_phylo$cv_975_wind= rep(-99, nrow(results_phylo))

#results_phylo$dip_chisq = rep(-99, nrow(results_phylo))
#results_phylo$dip_p = rep(-99, nrow(results_phylo))
#results_phylo$dip_lik = rep(-99, nrow(results_phylo))
#results_phylo$dip_AIC = rep(-99, nrow(results_phylo))

#results_phylo$dip_coef_insect= rep(-99, nrow(results_phylo))
#results_phylo$dip_025_insect= rep(-99, nrow(results_phylo))
#results_phylo$dip_975_insect= rep(-99, nrow(results_phylo))
#results_phylo$dip_coef_wind= rep(-99, nrow(results_phylo))
#results_phylo$dip_025_wind= rep(-99, nrow(results_phylo))
#results_phylo$dip_975_wind= rep(-99, nrow(results_phylo))

results_phylo$cv_ac_chisq = rep(-99, nrow(results_phylo))
results_phylo$cv_ac_p = rep(-99, nrow(results_phylo))
results_phylo$cv_ac_slope = rep(-99, nrow(results_phylo))
results_phylo$cv_ac_SE = rep(-99, nrow(results_phylo))
results_phylo$cv_ac_int = rep(-99, nrow(results_phylo))


 
#i=1 # this line is for testing the loop
for(i in 1:nrow(crop_repeats)){ # this line is for running the loop
  # tells us which row we are on
  print(paste("run", i, "out of", nrow(crop_repeats)))
  cur_outdat <- outdat_fiddle
  
  # inserts the appropriate Latin name from the dataframe of all combinations (crop_repeats) into the current iteration of outdat
  cur_outdat$sp_latin[cur_outdat$sp_latin == "lemons"] = crop_repeats$lemons[i]
  cur_outdat$sp_latin[cur_outdat$sp_latin == "plums"] = crop_repeats$plums[i]
  cur_outdat$sp_latin[cur_outdat$sp_latin == "walnuts"] = crop_repeats$walnuts[i]
  cur_outdat$sp_latin[cur_outdat$sp_latin == "mangos"] = crop_repeats$mangos[i]
  cur_outdat$sp_latin[cur_outdat$sp_latin == "persimmons"] = crop_repeats$persimmons[i]
  cur_outdat$sp_latin[cur_outdat$sp_latin == "oranges"] = crop_repeats$oranges[i]
  
  # make a matrix of the current outdat
  dat.mx = as.matrix(cur_outdat)
  # set the rownames as Latin names to check for coincidence between the data and the tree 
  rownames(dat.mx) = cur_outdat$sp_latin
  name.check(new_tree, dat.mx)
  
  # make an object to drop the extra species
  drop.sp = name.check(new_tree,dat.mx )$tree_not_data
  # drop the extra species
  pruned.tree<-drop.tip(new_tree,new_tree$tip.label[match(drop.sp, new_tree$tip.label)])
  # plot the tree
  plotTree(pruned.tree, fsize=1)
  
  # inspect tip labels
  pruned.tree$tip.label
  
  # make a Z matrix of the phylogeny 
  phyloZ_FAO = phylo.to.Z(pruned.tree)
  
  # make a variable called phylo for the models 
  cur_outdat$phylo = cur_outdat$sp_latin
  
  # run the linear mixed models with the current phylogeny
  
  # AC1
  m_ac1_phylo <- phylo_lmm(ac1~pollination_syndrome+(1|sp_latin)+(1|country),
                        data=cur_outdat,phylo=pruned.tree,
                        phyloZ=phyloZ_FAO, 
                        phylonm = "sp_latin", 
                        control=lmerControl(check.nobs.vs.nlev="ignore",
                                            check.nobs.vs.nRE="ignore"), REML = TRUE)
  
  # store p-value and chi-sq
  results_phylo$ac1_p[i]<- Anova(m_ac1_phylo)[1,3]
  results_phylo$ac1_chisq[i]<- Anova(m_ac1_phylo)[1,1]
  
  # run the "means" parameterization to calculate and store the coefficients
  m_ac1_phylo_mn <- phylo_lmm(ac1~0+pollination_syndrome+(1|sp_latin)+(1|country),
                           data=cur_outdat, phylo=pruned.tree,
                           phyloZ=phyloZ_FAO, 
                           phylonm = "sp_latin", 
                           control=lmerControl(check.nobs.vs.nlev="ignore",
                                               check.nobs.vs.nRE="ignore"), REML = TRUE)

  results_phylo$ac1_coef_insect[i]<- fixef(m_ac1_phylo_mn)[1]
  results_phylo$ac1_coef_wind[i]<- fixef(m_ac1_phylo_mn)[2]
  
  # calculate and store the confidence intervals
  temp=confint(m_ac1_phylo_mn)  
  temp[4,2]
  
  results_phylo$ac1_025_insect[i]<- temp[4,1]
  results_phylo$ac1_975_insect[i]<- temp[4,2]
  
  results_phylo$ac1_025_wind[i]<- temp[5,1]
  results_phylo$ac1_975_wind[i]<- temp[5,2]

  # calculate and store log likelihood and AIC
  results_phylo$ac_lik[i] = logLik(m_ac1_phylo)[1]
  results_phylo$ac_AIC[i] = AIC(m_ac1_phylo)[1]


  # CV
  m_cv_phylo <- phylo_lmm(cv~pollination_syndrome+(1|sp_latin)+(1|country),
                        data=cur_outdat,phylo=pruned.tree,
                        phyloZ=phyloZ_FAO, 
                        phylonm = "sp_latin", 
                        control=lmerControl(check.nobs.vs.nlev="ignore",
                                            check.nobs.vs.nRE="ignore"), REML = TRUE)
  
  # store p-value and chi-sq
  results_phylo$cv_p[i]<- Anova(m_cv_phylo)[1,3]
  results_phylo$cv_chisq[i]<- Anova(m_cv_phylo)[1,1]
  
  # run the "means" parameterization to calculate and store the coefficients
  m_cv_phylo_mn <- phylo_lmm(cv~0+pollination_syndrome+(1|sp_latin)+(1|country),
                          data=cur_outdat,phylo=pruned.tree,
                          phyloZ=phyloZ_FAO, 
                          phylonm = "sp_latin", 
                          control=lmerControl(check.nobs.vs.nlev="ignore",
                                              check.nobs.vs.nRE="ignore"), REML = TRUE)
  
  results_phylo$cv_coef_insect[i]<- fixef(m_cv_phylo_mn)[1]
  results_phylo$cv_coef_wind[i]<- fixef(m_cv_phylo_mn)[2]
  
  # caclulate and store confidence intervals
  temp_cv=confint(m_cv_phylo_mn)  
  
  results_phylo$cv_025_insect[i]<- temp_cv[4,1]
  results_phylo$cv_975_insect[i]<- temp_cv[4,2]
  
  results_phylo$cv_025_wind[i]<- temp_cv[5,1]
  results_phylo$cv_975_wind[i]<- temp_cv[5,2]
  
  
  # calculate and store log lik and AIC 

  results_phylo$cv_lik[i] = logLik(m_cv_phylo)[1]
  results_phylo$cv_AIC[i] = AIC(m_cv_phylo)[1]

  # DIP  

  m_d_phylo <- phylo_lmm(d~pollination_syndrome+(1|sp_latin)+(1|country),
                          data=cur_outdat,phylo=pruned.tree,
                          phyloZ=phyloZ_FAO, 
                          phylonm = "sp_latin", 
                          control=lmerControl(check.nobs.vs.nlev="ignore",
                                              check.nobs.vs.nRE="ignore"), REML = TRUE)
  results_phylo$dip_chisq[i] = Anova(m_d_phylo)[1,1]
  results_phylo$dip_p[i] = Anova(m_d_phylo)[1,3]

  
  # CV and AC1 (this is for figure 6)
  m_cv_ac = phylo_lmm(ac1 ~ cv + (1|sp_latin)+(1|country),
                      data=cur_outdat,phylo=pruned.tree,
                      phyloZ=phyloZ_FAO, 
                      phylonm = "sp_latin", 
                      control=lmerControl(check.nobs.vs.nlev="ignore",
                                          check.nobs.vs.nRE="ignore"), REML = TRUE)
  
  results_phylo$cv_ac_chisq[i] = Anova(m_cv_ac)[1,1]  
  results_phylo$cv_ac_p[i]= Anova(m_cv_ac)[1,3]

  
  
  sum<- summary(m_cv_ac)
  results_phylo$cv_ac_int[i] =   sum$coefficients[1,1]
  results_phylo$cv_ac_slope[i] =   sum$coefficients[2,1]
  results_phylo$cv_ac_SE[i] =   sum$coefficients[2,2]
  

  
}
head(results_phylo)
#summary(m1_phylo)


# calculating all the mean values and inspecting their range
# mean values chi-sq, p, and coefficient values are  reported in main text and plotted below

# AC1
mean(results_phylo$ac1_chisq)
range(results_phylo$ac1_chisq)

mean(results_phylo$ac1_p)
range(results_phylo$ac1_p)

### insect
mean(results_phylo$ac1_coef_insect) #plotted for graph 3b
std.error(results_phylo$ac1_coef_insect)

mean(results_phylo$ac1_025_insect) #error bars for graphs
mean(results_phylo$ac1_975_insect)

### wind
mean(results_phylo$ac1_coef_wind) #plotted for graph 3b
std.error(results_phylo$ac1_coef_wind)

mean(results_phylo$ac1_025_wind) #error bars for graphs
mean(results_phylo$ac1_975_wind)


# CV
mean(results_phylo$cv_chisq)
range(results_phylo$cv_chisq)

mean(results_phylo$cv_p)
range(results_phylo$cv_p)

### insect
mean(results_phylo$cv_coef_insect) #plot for fig 4b
std.error(results_phylo$cv_coef_insect)

mean(results_phylo$cv_025_insect) #error bars for graphs
mean(results_phylo$cv_975_insect)

### wind
mean(results_phylo$cv_coef_wind) #plot for fig 4b
std.error(results_phylo$cv_coef_wind)

mean(results_phylo$cv_025_wind) #error bars for graphs
mean(results_phylo$cv_975_wind)


#####supplement data######
# DIP

# mean(results_phylo$dip_chisq)
# range(results_phylo$dip_chisq)
# 
# mean(results_phylo$dip_p)
# range(results_phylo$dip_p)
# 
# ### insect
# mean(results_phylo$dip_coef_insect) #plot for supplemental
# std.error(results_phylo$dip_coef_insect)
# 
# mean(results_phylo$dip_025_insect) #error bars for graphs
# mean(results_phylo$dip_975_insect)
# 
# ### wind
# mean(results_phylo$dip_coef_wind) #plot for supplemental
# std.error(results_phylo$dip_coef_wind)
# 
# mean(results_phylo$dip_025_wind) #error bars for graphs
# mean(results_phylo$dip_975_wind)

# general syntax
plotCI(VALUES_FOR_THE_MEANS_IN_ONE_VECTOR, li = VALUES_FOR_LOWER_LIMIT_IN_ONE_VECTOR, ui = VALUES_FOR_UPPER_LIMIT_IN_ONE_VECTOR)

# for figure 3b and 4b
pdf("AB_GAM_long.pdf", height = 5, width = 7)
par(xpd=F)
par(mfrow = c(1,2))

# creating the vectors in advance for AC1
means = c(mean(results_phylo$ac1_coef_insect), mean(results_phylo$ac1_coef_wind))
lis = c(mean(results_phylo$ac1_025_insect), mean(results_phylo$ac1_025_wind))
uis = c(mean(results_phylo$ac1_975_insect), mean(results_phylo$ac1_975_wind))

plotCI(means, li = lis, ui = uis, xlab = "", xaxt = "n", ylab = "", pch = 21, pt.bg = c("black", "indian red"), cex = 1.4, xlim = c(0.5, 2.5))
points(c(0,3), c(0,0), type = "l", lty = "dotted")
mtext(side = 1, cex = 1.25, line = 0.5, at = 1:2, text = c("Insect", "Wind"))
mtext(side = 1, cex = 1.4, line = 2, "Pollination syndrome")
mtext(side = 2, cex = 1.4, line =2.5, "Lag-1 autocorrelation")

# creating the vectors in advance for CV
means_CV = c(mean(results_phylo$cv_coef_insect), mean(results_phylo$cv_coef_wind))
lis_CV = c(mean(results_phylo$cv_025_insect), mean(results_phylo$cv_025_wind))
uis_CV = c(mean(results_phylo$cv_975_insect), mean(results_phylo$cv_975_wind))

plotCI(means_CV, li = lis_CV, ui = uis_CV, xlab = "", xaxt = "n", ylab = "", pch = 21, pt.bg = c("black", "indian red"), cex = 1.4, xlim = c(0.5, 2.5), ylim = c(0, 0.4))
points(c(0,3), c(0,0), type = "l", lty = "dotted")
mtext(side = 1, cex = 1.25, line = 0.5, at = 1:2, text = c("Insect", "Wind"))
mtext(side = 1, cex = 1.4, line = 2, "Pollination syndrome")
mtext(side = 2, cex = 1.4, line =2.5, "Coefficient of variation")

# Supplement: creating the vectors in advance for dip statistic
#means_dip = c(mean(results_phylo$dip_coef_insect), mean(results_phylo$dip_coef_wind))
#lis_dip = c(mean(results_phylo$dip_025_insect), mean(results_phylo$dip_025_wind))
#uis_dip = c(mean(results_phylo$dip_975_insect), mean(results_phylo$dip_975_wind))

#plotCI(means_dip, li = lis_dip, ui = uis_dip, xlab = "", xaxt = "n", ylab = "", pch = 21, pt.bg = c("black", "indian red"), cex = 1.4, xlim = c(0.5, 2.5))
#points(c(0,3), c(0,0), type = "l", lty = "dotted")
#mtext(side = 1, cex = 1.25, line = 0.5, at = 1:2, text = c("Insect", "Wind"))
#mtext(side = 1, cex = 1.4, line = 2, "Pollination syndrome")
#mtext(side = 2, cex = 1.4, line =2.5, "Dip statistic")

dev.off()

# AC ~ CV #add this as a trendline to figure 6 using the mean slope and intercept
mean(results_phylo$cv_ac_chisq) 
range(results_phylo$cv_ac_chisq)

mean(results_phylo$cv_ac_p)
range(results_phylo$cv_ac_p)

mean(results_phylo$cv_ac_slope) 
range(results_phylo$cv_ac_slope)

mean(results_phylo$cv_ac_SE)
range(results_phylo$cv_ac_SE)



