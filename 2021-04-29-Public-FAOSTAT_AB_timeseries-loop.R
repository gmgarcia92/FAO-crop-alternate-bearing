# analysis of FAO data for alternate-bearing crops
# see data at http://www.fao.org/faostat/en/#data/QC


# This code is for compiling the dataset used in the phylogenetically corrected linear models in the "FAO_AB_phylogeny" R file


library(tidyverse)
library(mgcv)
library(lme4)
library(car)
library(emmeans)
library(plotrix)
library(moments)
library(bbmle)
library(diptest)

ntop = 10


setwd("/Users/gabrielagarcia/Documents/AB_Pollinator_FAO")
alldat = read.csv("GG-FAOSTAT_AB_crops.csv")
traits = read.csv("Polllination and Alternate Bearing in crops - GG_ABcrops.csv")

unique(alldat$Item)

summary(alldat)
alldat0 = alldat[alldat$Area != "China" & alldat$Area != "USSR" & 
                   alldat$Area != "Belgium-Luxembourg" & alldat$Area != "Belgium" & alldat$Area != "Luxembourg" & 
                   alldat$Area != "Serbia and Montenegro" & alldat$Area != "Serbia" & alldat$Area != "Montenegro" & 
                   alldat$Area != "Sudan (former)" &  alldat$Area != "Sudan" & alldat$Area != "South Sudan" &
                   alldat$Item != "Date",]
#removed China -> seperated into "China, mainland" and "China, Taiwan Province of"
#removed USSR -> all countries included in USSR began reporting in 1992 (Ukraine, Armenia, Geogria, Azerbijan, Belarus, Kazakhstan, Uzebekistan, Kyrgyzstan, Lativa, Turkmenestan, Lithuania, Tajikistan, Estonia)
#all other countries removed because their split resulted in less than 20yrs of data
#removed "Date" from anaylsis since they are largely hand pollinated and don't fit well into insect or wind pollinated categories

alldat0$Area[alldat0$Area =="Ethiopia PDR"] = "Ethiopia"
#included Ethiopia PDR as part of Ethopia

summary(traits)
unique(traits$crop)



(countries = unique(alldat0$Area)) # countries with at least some data
(crops = unique(alldat0$Item)) # focal alternate bearing crops
(responses = unique(alldat0$Element)) # two metrics - yield (per acre) and production (total amount per country)


ii = 1
output = tibble(country = character(), crop = character(), ac1 = double(), cv = double(), Ha = double(), d = double(), dip.p = double())
pdf("AB_gam_long series.pdf", height = 10, width = 7.5) #to generate supplement 2
par(mfrow = c(5,4))

# below loop results in long time series detrended using GAMs and autocorrelation and CV for each crop and country combination
# the outdat generated below is the input for the phylogenetically-corrected LMMs in the "FAO_AB_phylogeny" R file
# to generate supplement 2 and alternative outdats: sub line 73 for 74 (LOESS detrending), alternatively sub lines 78-80 for 73-77 (differencing); for short time series sub 1994 for 1960 in lines 63, 70 and 87
for(j in 1:length(crops)){
  my_crop = crops[j] 
  my_resp = "Production"
  usedatP = alldat0[alldat0$Item == my_crop & alldat0$Element == my_resp,]
  production = with(usedatP[usedatP$Year > 1960,], tapply(Value, Area, mean, na.rm = T))
  production = production[is.na(production) == F]
  countries_tmp = names(production[order(production, decreasing = T)][1:ntop])
  
  my_resp2 = "Yield"
  for(i in 1:length(countries_tmp)){
      my_country = countries_tmp[i] 
      usedat = alldat0[alldat0$Area == my_country & alldat0$Item == my_crop & alldat0$Element == my_resp2 & alldat0$Year > 1960,]
    
        if(sum(!is.na(usedat$Value)) > 20){
      m0 = gam(Value ~ s(Year), data = usedat) # smoothed trend of yield over time
      #m0 = loess(Value ~ Year, data = usedat) # smoothed trend of yield over time
      usedat$dev = usedat$Value - predict(m0) # calculate annual deviations from smoothed trend
      plot(usedat$Year, usedat$Value, pch = 19, cex = 0.75, xlab = "Year", ylab = "Yield", type = "o", main = paste(substr(my_crop, start = 1, stop =10), substr(my_country, start = 1, stop = 13)))
      points(usedat$Year, predict(m0), type = "l", col = "red", lwd = 2)
      #mylen = dim(usedat)[1]
      #usedat$dev = NA 
      #usedat$dev[1:(mylen-1)] = usedat$Value[2:mylen] - usedat$Value[1:(mylen-1)] # calculate annual deviations with differencing
      plot(usedat$Year, usedat$dev, pch = 19, cex = 0.75, type = "o", xlab = "Year",ylab = "Yield (detrended)", main = "")
      output[ii,1] = my_country
      output[ii,2] = my_crop
      output[ii,3] = 
        acf(usedat$dev, plot = T, main = "", na.action = na.pass)$acf[2] # lag-1 AC
      output[ii,4] = sd(usedat$dev, na.rm = T)/mean(usedat$Value, na.rm = T)
      usedat0 = alldat0[alldat0$Area == my_country & alldat0$Item == my_crop & alldat0$Element == "Area harvested" & alldat0$Year > 1960,]
      output[ii,5] = mean(usedat0$Value, na.rm = T)
      output[ii,6] = dip.test(usedat$dev)$statistic[1]
      output[ii,7] = dip.test(usedat$dev)$p.value[1]
      ii = ii+1   
      hist(scale(usedat$dev), xlab = "Scaled deviations", main = "")
      mtext(side = 3, line = 0.5, cex = 0.75, "Bimodal?")
      mtext(side = 3, line = -0.25, cex = 0.6, paste("D = ", round(dip.test(usedat$dev)$statistic[1],2), ", P = ", round(dip.test(usedat$dev)$p.value[1],3)))
      print(ii)
    }
  }
}
dev.off()

#number of dip statistics that are significant
sum(output$dip.p <= 0.05)

head(output)
output = data.frame(output)
summary(output)
outdat = merge(output, traits, by = "crop")
summary(outdat)
unique(outdat$self_pollination)
hist(outdat$ac1)

head(outdat)
cropname<- data.frame(unique(outdat$crop))

outdat[order(outdat$ac1),]

# merge with species Latin names
speciesname <- read.csv("crop_latin_names.csv", header=F)
outdat = merge(outdat, speciesname, by.x = "crop", by.y = "V2")
head(outdat_loess)

outdat = outdat %>% 
  rename("sp_latin" = V1)

#write_csv(outdat, "outdat.csv")


#########


#to plot figure 3a and 4a

pdf("AC1_Coef_species2.pdf", height = 5, width = 7)
par(xpd=F)
par(mfrow = c(1,2))

m0a = glm(ac1 ~ -1 + crop, data = outdat)
(means = coef(m0a)[order(coef(m0a))])
li = confint(m0a)[order(coef(m0a)),1]
ui = confint(m0a)[order(coef(m0a)),2]
mylab = substr(names(means), start = 5, stop = 15)

labs2 = substr(traits$crop, start = 1, stop = 10)

nums2 = array()
for(i in 1:length(mylab)){
  nums2[i] = (which(labs2 == mylab[i]))
}

## "indian red"
par(xpd=TRUE)
par(mfrow = c(1,1))
polltype = traits$pollination_syndrome[nums2]
plotCI(1:27, means, li = li, ui = ui, xaxt = "n", xlab = "", ylab = "", pch = 21, cex = 1.4, pt.bg = as.factor(polltype))
points(c(0,27), c(0,0), type = "l", lty = "dotted")
mtext(side = 2, cex = 1.4, line = 2.5, "Lag-1 autocorrelation")
mtext(side = 1, cex = 1.4, line = 3.5, "Species")
text(x = 1:27, y = par("usr")[3] - 0.05, labels = mylab, cex = 1, srt = 45, adj = c(1,.5))
legend("topleft", legend = c("Wind", "Insect"), pch = 21, cex = 1, pt.bg = c("indian red", "black"))

dev.off()



pdf("CV_Coef_species3.pdf", height = 5, width = 7)
par(xpd=F)
par(mfrow = c(1,1))

m1a = glm(cv ~ -1 + crop, data = outdat)
(means = coef(m1a)[order(coef(m1a))])
li = confint(m1a)[order(coef(m1a)),1]
ui = confint(m1a)[order(coef(m1a)),2]

mylab = substr(names(means), start = 5, stop = 15)

labs2 = substr(traits$crop, start = 1, stop = 10)

nums2 = array()
for(i in 1:length(mylab)){
  nums2[i] = (which(labs2 == mylab[i]))
}
par(xpd=TRUE)
par(mfrow = c(1,1))
polltype = traits$pollination_syndrome[nums2]
plotCI(1:27, means, li = li, ui = ui, xaxt = "n", xlab = "", ylab = "", pch = 21, cex = 1.4, pt.bg = as.factor(polltype))
mtext(side = 2, cex = 1.4, line = 2.5, "Coefficient of variation")
mtext(side = 1, cex = 1.4, line = 3.5, "Species")
text(x = 1:27, y = par("usr")[3] - 0.02, labels = mylab, cex = 1, srt = 45, adj = c(1,.5))
legend("topleft", legend = c("Wind", "Insect"), pch = 21, cex = 1, pt.bg = c("indian red", "black"))

dev.off()




#########

#to plot figure S6

# note: trendline values obtained from the "FAO_AB_phylogeny" R file

pdf("All crops coefvar2.pdf", height = 6, width = 7.5)
par(xpd=F)
crops2 = unique(outdat$crop)
mypch = order(crops2) #only works for 26 crops because symbols 26:31 are undefined
mypch[1:27] = c(0:26)
#mypch[27] = 34
mypch[26:27] = c(33:34)
mypch[8] = 126
mypch[16] = 20
mypch[17] = 38
mypch[20] = 37
mypch[21] = 35
mypch[which(crops2 == outdat$crop[1])]
outdat$pch = 1

for(i in 1:dim(outdat)[1]){outdat$pch[i] = mypch[which(crops2 == outdat$crop[i])]}
mycol = as.factor(traits$pollination_syndrome)[order(traits$crop)]

plot(outdat$cv, outdat$ac1, xlab = "", ylab = "", pch = outdat$pch, col = as.factor(outdat$pollination_syndrome), bg = "grey", xlim = c(0,1))
clip(0,0.88, -0.8, 0.6)
abline(a = 0.0492025, b = -0.6015263)
mtext(side = 1, cex = 1.4, line = 2.5, "Coefficient of variation")
mtext(side = 2, cex = 1.4, line = 2.5, "Lag-1 autocorrelation")
legend("topright", legend = substr(crops2, start = 1, stop = 10), pch = mypch, col = mycol, pt.bg = "grey", cex = 0.75)

dev.off()

##########

# to make figure 1 map 

library(maps)
library(ggmap)


register_google(key = "XXXXX") # in order to geocode countries, you need to obtain a Google key 
outdat99 = outdat
outdat99$country[outdat$country == "China, mainland"] = "China" #hack because google puts "China, mainland" near cuba.  Interesting...
outdat99$country[outdat$country == "Georgia"] = "Tbilisi" #hack because google puts "Georgia" in the SE USA
outdat99$country[outdat$country == "Burkina Faso"] = "Ouagadougou" #hack because google puts "Burkina Faso" in the SW USA
outdat99$country[outdat$country == "Togo"] = "Lome" #hack because google puts "Togo" at lon = -87

mylocs = geocode(outdat99$country)
mylocs2 = geocode(unique(outdat99$country)) # checking to make the geocoding work
cbind(unique(outdat99$country), mylocs2) # checking the geocoding 

maplon = jitter(mylocs$lon, amount = 5)
maplat = jitter(mylocs$lat, amount = 5)
maplon[outdat$country == "Bolivia"] = mylocs$lon[outdat$country == "Bolivia"] # removing jitter from Bolivia for chestnuts

# map with big legend
pdf("All crops map.pdf", height = 6, width = 10)
par(xpd=F)
maps::map(database = "world",regions = ".", lwd = 0.25, col = "black")
points(maplon, maplat,pch = outdat$pch, col = as.factor(outdat$pollination_syndrome), bg = "grey", cex = 0.75)
legend("left", legend = crops2, bg = "white", pch = mypch, col = mycol, pt.bg = "grey", cex = 0.7)
dev.off()



#############


# making a phylogeny 

#devtools::install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)
#install.packages("devtools")
library(devtools)
#install.packages("phytools")
library(phytools)
#install.packages("ape")
library(ape)
#install.packages("phylo_glm")

splist<- read.csv("crop_phylo_for_r.csv", header = T)
splist <- read.csv("/Users/gabrielagarcia/Google Drive/Graduate School/Masting and Pollinators Review/AB and Pollinators Team Files /R code for AB and Pollination Syndrome graphs/Revision R Code/crop_phylo_for_r.csv", header=T)
splist1 <- as.list(splist$V1)


phy3<-phylo.maker(splist, tree=GBOTB.extended, nodes=nodes.info.1, , output.tree = T, scenarios="S3")                    
phy3

# extracting the phylogeny 
tree3 = phy3$scenario.3

# plotting the tree
plotTree(tree3)

write.tree(tree3, "FAO_AB_tree_update.tre")



