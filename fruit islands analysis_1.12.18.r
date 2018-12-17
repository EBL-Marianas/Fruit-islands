#fruit islands analysis_1.12.18
library(Rcpp)
library(ggplot2)
library(MASS)
library(car)
library(pscl)
library(boot)
library(AER)
library(DHARMa)

setwd("C:/Users/Henry/Documents/Manuscripts/Guam/fruit islands paper")
fi<-read.table("fruit islands analysis_1.12.18.txt", header=TRUE, sep="\t", na.strings="")
fisf<-fi[(fi$location=="sf"),]
fiaafb<-fi[(fi$location=="aafb"),]
sd(fisf$dbh, na.rm = TRUE)
mean(fisf$dbh, na.rm=TRUE)
#question: is the sample of overstory trees we chose at each site similar in composition?
###need to do a tapply to tally the totals at each site, then a stacked bar graph for each site with species as fill
tree<-read.table("aafb vs sf_overstory trees.txt", sep="\t",header=TRUE, na.strings="")
ggplot(data=tree, aes(x=species, y= frequency, fill=site)) + geom_bar(stat="identity", position="dodge")

#chi-square test to see if proportions of focal trees differ between sites
tbl = table(fi$location, fi$tree_sp)
chi<-chisq.test(tbl)###proportions differ very significantly


#question: does AAFB have more MIST than SF?
#this will help confirm assumption that MIST are the ones dispersing seeds

#is there a correlation between # MIST and # seedlings?
#regular poisson
p<-glm(total_seedlings~location, family= "poisson",data = ficomplete)#no correlation - p = 0.59
dispersiontest(p)
sim_p <- simulateResiduals(p, refit=T)
testOverdispersion(sim_p)
plotSimulatedResiduals(sim_p)
1 - pchisq(summary(p)$deviance,summary(p)$df.residual)#poisson model does not fit the data

#zero-inflated poisson
z<- zeroinfl(total_seedlings~location, data = ficomplete)
vuong(p,z)

#negative binomial model
nb<-glm.nb(total_seedlings~location,data = ficomplete)
sim_nb <- simulateResiduals(nb, refit=T, n=99)
testOverdispersion(sim_nb)
plotSimulatedResiduals(sim_nb)
1 - pchisq(summary(nb)$deviance, summary(nb)$df.residual)

c(nb=AIC(nb), p=AIC(p))###the negative binomial models are wayyyy better!!!


#does the number of MIST differ between sites?
summary(lm(mist_number~location, data = fi))#yes!!! (p<0.0001)


##linear regression looks weird - OMIT
#ggplot(fi, aes(x=mist_number, y = total_vines.seedlings, shape=location))+geom_point(size=3)+
# scale_shape_manual(values=c(1,2))

#same as previous but as bar graph - OMIT - just use t-test
#ggplot(fi, aes(x=location, y = mist_number))+geom_boxplot()

#is there a relationship between size of the fruit island and # seedlings/species?
#remove outlier
fi2<-fi[c(0:32,34:206),]
summary(lm(total_vines.seedlings~dbh, data=fi2))#no; p =0.93
plot(total_vines.seedlings~dbh, data=fi2)# OMIT

####does average DBH differ between SF and AAFB?
summary(lm(dbh~location, data=fi))#no; p = 0.21
plot(dbh~location, data=fi)

#main question: are there differences in total seedlings/seedling diversity between sites??
# vines
summary(lm(total_vines~location,data=fi))#p<0.0001; more vines at aafb!
plot(total_vines~location, data=fi)
#seedlings
summary(lm(total_seedlings~location, data=fi))#p<0.0001; more seedlings at aafb!
plot(total_vines.seedlings~location, data=fi)
#seedlings+vines 
summary(lm(total_vines.seedlings~location, data=fi))
#species diversity
summary(aov(total_seedling.spp~location, data=fi))#p>0.0001 =  more seedlings at aafb!
plot(total_seedling.spp~location, data=fi)


#AIC analysis
library(MuMIn)
ficomplete<-read.table("fruit islands analysis_1.12.18_complete cases.txt", header=TRUE, sep="\t", na.strings="")
ficomplete<-ficomplete[ficomplete$location=="aafb",]
globalmodelseedlings <- aov(total_seedlings~location+tree_sp+dbh+mist_number, data=ficomplete,na.action = "na.fail")
combinations<- dredge(globalmodelseedlings)

mod1<-glm.nb(total_seedlings~location, data=ficomplete)
mod2<-glm.nb(total_seedlings~tree_sp, data=ficomplete)
mod3<-glm.nb(total_seedlings~dbh, data=ficomplete)
mod4<-glm.nb(total_seedlings~mist_number, data=ficomplete)
mod5<-glm.nb(total_seedlings~location*tree_sp, data=ficomplete)
mod6<-glm.nb(total_seedlings~location+tree_sp, data=ficomplete)
mod7<-glm.nb(total_seedlings~location+mist_number, data=ficomplete)
mod8<-glm.nb(total_seedlings~location+dbh, data=ficomplete)
mod9<-glm.nb(total_seedlings~location*dbh, data=ficomplete)
mod10<-glm.nb(total_seedlings~location*mist_number, data=ficomplete)


model.sel(mod1,mod2,mod3,mod4, mod6, mod7, mod8)#could be problematic because of artocarpus


globalmodeldiversity <- aov(total_seedling.spp~location+tree_sp+dbh+mist_number, data=ficomplete,na.action = "na.fail")
combinations<- dredge(globalmodeldiversity)
mod1<-glm.nb(total_seedling.spp~location, data=ficomplete)
mod2<-glm.nb(total_seedling.spp~tree_sp, data=ficomplete)
mod3<-glm.nb(total_seedling.spp~dbh, data=ficomplete)
mod4<-glm.nb(total_seedling.spp~mist_number, data=ficomplete)
mod5<-glm.nb(total_seedling.spp~location*tree_sp, data=ficomplete)
mod6<-glm.nb(total_seedling.spp~location+tree_sp, data=ficomplete)
mod7<-glm.nb(total_seedling.spp~location+mist_number, data=ficomplete)
mod8<-glm.nb(total_seedling.spp~location+dbh, data=ficomplete)
mod9<-glm.nb(total_seedling.spp~location*dbh, data=ficomplete)

model.sel(mod1,mod2,mod3,mod4,mod5, mod6, mod7, mod8, mod9)#could be problematic because of artocarpus



###do vines differ between sites? chi square test
tbl1<-table(fi$location,fi$Momordica)
tbl2<-table(fi$location,fi$Coccinia)
tbl3<-table(fi$location,fi$Passiflora.foetida)
tbl4<-table(fi$location,fi$Passiflora.suberosa)
vine<-cbind(tbl1,tbl2,tbl3,tbl4)
vine<-data.frame(vine)
names(vine)[1]<-"coccinia"
names(vine)[2]<-"momordica"
names(vine)[3]<-"pf"
names(vine)[4]<-"ps"
chisq.test(vine)

ggplot(fisf, aes(x=dbh, y = total_seedling.spp))+geom_point()+
  geom_smooth(method = "loess", size = 1.5)


##exploring correlations among predictor variables
summary(glm.nb(dbh~tree_sp, data=fi))#R^2 = 0.06
summary(glm.nb(mist_number~dbh, data=fi))#R^2 = 0.001
summary(glm.nb(mist_number~tree_sp, data=fi))#R^2 = 0.002

plot(dbh~tree_sp, data=fi)


