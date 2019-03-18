#Biotropica - Kastner et al. - Functional rescue of seed dispersal by a remnant frugivore population on a defaunated tropical island
library(Rcpp)
library(ggplot2)
library(MASS)
library(car)
library(pscl)
library(boot)
library(AER)
library(DHARMa)

fi<-read.table("fruit islands analysis_3.17.19.Biotropica.txt", header=TRUE, sep="\t", na.strings="")



##seedling abundance by location 
abund<-ggplot(data=fi, aes(y=total_seedlings, x = dbh, colour=location, shape=location))+geom_point(size=3)+
  scale_color_manual(values=c("black","red"), guide=FALSE)+
  scale_shape_manual(values=c(1,17), guide=FALSE)+
  geom_smooth(se=FALSE,size = 1, method="glm.nb")+
  theme(legend.text=element_blank(),
        legend.title =element_blank(),
        legend.background = element_blank())+
  xlab("DBH (cm)")+
  annotate("text", label="a)", x=0, y =120)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x =element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.title.y=element_text(size=14),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  ylab("Seedling Abundance")

#seedling richness by location
rich<-ggplot(data=fi, aes(y=total_seedling.spp, x = dbh, colour=location, shape=location))+geom_point(size=3)+
  scale_color_manual(values=c("black","red"))+
  scale_shape_manual(values=c(1,17))+
  geom_smooth(se=FALSE,size = 1, method="glm.nb")+
  theme(legend.position=c(0.92,0.92),
        legend.title =element_text(size=14),
        legend.background = element_blank())+
  xlab("DBH (cm)")+
  annotate("text", label="b)", x=0, y =5)+ 
  guides(fill=guide_legend(title="Location"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x =element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.title.y=element_text(size=14),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  ylab("Seedling Richness")


#Supplementary Information 4
library(gridExtra)
grid.arrange(abund,rich,ncol=2,nrow=1, widths=c(1.05,1))


#Is the sample of overstory trees we chose at each site similar in composition?
tbl = table(fi$location, fi$tree_sp)
chi<-chisq.test(tbl)###proportions differ very significantly

tree<-read.table("aafb vs sf_overstory trees.txt", sep="\t",header=TRUE, na.strings="")
tree$species<-factor(tree$species, levels=c('Casuarina','Vitex','Other','Calophyllum','Tabebuia'))
ggplot(data=tree, aes(x=species, y= frequency, fill=site)) + geom_bar(stat="identity", position="dodge")+
  scale_fill_discrete(name="Site", labels=c('AAFB','SF'))+
  theme(legend.position=c(0.92,0.88),
        legend.title =element_text(size=14),
        legend.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x =element_text(size=12, face='italic'),
        axis.title.x=element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y=element_text(size=14),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
ylab("Frequency")


#Does the number of MIST differ between sites?
summary(lm(mist_number~location, data = fi))#yes!!! (p<0.0001)

#Does average DBH differ between SF and AAFB?
summary(lm(dbh~location, data=fi))#no; p = 0.21
plot(dbh~location, data=fi)

#Are there differences in total seedlings/seedling diversity between sites??
# vine abundance
summary(lm(total_vines~location,data=fi))#p<0.0001; more vines at aafb
plot(total_vines~location, data=fi)
#seedling abundance
summary(lm(total_seedlings~location, data=fi))#p<0.0001; more seedlings at aafb
plot(total_vines.seedlings~location, data=fi)
#seedlings+vine abundance
summary(lm(total_vines.seedlings~location, data=fi))
#species richness
summary(aov(total_seedling.spp~location, data=fi))#p>0.0001 =  more seedlings at aafb
plot(total_seedling.spp~location, data=fi)


#AIC analysis
library(MuMIn)

#seedling abundance
mod1<-glm.nb(total_seedlings~location, data=fi)
mod2<-glm.nb(total_seedlings~tree_sp, data=fi)
mod3<-glm.nb(total_seedlings~dbh, data=fi)
mod4<-glm.nb(total_seedlings~location*tree_sp, data=fi)
mod5<-glm.nb(total_seedlings~location+tree_sp, data=fi)
mod6<-glm.nb(total_seedlings~location+dbh, data=fi)
mod7<-glm.nb(total_seedlings~location*dbh, data=fi)

#compare model AICs
model.sel(mod1,mod2,mod3,mod4, mod5, mod6,mod7)

#seedling richness
mod1<-glm.nb(total_seedling.spp~location, data=fi)
mod2<-glm.nb(total_seedling.spp~tree_sp, data=fi)
mod3<-glm.nb(total_seedling.spp~dbh, data=fi)
mod4<-glm.nb(total_seedling.spp~location*tree_sp, data=fi)
mod5<-glm.nb(total_seedling.spp~location+tree_sp, data=fi)
mod6<-glm.nb(total_seedling.spp~location+dbh, data=fi)
mod7<-glm.nb(total_seedling.spp~location*dbh, data=fi)

#compare model AICs
model.sel(mod1,mod2,mod3,mod4, mod5, mod6,mod7)


###Do vines differ between sites? chi square test
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
chisq.test(vine)# yes, they are more abundant at AAFB



##Exploring correlations among predictor variables
summary(glm.nb(dbh~tree_sp, data=fi))#R^2 = 0.06
summary(glm.nb(mist_number~dbh, data=fi))#R^2 = 0.001
summary(glm.nb(mist_number~tree_sp, data=fi))#R^2 = 0.002
summary(lm(mist_number~location, data=fi))#R^2 = 0.002