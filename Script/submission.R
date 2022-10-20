#### Script

#working directory 
setwd("../Data")

###################################################################
## trait analysis and comparisons
###################################################################

trait <- read.csv("traits.csv")

# Summary of fruits
unique(trait$SpecName[which(trait$organism=="Fruit")])
summary(trait$diam[which(trait$SpecName == "Acrocomia crispa")])
summary(trait$diam[which(trait$SpecName == "Catesbaea spinosa")])
summary(trait$diam[which(trait$SpecName == "Chrysobalanus icaco")])
summary(trait$diam[which(trait$SpecName == "Coccoloba uvifera")])
summary(trait$diam[which(trait$SpecName == "Coccotherinax barbadensis")])
summary(trait$diam[which(trait$SpecName == "Cordia sebestena")])
summary(trait$diam[which(trait$SpecName == "Goetzea elegans")])
summary(trait$diam[which(trait$SpecName == "Pseudophoenix vinifera")])
summary(trait$diam[which(trait$SpecName == "Syagrus amara")])
summary(trait$diam[which(trait$SpecName == "Theophrasta jussieui")])

fruit_dat <- trait[which(trait$organism == "Fruit"),]
boxplot(fruit_dat$diam ~ fruit_dat$SpecName, las=2, ylab="diameter (mm)", xlab="", ylim=c(0,60))  #overview of fruit traits

trait$organism2 <- ifelse(trait$organism == "Tortoise", "large reptile", trait$organism)
trait$organism2 <- ifelse(trait$organism2 == "Iguana", "large reptile", trait$organism2)

### We are using the averages of each species
# Creating a new dataset for just the averages
data_average <- data.frame(matrix(ncol=3, nrow=length(unique(trait$SpecName))))
colnames(data_average) <- c('Species', 'diam', 'organism')

species <- unique(trait$SpecName)

for(i in c(1:length(species))){
  data_average$Species[i] <- species[i]
  data_average$diam[i] <- mean(trait$diam[which(trait$SpecName == species[i])])
  data_average$organism[i] <- unique(trait$organism2[which(trait$SpecName == species[i])])
}

data_average$organism <- as.factor(data_average$organism)

## Kruskal wallis
# Total
m2<-kruskal.test(diam ~ organism, data=data_average)
print(m2)   #p < 0.001

library("PMCMRplus")
posthocs2<-dscfAllPairsTest(diam ~ organism, data=data_average)
print(posthocs2) 

## Figure 1

data_average$code <- NA
small_lizard <- data_average[which(data_average$organism=="Small_lizards"),]
large_reptile <- data_average[which(data_average$organism=="large reptile"),]
fruit <- data_average[which(data_average$organism=="Fruit"),]

fruit$code <- "A"
small_lizard$code <- "B"
large_reptile$code <- "C"

dat_average <- rbind(small_lizard, large_reptile, fruit)
dat_average$diam <- as.numeric(dat_average$diam)
dat_average$organism <- as.factor(dat_average$organism)

library(ggridges)
library(ggplot2)
library(plyr)

## Comparing small lizard vs cyclura vs tortoise
ggplot(dat_average, aes(x = diam, y = code, fill = code)) +
  geom_density_ridges2(
    aes(point_color = code, point_fill = code, point_shape = code),
    alpha = .2, point_alpha = 1, jittered_points = TRUE, scale=0.9) +
  scale_point_color_hue(l = 40) +
  scale_y_discrete(labels=c("A" = "Fruit", "B" = "Small Lizard", "C" = "Large Reptile")) +
  scale_discrete_manual(aesthetics = "point_shape",values = c(21, 22, 23))+
  theme_ridges()+
  xlab("diam (mm)")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6))+
  theme(legend.position="none") #removing legend

###################################################################
## Germination trial ######################################
###################################################################

library(tidyverse)

### Calculating PV and GV ############ (Table 1)
data.germ.jab.analyse <- read_csv2("dat_germ.csv")

name <- unique(data.germ.jab.analyse$Species)

## Need to input manually for each species. 
temp <- data.germ.jab.analyse[which(data.germ.jab.analyse$Species==name[1]),]    #change name[1] to other integers for other species

### Input integers using values from the "no_alive" column in Supp. Table 3. Example below is for Acrocomia crispa, which had 14, 7, and 13 potentially alive seeds at the end of the greenhouse trials for Control, Inguana-ingested, and Tortoise-ingested respectively.
library(germinationmetrics)
PV <- PeakValue(germ.counts = temp$Control, intervals = temp$Days, total.seeds = 14)
PV
GV <- GermValue(germ.counts = temp$Control, intervals = temp$Days, total.seeds = 14,
                method = "czabator")
GV$`Germination Value`

PV <- PeakValue(germ.counts = temp$Iguana, intervals = temp$Days, total.seeds = 7)
PV
GV <- GermValue(germ.counts = temp$Iguana, intervals = temp$Days, total.seeds = 7,
                method = "czabator")
GV$`Germination Value`

PV <- PeakValue(germ.counts = temp$Tortoise, intervals = temp$Days, total.seeds = 13)
PV
GV <- GermValue(germ.counts = temp$Tortoise, intervals = temp$Days, total.seeds = 13,
                method = "czabator")
GV$`Germination Value`

##### Survival curves ####### (Supp. Figure 1; stats)
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condsurv)
library(dplyr)
library(ggplot2)

dat_temp <- read.csv("germination_survival_dat.csv")
dat <- dat_temp[which(dat_temp$Seed.fate == 'Alive'),]

dat <- 
  dat %>% 
  mutate(status = recode(Germ.fate, 'No' = 0, 'Yes' = 1))

dat$Sown <- as.POSIXct(strptime(dat$Sown, "%m/%d/%Y"))
dat$Germinated <- as.POSIXct(strptime(dat$Germinated, "%m/%d/%Y"))
dat$time <- dat$Germinated - dat$Sown
dat$time <- as.numeric(dat$time)

## Total
survfit2(Surv(time, status) ~ 1, data = dat) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + xlim(0,200)+
  add_confidence_interval() +
  add_risktable()

survdiff(Surv(time, status) ~ Treatment + (1 | ID), data = dat) #p < 0.001
coxph(Surv(time, status) ~ Treatment + (1|ID), data = dat) %>%  
  tbl_regression(exp = TRUE)   

#### Acrocomia crispa
dat_ac <- dat[which(dat$Species == "Acrocomia crispa"),]
survdiff(Surv(time, status) ~ Treatment, data = dat_ac)   #p=0.09

survfit2(Surv(time, status) ~ Treatment, data = dat_ac) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + xlim(0,200)+
  add_confidence_interval()

survdiff(Surv(time, status) ~ Treatment + (1 | ID), data = dat_ac) #p = 0.09
coxph(Surv(time, status) ~ Treatment + (1|ID), data = dat_ac) %>%  
  tbl_regression(exp = TRUE)   

#### Catesbae_spinosa
dat_cs <- dat[which(dat$Species == "C. spinosa"),]

survfit2(Surv(time, status) ~ Treatment, data = dat_cs) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + xlim(0,40)+
  add_confidence_interval()

survdiff(Surv(time, status) ~ Treatment + (1 | ID), data = dat_cs) #p =0.2
coxph(Surv(time, status) ~ Treatment + (1|ID), data = dat_cs) %>%  
  tbl_regression(exp = TRUE)  

#### Chrysobalanus_icaco
dat_ci <- dat[which(dat$Species == "Chrysobalanus icaco"),]

survfit2(Surv(time, status) ~ Treatment, data = dat_ci) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + xlim(0,200)+
  add_confidence_interval()

survdiff(Surv(time, status) ~ Treatment + (1 | ID), data = dat_ci) #p < 0.001
coxph(Surv(time, status) ~ Treatment + (1|ID), data = dat_ci) %>%  
  tbl_regression(exp = TRUE)  

#### Cocoloba_uvifera
dat_cu <- dat[which(dat$Species == "Cocoloba uvifera"),]

survfit2(Surv(time, status) ~ Treatment, data = dat_cu) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + xlim(0,75)+
  add_confidence_interval()

survdiff(Surv(time, status) ~ Treatment + (1 | ID), data = dat_cu) #p =0.008
coxph(Surv(time, status) ~ Treatment + (1|ID), data = dat_cu) %>%  
  tbl_regression(exp = TRUE) 

#### Cordia_sebestena
dat_cos <- dat[which(dat$Species == "Cordia sebestena"),]
survfit2(Surv(time, status) ~ Treatment, data = dat_cos) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + xlim(0,75)+
  add_confidence_interval()

survdiff(Surv(time, status) ~ Treatment + (1 | ID), data = dat_cos) #p < 0.001
coxph(Surv(time, status) ~ Treatment + (1|ID), data = dat_cos) %>%  
  tbl_regression(exp = TRUE)  

#### Pseudophoenix vinifera   ## Unable to run due to lack of germinated iguana-ingested individuals

##### Figure 2 ######### (and SUpp. Figure 2)

data.germ.jab.graph_temp <- read_csv("dat_germ_perc_a.csv")

data.germ.jab.graph <- subset(data.germ.jab.graph_temp, Species != "Cordia_sebestena" & Species != "Pseudophoenix_vinifera")
data.germ.jab.graph1 <- subset(data.germ.jab.graph_temp, Species == "Cordia_sebestena" | Species == "Pseudophoenix_vinifera")
data.germ.jab.graph1$Tortoise <- NULL

data.germ.jab.graph2 <- subset(data.germ.jab.graph_temp, Species == "Cocoloba uvifera")

glimpse(data.germ.jab.graph)

data.germ.jab.graph <- data.germ.jab.graph_temp %>% 
  gather(key = "treatment", value = "germination", 3:5) 

data.germ.jab.graph2 <- data.germ.jab.graph[which(data.germ.jab.graph$Species =="Cocoloba_uvifera" ),]

### Figure 2   #Cocoloba_uvifera only

ggplot(data = data.germ.jab.graph2, aes(y = germination, x = Days, fill = treatment, color = treatment)) +
  geom_line(size = 1) +
  geom_point(shape=21, color="black", size=3) +
  scale_fill_manual(values=c("blue","cyan","orange"), 
                    name = "", labels = c("Control", "Iguana","Tortoise")) +
  scale_color_manual(values=c("blue","cyan","orange"
  ),
  name = "", labels = c("Control", "Iguana","Tortoise")) +
  labs(title = "",
       y = "Proportion germinated (%)", x = "Time (days)") + 
  theme_bw() +
  theme(legend.position = c(0.92, .35),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        strip.text = element_text(face = "italic"),
        legend.background=element_blank(),
        legend.key=element_blank()) +
  facet_wrap(.~ Species, labeller = labeller
             (Species = c(
                          "Cocoloba_uvifera" = "Cocoloba uvifera"
             )))

## Supp. Figure 2
ggplot(data = data.germ.jab.graph, aes(y = germination, x = Days, fill = treatment, color = treatment)) +
  geom_line(size = 1) +
  geom_point(shape=21, color="black", size=3) +
  scale_fill_manual(values=c("blue","cyan","orange"), 
                    name = "", labels = c("Control", "Iguana","Tortoise")) +
  scale_color_manual(values=c("blue","cyan","orange"
  ),
  name = "", labels = c("Control", "Iguana","Tortoise")) +
  labs(title = "",
       y = "Germination rate (%)", x = "Time (days)") + 
  theme_bw() +
  theme(legend.position = c(0.92, .35),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        strip.text = element_text(face = "italic"),
        legend.background=element_blank(),
        legend.key=element_blank()) +
  facet_wrap(.~ Species, labeller = labeller
             (Species = c("Catesbae_spinosa" = "Catesbae spinosa",
                          "Chrysobalanus_icaco" = "Chrysobalanus icaco",
                          "Cocoloba_uvifera" = "Cocoloba uvifera",
                          "Theophrasta" = "Theophrasta",
                          "Acrocomia_crispa" = "Acrocomia crispa"
                          , "Cordia_sebestena" = "Cordia sebestena",
                          "Pseudophoenix_vinifera" = "Pseudophoenix vinifera"
             )))


