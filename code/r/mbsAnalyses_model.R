
####################################################################
################# Computational model ##############################
####################################################################
################# Katyal, Huys, Dolan, Fleming #####################
### How underconfidence is maintained in anxiety and depression ####
####################################################################
##  This file reproduces all figures for the computational model ###
####################################################################

library(R.matlab)
library(tidyverse)
library(ggplot2)
library(see)
library(ggpubr)
library(ggbeeswarm) # for geom_quasirandom()
library(bayestestR)

font_size <- 20
tag_size <- 22
pval_size <- 6
alpha <- .2

percmemColours = c("cyan3", "tomato2")
taskNames = c('Perception', 'Memory')

# update this folder with the folder location with data
data.folder <- "~/OneDrive - University of Copenhagen/Projects/Experiments/metaBiasShift/data/exp"

expNum <- 1
setwd(paste(data.folder, expNum, sep = ''))
setwd('processed/')
mbsDataExp1 = readMat(paste('mbsDataExp',expNum,'.mat',sep='')) # load data
# mbsDataExp1 = readMat(paste('mbsDataExp',expNum,'_excPR.mat',sep='')) # load data with more stringent exclusion criterion
mbsDataExp1 <- mbsDataExp1[[paste('mbsDataExp',expNum,sep='')]]

expNum <- 2
setwd(paste(data.folder, expNum, sep = ''))
setwd('processed/')
mbsDataExp2 = readMat(paste('mbsDataExp',expNum,'.mat',sep='')) # load data
# mbsDataExp2 = readMat(paste('mbsDataExp',expNum,'_excPR.mat',sep='')) # load data with more stringent exclusion criterion
mbsDataExp2 <- mbsDataExp2[[paste('mbsDataExp',expNum,sep='')]]

####################
#########################
####################
##########
#####
##### Plot model fits - basic no distortion model

dic <- array(dim = c(2,5))

mbsDataFit <- mbsDataExp1[,,1]$fitZ
dic[1,1] <-  mbsDataFit[,,1]$model.0[,,1]$dic
dic[1,2] <-  mbsDataFit[,,1]$model.1[,,1]$dic
dic[1,3] <-  mbsDataFit[,,1]$model.2[,,1]$dic
dic[1,4] <-  mbsDataFit[,,1]$model.3[,,1]$dic
dic[1,5] <-  mbsDataFit[,,1]$model.4[,,1]$dic
dic[1,] <- 1*(dic[1,] - min(dic[1,]))/1#min(dic[1,])

mbsDataFit <- mbsDataExp2[,,1]$fitZ # select this one to plot exp 2
dic[2,1] <-  mbsDataFit[,,1]$model.0[,,1]$dic
dic[2,2] <-  mbsDataFit[,,1]$model.1[,,1]$dic
dic[2,3] <-  mbsDataFit[,,1]$model.2[,,1]$dic
dic[2,4] <-  mbsDataFit[,,1]$model.3[,,1]$dic
dic[2,5] <-  mbsDataFit[,,1]$model.4[,,1]$dic
dic[2,] <- 1*(dic[2,] - min(dic[2,]))/1#min(dic[2,])

dfdic <- data.frame(c(t(dic)), c(rep('Exp 1',5), rep('Exp 2', 5)), c(seq(5)-1, seq(5)-1))
colnames(dfdic) <- c('dic', 'expNum', 'modelNum')
dfdic <- dfdic %>%
  mutate_at(c('expNum', 'modelNum'), as.factor) %>%
  mutate(modelNum = recode_factor(modelNum, "0" = 'M0: No asymm.', 
                                "1" = 'M1: Dom. gen.  \n   FB-Conf Same',
                                "2" = 'M2: Dom. gen.  \n   FB-Conf Diff.',
                                "3" = 'M3: Dom. spec.\n   FB-Conf Same',
                                "4" = 'M4: Dom. spec.\n   FB-Conf Diff.') )

####  Supplementary Figure 12 - Model comparison of basic no-distortion models
ggplot(dfdic, aes(y=modelNum,x=dic,fill=expNum)) +
  scale_fill_manual(breaks = c('Exp 1', 'Exp 2'),
                    values=c('gray75', 'gray15')) +
  geom_vline(xintercept = 0, color = 'grey20') +
  geom_bar(stat="identity", width=.7, position = "dodge") +
  theme_pubclean() +
  xlab("DIC") +
  ylab("") +
  # labs(tag = "A") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        legend.position = c(.3,1.07),
        plot.margin = margin(t=35),
        # legend.justification = c('left', 'top'),
        # legend.box.background = element_rect(color = 'black'),
        # legend.margin = margin(1,5,5,5),
        legend.title = element_blank(),
        text = element_text(size=font_size-4), 
        # legend.box.just = 'left',
        plot.caption = element_text(size = font_size-4),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0))+ 
  guides(fill = guide_legend(nrow = 1))


#######################
####################
###############
##### Model fits for mental health distortions models

n.models <- 7
dic <- array(dim = c(2,n.models))

mbsDataFit <- mbsDataExp1[,,1]$fitRegZ
dic[1,1] <-  mbsDataExp1[,,1]$fitZ[,,1]$model.4[,,1]$dic
dic[1,2] <-  mbsDataFit[,,1]$model.6.q2[,,1]$dic
dic[1,3] <-  mbsDataFit[,,1]$model.5.q2[,,1]$dic
dic[1,4] <-  mbsDataFit[,,1]$model.7.q2[,,1]$dic
dic[1,5] <-  mbsDataFit[,,1]$model.8.q2[,,1]$dic
dic[1,6] <-  mbsDataFit[,,1]$model.10.q2[,,1]$dic
dic[1,7] <-  mbsDataFit[,,1]$model.9.q2[,,1]$dic
# dic[1,8] <-  mbsDataFit[,,1]$model.12.q1[,,1]$dic
dic[1,] <- 1*(dic[1,] - min(dic[1,]))/1#min(dic[1,])

mbsDataFit <- mbsDataExp2[,,1]$fitRegZ # select this one to plot exp 2
dic[2,1] <-  mbsDataExp2[,,1]$fitZ[,,1]$model.4[,,1]$dic
dic[2,2] <-  mbsDataFit[,,1]$model.6.q1[,,1]$dic
dic[2,3] <-  mbsDataFit[,,1]$model.5.q1[,,1]$dic
dic[2,4] <-  mbsDataFit[,,1]$model.7.q1[,,1]$dic
dic[2,5] <-  mbsDataFit[,,1]$model.8.q1[,,1]$dic
dic[2,6] <-  mbsDataFit[,,1]$model.10.q1[,,1]$dic
dic[2,7] <-  mbsDataFit[,,1]$model.9.q1[,,1]$dic
# dic[2,8] <-  mbsDataFit[,,1]$model.12.q1[,,1]$dic
dic[2,] <- 1*(dic[2,] - min(dic[2,]))/1#min(dic[2,])

dfdic <- data.frame(c(t(dic)), c(rep('Exp 1',n.models), rep('Exp 2', n.models)), c(seq(n.models)-1, seq(n.models)-1))
colnames(dfdic) <- c('dic', 'expNum', 'modelNum')
dfdic <- dfdic %>%
  mutate_at(c('expNum', 'modelNum'), as.factor) %>%
  mutate(modelNum = recode_factor(modelNum, "0" = 'D0: No distortion', 
                                            "1" = 'D1: Confidence   ',
                                            "2" = 'D2: Feedback     ',
                                            "3" = 'D3: Resp.    bias',
                                            "4" = 'D4: Conf. + FB   ',
                                            "5" = 'D5: Conf. +  Bias',
                                            "6" = 'D6: FB    +  Bias',
                                            "7" = 'D7: FB + Conf. + Bias') ) %>%
  mutate(modelNum = factor(modelNum, levels = c("0" = 'D0: No distortion', 
                                  "1" = 'D1: Confidence   ',
                                  "2" = 'D2: Feedback     ',
                                  "3" = 'D3: Resp.    bias',
                                  "4" = 'D4: Conf. + FB   ',
                                  "5" = 'D5: Conf. +  Bias',
                                  "6" = 'D6: FB    +  Bias',
                                  "7" = 'D7: FB + Conf. + Bias') ))

#### Figure 4E - Model comparison of distortion models

ggplot(filter(dfdic, expNum=='Exp 1'), aes(y=modelNum,x=dic)) +
  # scale_fill_manual(breaks = c('Exp 1', 'Exp 2'),
  #                   values=c('gray75', 'gray15')) +
  # values=c('darkorange1', 'purple4')) +
  geom_vline(xintercept = 0, color = 'grey20') +
  geom_bar(stat="identity", width=.6, position = "dodge") +
  theme_pubclean() +
  xlab("DIC") +
  ylab("") +
  scale_x_continuous(limits=c(0,380), breaks=c(0,150,300))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        legend.position = c(.7,.95),
        plot.margin = margin(t=35, r = 15),
        legend.title = element_blank(),
        text = element_text(size=font_size-2), 
        plot.caption = element_text(size = font_size-2),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0))+ 
  guides(fill = guide_legend(nrow = 2))

ggplot(filter(dfdic, expNum=='Exp 2'), aes(y=modelNum,x=dic)) +
  # scale_fill_manual(breaks = c('Exp 1', 'Exp 2'),
  #                   values=c('gray75', 'gray15')) +
  # values=c('darkorange1', 'purple4')) +
  geom_vline(xintercept = 0, color = 'grey20') +
  geom_bar(stat="identity", width=.6, position = "dodge") +
  theme_pubclean() +
  xlab("DIC") +
  ylab("") +
  scale_x_continuous(limits=c(0,45), breaks=c(0,20,40))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        legend.position = c(.7,.95),
        plot.margin = margin(t=35, r = 15),
        legend.title = element_blank(),
        text = element_text(size=font_size-2), 
        plot.caption = element_text(size = font_size-2),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0))+ 
  guides(fill = guide_legend(nrow = 2))

##################################################################
######## Read in the fitted models to plot posteriors

ci.lvl = .99

mbsDataFit1 <- mbsDataExp1[,,1]$fitRegZ
mbsDataFit2 <- mbsDataExp2[,,1]$fitRegZ

## for reporting stats
phq <- mbsDataFit1[,,1]$model.8.q1 # conf + add distortions 
phq <- phq[,,1]$samples
phq <- phq[,,1]

hphq <- hdi(phq$beta.fb.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(phq$beta.conf.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)

gad <- mbsDataFit1[,,1]$model.8.q2
gad <- gad[,,1]$samples
gad <- gad[,,1]
hphq <- hdi(gad$beta.fb.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(gad$beta.conf.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)

spin <- mbsDataFit1[,,1]$model.8.q3
spin <- spin[,,1]$samples
spin <- spin[,,1]
hphq <- hdi(spin$beta.fb.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(spin$beta.conf.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)

td <- mbsDataFit2[,,1]$model.8.q1[,,1]$samples[,,1]
# ad <- ad[,,1]$samples
# ad <- ad[,,1]
ad <- c()
ad$beta.fb.lr <- td$beta.fb.lr
ad$beta.conf.lr <- td$beta.conf.lr
ad$beta.post.bias <- td$beta.post.bias
hphq <- hdi(ad$beta.fb.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(ad$beta.conf.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)


phq <- mbsDataFit1[,,1]$model.10.q1 # conf + add distortions 
phq <- phq[,,1]$samples
phq <- phq[,,1]

hphq <- hdi(phq$beta.fb.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(phq$beta.conf.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(phq$beta.post.bias,ci.lvl)
c(hphq$CI_low, hphq$CI_high)

gad <- mbsDataFit1[,,1]$model.10.q2
gad <- gad[,,1]$samples
gad <- gad[,,1]
hphq <- hdi(gad$beta.fb.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(gad$beta.conf.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
# 
# spin <- mbsDataFit2[,,1]$model.6.q3
# spin <- spin[,,1]$samples
# spin <- spin[,,1]
# hphq <- hdi(spin$beta.fb.lr,ci.lvl)
# c(hphq$CI_low, hphq$CI_high)
# hphq <- hdi(spin$beta.conf.lr,ci.lvl)
# c(hphq$CI_low, hphq$CI_high)

td <- mbsDataFit2[,,1]$model.10.q1[,,1]$samples[,,1]
# ad <- ad[,,1]$samples
# ad <- ad[,,1]
ad <- c()
ad$beta.fb.lr <- td$beta.fb.lr
ad$beta.conf.lr <- td$beta.conf.lr
ad$beta.post.bias <- td$beta.post.bias
hphq <- hdi(ad$beta.fb.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(ad$beta.conf.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(ad$beta.post.bias,ci.lvl)
c(hphq$CI_low, hphq$CI_high)

posneg = data.frame(c(gad$beta.fb.lr),
                    c(gad$beta.conf.lr),
                    c(gad$beta.post.bias),
                    c(ad$beta.fb.lr)/5,
                    c(ad$beta.conf.lr)/5,
                    c(ad$beta.post.bias)/5)
colnames(posneg) <- c('gad_fb', 'gad_conf', 'gad_add',
                      'ad_fb', 'ad_conf', 'ad_add')

posneg2 <- posneg %>% pivot_longer(
  cols = c(1:6),
  names_sep = '_',
  values_to = 'beta',
  names_to = c('questname', 'fbconf')) %>% 
  mutate_at(c('fbconf', 'questname'), as.factor) %>%
  mutate(fbconf = recode_factor(fbconf, "fb" = 'Feedback', 
                                "conf" = 'Confidence',
                                'add' = 'Additive') ) %>%
  mutate(expNum = recode_factor(questname, "phq" = 'Exp 1', 
                                "gad" = 'Exp 1', "spin" = 'Exp 1', 'ad' = 'Exp 2'))%>%
  mutate(questname = recode_factor(questname, "phq" = 'PHQ', "ad" = 'AD')) %>%
  mutate(questname = factor(questname, levels = c('PHQ', 'AD')))

#### #### #### 
#### Supplementary Figure 6A - Posterior distributions of Beta-confidence and Beta-additive for model 10

ggplot(filter(posneg2, fbconf=='Confidence'), aes(x = expNum, y = beta, fill = expNum))+
  scale_fill_manual(breaks = c('Exp 1', 'Exp 2'),
                    values=c('darkorange1', 'darkorange1')) +
  geom_hline(yintercept = 0, color = 'grey50', size=1) +
  geom_violinhalf(alpha=.7)+
  theme_pubclean() +
  xlab("Mental health score") +
  ylab("β-confidence") +
  labs(fill = 'Dataset') +
  coord_cartesian(ylim=c(-.04,0.02))+
  # labs(tag = "E") +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_blank()  ,
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none',
        # legend.direction = "horizontal"
  ) 

ggplot(filter(posneg2, fbconf=='Additive'), aes(x = expNum, y = beta, fill = expNum))+
  geom_hline(yintercept = 0, color = 'grey50', size=1) +
  scale_fill_manual(breaks = c('Exp 1', 'Exp 2'),
                    values=c('maroon4', 'maroon4')) +
  geom_violinhalf(alpha=.7)+
  theme_pubclean() +
  ylab("β-additive") +
  xlab("Mental health score") +
  labs(fill = 'Dataset') +
  coord_cartesian(ylim=c(-.0001,0.0002))+
  # labs(tag = "E") +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        # axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none',
        axis.title.x=element_blank()  ) 
# + 
#   guides(fill = guide_legend(nrow = 1))

#################
##### Plot model where beta(confidence) and beta(feedback) are modelled simultaneously
#####

gad <- mbsDataFit1[,,1]$model.8.q2
gad <- gad[,,1]$samples
gad <- gad[,,1]

td <- mbsDataFit2[,,1]$model.8.q1[,,1]$samples[,,1]
ad <- c()
ad$beta.fb.lr <- td$beta.fb.lr
ad$beta.conf.lr <- td$beta.conf.lr
ad$beta.post.bias <- td$beta.post.bias

posneg = data.frame(c(gad$beta.fb.lr),
                    c(gad$beta.conf.lr),
                    c(gad$beta.post.bias),
                    c(ad$beta.fb.lr)/5,
                    c(ad$beta.conf.lr)/5,
                    c(ad$beta.post.bias)/5)
colnames(posneg) <- c('gad_fb', 'gad_conf', 'gad_add',
                      'ad_fb', 'ad_conf', 'ad_add')

posneg2 <- posneg %>% pivot_longer(
  cols = c(1:6),
  names_sep = '_',
  values_to = 'beta',
  names_to = c('questname', 'fbconf')) %>% 
  mutate_at(c('fbconf', 'questname'), as.factor) %>%
  mutate(fbconf = recode_factor(fbconf, "fb" = 'Feedback', 
                                "conf" = 'Confidence',
                                'add' = 'Additive') ) %>%
  mutate(expNum = recode_factor(questname, "phq" = 'Exp 1', 
                                "gad" = 'Exp 1', "spin" = 'Exp 1', 'ad' = 'Exp 2'))%>%
  mutate(questname = recode_factor(questname, "phq" = 'PHQ', "ad" = 'AD')) %>%
  mutate(questname = factor(questname, levels = c('PHQ', 'AD')))

#### #### #### 
#### Figure 4F - Posterior distributions of Beta-confidence and Beta-feedback for model 8

ggplot(filter(posneg2, fbconf=='Confidence'), aes(x = expNum, y = beta, fill = expNum))+
  scale_fill_manual(breaks = c('Exp 1', 'Exp 2'),
                    values=c('darkorange1', 'darkorange1')) +
  geom_hline(yintercept = 0, color = 'grey50', size=1) +
  geom_violinhalf(alpha=.7)+
  theme_pubclean() +
  xlab("Mental health score") +
  ylab("β-confidence") +
  labs(fill = 'Dataset') +
  coord_cartesian(ylim=c(-.05,0.05))+
  scale_y_continuous(breaks=c(-.04,0.04))+
  theme(text = element_text(size=font_size-2), 
        plot.caption = element_text(size = font_size-2),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_blank()  ,
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none') 

ggplot(filter(posneg2, fbconf=='Feedback'), aes(x = expNum, y = beta, fill = expNum))+
  geom_hline(yintercept = 0, color = 'grey50', size=1) +
  scale_fill_manual(breaks = c('Exp 1', 'Exp 2'),
                    values=c('purple4', 'purple4')) +
  geom_violinhalf(alpha=.7)+
  theme_pubclean() +
  ylab("β-feedback") +
  xlab("Mental health score") +
  labs(fill = 'Dataset') +
  coord_cartesian(ylim=c(-.1,0.1))+
  scale_y_continuous(breaks = c(-.1, 0, .1))+
  # labs(tag = "E") +
  theme(text = element_text(size=font_size-2), 
        plot.caption = element_text(size = font_size-2),
        plot.tag = element_text(size = tag_size, face = "bold"),
        # axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none',
        axis.title.x=element_blank()  ) 


################
#################
#### Supplementary Analysis of all 3 transdiag axes

td <- mbsDataFit2[,,1]$model.12.q2[,,1]$samples[,,1]
ad <- c()
ad$beta.fb.lr <- td$beta.fb.lr[,,1]
ad$beta.conf.lr <- td$beta.conf.lr[,,1]
hphq <- hdi(ad$beta.conf.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)

cit <- c()
cit$beta.fb.lr <- td$beta.fb.lr[,,2]
cit$beta.conf.lr <- td$beta.conf.lr[,,2]
hphq <- hdi(cit$beta.fb.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(cit$beta.conf.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)

sw <- c()
sw$beta.fb.lr <- td$beta.fb.lr[,,3]
sw$beta.conf.lr <- td$beta.conf.lr[,,3]
hphq <- hdi(sw$beta.fb.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(sw$beta.conf.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)

posneg = data.frame(c(ad$beta.conf.lr),
                    c(cit$beta.conf.lr),
                    c(sw$beta.conf.lr))
colnames(posneg) <- c('ad_conf', 'cit_conf', 'sw_conf')

posneg2 <- posneg %>% pivot_longer(
  cols = c(1:3),
  names_sep = '_',
  values_to = 'beta',
  names_to = c('questname', 'fbconf')) %>% 
  mutate_at(c('fbconf', 'questname'), as.factor) %>%
  mutate(fbconf = recode_factor(fbconf, "fb" = 'Feedback', 
                                "conf" = 'Confidence',
                                'add' = 'Additive') ) %>%
  mutate(questname = recode_factor(questname, "ad" = 'AD', "cit" = 'CIT', 'sw' = 'SW')) %>%
  mutate(questname = factor(questname, levels = c('AD', 'CIT', 'SW')))

### Supplementary Figure 14A
ggplot(filter(posneg2, fbconf=='Confidence'), aes(x = questname, y = beta, fill = questname))+
  geom_hline(yintercept = 0, color = 'grey50', size=1) +
  geom_violinhalf(alpha=.7)+
  theme_pubclean() +
  xlab("Transdiagnostic axis") +
  ylab("β-confidence") +
  labs(fill = 'Dataset') +
  # coord_cartesian(ylim=c(-.05,0.025))+
  # labs(tag = "E") +
  theme(text = element_text(size=font_size-4), 
        plot.caption = element_text(size = font_size-4),
        plot.tag = element_text(size = tag_size, face = "bold"),
        # axis.title.x = element_blank()  ,
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none',
        # legend.direction = "horizontal"
) 


#####################
#####################
#### Plot symptom betas for all models and 3 distortions 
# conf
  gad <- mbsDataFit1[,,1]$model.6.q2
  gad <- gad[,,1]$samples
  gad <- gad[,,1]
  gad.fb.d1 <- gad$beta.fb.lr
  gad.conf.d1 <- gad$beta.conf.lr
  gad.bias.d1 <- gad$beta.post.bias
  
  gad <- mbsDataFit1[,,1]$model.5.q2
  gad <- gad[,,1]$samples
  gad <- gad[,,1]
  gad.fb.d2 <- gad$beta.fb.lr
  gad.conf.d2 <- gad$beta.conf.lr
  gad.bias.d2 <- gad$beta.post.bias
  
  gad <- mbsDataFit1[,,1]$model.7.q2
  gad <- gad[,,1]$samples
  gad <- gad[,,1]
  gad.fb.d3 <- gad$beta.fb.lr
  gad.conf.d3 <- gad$beta.conf.lr
  gad.bias.d3 <- gad$beta.post.bias
  
  gad <- mbsDataFit1[,,1]$model.8.q2
  gad <- gad[,,1]$samples
  gad <- gad[,,1]
  gad.fb.d4 <- gad$beta.fb.lr
  gad.conf.d4 <- gad$beta.conf.lr
  gad.bias.d4 <- gad$beta.post.bias
  
  gad <- mbsDataFit1[,,1]$model.10.q2
  gad <- gad[,,1]$samples
  gad <- gad[,,1]
  gad.fb.d5 <- gad$beta.fb.lr
  gad.conf.d5 <- gad$beta.conf.lr
  gad.bias.d5 <- gad$beta.post.bias
  
  gad <- mbsDataFit1[,,1]$model.9.q2
  gad <- gad[,,1]$samples
  gad <- gad[,,1]
  gad.fb.d6 <- gad$beta.fb.lr
  gad.conf.d6 <- gad$beta.conf.lr
  gad.bias.d6 <- gad$beta.post.bias
  
  posneg = data.frame(c(gad.fb.d1), c(gad.conf.d1), c(gad.bias.d1),
                      c(gad.fb.d2), c(gad.conf.d2), c(gad.bias.d2),
                      c(gad.fb.d3), c(gad.conf.d3), c(gad.bias.d3),
                      c(gad.fb.d4), c(gad.conf.d4), c(gad.bias.d4),
                      c(gad.fb.d5), c(gad.conf.d5), c(gad.bias.d5),
                      c(gad.fb.d6), c(gad.conf.d6), c(gad.bias.d6))
  colnames(posneg) <- c('gad_fb_d1', 'gad_conf_d1', 'gad_bias_d1',
                        'gad_fb_d2', 'gad_conf_d2', 'gad_bias_d2',
                        'gad_fb_d3', 'gad_conf_d3', 'gad_bias_d3',
                        'gad_fb_d4', 'gad_conf_d4', 'gad_bias_d4',
                        'gad_fb_d5', 'gad_conf_d5', 'gad_bias_d5',
                        'gad_fb_d6', 'gad_conf_d6', 'gad_bias_d6')

posneg2 <- posneg %>% pivot_longer(cols = c(1:18),
                                  names_sep = '_',
                                  values_to = 'beta',
                                  names_to = c('AD', 'distortion', 'model')) %>%
  mutate_at(c('distortion', 'AD', 'model'), as.factor)  %>%
  mutate(distortion = recode_factor(distortion, "fb" = 'Feedback',
                                "conf" = 'Confidence',
                                'bias' = 'Response bias') ) %>%
  mutate(model = recode_factor(model, "0" = 'D0: No distortion', 
                                  "d1" = 'D1: Confidence ',
                                  "d2" = 'D2: Feedback ',
                                  "d3" = 'D3: Resp. bias',
                                  "d4" = 'D4: Conf. + FB   ',
                                  "d5" = 'D5: Conf. +  Bias',
                                  "d6" = 'D6: FB + Bias') ) %>%
  mutate(model = factor(model, levels = c("0" = 'D0: No distortion', 
                                                "d1" = 'D1: Confidence ',
                                                "d2" = 'D2: Feedback ',
                                                "d3" = 'D3: Resp. bias',
                                                "d4" = 'D4: Conf. + FB   ',
                                                "d5" = 'D5: Conf. +  Bias',
                                                "d6" = 'D6: FB + Bias') ))
# %>%
#   mutate(expNum = recode_factor(questname, "phq" = 'Exp 1',
#                                 "gad" = 'Exp 1', "spin" = 'Exp 1', 'ad' = 'Exp 2'))%>%
#   mutate(questname = recode_factor(questname, "phq" = 'PHQ',
#                                    "gad" = 'GAD', "spin" = 'SPIN', "ad" = 'AD')) %>%
#   mutate(questname = factor(questname, levels = c('PHQ', 'GAD', 'SPIN', 'AD')))

posneg2$beta[posneg2$beta==0] <- NA


ggplot(posneg2 %>% filter(distortion=='Confidence'), 
       aes(x = model, y = beta, fill = distortion))+
  scale_fill_manual(breaks = c('Confidence'), values=c('darkorange1')) +
  geom_hline(yintercept = 0, color = 'grey50', size=1.5) +
  geom_violinhalf(alpha=.7)+
  theme_pubclean() +
  ylab("Beta") +
  coord_cartesian(ylim = c(-.04,.02))+
  xlab("") +
  labs(fill = 'Distortion type') +
  theme(text = element_text(size=font_size),
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(angle = 300, hjust = 0, size = font_size-6)) 

ggplot(posneg2 %>% filter(distortion=='Feedback'), 
       aes(x = model, y = beta, fill = distortion))+
  scale_fill_manual(breaks = c('Feedback'), values=c('purple4')) +
  geom_hline(yintercept = 0, color = 'grey50', size=1.5) +
  geom_violinhalf(alpha=.7)+
  theme_pubclean() +
  ylab("Beta") +
  # coord_cartesian(ylim = c(-2e-5,4e-5))+
  xlab("") +
  labs(fill = 'Distortion type') +
  theme(text = element_text(size=font_size),
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(angle = 300, hjust = 0, size = font_size-6))

ggplot(posneg2 %>% filter(distortion=='Response bias'), 
       aes(x = model, y = beta, fill = distortion))+
  scale_fill_manual(breaks = c('Response bias'), values=c('maroon4')) +
  geom_hline(yintercept = 0, color = 'grey50', size=1.5) +
  geom_violinhalf(alpha=.7)+
  theme_pubclean() +
  ylab("Beta") +
  coord_cartesian(ylim = c(-1e-5,5e-5))+
  xlab("") +
  labs(fill = 'Distortion type') +
  theme(text = element_text(size=font_size-2),
        plot.caption = element_text(size = font_size-2),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(angle = 300, hjust = 0, size = font_size-6)) 


# conf
ad <- mbsDataFit2[,,1]$model.6.q1
ad <- ad[,,1]$samples
ad <- ad[,,1]
ad.fb.d1 <- ad$beta.fb.lr
ad.conf.d1 <- ad$beta.conf.lr
ad.bias.d1 <- ad$beta.post.bias

ad <- mbsDataFit2[,,1]$model.5.q1
ad <- ad[,,1]$samples
ad <- ad[,,1]
ad.fb.d2 <- ad$beta.fb.lr
ad.conf.d2 <- ad$beta.conf.lr
ad.bias.d2 <- ad$beta.post.bias

ad <- mbsDataFit2[,,1]$model.7.q1
ad <- ad[,,1]$samples
ad <- ad[,,1]
ad.fb.d3 <- ad$beta.fb.lr
ad.conf.d3 <- ad$beta.conf.lr
ad.bias.d3 <- ad$beta.post.bias

ad <- mbsDataFit2[,,1]$model.8.q1
ad <- ad[,,1]$samples
ad <- ad[,,1]
ad.fb.d4 <- ad$beta.fb.lr
ad.conf.d4 <- ad$beta.conf.lr
ad.bias.d4 <- ad$beta.post.bias

ad <- mbsDataFit2[,,1]$model.10.q1
ad <- ad[,,1]$samples
ad <- ad[,,1]
ad.fb.d5 <- ad$beta.fb.lr
ad.conf.d5 <- ad$beta.conf.lr
ad.bias.d5 <- ad$beta.post.bias

ad <- mbsDataFit2[,,1]$model.9.q1
ad <- ad[,,1]$samples
ad <- ad[,,1]
ad.fb.d6 <- ad$beta.fb.lr
ad.conf.d6 <- ad$beta.conf.lr
ad.bias.d6 <- ad$beta.post.bias

posneg = data.frame(c(ad.fb.d1), c(ad.conf.d1), c(ad.bias.d1),
                    c(ad.fb.d2), c(ad.conf.d2), c(ad.bias.d2),
                    c(ad.fb.d3), c(ad.conf.d3), c(ad.bias.d3),
                    c(ad.fb.d4), c(ad.conf.d4), c(ad.bias.d4),
                    c(ad.fb.d5), c(ad.conf.d5), c(ad.bias.d5),
                    c(ad.fb.d6), c(ad.conf.d6), c(ad.bias.d6))
colnames(posneg) <- c('ad_fb_d1', 'ad_conf_d1', 'ad_bias_d1',
                      'ad_fb_d2', 'ad_conf_d2', 'ad_bias_d2',
                      'ad_fb_d3', 'ad_conf_d3', 'ad_bias_d3',
                      'ad_fb_d4', 'ad_conf_d4', 'ad_bias_d4',
                      'ad_fb_d5', 'ad_conf_d5', 'ad_bias_d5',
                      'ad_fb_d6', 'ad_conf_d6', 'ad_bias_d6')

posneg2 <- posneg %>% pivot_longer(cols = c(1:18),
                                   names_sep = '_',
                                   values_to = 'beta',
                                   names_to = c('AD', 'distortion', 'model')) %>%
  mutate_at(c('distortion', 'AD', 'model'), as.factor)  %>%
  mutate(distortion = recode_factor(distortion, "fb" = 'Feedback',
                                    "conf" = 'Confidence',
                                    'bias' = 'Response bias') ) %>%
  mutate(model = recode_factor(model, "0" = 'D0: No distortion', 
                               "d1" = 'D1: Confidence ',
                               "d2" = 'D2: Feedback ',
                               "d3" = 'D3: Resp. bias',
                               "d4" = 'D4: Conf. + FB',
                               "d5" = 'D5: Conf. + Bias',
                               "d6" = 'D6: FB + Bias') ) %>%
  mutate(model = factor(model, levels = c("0" = 'D0: No distortion', 
                                          "d1" = 'D1: Confidence ',
                                          "d2" = 'D2: Feedback ',
                                          "d3" = 'D3: Resp. bias',
                                          "d4" = 'D4: Conf. + FB',
                                          "d5" = 'D5: Conf. + Bias',
                                          "d6" = 'D6: FB + Bias') ))
# %>%
#   mutate(expNum = recode_factor(questname, "ad" = 'Exp 1',
#                                 "gad" = 'Exp 1', "spin" = 'Exp 1', 'ad' = 'Exp 2'))%>%
#   mutate(questname = recode_factor(questname, "ad" = 'ad',
#                                    "gad" = 'GAD', "spin" = 'SPIN', "ad" = 'AD')) %>%
#   mutate(questname = factor(questname, levels = c('ad', 'GAD', 'SPIN', 'AD')))

posneg2$beta[posneg2$beta==0] <- NA


ggplot(posneg2 %>% filter(distortion=='Confidence'), 
       aes(x = model, y = beta, fill = distortion))+
  scale_fill_manual(breaks = c('Confidence'), values=c('darkorange1')) +
  geom_hline(yintercept = 0, color = 'grey50', size=1.5) +
  geom_violinhalf(alpha=.7)+
  theme_pubclean() +
  ylab("Beta") +
  coord_cartesian(ylim = c(-.1,.03))+
  xlab("") +
  labs(fill = 'Distortion type') +
  theme(text = element_text(size=font_size),
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(angle = 300, hjust = 0, size = font_size-6)) 

ggplot(posneg2 %>% filter(distortion=='Feedback'), 
       aes(x = model, y = beta, fill = distortion))+
  scale_fill_manual(breaks = c('Feedback'), values=c('purple4')) +
  geom_hline(yintercept = 0, color = 'grey50', size=1.5) +
  geom_violinhalf(alpha=.7)+
  theme_pubclean() +
  ylab("Beta") +
  # coord_cartesian(ylim = c(-2e-5,4e-5))+
  xlab("") +
  labs(fill = 'Distortion type') +
  theme(text = element_text(size=font_size),
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(angle = 300, hjust = 0, size = font_size-6))

ggplot(posneg2 %>% filter(distortion=='Response bias'), 
       aes(x = model, y = beta, fill = distortion))+
  scale_fill_manual(breaks = c('Response bias'), values=c('maroon4')) +
  geom_hline(yintercept = 0, color = 'grey50', size=1.5) +
  geom_violinhalf(alpha=.7)+
  theme_pubclean() +
  ylab("Beta") +
  coord_cartesian(ylim = c(-2e-4,8e-4))+
  xlab("") +
  labs(fill = 'Distortion type') +
  theme(text = element_text(size=font_size-2),
        plot.caption = element_text(size = font_size-2),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(angle = 300, hjust = 0, size = font_size-6)) 

#########
### plot model fits to empirical data

spe <- mbsDataExp1[,,1]$spe
nsubjrun <- dim(spe)
# speFit <- mbsDataExp1[,,1]$fitZ[,,1]$model.4[,,1]$spe.est
speFit <- mbsDataExp1[,,1]$fitRegZ[,,1]$model.6[,,1]$spe.est

groups <- t(mbsDataExp1[,,1]$groups)
groups <- kronecker(matrix(1,1,nsubjrun[2]),groups)

subj <- kronecker(matrix(1,1,nsubjrun[2]), c(1:nsubjrun[1]))
runnum <- t(kronecker(matrix(1,1,nsubjrun[1]), c(1:nsubjrun[2])))

df <- data.frame(c(spe), c(speFit), c(groups), c(subj), c(runnum))
colnames(df) <- c('spe.emp', 'spe.fit', 'group', 'subj', 'runnum')
df <- df %>% pivot_longer(cols = starts_with("spe."),
                          names_to = "emporfit",
                          names_prefix = "spe.",
                          values_to = "spe") %>%
  mutate_at(c("group", "subj", "runnum"), as.factor) %>%
  mutate(firstFeedback = recode_factor(group, "1" = 'Positive feedback first', 
                                       "2" = 'Negative feedback first', "3" = 'Positive feedback first', "4" = 'Negative feedback first', 
                                       "5" = 'Positive feedback first', "6" = 'Negative feedback first', "7" = 'Positive feedback first', 
                                       "8" = 'Negative feedback first')) %>%
  mutate(group = recode_factor(group, "1" = 'Group 1', 
                               "2" = 'Group 2', "3" = 'Group 3', "4" = 'Group 4', 
                               "5" = 'Group 5', "6" = 'Group 6', "7" = 'Group 7', 
                               "8" = 'Group 8')) 


### Supplementary Figure 13 - fits for each group separately
ggplot(df, aes(x=runnum, y=spe, colour=emporfit, group=emporfit)) + 
  scale_color_manual(breaks = c('emp', 'fit'),
                     values=c('black', 'darkorchid1'),
                     labels = c('Observed data', 'Fit')) +
  facet_wrap(~factor(group), nrow = 4) + theme_pubclean() +
  stat_summary(fun.data = mean_se,na.rm = T) +
  stat_summary(fun.data = mean_se,na.rm = T, geom = c('line'),size = 1.25,
               aes(linetype = emporfit), show.legend = F) +
  scale_linetype_manual(breaks = c('emp', 'fit'), values = c('solid', 'twodash'))  +
  ylab("SPE") +
  xlab("Block number") +
  labs(colour = '') + ggtitle('Exp 1') +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        plot.title = element_text(hjust = .5),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

### fits separated by whether pos or neg feedback came first
ggplot(df, aes(x=runnum, y=spe, colour=emporfit, group=emporfit)) + 
  scale_color_manual(breaks = c('emp', 'fit'),
                     values=c('black', 'darkorchid1'),
                     labels = c('Observed', 'Fitted')) +
  facet_wrap(~factor(firstFeedback), nrow = 2) + theme_pubclean() +
  stat_summary(fun.data = mean_se,na.rm = T, 
               position=position_dodge(width = .1),
               geom = c('errorbar'), width = .3, size = 1.2) +
  stat_summary(fun.data = mean_se,na.rm = T, 
               position=position_dodge(width = .1),
               geom = c('line'), size = 1.25,
               aes(linetype = emporfit), show.legend = F) +
  scale_linetype_manual(breaks = c('emp', 'fit'), values = c('solid', 'twodash'))  +
  ylab("SPE") +
  xlab("Block number") +
  labs(colour = 'Data') + 
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        # legend.position = "right",
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

####################################
####################################
####################################
#### Plot model recovery

setwd('~/OneDrive - University of Copenhagen/Projects/Experiments/metaBiasShift/data/')
mbsRec1 <- readMat('paramRecovery_fbconf.mat') # load questionnaire data
mbsRec2 <- readMat('paramRecovery_postbias.mat') # load questionnaire data

nsim <- length(c(mbsRec1$fit.lr.conf))
sim.df <- data.frame(c(mbsRec1$fit.blr.fb, mbsRec1$fit.blr.conf), #
                     c(mbsRec1$sim.blr.fb,mbsRec1$sim.blr.conf),
                     c(array('β-feedback', nsim), array('β-confidence', nsim)))
colnames(sim.df) <- c('Fitted', 'Simulated', 'fbconf')

####################################
### Supplementary Figure 5
ggplot(sim.df, aes(x = Simulated, y = Fitted)) +
  facet_wrap(~factor(fbconf, levels = c('β-feedback', 'β-confidence')), ncol=2) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  theme_pubclean()+
  ylab("Fitted") +
  xlab("Simulated") +
  labs(colour = 'Block number') + 
  # labs(tag = "B") +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        plot.title = element_text(hjust = .5),
        panel.spacing.x = unit(1.5, "lines"))

nsim <- length(mbsRec2$fit.bpostbias)

sim.df <- data.frame(c(mbsRec2$fit.bpostbias[,,2,2,]), # recovered using model 7 with only additive dist. so de-select the other simulated beta parameters
                     c(mbsRec2$sim.bpostbias[,,2,2,]),
                     c(array('β-additive', nsim)))
colnames(sim.df) <- c('Fitted', 'Simulated', 'fbconf')

ggplot(sim.df, aes(x = Simulated, y = Fitted)) +
  facet_wrap(~factor(fbconf), ncol=2) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  theme_pubclean()+
  ylab("Fitted") +
  xlab("Simulated") +
  labs(colour = 'Block number') + 
  # labs(tag = "B") +
  scale_y_continuous(breaks = seq(-.005,.005,.0025)) +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        plot.title = element_text(hjust = .5),
        panel.spacing.x = unit(1.5, "lines"),
        plot.margin = margin(r=15)
        # axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
        )



########################################
##############################
##### Plot model simulations
##### The following code is used to generate Supplementary Figure 8. 
##### To generate the figure, the code needs to be run once each by setting the model2Plot
##### variable to values between 1 and 4
###############

setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/metaBiasShift/data/simulation/")

### model 1 - feedback distortion (for Supp Figure 8B)
### model 2 - confidence distortion (for Supp Figure 8A)
### model 3 - no distortion (for Supp Figure 8D)
### model 4 - additive distortion (for Supp Figure 8c)

### model 5 - feedback + confidence distortions
### model 6 - additive distortion with empirical values (for Supplementary Figure 11C)
### model 7 - additive distortion with empirical values (for Supplementary Figure 11C)

model2Plot <- 4

mod <- switch(model2Plot, 3,4,5,6,7,9,17)

simExp = readMat(paste('simExp_m',mod,'.mat', sep='')) # load questionnaire data
spe <- simExp$spe
conf <- simExp$conf
phq <- simExp$phq
nsubjrun <- dim(conf)
fbblock <- simExp$fbblock
task <- simExp$task
fb <- simExp$feedback
group <- simExp$group
accu <- simExp$correct

subj <- array(c(1:nsubjrun[1]), dim = nsubjrun)
runnum <- aperm(array(c(1:nsubjrun[2]), dim = c(nsubjrun[2], nsubjrun[1], nsubjrun[3])), c(2,1,3))
phq <- array(phq, dim = nsubjrun)
spe <- array(spe, dim = nsubjrun)
fbblock <- array(fbblock, dim = nsubjrun)
task <- array(task, dim = nsubjrun)
trialnum <- aperm(array(c(1:nsubjrun[3]), dim = c(nsubjrun[3], nsubjrun[1], nsubjrun[2])), c(2,3,1))
firstFeedback <- array(c(fbblock[,3,1]), dim = nsubjrun)

exp1.sim <- data.frame(c(spe), c(phq), c(runnum), c(subj), c(fbblock), c(task), 
                       c(conf), c(trialnum), c(fb), c(firstFeedback), c(accu))
colnames(exp1.sim) <- c('spe', 'phq', 'runnum', 'subj', 'fbblock', 'task', 
                        'conf', 'trialnum', 'feedback', 'firstFeedback', 'accu')

exp1.sim <- exp1.sim %>%
  mutate_at(c('subj', 'fbblock', 'task', 'firstFeedback'), as.factor) %>%
  mutate(firstFeedback = recode_factor(firstFeedback, '1' = 'Positive',
                                       '2' = 'Negative'))

##########################################################
#### regression model predictions from computational models

exp1.sim$confZ <- exp1.sim$conf
subj <- unique(exp1.sim$subj)
for (ns in subj){
  inds <- exp1.sim$subj==ns
  exp1.sim$confZ[inds] <- scale(exp1.sim$conf[inds])
}

sim.model <- exp1.sim %>% 
  dplyr::select(c(subj, runnum, trialnum, spe, confZ, phq)) %>%
  group_by(subj,trialnum) %>%
  mutate(spe.tile = ntile(spe, 6))%>%
  mutate(conf.tile = ntile(confZ, 4)) %>%
  mutate(confhilo = confZ > 0) %>%
  group_by(runnum) %>%
  mutate(AD.tile = ntile(phq, 3)) %>%
  mutate_at(c('AD.tile'), as.factor)
levels(sim.model$AD.tile) <- c(levels(sim.model$AD.tile), 'High', 'Low')
sim.model <- sim.model %>%
  mutate(AD.tile = recode_factor(AD.tile, '1' = 'Low', '3' = 'High'))
sim.model$conftile.hilo <- as.factor(ceiling(sim.model$conf.tile/2))

## test

sim.obs.lo <- lmer(spe ~ phq*confZ + (1|subj) + (1|runnum), 
                        filter(sim.model, confhilo==F))
summary(sim.obs.lo)
# plot_model(sim.obs.lo,type = 'pred', terms= c('confZ', 'phq'))

sim.obs.hi <- lmer(spe ~ phq*confZ + (1|subj) + (1|runnum), 
                     filter(sim.model, confhilo==T))
summary(sim.obs.hi)
# plot_model(sim.obs.hi,type = 'pred', terms= c('confZ', 'phq'))


AD.tile.cents  <- c(mean(sim.model$phq[sim.model$AD.tile=='Low']),
                    mean(sim.model$phq[sim.model$AD.tile=='High']))

conf.tile.cents <- c(mean(sim.model$confZ[sim.model$conf.tile==1]),
                     mean(sim.model$confZ[sim.model$conf.tile==2]))

emm.obs.lo <- summary(emmeans(sim.obs.lo, ~ confZ|phq, 
                              at = list(confZ = conf.tile.cents, phq = AD.tile.cents)))
emm.obs.lo

conf.tile.cents <- c(mean(exp1.model$confZ[exp1.model$conf.tile==3]),
                     mean(exp1.model$confZ[exp1.model$conf.tile==4]))

emm.obs.hi <- summary(emmeans(sim.obs.hi, ~ confZ|phq, 
                           at = list(confZ = conf.tile.cents, phq = AD.tile.cents)))
emm.obs.hi

spe <- c(emm.obs.lo$emmean,   emm.obs.hi$emmean)
if (is.null(emm.obs.hi$asymp.LCL)){
  emm.obs.hi$asymp.LCL <- emm.obs.hi$lower.CL
  emm.obs.hi$asymp.UCL <- emm.obs.hi$upper.CL
}
if (is.null(emm.obs.lo$asymp.LCL)){
  emm.obs.lo$asymp.LCL <- emm.obs.lo$lower.CL
  emm.obs.lo$asymp.UCL <- emm.obs.lo$upper.CL
}
lci <- c(emm.obs.lo$asymp.LCL, emm.obs.hi$asymp.LCL)
uci <- c(emm.obs.lo$asymp.UCL, emm.obs.hi$asymp.UCL)

conf.tile <- c(c(1:2,1:2), c(3:4,3:4))
AD.tile <- rep(c(rep('Low',2), rep('High',2)),2)
conftile.hilo <- c(rep(1,4), rep(2,4))
sim.model.df <- data.frame(spe, lci, uci, 
                            conf.tile, AD.tile, conftile.hilo)
sim.model.df$conftile.hilo <- as.factor(sim.model.df$conftile.hilo)

ggplot(filter(sim.model.df),
       aes(x = conf.tile, y = spe, colour = AD.tile,
           linetype = conftile.hilo)) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', 'royalblue2')) +
  scale_shape_manual(values = c(19,2)) +
  geom_point(position=position_dodge(.4), size=3.2, stroke=.5) +
  geom_errorbar(position=position_dodge(.4), aes(ymin = lci, ymax = uci), 
                linetype = 1, linewidth=.6, width=.3) +
  stat_summary(position=position_dodge(.4), geom="line") +
  scale_linetype_manual(values = c('dashed', 'solid')) +
  theme_classic2() +
  coord_cartesian(ylim = c(.1,.4)) +
  labs(color = 'AD score', linetype = 'Confidence level')+
  xlab('Confidence (tiled)') +
  ylab('Global SPE') +
  theme(text = element_text(size=font_size-2),
        legend.direction = 'vertical',
        # axis.title.y = element_text(hjust=.9),
        legend.text=element_text(size=font_size-6),
        plot.caption = element_text(size = font_size-2),
        legend.position = 'none') +
  geom_segment(aes(x = 3.3, y = .29, xend = 3.3, yend = .35), color = 'grey60') +
  # geom_segment(aes(x = 3.3, y = .3, xend = 3.5, yend = .3), color = 'grey60') +
  # geom_segment(aes(x = 3.3, y = .29, xend = 3.3, yend = .35), color = 'grey60') +
  geom_segment(aes(x = 3.3, y = .32, xend = 3.5, yend = .32), color = 'grey60') +
  annotate('text', label = 'n.s.', x = 3.78, y = .33, size=4.5)

#### orig
sim.reg <- lmer(confZ ~ spe*phq + (1|subj), filter(sim.model, confhilo==T))
summary(sim.reg)
# plot_model(sim.reg,type = 'pred', terms= c('spe', 'phq'))


# ggplot(filter(sim.model, AD.tile!=2),
#        aes(x = conf.tile, y = spe, colour = AD.tile)) +
#   scale_color_manual(breaks = c('Low', 'High'),
#                      values = c('tan3', 'royalblue2')) +
#   # stat_summary(geom='line', aes(linetype = conftile.hilo), size=1) +
#   geom_smooth(method='lm', aes(linetype = conftile.hilo), se = F) +
#   scale_linetype_manual(values = c('dotted', 'solid')) +
#   stat_summary(fun.data = mean_cl_boot, size=.5, shape = 1) +
#   stat_summary(fun.y = mean, size=.5, alpha = .5) +
#   labs(color = 'AD score', linetype = 'Confidence level')+
#   theme_classic2() +
#   xlab('Confidence (tiled)') +
#   ylab('Global SPE') +
#   coord_cartesian(ylim = c(.1,.4)) +
#   theme(text = element_text(size=font_size-2),
#         legend.direction = 'vertical',
#         # axis.title.y = element_text(hjust=.9),
#         legend.text=element_text(size=font_size-6),
#         plot.caption = element_text(size = font_size-2),
#         legend.position = 'none') #+
#   # geom_segment(aes(x = 3.3, y = .25, xend = 3.3, yend = .35), color = 'grey60') +
#   # geom_segment(aes(x = 3.3, y = .3, xend = 3.5, yend = .3), color = 'grey60') +
#   # geom_segment(aes(x = 3.3, y = .29, xend = 3.3, yend = .35), color = 'grey60') +
#   # geom_segment(aes(x = 3.3, y = .32, xend = 3.5, yend = .32), color = 'grey60') +
#   # annotate('text', label = 'n.s.', x = 3.75, y = .32, size=4.5)



sim.model <- exp1.sim %>% group_by(subj, fbblock, task, runnum) %>%
  summarise(spe = mean(spe),
            phq = mean(phq)) %>%
  # mutate(speZ = c(scale(spe))) %>%
  group_by(runnum) %>%
  mutate(AD.tile = ntile(phq, 3)) %>%
  mutate_at(c('AD.tile'), as.factor)

levels(sim.model$fbblock) <- c(levels(sim.model$fbblock), 'None', 'Positive', 'Negative')
levels(sim.model$AD.tile) <- c(levels(sim.model$AD.tile), 'High', 'Low')
sim.model <- sim.model %>%
  mutate(AD.tile = recode_factor(AD.tile, '1' = 'Low', '3' = 'High'))  %>%
  mutate(fbblock = recode_factor(fbblock, '0' = 'None', '1' = 'Positive',
                                 '2' = 'Negative'))

AD.tile.cents  <- c(mean(sim.model$phq[sim.model$AD.tile=='Low']),
                    mean(sim.model$phq[sim.model$AD.tile=='High']))

exp1.reg <- lmer(spe ~ fbblock*phq + (1|subj) + (1|runnum), sim.model)
summary(exp1.reg)
# plot(exp1.reg)

ann1 <- c('sig.', 'sig.', 'n.s.', 'n.s.')
ann2 <- c('sig.', 'sig.', 'n.s.', 'n.s.')
# ggplot(filter(sim.model %>% mutate(fbblock = factor(
#   fbblock, levels = c ( 'Negative','None', 'Positive'))), AD.tile!=2),
#   aes(x = fbblock, y = spe, colour = AD.tile)) +
#   scale_color_manual(breaks = c('Low', 'High'),
#                      values = c('tan3', '#0066cc')) +
#   stat_summary(fun.y = mean, size=.8, alpha = .3,
#                position = position_dodge(.5)) +
#   stat_summary(fun.data = mean_cl_boot, size=.8, shape = 1,
#                position = position_dodge(.5)) +
#   xlab('Feedback type') +
#   ylab('Global SPE') +
#   theme_classic2() +
#   coord_cartesian(ylim = c(.2, .5)) +
#   theme(text = element_text(size=font_size-2),
#         # axis.title.y = element_text(hjust=.9),
#         plot.caption = element_text(size = font_size-2),
#         legend.position = 'none') +
#   # geom_segment(aes(x = 2, y = .43, xend = 3, yend = .43), color = 'grey60') +
#   # geom_segment(aes(x = 2, y = .4, xend = 2, yend = .43), color = 'grey60') +
#   # geom_segment(aes(x = 3, y = .43, xend = 3, yend = .45), color = 'grey60') +
#   # annotate('text', label = ann2[model2Plot], x = 2.5, y = .46, size=5) +
#   # geom_segment(aes(x = 1, y = .32, xend = 2, yend = .32), color = 'grey60') +
#   # geom_segment(aes(x = 2, y = .32, xend = 2, yend = .35), color = 'grey60') +
#   # geom_segment(aes(x = 1, y = .3, xend = 1, yend = .32), color = 'grey60') +
#   # annotate('text', label = ann1[model2Plot], x = 1.5, y = .345, size=5)+
#   geom_segment(aes(x = 2, y = .37, xend = 3, yend = .37), color = 'grey60') +
#   geom_segment(aes(x = 2, y = .32, xend = 2, yend = .37), color = 'grey60') +
#   geom_segment(aes(x = 3, y = .37, xend = 3, yend = .4), color = 'grey60') +
#   annotate('text', label = ann2[model2Plot], x = 2.5, y = .39, size=5) +
#   geom_segment(aes(x = 1, y = .27, xend = 2, yend = .27), color = 'grey60') +
#   geom_segment(aes(x = 2, y = .27, xend = 2, yend = .29), color = 'grey60') +
#   geom_segment(aes(x = 1, y = .24, xend = 1, yend = .27), color = 'grey60') +
#   annotate('text', label = ann1[model2Plot], x = 1.5, y = .29, size=5)

sim.fb <- lmer(spe ~ fbblock*phq + (1|subj) + (1|runnum), sim.model)
summary(sim.fb)
emm.fb <- summary(emmeans(sim.fb, ~ fbblock|phq, at = list(phq = AD.tile.cents)))

spe <- emm.fb$emmean
lci <- emm.fb$lower.CL
uci <- emm.fb$upper.CL

fb.type <- rep(c('None', 'Positive', 'Negative'),2)
AD.tile <- c(rep('Low',3), rep('High',3))
sim.df <- data.frame(spe, lci, uci, fb.type, AD.tile)
sim.df$AD.tile <- as.factor(sim.df$AD.tile)
sim.df <- sim.df %>% mutate(AD.tile = factor(AD.tile, levels = c('Low','High'))) 


ggplot(sim.df, aes(x = fb.type, y = spe, colour = AD.tile)) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', 'royalblue2')) +
  geom_point(position=position_dodge(.5), size=3.2, stroke=.5) +
  geom_errorbar(aes(ymin = lci, ymax = uci), linetype = 1,
                position=position_dodge(.5), linewidth=.6, width=.3) +
  theme_classic2() +
  labs(color = 'AD score', linetype = 'Confidence level')+
  xlab('Confidence (tiled)') +
  ylab('Global SPE') +
  scale_y_continuous(breaks = seq(.1,.5,.1)) +
  coord_cartesian(ylim = c(.15, .5)) +
  theme(text = element_text(size=font_size-4),
        legend.direction = 'vertical',
        legend.text=element_text(size=font_size-6),
        plot.caption = element_text(size = font_size-2),
        legend.position = 'none')  +
  # geom_segment(aes(x = 2, y = .43, xend = 3, yend = .43), color = 'grey60') +
  # geom_segment(aes(x = 2, y = .4, xend = 2, yend = .43), color = 'grey60') +
  # geom_segment(aes(x = 3, y = .43, xend = 3, yend = .45), color = 'grey60') +
  # annotate('text', label = ann2[model2Plot], x = 2.5, y = .46, size=5) +
  # geom_segment(aes(x = 1, y = .32, xend = 2, yend = .32), color = 'grey60') +
  # geom_segment(aes(x = 2, y = .32, xend = 2, yend = .35), color = 'grey60') +
  # geom_segment(aes(x = 1, y = .3, xend = 1, yend = .32), color = 'grey60') +
  # annotate('text', label = ann1[model2Plot], x = 1.5, y = .345, size=5)
  geom_segment(aes(x = 2, y = .36, xend = 3, yend = .36), color = 'grey60') +
  geom_segment(aes(x = 2, y = .32, xend = 2, yend = .36), color = 'grey60') +
  geom_segment(aes(x = 3, y = .36, xend = 3, yend = .4), color = 'grey60') +
  annotate('text', label = ann2[model2Plot], x = 2.5, y = .39, size=5) +
  geom_segment(aes(x = 1, y = .27, xend = 2, yend = .27), color = 'grey60') +
  geom_segment(aes(x = 2, y = .27, xend = 2, yend = .29), color = 'grey60') +
  geom_segment(aes(x = 1, y = .24, xend = 1, yend = .27), color = 'grey60') +
  annotate('text', label = ann1[model2Plot], x = 1.5, y = .3, size=5)
  # geom_segment(aes(x = 2.1, y = .59, xend = 2.9, yend = .59), color = 'grey60') +
  # geom_segment(aes(x = 2.1, y = .55, xend = 2.1, yend = .59), color = 'grey60') +
  # geom_segment(aes(x = 2.9, y = .59, xend = 2.9, yend = .63), color = 'grey60') +
  # annotate('text', label = ann2[model2Plot], x = 2.5, y = .37, size=5) +
  # geom_segment(aes(x = 1.1, y = .48, xend = 1.9, yend = .48), color = 'grey60') +
  # geom_segment(aes(x = 1.9, y = .48, xend = 1.9, yend = .52), color = 'grey60') +
  # geom_segment(aes(x = 1.1, y = .44, xend = 1.1, yend = .48), color = 'grey60') +
  # annotate('text', label = ann1[model2Plot], x = 1.5, y = .27, size=5)
