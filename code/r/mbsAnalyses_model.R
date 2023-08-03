
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
data.folder <- "/Users/skatyal/OneDrive - University College London/Projects/Experiments/metaBiasShift/data/exp"

expNum <- 1
setwd(paste(data.folder, expNum, sep = ''))
mbsDataExp1 = readMat(paste('mbsDataExp',expNum,'.mat',sep='')) # load questionnaire data
mbsDataExp1 <- mbsDataExp1[[paste('mbsDataExp',expNum,sep='')]]

expNum <- 2
setwd(paste(data.folder, expNum, sep = ''))
mbsDataExp2 = readMat(paste('mbsDataExp',expNum,'.mat',sep='')) # load questionnaire data
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
dic[1,] <- dic[1,] - min(dic[1,])

mbsDataFit <- mbsDataExp2[,,1]$fitZ # select this one to plot exp 2
dic[2,1] <-  mbsDataFit[,,1]$model.0[,,1]$dic
dic[2,2] <-  mbsDataFit[,,1]$model.1[,,1]$dic
dic[2,3] <-  mbsDataFit[,,1]$model.2[,,1]$dic
dic[2,4] <-  mbsDataFit[,,1]$model.3[,,1]$dic
dic[2,5] <-  mbsDataFit[,,1]$model.4[,,1]$dic
dic[2,] <- dic[2,] - min(dic[2,])

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
        text = element_text(size=font_size), 
        # legend.box.just = 'left',
        plot.caption = element_text(size = font_size),
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
dic[1,2] <-  mbsDataFit[,,1]$model.5.q1[,,1]$dic
dic[1,3] <-  mbsDataFit[,,1]$model.6.q1[,,1]$dic
dic[1,4] <-  mbsDataFit[,,1]$model.7.q1[,,1]$dic
dic[1,5] <-  mbsDataFit[,,1]$model.8.q1[,,1]$dic
dic[1,6] <-  mbsDataFit[,,1]$model.9.q1[,,1]$dic
dic[1,7] <-  mbsDataFit[,,1]$model.10.q1[,,1]$dic
dic[1,] <- dic[1,] - min(dic[1,])

mbsDataFit <- mbsDataExp2[,,1]$fitRegZ # select this one to plot exp 2
dic[2,1] <-  mbsDataExp2[,,1]$fitZ[,,1]$model.4[,,1]$dic
dic[2,2] <-  mbsDataFit[,,1]$model.5.q1[,,1]$dic
dic[2,3] <-  mbsDataFit[,,1]$model.6.q1[,,1]$dic
dic[2,4] <-  mbsDataFit[,,1]$model.7.q1[,,1]$dic
dic[2,5] <-  mbsDataFit[,,1]$model.8.q1[,,1]$dic
dic[2,6] <-  mbsDataFit[,,1]$model.9.q1[,,1]$dic
dic[2,7] <-  mbsDataFit[,,1]$model.10.q1[,,1]$dic
dic[2,] <- dic[2,] - min(dic[2,])

dfdic <- data.frame(c(t(dic)), c(rep('Exp 1',n.models), rep('Exp 2', n.models)), c(seq(n.models)-1, seq(n.models)-1))
colnames(dfdic) <- c('dic', 'expNum', 'modelNum')
dfdic <- dfdic %>%
  mutate_at(c('expNum', 'modelNum'), as.factor) %>%
  mutate(modelNum = recode_factor(modelNum, "0" = 'D0: No dist.', 
                                  "1" = 'D1: FB  ',
                                  "2" = 'D2: Conf.  ',
                                  "3" = 'D3: Add.   ',
                                  "4" = 'D4: FB +\n    Conf.',
                                  "5" = 'D5: FB +\n    Add. ',
                                  "6" = 'D6: Conf. +\n    Add.',
                                  "7" = 'D7: FB + Conf.\n  + Add.') )

#### Figure 4E - Model comparison of distortion models

ggplot(dfdic, aes(y=modelNum,x=dic,fill=expNum)) +
  scale_fill_manual(breaks = c('Exp 1', 'Exp 2'),
                    values=c('gray75', 'gray15')) +
  # values=c('darkorange1', 'purple4')) +
  geom_vline(xintercept = 0, color = 'grey20') +
  geom_bar(stat="identity", width=.7, position = "dodge") +
  theme_pubclean() +
  xlab("DIC") +
  ylab("") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        legend.position = c(.7,.95),
        plot.margin = margin(t=35, r = 15),
        legend.title = element_blank(),
        text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0))+ 
  guides(fill = guide_legend(nrow = 2))

##################################################################
######## Read in the fitted models to plot posteriors

ci.lvl = .99

mbsDataFit1 <- mbsDataExp1[,,1]$fitRegZ
mbsDataFit2 <- mbsDataExp2[,,1]$fitRegZ

phq <- mbsDataFit1[,,1]$model.10.q1 # conf + add distortions 
phq <- phq[,,1]$samples
phq <- phq[,,1]

hphq <- hdi(phq$beta.fb.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(phq$beta.conf.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(phq$beta.post.bias,ci.lvl)
c(hphq$CI_low, hphq$CI_high)

gad <- mbsDataFit1[,,1]$model.6.q2
gad <- gad[,,1]$samples
gad <- gad[,,1]
hphq <- hdi(gad$beta.fb.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(gad$beta.conf.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)

spin <- mbsDataFit1[,,1]$model.6.q3
spin <- spin[,,1]$samples
spin <- spin[,,1]
hphq <- hdi(spin$beta.fb.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)
hphq <- hdi(spin$beta.conf.lr,ci.lvl)
c(hphq$CI_low, hphq$CI_high)

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

posneg = data.frame(c(phq$beta.fb.lr),
                    c(phq$beta.conf.lr),
                    c(phq$beta.post.bias),
                    c(ad$beta.fb.lr)/5,
                    c(ad$beta.conf.lr)/5,
                    c(ad$beta.post.bias)/5)
colnames(posneg) <- c('phq_fb', 'phq_conf', 'phq_add',
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
  coord_cartesian(ylim=c(-.05,0.025))+
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

phq <- mbsDataFit1[,,1]$model.8.q1
phq <- phq[,,1]$samples
phq <- phq[,,1]

td <- mbsDataFit2[,,1]$model.8.q1[,,1]$samples[,,1]
ad <- c()
ad$beta.fb.lr <- td$beta.fb.lr
ad$beta.conf.lr <- td$beta.conf.lr
ad$beta.post.bias <- td$beta.post.bias

posneg = data.frame(c(phq$beta.fb.lr),
                    c(phq$beta.conf.lr),
                    c(phq$beta.post.bias),
                    c(ad$beta.fb.lr)/5,
                    c(ad$beta.conf.lr)/5,
                    c(ad$beta.post.bias)/5)
colnames(posneg) <- c('phq_fb', 'phq_conf', 'phq_add',
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
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
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
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
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
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        # axis.title.x = element_blank()  ,
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none',
        # legend.direction = "horizontal"
) 

# 
# posneg = data.frame(c(phq$beta.fb.lr),
#                     c(phq$beta.conf.lr),
#                     c(gad$beta.fb.lr),
#                     c(gad$beta.conf.lr),
#                     c(spin$beta.fb.lr),
#                     c(spin$beta.conf.lr),
#                     c(ad$beta.fb.lr)*.25,
#                     c(ad$beta.conf.lr)*.25)
# colnames(posneg) <- c('phq_fb', 'phq_conf', 'gad_fb', 'gad_conf',
#                       'spin_fb', 'spin_conf', 'ad_fb', 'ad_conf')
# 
# posneg2 <- posneg %>% pivot_longer(
#   cols = c(1:8),
#   names_sep = '_',
#   values_to = 'beta',
#   names_to = c('questname', 'fbconf')) %>% 
#   mutate_at(c('fbconf', 'questname'), as.factor) %>%
#   mutate(fbconf = recode_factor(fbconf, "fb" = 'Feedback', 
#                                 "conf" = 'Confidence') ) %>%
#   mutate(expNum = recode_factor(questname, "phq" = 'Exp 1', 
#                                 "gad" = 'Exp 1', "spin" = 'Exp 1', 'ad' = 'Exp 2'))%>%
#   mutate(questname = recode_factor(questname, "phq" = 'PHQ', 
#                                    "gad" = 'GAD', "spin" = 'SPIN', "ad" = 'AD')) %>%
#   mutate(questname = factor(questname, levels = c('PHQ', 'GAD', 'SPIN', 'AD')))
# 
# 
# ggplot(filter(posneg2, fbconf=='Feedback'), aes(x = questname, y = beta, fill = expNum))+
#   geom_hline(yintercept = 0, color = 'grey50', size=1.5) +
#   geom_violinhalf(alpha=.7)+
#   theme_pubclean() +
#   ylab("Beta (feedback)") +
#   xlab("Mental health score") +
#   labs(fill = 'Dataset') + 
#   # labs(tag = "E") +
#   theme(text = element_text(size=font_size), 
#         plot.caption = element_text(size = font_size),
#         plot.tag = element_text(size = tag_size, face = "bold"),
#         axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
#         axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
#         # legend.position = c(.35,.92),
#         # legend.direction = "horizontal"
#   ) +
#   coord_cartesian(ylim=c(-.15,0.15))
# # + 
# #   guides(fill = guide_legend(nrow = 1))
# 
# ggplot(filter(posneg2, fbconf=='Confidence'), aes(x = questname, y = beta, fill = expNum))+
#   geom_hline(yintercept = 0, color = 'grey50', size=1.5) +
#   geom_violinhalf(alpha=.7)+
#   theme_pubclean() +
#   ylab("Beta (confidence)") +
#   xlab("Mental health score") +
#   labs(fill = 'Dataset') + 
#   theme(text = element_text(size=font_size), 
#         plot.caption = element_text(size = font_size),
#         plot.tag = element_text(size = tag_size, face = "bold"),
#         # legend.position = "none",
#         # panel.border = element_rect(color = "black", fill = NA, size = .5),
#         axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
#         axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
#   coord_cartesian(ylim=c(-.05,0.05)) #+ guides(fill="none")

# 
# #####
# ## CIT model fits
# posneg = data.frame(c(cit$beta.fb.lr*.25),
#                     c(cit$beta.conf.lr),
#                     c(sw$beta.fb.lr*.25),
#                     c(sw$beta.conf.lr))
# colnames(posneg) <- c('cit_fb', 'cit_conf', 'sw_fb', 'sw_conf')
# 
# posneg2 <- posneg %>% pivot_longer(
#   cols = c(1:4),
#   names_sep = '_',
#   values_to = 'beta',
#   names_to = c('questname', 'fbconf')) %>% 
#   mutate_at(c('fbconf', 'questname'), as.factor) %>%
#   mutate(fbconf = recode_factor(fbconf, "fb" = 'Feedback', 
#                                 "conf" = 'Confidence') ) %>%
#   mutate(questname = recode_factor(questname, "cit" = 'CIT', 
#                                    "sw" = 'SW')) %>%
#   mutate(questname = factor(questname, levels = c('CIT', 'SW')))
# 
# ggplot(filter(posneg2), aes(x = questname, y = beta, fill = fbconf))+
#   geom_hline(yintercept = 0, color = 'grey50', size=1.5) +
#   geom_violinhalf(alpha=.7)+
#   theme_pubclean() +
#   ylab("Beta") +
#   xlab("Mental health score") +
#   labs(fill = '') + 
#   theme(text = element_text(size=font_size), 
#         plot.caption = element_text(size = font_size),
#         plot.tag = element_text(size = tag_size, face = "bold"),
#         # legend.position = "none",
#         # panel.border = element_rect(color = "black", fill = NA, size = .5),
#         axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
#         axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) 
# 
# 
#########
### plot model fits to empirical data

spe <- mbsDataExp1[,,1]$spe
nsubjrun <- dim(spe)
speFit <- mbsDataExp1[,,1]$fitZ[,,1]$model.4[,,1]$spe.est

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

setwd('~/OneDrive - University College London/Projects/Experiments/metaBiasShift/data/')
mbsRec1 <- readMat('paramRecovery_fbconf.mat') # load questionnaire data
mbsRec2 <- readMat('paramRecovery_postbias.mat') # load questionnaire data

nsim <- length(c(mbsRec1$fit.lr.conf))
sim.df <- data.frame(c(mbsRec1$fit.lr.fb, mbsRec1$fit.lr.conf), 
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


##########
##### Plot model simulations
#####
setwd("/Users/skatyal/OneDrive - University College London/Projects/Experiments/metaBiasShift/data/simulation/")

### model 1 - feedback distortion (for Figure 3C)
### model 2 - confidence distortion (for Figure 3D)
### model 3 - no distortion (for Figure 3B)
### model 4 - additive distortion (for Figure 3E)

### model 5 - additive distortion with empirical values (for Supplementary Figure 6B)

model2Plot <- 5

mod <- switch(model2Plot, 3,4,5,6,9)

simExp = readMat(paste('simExp_m',mod,'.mat', sep='')) # load questionnaire data
spe <- simExp$spe
confZ <- simExp$conf
phq <- simExp$phq
nsubjrun <- dim(confZ)
fbblock <- simExp$fbblock
task <- simExp$task
fb <- simExp$feedback
group <- simExp$group

subj <- array(c(1:nsubjrun[1]), dim = nsubjrun)
runnum <- aperm(array(c(1:nsubjrun[2]), dim = c(nsubjrun[2], nsubjrun[1], nsubjrun[3])), c(2,1,3))
phq <- array(phq, dim = nsubjrun)
spe <- array(spe, dim = nsubjrun)
fbblock <- array(fbblock, dim = nsubjrun)
task <- array(task, dim = nsubjrun)
trialnum <- aperm(array(c(1:nsubjrun[3]), dim = c(nsubjrun[3], nsubjrun[1], nsubjrun[2])), c(2,3,1))
firstFeedback <- array(c(fbblock[,3,1]), dim = nsubjrun)

exp1.sim <- data.frame(c(spe), c(phq), c(runnum), c(subj), c(fbblock), c(task), 
                       c(confZ), c(trialnum), c(fb), c(firstFeedback))
colnames(exp1.sim) <- c('spe', 'phq', 'runnum', 'subj', 'fbblock', 'task', 
                        'confZ', 'trialnum', 'feedback', 'firstFeedback')

exp1.sim <- exp1.sim %>%
  mutate_at(c('subj', 'fbblock', 'task', 'firstFeedback'), as.factor) %>%
  mutate(firstFeedback = recode_factor(firstFeedback, '1' = 'Positive',
                                       '2' = 'Negative'))

exp1.sim.bas <- filter(exp1.sim, runnum %in% c(1,2)) %>%
  mutate(phqTile = ntile(phq, 2)) %>%
  group_by(subj) %>%
  summarise(spe.bas = mean(spe),
            phqTile = median(phqTile),
            phq = mean(phq)) %>%
  mutate_at(c('phqTile'), as.factor) %>%
  mutate(phqTile = recode_factor(phqTile, '1' = 'Low', '2' = 'High', '3' = 'High'))

exp1.sim.all <- exp1.sim  %>% 
  left_join(exp1.sim.bas) %>%
  # mutate_at('feedback', as.numeric) %>%
  group_by(subj, runnum, fbblock, task, phqTile, firstFeedback) %>%
  summarise(spe = mean(spe),
            spe_b = mean(spe) - mean(spe.bas),
            phq = mean(phq),
            confZ = mean(confZ),
            fbdif = sum(feedback)
  )  %>%
  mutate_at(c('runnum'), as.factor) %>% 
  mutate(isfb = ifelse(fbblock==0, 0, 1)) %>%
  mutate_at('isfb', as.factor)

exp1.sim.fb <- exp1.sim  %>% filter(runnum %in% c(3,5))%>% 
  left_join(exp1.sim.bas) %>%
  group_by(subj, runnum, fbblock, task, phqTile) %>%
  summarise(spe = mean(spe),
            spe_b = mean(spe) - mean(spe.bas)
  )  %>%
  mutate_at(c('runnum'), as.numeric) %>%
   mutate_at(c('fbblock'), as.factor) 
 
exp1.sim.fb$fbblock = recode_factor(exp1.sim.fb$fbblock, "1" = "Positive", "2" = "Negative")

# baseline subtracted by feedback
ggplot(exp1.sim.fb , 
       aes(x = fbblock, y = spe_b, color = phqTile)) +
  scale_color_manual(breaks = c('Low', 'High'), 
                     values = c('tan3', 'royalblue3')) +
  stat_summary(fun.data = mean_se, size = .6, 
               position=position_dodge(width = .1), alpha = .8) +
  stat_summary(fun = mean, linewidth = 1.2, geom = 'line', aes(group=phqTile),
               position=position_dodge(width = .1), alpha = .8) +
  theme_pubclean() +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none',
        legend.title = element_text(size = 26),
        legend.text = element_text(size=22)
  ) +
  coord_cartesian(ylim=c(-.2,.15)) +
  ylab("SPE-b") +
  xlab("Feedback") + labs(color = 'AD score')


exp1.sim.conf <- exp1.sim %>%
  mutate(phqTile = ntile(phq, 2)) %>%
  group_by(subj, runnum) %>%
  summarise(spe = mean(spe),
            confZ = mean(confZ)) %>% 
  mutate(confTile = ntile(confZ, 6)) %>%
  left_join(exp1.sim.bas)

##############################
########## Figure 3B-E

ggplot(exp1.sim.conf %>% filter(runnum %in% c(1,2,4,6))
         , aes(x = runnum, y = spe, color = phqTile)) +
  stat_summary(fun.data = mean_se, size = .5,
               position=position_dodge(width = .1), alpha=.8) +
  scale_color_manual(breaks = c('Low', 'High'), 
                     values = c('tan3', 'royalblue3')) +
  stat_summary(fun = mean, geom = 'line', size=1.2,
               position=position_dodge(width = .1), alpha=1) +
  theme_pubclean() +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        # legend.position = 'none'
        ) +
  scale_y_continuous(breaks = seq(.55,.7,.05))+
  scale_x_continuous(breaks = c(1,2,4,6))+
  coord_cartesian(ylim=c(.53,.73)) +
  ylab("SPE") +
  xlab("Block #") +
  labs(colour = 'AD score', linetype = 'Conf level') 

