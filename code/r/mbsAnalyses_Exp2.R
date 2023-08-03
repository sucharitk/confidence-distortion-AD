
####################################################################
################# Analysis and Figure generation for Exp 2 #########
################# Katyal, Huys, Dolan, Fleming #####################
####################################################################
### How underconfidence is maintained in anxiety and depression ####
####################################################################
####  This file reproduces all analyses and figures for Exp 2 ######
####################################################################

library(R.matlab)
library(lme4) #Linear mixed effect model
library(emmeans) # Least squares means
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm) # for geom_quasirandom()
library(ggpp)
library(sjPlot)
library(see)

## mediation analysis on mixed models does not work with lmerTest package
## so use only one or the other 

doMediation = F

if (doMediation){
  library(mediation) # mediation analysis
} else {
  library(lmerTest) # anova and summary with p-values
  
}


emm_options(lmer.df = "satterthwaite")
emm_options(lmerTest.limit = 35000)
font_size <- 22
tag_size <- 24
pval_size <- 6
alpha <- .2
nLags = 0
subtractBaselineConfSpe = T
taskNames = c('Perception', 'Memory')
transfertaskNames = c('Same', 'Opposite')
posnegNames = c('Positive', 'Negative')
posnegColours = c("green4", "firebrick2")
percmemColours = c("cyan3", "tomato2")

if (Sys.info()["sysname"]=='Darwin'){
  setwd("~/OneDrive - University College London/Projects/Experiments/metaBiasShift/")
} else{
  setwd("C:/Users/skatyal/OneDrive - University College London/Projects/Experiments/metaBiasShift/")
}
## read data file
exp2.mbs = read.csv(paste('data/exp2/mbsExp2.csv', sep=''),
               header = TRUE, sep = ',')
## recode columns as factors
exp2.mbs$subj <- as.factor(exp2.mbs$subj)
exp2.mbs$task <- as.factor(exp2.mbs$task)
exp2.mbs$awarepos[exp2.mbs$awarepos==0] <- NA
exp2.mbs$awarepos[is.nan(exp2.mbs$awarepos)] <- NA
exp2.mbs$awareneg[exp2.mbs$awareneg==0] <- NA
exp2.mbs$awareneg[is.nan(exp2.mbs$awareneg)] <- NA
exp2.mbs$awarepos <- as.factor(exp2.mbs$awarepos)
exp2.mbs$awareneg <- as.factor(exp2.mbs$awareneg)
exp2.mbs$affectpos[is.nan(exp2.mbs$affectpos)] <- NA
exp2.mbs$affectneg[is.nan(exp2.mbs$affectneg)] <- NA
exp2.mbs$affectpos[exp2.mbs$affectpos!=1] <- 0
exp2.mbs$affectneg[exp2.mbs$affectneg!=2] <- 0
exp2.mbs$affectpos <- as.factor(exp2.mbs$affectpos)
exp2.mbs$affectneg <- as.factor(exp2.mbs$affectneg)
exp2.mbs$gender[is.nan(exp2.mbs$gender)] <- NA
exp2.mbs$gender <- as.factor(exp2.mbs$gender)

exp2.mbs <- exp2.mbs %>%
  mutate_at(vars(contains("fb")), as.factor) %>%
  rename(accu = corr) %>%
  rename(stair = incdec)

exp2.mbs <- exp2.mbs %>%
  mutate(fbblock = recode_factor(fbblock, "1" = posnegNames[1], 
                                 "2" = posnegNames[2], "0" = "None"),
         task = recode_factor(task, "0" = taskNames[1], "1" = taskNames[2]),
         awarepos = recode_factor(awarepos, '1' = 'Yes', '2'='No'),
         awareneg = recode_factor(awareneg, '1' = 'Yes', '2'='No'),
         affectpos = recode_factor(affectpos, '1' = 'Yes', '0'='No'),
         affectneg = recode_factor(affectneg, '2' = 'Yes', '0'='No')
  )

## code highly deviant RTs as NA
max_RT_deviation = 3
rt1_max = c(median(exp2.mbs$rt1[exp2.mbs$task==taskNames[1]], na.rm=T)+ 
              max_RT_deviation*IQR(exp2.mbs$rt1[exp2.mbs$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(exp2.mbs$rt1[exp2.mbs$task==taskNames[2]], na.rm=T)+ 
              max_RT_deviation*IQR(exp2.mbs$rt1[exp2.mbs$task==taskNames[2]], na.rm=T ))
exp2.mbs <- filter(exp2.mbs, trialnum > nLags) %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[1] & .>rt1_max[1],NA)) %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[2] & .>rt1_max[2],NA))

## scale RTs 
exp2.mbs$scale_rt <- exp2.mbs$rt1
exp2.mbs$scale_rt[exp2.mbs$task==taskNames[1]] <- scale(exp2.mbs$rt1[exp2.mbs$task==taskNames[1]])
exp2.mbs$scale_rt[exp2.mbs$task==taskNames[2]] <- scale(exp2.mbs$rt1[exp2.mbs$task==taskNames[2]])

rt2_max = c(median(exp2.mbs$rt2[exp2.mbs$task==taskNames[1]], na.rm=T)+ 
              max_RT_deviation*IQR(exp2.mbs$rt2[exp2.mbs$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(exp2.mbs$rt2[exp2.mbs$task==taskNames[2]], na.rm=T)+ 
              max_RT_deviation*IQR(exp2.mbs$rt2[exp2.mbs$task==taskNames[2]], na.rm=T ))
exp2.mbs <- filter(exp2.mbs, trialnum > nLags) %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[1] & .>rt2_max[1],NA)) %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[2] & .>rt2_max[2],NA))

exp2.mbs$scale_rt2 <- exp2.mbs$rt2
exp2.mbs$scale_rt2[exp2.mbs$task==taskNames[1]] <- scale(exp2.mbs$rt2[exp2.mbs$task==taskNames[1]])
exp2.mbs$scale_rt2[exp2.mbs$task==taskNames[2]] <- scale(exp2.mbs$rt2[exp2.mbs$task==taskNames[2]])

exp2.mbs1 <- exp2.mbs

## if to subtract the baseline blocks confidence and spe from the rest of the blocks
  # get the baseline values for the two tasks from runs 1 and 2
  exp2.mbs2 <- filter(exp2.mbs, runnum<3) %>%
    group_by(subj, task, group) %>%
    dplyr::summarise(
      conf.bas = mean(conf),
      spe.bas = mean(spe),
      accu.bas = mean(accu)
    )
  
  # save the baseline data separately
  exp2.mbs.bas <- filter(exp2.mbs, runnum<3) %>%
    group_by(subj, task, runnum, group, trialnum, gender) %>%
    dplyr::summarise(
      conf = mean(conf),
      accu = mean(accu),
      stair = mean(stair),
      scale_rt = mean(scale_rt),
      rt1 = mean(rt1),
      spe = mean(spe),
      spe0_perc = mean(spe0_perc),
      spe0_mem = mean(spe0_mem),
      age = mean(age),
      SDS = mean(SDS),
      STAI = mean(STAI),
      LSAS = mean(LSAS)
    )
  exp2.mbs <- exp2.mbs %>% left_join(exp2.mbs2)
  
  if (subtractBaselineConfSpe){
    # left join the baseline values an subtract them from confidence and spe
  exp2.mbs <- exp2.mbs %>%
    mutate(accu_b = accu - accu.bas) %>%
    mutate(conf = conf-conf.bas)%>%
    mutate(spe = spe-spe.bas) 
}

## take away the baseline (first two) blocks from the dataset
exp2.mbs <- exp2.mbs %>%
  filter(runnum>2) #%>%
  # mutate(runnum = runnum-2)
exp2.mbs <- exp2.mbs %>%
  mutate(blockType = factor(runnum%%2)) %>%
  mutate(blockType = recode_factor(blockType, "0" = "transf", "1" = "interv"))
exp2.mbs <- exp2.mbs %>%
  mutate_at(vars(contains("runnum")), as.factor)
exp2.mbs$fbblock <- droplevels(exp2.mbs$fbblock)


## create new columns for trial-by-trial positive and negative feedback
exp2.mbs <- exp2.mbs %>%
  mutate(fbpos=feedback, fbneg=feedback) %>%
  mutate_at("fbpos", ~replace(., feedback==1 & accu==1, 1)) %>%
  mutate_at("fbpos", ~replace(., feedback==1 & accu==0, 0)) %>%
  mutate_at("fbneg", ~replace(., feedback==1 & accu==1, 0)) %>%
  mutate_at("fbneg", ~replace(., feedback==1 & accu==0, 1)) 

# give both interv and transf run the same value of fbblock (feedback block type)
exp2.mbs$fbblock[exp2.mbs$runnum==2] = exp2.mbs$fbblock[exp2.mbs$runnum==1 & 
                                                          exp2.mbs$trialnum<=20]
exp2.mbs$fbblock[exp2.mbs$runnum==4] = exp2.mbs$fbblock[exp2.mbs$runnum==3
                                                        & exp2.mbs$trialnum<=20]

## code intervention task by group
exp2.mbs$intervTask = exp2.mbs$group
exp2.mbs$intervTask[exp2.mbs$group %in% c(1,2,3,4)] = taskNames[1]
exp2.mbs$intervTask[exp2.mbs$group %in% c(5,6,7,8)] = taskNames[2]

# transfer task by group (same as task column)
exp2.mbs$testTask = exp2.mbs$group
exp2.mbs$testTask[exp2.mbs$group %in% c(1,2,3,4)] = taskNames[2]
exp2.mbs$testTask[exp2.mbs$group %in% c(5,6,7,8)] = taskNames[1]

## code whether pos or neg is the first block of feedback
exp2.mbs$firstFeedback = exp2.mbs$group
exp2.mbs$firstFeedback[exp2.mbs$group %in% c(1,3,5,7)] = posnegNames[1]
exp2.mbs$firstFeedback[exp2.mbs$group %in% c(2,4,6,8)] = posnegNames[2]


# exp2.mbs$group <- as.factor(exp2.mbs$group)
exp2.mbs <- exp2.mbs %>%
  mutate_at(c("group", "intervTask", "firstFeedback", 'testTask'), as.factor) %>%
  mutate(group = recode_factor(group, "1" = 'Group 1', 
                               "2" = 'Group 2', "3" = 'Group 3', "4" = 'Group 4', 
                               "5" = 'Group 5', "6" = 'Group 6', "7" = 'Group 7', 
                               "8" = 'Group 8')) 

# #####################
# #####################
# ## read in the psychiatric scores
# ## read data file
psych = read.csv(paste('data/exp2/factor_scores.csv', sep=''),
                 header = TRUE, sep = ',')
psych <- psych %>% rename(subj = subjIDs) %>% dplyr::select(-X) %>%
  mutate_at("subj", as.factor)

exp2.mbs <- exp2.mbs %>% left_join(psych)
exp2.mbs.bas <- exp2.mbs.bas %>% left_join(psych)

exp2.mbs3 <- exp2.mbs # save the non NA removed data for sret analysis which does not depend on trials

exp2.mbs <- exp2.mbs %>% drop_na(scale_rt, scale_rt2)

exp2.mbs$scale_rt[as.logical(exp2.mbs$invalid_rt1)] = NA # for trials where participants changed their response, mark rt1 as NA


#####################
#####################


exp2.df <- exp2.mbs %>%
  group_by(subj, task, runnum, fbblock, group, 
           blockType, intervTask, testTask,sretset,
           firstFeedback, awareneg, awarepos,
           affectpos, affectneg, gender) %>%
  dplyr::summarise(
    conf = mean(conf),
    conf.bas = mean(conf),
    accu_b = mean(accu)-mean(accu.bas),
    accu = mean(accu),
    scale_rt = mean(scale_rt, na.rm=T),
    rt1 = mean(rt1, na.rm=T),
    fbneg = mean(fbneg),
    fbpos = mean(fbpos),
    spe = mean(spe),
    spe0_perc = mean(spe0_perc),
    spe0_mem = mean(spe0_mem),
    spe.bas = mean(spe.bas),
    age = mean(age),
    stair = mean(stair),
    AD = mean(AD),
    Compul = mean(Compul),
    SW = mean(SW))

#####################
exp2.mbs1 <- exp2.mbs1 %>% left_join(psych)


######################################
##### Supp Figure 4. SPE of intervention block

exp2.df.2 <- filter(exp2.df, blockType=="interv")

sf4.left <- ggplot(exp2.df.2, aes(x = fbblock, y = spe, color = fbblock)) +
  geom_hline(yintercept = 0, color = 'grey30') +
  scale_color_manual(breaks = posnegNames, values=posnegColours) +
  theme_pubclean()  +
  # scale_x_discrete(limits = taskNames) +
  geom_quasirandom(dodge.width=.7, alpha = .7) +
  geom_violin(alpha=.6, position = position_dodge(.7)) +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .7),
               geom = 'errorbar',
               size = .9, aes(width = .2)) +
  ylab("Global confidence:\nSelf-estimated performance") +
  xlab("Feedback") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) + 
  coord_cartesian(ylim = c(-.65,.65)) +
  labs(color = 'Feedback') + 
  # labs(tag = "A") +
  theme(text = element_text(size=font_size), 
        legend.position = 'none', 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  geom_signif(y_position = c(.56), vjust = -.3, hjust=.4,
              xmin = c(1.1), color = 'black',
              xmax = c(1.9),
              annotation = c("p < .0001"), tip_length = .04,
              textsize = pval_size, size = .8)
sf4.left

sf6.right <- ggplot(exp2.df.2, aes(x = fbblock, y = accu_b, color = fbblock)) +
  geom_hline(yintercept = 0, color = 'grey30') +
  scale_color_manual(breaks = posnegNames, values=posnegColours) +
  theme_pubclean()  +
  # scale_x_discrete(limits = taskNames) +
  geom_quasirandom(dodge.width=.7, alpha = .7) +
  geom_violin(alpha=.6, position = position_dodge(.7)) +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .7),
               geom = 'errorbar',
               size = .9, aes(width = .2)) +
  ylab("Actual performance") +
  xlab("Feedback") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.position = 'none') + 
  coord_cartesian(ylim = c(-.2,.2)) +
  labs(color = 'Feedback') + 
  # labs(tag = "A") +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  geom_signif(y_position = c(.19), vjust = -.3,
              xmin = c(1.1), color = 'black',
              xmax = c(1.9),
              annotation = c("n.s."), tip_length = .04,
              textsize = pval_size, size = .8)
sf6.right


m.reg <- lm(accu ~ fbblock*task,
            exp2.df.2)
summary(m.reg)
m.reg <- lm(accu ~ fbblock ,
            exp2.df.2)
summary(m.reg)
m.reg <- lmer(stair ~ fbblock*task + (1|subj),
              exp2.df.2)
summary(m.reg)
m.reg <- lmer(stair ~ fbblock + (1|subj),
              exp2.df.2)
summary(m.reg)

######################################
##### Supp Figure 5 (lower). Task difficulty during intervention blocks

sf5.bottom <- ggplot(exp2.df.2, aes(x = task, y = stair, color = fbblock)) +
  # geom_hline(yintercept = .74, color = 'grey30') +
  scale_color_manual(breaks = posnegNames, values=posnegColours) +
  theme_pubclean()  +
  scale_x_discrete(limits = taskNames) +
  geom_quasirandom(dodge.width=.7, alpha = .7) +
  geom_violin(alpha=.6, position = position_dodge(.7)) +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .7),
               geom = 'errorbar',
               size = .9, aes(width = .2)) +
  ylab("Average difficulty achieved") +
  xlab("Intervention task") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ggtitle('Exp 2') +
  labs(color = 'Feedback') + 
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  coord_cartesian(ylim = c(0,11))  +
  geom_signif(y_position = c(9, 9), vjust = -.3,
              xmin = c(0.9, 1.9), color = 'black',
              xmax = c(1.1, 2.1),
              annotation = c("n.s.", "n.s."), tip_length = .04,
              textsize = pval_size, size = .8)
sf5.bottom

######################################
##### Supp Figure 9. Transfer effect of feedback upon test block confidence

exp2.df.2 <- filter(exp2.mbs, blockType=="transf") %>% 
  mutate(intervTask = relevel(intervTask, "Perception")) %>%
  mutate(testTask = relevel(testTask, "Perception"))

sf9a <- ggplot(exp2.df.2, aes(x = trialnum, y = conf, colour = fbblock)) +
  geom_hline(yintercept = 0, color = 'gray') +
  scale_color_manual(breaks = posnegNames, values=posnegColours) +
  # scale_linetype_manual(values = c('twodash')) +
  facet_grid(testTask~intervTask,
             labeller = labeller(intervTask = interv.labs, testTask = test.labs)) + 
  theme_pubclean() +
  geom_smooth(method='loess', span=.6,se=F, size=1.4,
              linetype = 'twodash')+
  ylab("Confidence (test)") +
  xlab("Trial number") +
  labs(colour = 'Interv. block feedback') +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = .5),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15)) +
  guides(color = guide_legend(nrow = 1))

sf9a

exp2.df.2 <- exp2.mbs %>% 
  filter(blockType=="transf") %>% 
  group_by(intervTask, testTask, fbblock, subj) %>%
  summarise(conf = mean(conf)) %>%
  pivot_wider(
    names_from = fbblock,
    values_from = conf
  ) %>% group_by(subj, intervTask, testTask) %>%
  summarise(
    confPosNeg = mean(Positive) - mean(Negative)
  )%>% mutate(intervTask = relevel(intervTask, "Perception")) %>%
  mutate(testTask = relevel(testTask, "Perception"))

sf9b <- ggplot(exp2.df.2, aes(x = intervTask, y = confPosNeg, 
                              color = forcats::fct_infreq(testTask))) +
  scale_color_manual(values=percmemColours, breaks = taskNames) +
  geom_hline(yintercept = 0, color = 'grey30') +
  geom_quasirandom(dodge.width=.8, alpha = .7) +
  geom_violin(alpha=.6, position = position_dodge(.8)) +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .8),
               geom = 'errorbar',
               size = .9, aes(width = .2)) +
  scale_x_discrete(limits = taskNames) +
  ylab("Confidence (test)\nPositive - Negative") +
  xlab("Intervention task") +
  labs(color = 'Test task') +
  coord_cartesian( ylim = (c(-.3,.45))) +
  theme_pubclean() +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  geom_signif(map_signif_level = TRUE,
              annotations = c("p = .004"),vjust = -.3,
              y_position = c(.42),
              xmin = 1.03, xmax = 1.97,color= 'black',
              textsize = pval_size, size = .8,
              tip_length = .03) +
  geom_signif(y_position = c(.4, .4),
              xmin = c(.85, 1.85),
              xmax = c(1.15, 2.15),
              annotation = c("", ""),  color= 'black',
              textsize = pval_size, size = .8,
              tip_length = .0)
sf9b

######################################
##### Supp Figure 14b. transfer of feedback to test block SPE

exp2.df.2 <- exp2.mbs %>% 
  filter(blockType=="transf") %>% 
  group_by(fbblock, subj, intervTask, testTask) %>%
  summarise(spe = mean(spe)) %>%
  pivot_wider(
    names_from = fbblock,
    values_from = spe
  ) %>% group_by(subj, intervTask, testTask) %>%
  summarise(
    spe = mean(Positive) - mean(Negative)
  )%>% mutate(intervTask = relevel(intervTask, "Perception")) %>%
  mutate(testTask = relevel(testTask, "Perception"))

sf14b <- ggplot(exp2.df.2, aes(x = intervTask, y = spe, color = testTask)) +
  geom_hline(yintercept = 0, color = 'grey30') +
  scale_color_manual(values=percmemColours, breaks = taskNames) +
  geom_violin(alpha=.8, position = position_dodge(.8)) +
  stat_summary(fun.data = mean_se,
               position=position_dodge(width = .8),
               geom = 'errorbar',
               size = 1, aes(width = .1))  +
  geom_quasirandom(dodge.width=.8,
                   alpha = .15) +
  scale_x_discrete(limits = taskNames) +
  ylab("SPE-b (test)\nPositive - Negative") +
  xlab("Intervention task") +
  labs(color = 'Test task') +
  coord_cartesian(ylim = c(-.5,.75)) +
  theme_pubclean() +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  geom_signif(y_position = c(.73),
              xmin = c(.87), vjust = -.3,
              xmax = c(1.17), color = 'black',
              annotation = c("p = .048"), 
              tip_length = .0, textsize = pval_size, size = .8)
sf14b



############
############
###############################################
### SPE (intervention) analyses
################################################

exp2.df.2 <- filter(exp2.df, blockType=="interv")


m.reg <- lmer(spe ~ fbblock*task*firstFeedback + 
                accu + conf + stair + scale_rt +
                (1|subj) + (1|group) + (1|runnum), REML = F, exp2.df.2)
summary(m.reg)
m.reg <- lmer(spe ~ fbblock*task*firstFeedback + 
                accu + conf + stair + scale_rt +
                (1|subj), REML = F, exp2.df.2)
summary(m.reg)
plot(m.reg)

m2 <- update(m.reg, ~.-fbblock:task:firstFeedback)
anova(m.reg,m2)
m.reg <- m2 # accept the reduced model
summary(m.reg)

m2 <- update(m.reg, ~.-firstFeedback:task)
anova(m.reg,m2)
m.reg <- m2 # accept the reduced model
summary(m.reg)

m2 <- update(m.reg, ~.-fbblock:firstFeedback)
anova(m.reg,m2)
m.reg <- m2 # accept the reduced model
summary(m.reg)

m.reg <- update(m.reg, REML=T)
summary(m.reg)
plot(m.reg)

emmeans(m.reg, pairwise ~ fbblock|task)$contrasts
emmeans(m.reg, pairwise ~ task|fbblock)$contrasts
emm <- emmeans(m.reg, specs = c("fbblock"))
test(emm)
emmeans(m.reg, pairwise ~ fbblock)$contrasts


plot_model(m.reg, type = c("pred"),
           terms = c("fbblock", 'task'), 
           title = "") + 
  theme_pubclean() +
  xlab('Feedback type') + 
  ylab('SPE (intervention)') +
  theme_pubclean() + labs(tag = 'B') +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"))
# highly sig eff of fb block



############
############
############
### Effect of feedback on local conf (test block)

exp2.df3 <- filter(exp2.mbs, blockType=="transf")
exp2.mbs.reg = lmer(conf ~ fbblock*task + accu + stair + 
                 (1|subj) + (1|trialnum)  , exp2.df3)
summary(exp2.mbs.reg)
emm <- emmeans(exp2.mbs.reg, ~fbblock|task)
contrast(emm)
m2 <- update(exp2.mbs.reg, ~.-task:fbblock)
anova(exp2.mbs.reg,m2)
exp2.mbs.reg <- m2 # accept the reduced model
summary(exp2.mbs.reg)


exp2.mbs.reg <- update(exp2.mbs.reg, REML=T)
plot(exp2.mbs.reg)
summary(exp2.mbs.reg)
summary(exp2.mbs.reg)$coefficients

emm <- emmeans(exp2.mbs.reg, ~fbblock|task)
contrast(emm)

plot_model(exp2.mbs.reg, type=c("pred"),
           terms = c( "task", "fbblock"))


############
## Mediation analysis to test if spe(intervention) mediates the effect
## of pos and neg fb on conf(test)

first_n_of_trans <- c(1:20)

exp2.df7 <- filter(exp2.mbs, blockType=="transf",
              trialnum %in% first_n_of_trans) %>%
  group_by(subj, runnum, fbblock, task, group) %>%
  dplyr::summarise(
    spe = mean(spe),
    accu = mean(accu),
    conf = mean(conf),
    stair = mean(stair),
    scale_rt = mean(scale_rt, na.rm = T)
  )

# num of divisions of intervention block
ndivs_per_run <- 1

exp2.df.temp <- filter(exp2.mbs, blockType=="interv")

# create a new factor by combining run and trialnumber and then divide it into equal divisions
exp2.df.temp <- exp2.df.temp %>%
  mutate(rundiv = factor(ceiling(trialnum/(40/ndivs_per_run))))

exp2.df.temp <- exp2.df.temp %>%
  group_by(rundiv, subj, task, runnum, fbblock) %>%
  dplyr::summarise(
    conf = mean(conf),
    accu = mean(accu),
    scale_rt = mean(scale_rt, na.rm=T),
    fbneg = mean(fbneg),
    fbpos = mean(fbpos),
    spe = mean(spe),
    stair = mean(stair)
  ) %>%
  ungroup() %>%
  mutate_at("rundiv", as.factor)

exp2.df.temp <- exp2.df.temp %>% dplyr::select(
  subj, rundiv, conf, spe, accu, stair, scale_rt,
  fbneg, fbpos,fbblock
) %>%
  pivot_wider(
    names_from = c( rundiv),
    values_from = c(conf, fbneg, fbpos, spe, accu, stair, scale_rt)
  )

exp2.df8 <- exp2.df7 %>% left_join(exp2.df.temp)


## mediation combined
med.m1 <- lmer(conf ~ fbpos_1 + fbneg_1+  (1|subj), exp2.df8)
summary(med.m1)

med.m2 <- lmer(spe_1 ~ fbpos_1 + fbneg_1+  (1|subj), exp2.df8)
summary(med.m2)

med.m3 <- lmer(conf ~ fbpos_1 + fbneg_1+spe_1 + (1|subj), exp2.df8)
summary(med.m3)

med.pos = mediate(med.m2, med.m3, treat='fbpos_1', mediator='spe_1',
                  boot=F, sims = 1000)
summary(med.pos)

med.pos = mediate(med.m2, med.m3, treat='fbneg_1', mediator='spe_1',
                  boot=F, sims = 1000)
summary(med.pos)


############
############
### Effect of feedback on SPE (test block)

exp2.df9 <- filter(exp2.mbs, blockType=="transf") %>%
  group_by(subj, runnum, fbblock, task, group) %>%
  dplyr::summarise(
    spe = mean(spe),
    accu = mean(accu),
    conf = mean(conf),
    scale_rt = mean(scale_rt, na.rm = T),
    stair = mean(stair)
  )

exp2.mbs.reg = lmer(spe ~ fbblock*task +
                 accu*task + scale_rt*task + stair*task + 
                 (1|subj), REML = F,
               exp2.df9)
summary(exp2.mbs.reg)

emm <- emmeans(exp2.mbs.reg, pairwise~fbblock|task)
test(emm)

m2 <- update(exp2.mbs.reg, ~.-fbblock:task)
anova(exp2.mbs.reg,m2)
exp2.mbs.reg <- m2 # accept the reduced model
summary(exp2.mbs.reg)

m2 <- update(exp2.mbs.reg, ~.-task:accu)
anova(exp2.mbs.reg,m2)
exp2.mbs.reg <- m2 # accept the reduced model
summary(exp2.mbs.reg)
m2 <- update(exp2.mbs.reg, ~.-scale_rt:task)
anova(exp2.mbs.reg,m2)
exp2.mbs.reg <- m2 # accept the reduced model
summary(exp2.mbs.reg)

exp2.mbs.reg <- update(exp2.mbs.reg, REML=T)
summary(exp2.mbs.reg)

plot_model(exp2.mbs.reg,
           title="Perception task - interaction",
           show.p=TRUE, show.values = T,
           value.offset = .4, value.size = 3.8) +
  coord_cartesian( ylim = (c(-.65,.65)))

emm <- emmeans(exp2.mbs.reg, pairwise~fbblock)
test(emm)


############
## mediation

# num of divisions of intervention block
ndivs_per_run <- 1

exp2.df.temp <- filter(exp2.mbs, blockType=="interv")

# create a new factor by combining run and trialnumber and then divide it into equal divisions
exp2.df.temp <- exp2.df.temp %>%
  mutate(rundiv = factor(ceiling(trialnum/(40/ndivs_per_run))))

exp2.df.temp <- exp2.df.temp %>%
  group_by(rundiv, subj, task, runnum, fbblock) %>%
  dplyr::summarise(
    conf = mean(conf),
    accu = mean(accu),
    scale_rt = mean(scale_rt, na.rm=T),
    fbneg = mean(fbneg),
    fbpos = mean(fbpos),
    spe = mean(spe),
    stair = mean(stair)
  ) %>% 
  ungroup() %>%
  mutate_at("rundiv", as.factor)


exp2.df.temp <- exp2.df.temp %>% dplyr::select(
  subj, rundiv, conf, accu, scale_rt, spe, 
  fbneg, fbpos,stair, fbblock
) %>%
  pivot_wider(
    names_from = c( rundiv),
    values_from = c(conf, accu, scale_rt, fbneg, fbpos, stair, spe)
  ) 

exp2.df10 <- exp2.df9 %>% left_join(exp2.df.temp)


## mediation analyses

med.m1 <- lmer(spe ~ fbpos_1 + fbneg_1 + (1|subj), exp2.df10)
summary(med.m1)

med.m2 <- lmer(spe_1 ~ fbpos_1 + fbneg_1 + (1|subj), exp2.df10)
summary(med.m2)

med.m3 <- lmer(spe ~ fbpos_1 + fbneg_1 + spe_1 + (1|subj), exp2.df10)
summary(med.m3)

med.pos = mediate(med.m2, med.m3, treat='fbpos_1', mediator='spe_1',
                  boot=F, sims = 1000)
summary(med.pos)

med.pos = mediate(med.m2, med.m3, treat='fbneg_1', mediator='spe_1',
                  boot=F, sims = 1000)
summary(med.pos)



#######################################
#####3##############################
#### Baseline correlation with mhq

exp2.df11 <- exp2.mbs.bas %>% group_by(task, runnum, group, subj, gender, trialnum) %>%
  summarise(conf = mean(conf),
            accu = mean(accu),
            stair = mean(stair),
            scale_rt = mean(scale_rt, na.rm=T),
            # rt1 = mean(rt1, na.rm=T),
            spe = mean(spe),
            AD = mean(AD),
            Compul = mean(Compul),
            SW = mean(SW),
            age = mean(age))

exp2.df4 <- exp2.mbs.bas %>% group_by(task, runnum, group, subj, gender) %>%
  summarise(conf = mean(conf),
            accu = mean(accu),
            stair = mean(stair),
            scale_rt = mean(scale_rt, na.rm=T),
            rt1 = mean(rt1, na.rm=T),
            spe = mean(spe),
            AD = mean(AD),
            Compul = mean(Compul),
            SW = mean(SW),
            sds = mean(SDS),
            stai = mean(STAI),
            lsas = mean(LSAS),
            age = mean(age))

exp2.spe0 <- exp2.mbs.bas %>% group_by(group, subj)%>%
  summarise(
    spe0_perc = mean(spe0_perc),
    spe0_mem = mean(spe0_mem),
  ) %>%
  pivot_longer(
    cols = starts_with("spe0_"),
    names_to = "task",
    names_prefix = "spe0_",
    values_to = "spe0"
  ) %>%  mutate(
    task = recode_factor(task, "perc" = taskNames[1], "mem" = taskNames[2])
  )

exp2.df4 <- exp2.df4 %>% left_join(exp2.spe0)

### 
#
# for correlation of transdiag axes with confidence use individual trial level data
m.reg = lmer(conf ~ AD + Compul + SW +
               age + gender +
               accu + stair + scale_rt +
               (1|task) + (1|runnum) + (1|trialnum) + (1|group),exp2.df11)
summary(m.reg)
m.reg = lmer(conf ~ AD + Compul + SW +
               age + gender +
               accu + stair + scale_rt +
               (1|task) + (1|runnum) + (1|group),exp2.df11)
summary(m.reg)

lconf.beta <- summary(m.reg)$coefficients[2:4,1]
lconf.ci <- confint(m.reg)
lconf.ci <- lconf.ci[6:8,]

# for correlation of transdiag axes with SPEs use at block level
m.reg = lmer(spe ~ AD + Compul + SW + 
               age + gender + 
               accu + stair + scale_rt +
               (1|task),exp2.df4)
summary(m.reg)

spe.beta <- summary(m.reg)$coefficients[2:4,1]
spe.ci <- confint(m.reg)
spe.ci <- spe.ci[4:6,]


m.reg = lmer(spe0 ~ AD + Compul + SW + 
               age + gender + 
               accu + stair + scale_rt +
               (1|task),exp2.df4)
summary(m.reg)

pspe.beta <- summary(m.reg)$coefficients[2:4,1]
pspe.ci <- confint(m.reg)
pspe.ci <- pspe.ci[4:6,]



betas <- c(lconf.beta, spe.beta, pspe.beta)
cilo <- c(lconf.ci[,1], spe.ci[,1], pspe.ci[,1])
cihi <- c(lconf.ci[,2], spe.ci[,2], pspe.ci[,2])
tdims <- c('AD', 'CIT', 'SW')
tdims <- c(tdims, tdims, tdims)
# ctype <- c(array('Local',3), 
#            array('Global (retrosp.)',3), 
#            array('Global (prosp.)',3))
ctype <- c(array('Local confidence',3), 
           array('Global SPE (retrosp.)',3), 
           array('Global SPE (prosp.)',3))

exp2.mhqbeta <- data.frame(betas, tdims, ctype, cilo, cihi)
colnames(exp2.mhqbeta) <- c('betas', 'tdims',
                            'ctype', 'cilo', 'cihi')
exp2.mhqbeta <- exp2.mhqbeta %>%
  mutate(ctype = factor(ctype, levels = c( 'Local confidence',
                        'Global SPE (retrosp.)',
                        'Global SPE (prosp.)'))) %>%
  mutate(tdims = recode_factor(tdims, 'AD' = 'Anxious-\nDepression',
                               'CIT' = 'Compulsivity &\nIntrusive Thought',
                               'SW' = 'Social\nWithdrawal'))

##################################
##### Figure 2B. Transdiagnostic axes and local-global confidence 
f2b <- ggplot(exp2.mhqbeta, aes(x = tdims, y = betas, fill =ctype)) +
  scale_fill_manual(breaks = c( 'Local confidence', 'Global SPE (retrosp.)', 'Global SPE (prosp.)'), 
                     values= c('orchid', 'deepskyblue' , 'khaki3')) +
  geom_bar(position=position_dodge(.9), stat='identity',
           width=.8) +
  geom_errorbar(aes(ymin=cilo, ymax=cihi),
                width=.2, size=1,                   # Width of the error bars
                position=position_dodge(.9)) +
  theme_pubclean() +
  labs(fill = '') +
  coord_cartesian(ylim = c(-.09,.058)) +
  ylab('Regression slopes') +
  xlab('') +
  theme(text = element_text(size=font_size),
        legend.direction = 'horizontal',
        legend.position = c(.42,1.05),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        plot.margin = margin(t = 30, r = 0, b = 0, l = 5, unit = "pt"),
        plot.caption = element_text(size = font_size)) +
  geom_signif(y_position = c(-.087, -.092, -.085, .037),
              xmin = c(.7, 1, 1.3,1.7),
              xmax = c(.7, 1, 1.3, 1.7),
              annotation = c('****', '****', '****', '****'),
              tip_length = 0, color = 'black', textsize = pval_size, size = 0)
f2b

####################################
####################################
###### model predictions with mhq

exp2.df2 <- exp2.mbs1 %>%
  select(c(subj, runnum, AD, Compul, SW, fbblock, task, spe, feedback, accu, group)) %>%
  mutate(firstFeedback = ifelse(group %in% c(1,3,5,7), 'Positive', 'Negative')) %>%
  group_by(subj, runnum, fbblock, task, firstFeedback) %>%
  summarise(spe = mean(spe),
            fbpos = sum(feedback==1 & accu==1),
            fbneg = sum(feedback==1 & accu==0),
            fbboth = sum(feedback!=0),
            AD = mean(AD),
            CIT = mean(Compul),
            SW = mean(SW)) %>%
  ungroup() %>%
  mutate(ADTile = ntile(AD,2)) %>%
  mutate(CTile = ntile(CIT, 2)) %>%
  mutate(STile = ntile(SW, 2)) %>%
  mutate_at(c('ADTile', 'CTile', 'STile'), as.factor) %>%
  mutate(ADTile = recode_factor(ADTile, '1' = 'Low', '2' = 'High')) %>%
  mutate(CTile = recode_factor(CTile, '1' = 'Low', '2' = 'High')) %>%
  mutate(STile = recode_factor(STile, '1' = 'Low', '2' = 'High')) %>%
  mutate_at(c('runnum'), as.numeric) #%>%# left_join(exp2.speFit)

exp2.df2$fbblock[exp2.df2$runnum==4|exp2.df2$runnum==6] <- 'None'

exp2.df2.bas <- filter(exp2.df2, runnum %in% c(1,2)) %>%
  group_by(subj) %>%
  summarise(spe.bas = mean(spe),
            # speFit.bas = mean(speFit)
            )

exp2.df3 <- exp2.df2 %>% left_join(exp2.df2.bas) %>%
  group_by(subj, runnum, fbblock, ADTile, CTile, STile,
           firstFeedback) %>%
  summarise(
    spe_b = mean(spe) - mean(spe.bas),
    # speFit_b = mean(speFit) - mean(speFit.bas),
    spe = mean(spe),
    # speFit = mean(speFit),
    AD = mean(AD),
    # STile = mean(STile),
    fbpos = mean(fbpos),
    fbdif = mean(fbpos) - mean(fbneg),
  ) 

##################################
##### Figure 3B. Model predictions for Exp 2

f3b <- ggplot(exp2.df3%>% filter(runnum %in% c(3,5)), 
       aes(x = fbblock, y = spe_b, color = ADTile)) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', 'royalblue3')) +
  stat_summary(fun.data = mean_se, size = .6,
               position=position_dodge(width = .04), alpha =.8) +
  stat_summary(fun = mean, size = 1.5, geom = 'line', aes(group=ADTile),
               position=position_dodge(width = .04), alpha =.8) +
  theme_pubclean() +
  theme(text = element_text(size=font_size),
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        # legend.position = 'none',
        legend.text = element_text(size=25), legend.title = element_text(size=25)
        ) +
  coord_cartesian(ylim=c(-.2,.15)) +
  ylab("SPE-b") +
  xlab("Feedback") + labs(color = 'AD score')
f3b

##################################
##### Supp Figure 12B. Regression of individual SRET words on baseline SPE

ggplot(exp2.df3%>% filter(runnum %in% c(3,5)), aes(x = fbblock, y = spe_b, color = CTile)) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', 'royalblue3')) +
  stat_summary(fun.data = mean_se, size = .6,
               position=position_dodge(width = .1), alpha =.8) +
  stat_summary(fun = mean, size = 1.5, geom = 'line', aes(group=CTile),
               position=position_dodge(width = .1), alpha =.8) +
  theme_pubclean() +
  theme(text = element_text(size=font_size),
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = c(.8,.8)) +
  ylab("SPE-b") +
  xlab("Feedback") + labs(color = 'CIT score')

ggplot(exp2.df3%>% filter(runnum %in% c(3,5)), aes(x = fbblock, y = spe_b, color = STile)) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', 'royalblue3')) +
  stat_summary(fun.data = mean_se, size = .6,
               position=position_dodge(width = .1), alpha =.8) +
  stat_summary(fun = mean, size = 1.5, geom = 'line', aes(group=STile),
               position=position_dodge(width = .1), alpha =.8) +
  theme_pubclean() +
  theme(text = element_text(size=font_size),
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = c(.8,.8)) +
  ylab("SPE-b") +
  xlab("Feedback") + labs(color = 'SW score')


exp2.df2 <- exp2.mbs1 %>%
  mutate(ADTile = ntile(AD, 2)) %>%
  mutate(CTile = ntile(Compul, 2)) %>%
  mutate(STile = ntile(SW, 2)) %>%
  drop_na(scale_rt) %>%
  group_by(subj, runnum) %>%
  summarise(spe = mean(spe),
            confZ = mean(confZ),
            ADTile = median(ADTile),
            CTile = median(CTile),
            STile = median(STile)
  ) %>%
  mutate(confTile = ntile(confZ, 6)) %>%
  mutate(hiloConf = ifelse(confTile<=2, 0, 1)) %>%
  mutate_at(c('ADTile', 'hiloConf', 'CTile'), as.factor) %>%
  mutate(hiloConf = recode_factor(hiloConf, '0' = 'lo', '1' = 'hi')) %>%
  mutate(ADTile = recode_factor(ADTile, '1' = 'Low', '2' ='High')) %>%
  mutate(CTile = recode_factor(CTile, '1' = 'Low', '2' = 'High')) %>%
  mutate(STile = recode_factor(STile, '1' = 'Low', '2' = 'High')) #%>% left_join(exp2.speFit)


##################################
##### Figure 3D. Model predictions for Exp 2

f3d <- ggplot(exp2.df2 %>% filter(runnum %in% c(1,2,4,6)), 
       aes(x = runnum, y = spe, color = ADTile)) +
  stat_summary(fun.data = mean_se, size = .6,
               position=position_dodge(width = .04), alpha=.8) +
  scale_color_manual(breaks = c('Low', 'High'), 
                     values = c('tan3', 'royalblue3')) +
  stat_summary(fun = mean, geom = 'line', size=1.5,
               position=position_dodge(width = .04), alpha=.8) +
  theme_pubclean() +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none') +
  ylab("SPE") +
  coord_cartesian(ylim=c(.47,.65)) +
  scale_x_continuous(breaks = c(1,2,4,6))+
  xlab("Block #") +
  labs(colour = 'AD score', linetype = 'Conf level') 
f3d

ggplot(exp2.df2 %>% filter(runnum %in% c(1,2,4,6)), 
       aes(x = runnum, y = spe, color = CTile)) +
  stat_summary(fun.data = mean_se, size = .6,
               position=position_dodge(width = .1), alpha=.8) +
  scale_color_manual(breaks = c('Low', 'High'), 
                     values = c('tan3', 'royalblue3')) +
  stat_summary(fun = mean, geom = 'line', size=1.5,
               position=position_dodge(width = .1), alpha=.8) +
  theme_pubclean() +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none') +
  ylab("SPE") +
  coord_cartesian(ylim=c(.47,.65)) +
  scale_x_continuous(breaks = c(1,2,4,6)) +
  xlab("Block #") +
  labs(colour = 'CIT score', linetype = 'Conf level')

ggplot(exp2.df2 %>% filter(runnum %in% c(1,2,4,6)), 
       aes(x = runnum, y = spe, color = STile)) +
  stat_summary(fun.data = mean_se, size = .6,
               position=position_dodge(width = .04), alpha=.8) +
  scale_color_manual(breaks = c('Low', 'High'), 
                     values = c('tan3', 'royalblue3')) +
  stat_summary(fun = mean, geom = 'line', size=1.5,
               position=position_dodge(width = .04), alpha=.8) +
  theme_pubclean() +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none') +
  ylab("SPE") +
  coord_cartesian(ylim=c(.47,.65)) +
  scale_x_continuous(breaks = c(1,2,4,6)) +
  xlab("Block #") +
  labs(colour = 'AD score', linetype = 'Conf level')


##
##
library(scales)
if (Sys.info()["sysname"]=='Darwin'){
  setwd("~/OneDrive - University College London/Projects/Experiments/metaBiasShift/")
} else{
  setwd("C:/Users/skatyal/OneDrive - University College London/Projects/Experiments/metaBiasShift/")
}
## read data file
exp2.mbs3 = read.csv(paste('data/exp2/mbsExp2_noperfexcl.csv', sep=''),
                    header = TRUE, sep = ',')
## recode columns as factors
exp2.mbs3$subj <- as.factor(exp2.mbs3$subj)
exp2.mbs3$task <- as.factor(exp2.mbs3$task)
exp2.mbs3$gender[is.nan(exp2.mbs3$gender)] <- NA
exp2.mbs3$gender <- as.factor(exp2.mbs3$gender)

exp2.mbs3 <- exp2.mbs3 %>%
  mutate_at(vars(contains("fb")), as.factor) %>%
  rename(accu = corr) %>%
  rename(stair = incdec)

exp2.mbs3 <- exp2.mbs3 %>%
  mutate(fbblock = recode_factor(fbblock, "1" = posnegNames[1], 
                                 "2" = posnegNames[2], "0" = "None"),
         task = recode_factor(task, "0" = taskNames[1], "1" = taskNames[2])
  )

## code highly deviant RTs as NA
max_RT_deviation = 3
rt1_max = c(median(exp2.mbs3$rt1[exp2.mbs3$task==taskNames[1]], na.rm=T)+ 
              max_RT_deviation*IQR(exp2.mbs3$rt1[exp2.mbs3$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(exp2.mbs3$rt1[exp2.mbs3$task==taskNames[2]], na.rm=T)+ 
              max_RT_deviation*IQR(exp2.mbs3$rt1[exp2.mbs3$task==taskNames[2]], na.rm=T ))
exp2.mbs3 <- filter(exp2.mbs3, trialnum > nLags) %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[1] & .>rt1_max[1],NA)) %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[2] & .>rt1_max[2],NA))

## scale RTs 
exp2.mbs3$scale_rt <- exp2.mbs3$rt1
exp2.mbs3$scale_rt[exp2.mbs3$task==taskNames[1]] <- scale(exp2.mbs3$rt1[exp2.mbs3$task==taskNames[1]])
exp2.mbs3$scale_rt[exp2.mbs3$task==taskNames[2]] <- scale(exp2.mbs3$rt1[exp2.mbs3$task==taskNames[2]])

rt2_max = c(median(exp2.mbs3$rt2[exp2.mbs3$task==taskNames[1]], na.rm=T)+ 
              max_RT_deviation*IQR(exp2.mbs3$rt2[exp2.mbs3$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(exp2.mbs3$rt2[exp2.mbs3$task==taskNames[2]], na.rm=T)+ 
              max_RT_deviation*IQR(exp2.mbs3$rt2[exp2.mbs3$task==taskNames[2]], na.rm=T ))
exp2.mbs3 <- filter(exp2.mbs3, trialnum > nLags) %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[1] & .>rt2_max[1],NA)) %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[2] & .>rt2_max[2],NA))

exp2.mbs3$scale_rt2 <- exp2.mbs3$rt2
exp2.mbs3$scale_rt2[exp2.mbs3$task==taskNames[1]] <- scale(exp2.mbs3$rt2[exp2.mbs3$task==taskNames[1]])
exp2.mbs3$scale_rt2[exp2.mbs3$task==taskNames[2]] <- scale(exp2.mbs3$rt2[exp2.mbs3$task==taskNames[2]])

# exp2.mbs31 <- exp2.mbs3

## if to subtract the baseline blocks confidence and spe from the rest of the blocks
# get the baseline values for the two tasks from runs 1 and 2
exp2.mbs3.bas <- filter(exp2.mbs3, runnum<3) %>%
  group_by(subj, task, group) %>%
  dplyr::summarise(
    conf.bas = mean(conf),
    spe.bas = mean(spe),
    accu.bas = mean(accu)
  )

exp2.mbs3 <- exp2.mbs3 %>% left_join(exp2.mbs3.bas)

if (subtractBaselineConfSpe){
  # left join the baseline values an subtract them from confidence and spe
  exp2.mbs3 <- exp2.mbs3 %>%
    mutate(accu_b = accu - accu.bas) %>%
    mutate(conf = conf-conf.bas)%>%
    mutate(spe = spe-spe.bas) 
}

## take away the baseline (first two) blocks from the dataset
exp2.mbs3 <- exp2.mbs3 %>%
  filter(runnum>2) #%>%
# mutate(runnum = runnum-2)
exp2.mbs3 <- exp2.mbs3 %>%
  mutate(blockType = factor(runnum%%2)) %>%
  mutate(blockType = recode_factor(blockType, "0" = "transf", "1" = "interv"))
exp2.mbs3 <- exp2.mbs3 %>%
  mutate_at(vars(contains("runnum")), as.factor)
exp2.mbs3$fbblock <- droplevels(exp2.mbs3$fbblock)


## create new columns for trial-by-trial positive and negative feedback
exp2.mbs3 <- exp2.mbs3 %>%
  mutate(fbpos=feedback, fbneg=feedback) %>%
  mutate_at("fbpos", ~replace(., feedback==1 & accu==1, 1)) %>%
  mutate_at("fbpos", ~replace(., feedback==1 & accu==0, 0)) %>%
  mutate_at("fbneg", ~replace(., feedback==1 & accu==1, 0)) %>%
  mutate_at("fbneg", ~replace(., feedback==1 & accu==0, 1)) 

# give both interv and transf run the same value of fbblock (feedback block type)
exp2.mbs3$fbblock[exp2.mbs3$runnum==2] = exp2.mbs3$fbblock[exp2.mbs3$runnum==1 & 
                                                          exp2.mbs3$trialnum<=20]
exp2.mbs3$fbblock[exp2.mbs3$runnum==4] = exp2.mbs3$fbblock[exp2.mbs3$runnum==3
                                                        & exp2.mbs3$trialnum<=20]

## code intervention task by group
exp2.mbs3$intervTask = exp2.mbs3$group
exp2.mbs3$intervTask[exp2.mbs3$group %in% c(1,2,3,4)] = taskNames[1]
exp2.mbs3$intervTask[exp2.mbs3$group %in% c(5,6,7,8)] = taskNames[2]

# transfer task by group (same as task column)
exp2.mbs3$testTask = exp2.mbs3$group
exp2.mbs3$testTask[exp2.mbs3$group %in% c(1,2,3,4)] = taskNames[2]
exp2.mbs3$testTask[exp2.mbs3$group %in% c(5,6,7,8)] = taskNames[1]

## code whether pos or neg is the first block of feedback
exp2.mbs3$firstFeedback = exp2.mbs3$group
exp2.mbs3$firstFeedback[exp2.mbs3$group %in% c(1,3,5,7)] = posnegNames[1]
exp2.mbs3$firstFeedback[exp2.mbs3$group %in% c(2,4,6,8)] = posnegNames[2]


# exp2.mbs3$group <- as.factor(exp2.mbs3$group)
exp2.mbs3 <- exp2.mbs3 %>%
  mutate_at(c("group", "intervTask", "firstFeedback", 'testTask'), as.factor) %>%
  mutate(group = recode_factor(group, "1" = 'Group 1', 
                               "2" = 'Group 2', "3" = 'Group 3', "4" = 'Group 4', 
                               "5" = 'Group 5', "6" = 'Group 6', "7" = 'Group 7', 
                               "8" = 'Group 8')) 

# #####################
# #####################
# ## read in the psychiatric scores
# ## read data file
psych = read.csv(paste('data/exp2/factor_scores_noperfexc.csv', sep=''),
                 header = TRUE, sep = ',')
psych <- psych %>% rename(subj = subjIDs) %>% dplyr::select(-X) %>%
  mutate_at("subj", as.factor)

exp2.mbs3 <- exp2.mbs3 %>% left_join(psych)

####################
####################
#### A5.2 sret analysis
posnegWordColours.exp2 =  c('springgreen3','maroon2')

exp2.df.2 <- filter(exp2.mbs3, runnum==3) %>%
  pivot_longer(cols = starts_with("wordRT_"),
               names_to = "wordnum",
               names_prefix = "wordRT_",
               values_to = "sretRT") %>%
  dplyr::select(-starts_with("word_"))%>%
  dplyr::select(-starts_with("qn_")) %>%
  group_by(subj, wordnum, sretset, intervTask, task, fbblock, group) %>%
  summarise(
    accu = mean(accu),
    stair = mean(stair),
    sretRT = mean(sretRT),
    AD = mean(AD),
    Compul = mean(Compul),
    SW = mean(SW),
    spe.bas = mean(spe.bas),
    conf.bas = mean(conf.bas),
    fbneg = mean(fbneg),
    fbpos = mean(fbpos),
    spe = mean(spe)
  )
  
exp2.temp <- filter(exp2.mbs3, runnum==3) %>%
  pivot_longer(cols = starts_with("word_"),
               names_to = "wordnum",
               names_prefix = "word_",
               values_to = "sretVal")  %>%
  dplyr::select(-starts_with("wordRT_"))%>%
  dplyr::select(-starts_with("qn_")) %>%
  group_by(subj, wordnum, sretset) %>%
  summarise(
    accu = mean(accu),
    sretVal = mean(sretVal)
  )

exp2.df.2 <- exp2.df.2 %>% left_join(exp2.temp) 

exp2.df.2$wordvalence = exp2.df.2$wordnum
exp2.df.2$wordvalence[exp2.df.2$wordnum %in% c(1:10, 21:30)] <- 'Positive'
exp2.df.2$wordvalence[exp2.df.2$wordnum %in% c(11:20, 31:40)] <- 'Negative'
exp2.df.2$sretTimepoint = exp2.df.2$wordnum
exp2.df.2$sretTimepoint[exp2.df.2$sretset==1 & exp2.df.2$wordnum %in% c(1:20)] <- 1
exp2.df.2$sretTimepoint[exp2.df.2$sretset==1 & exp2.df.2$wordnum %in% c(21:40)] <- 2
exp2.df.2$sretTimepoint[exp2.df.2$sretset==2 & exp2.df.2$wordnum %in% c(1:20)] <- 2
exp2.df.2$sretTimepoint[exp2.df.2$sretset==2 & exp2.df.2$wordnum %in% c(21:40)] <- 1

outlierRTCutoff = median(exp2.df.2$sretRT, na.rm = T) + 5*IQR(exp2.df.2$sretRT, na.rm = T)
exp2.df.2$sretRT[exp2.df.2$sretRT > outlierRTCutoff] = NA

exp2.df.2 <- exp2.df.2 %>% 
  mutate(wordvalence = ordered(wordvalence, levels = c("Positive", "Negative")))

include.words <- c(1:40)
exp2.df.3 <- filter(exp2.df.2, wordnum %in% include.words) %>% 
  group_by(subj, wordvalence, sretTimepoint,
                                    intervTask, task, fbblock, group,
                                    sretset) %>%
  summarise(sretRT = median(sretRT),
            sretVal = mean(sretVal),
            AD = mean(AD),
            accu = mean(accu),
            stair = mean(stair),
            fbneg = mean(fbneg),
            fbpos = mean(fbpos),
            spe = mean(spe)
            )


exp2.df.4 <- exp2.df.3 %>% pivot_wider(
  names_from = sretTimepoint,
  values_from = c(sretVal, sretRT)
) %>% mutate(sretValDiff = sretVal_2 - sretVal_1) %>%
  mutate(sretRTDiff = sretRT_2 - sretRT_1)

exp2.df.4 <- exp2.df.4 %>% drop_na(sretRTDiff)

m.reg <- glmer(sretVal/4 ~ AD*wordvalence + 
                Compul*wordvalence +
                SW*wordvalence + (1|wordnum), 
               family = binomial, control=glmerControl(optimizer = 'bobyqa'),
              filter(exp2.df.2, sretTimepoint==1))
summary(m.reg)

##################################
##### Figure 5A. SRET Figure

f5a <- plot_model(m.reg, type = ("pred"),
           terms = c("AD", "wordvalence"), 
           line.size = 1.5, ci.lvl = .95,
           title = '', colors = posnegWordColours.exp2) +
  ylab("Self-endorsement") +
  xlab("AD axis") +
  labs(colour = 'Valence') +
  theme_pubclean() +
  coord_cartesian(ylim = c(.0,.9))   +
  scale_y_continuous(breaks = seq(0, .9, len = 4), labels = percent)+
  annotate('text', x=0, y=.5, label="p < .0001", size=6) +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.position = 'none')

f5a

emm <- emtrends(m.reg, ~wordvalence, "AD")
test(emm)

##################################
##### Supp Figure 11. SRET
sf11a <- plot_model(m.reg, type = ("pred"),
           terms = c("Compul", "wordvalence"), 
           line.size = 1.5,ci.lvl = .95,
           title = '', colors = posnegWordColours.exp2) +
  ylab("Self-endorsement") +
  xlab("CIT axis") +
  labs(colour = 'Valence') +
  theme_pubclean() +
  coord_cartesian(ylim = c(.0,.9))   +
  scale_y_continuous(breaks = seq(0, .9, len = 4), labels = percent)+
  annotate('text', x=0, y=.5, label="p < .0001", size=6) +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.position = 'none')
sf11a
emm <- emtrends(m.reg, ~wordvalence, "Compul")
test(emm)

sf11b <- plot_model(m.reg, type = ("pred"),
           terms = c("SW", "wordvalence"), 
           line.size = 1.5,ci.lvl = .95,
           title = '', colors = posnegWordColours.exp2) +
  ylab("Self-endorsement") +
  xlab("SW axis") +
  labs(colour = 'Valence') +
  theme_pubclean() +
  coord_cartesian(ylim = c(.0,.9))   +
  scale_y_continuous(breaks = seq(0, .9, len = 4), labels = percent)+
  annotate('text', x=0, y=.5, label="p < .0001", size=6) +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.position = 'none')
sf11b

emm <- emtrends(m.reg, ~wordvalence, "SW")
test(emm)

m.reg <- glmer(sretVal/4 ~ spe.bas*wordvalence + conf.bas*wordvalence + 
                 (1|wordnum) , 
               family = binomial(link = 'logit'), control = glmerControl(optimizer = c('bobyqa')),
               filter(exp2.df.2, sretTimepoint=='1'))
summary(m.reg)
# 

sf11c<-plot_model(m.reg, type = ("pred"),
           terms = c("spe.bas", "wordvalence"), 
           line.size = 1.5, ci.lvl = .95,
           title = '', colors = posnegWordColours.exp2) +
  theme_pubclean() +
  ylab("Self-endorsement") +
  xlab("SPE (baseline)") +
  labs(colour = 'Word valence') +
  coord_cartesian(ylim = c(.0,.9))   +
  annotate('text', x=0.5, y=.5, label="p = .018", size=6) +
  scale_y_continuous(breaks = seq(0, .9, len = 4), labels = percent)+
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.position = 'none')
sf11c

emm <- emtrends(m.reg, ~wordvalence, "spe.bas")
test(emm)

f5b<-plot_model(m.reg, type = ("pred"),
           terms = c("conf.bas", "wordvalence"), 
           line.size = 1.5, ci.lvl = .95,
           title = '', colors = posnegWordColours.exp2) +
  ylab("Self-endorsement") +
  xlab("Local confidence (baseline)") +
  labs(colour = 'Valence') +
  coord_cartesian(ylim = c(.0,.9))   +
  scale_y_continuous(breaks = seq(0, .9, len = 4), labels = percent)+
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  theme_pubclean() +
  annotate('text', x=0.5, y=.5, label="p < .0001", size=6) +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.position = 'none')
f5b
emm <- emtrends(m.reg, ~wordvalence, "conf.bas")
test(emm)



f5c <- ggplot(exp2.df.4, aes(x = fbblock, y = sretValDiff, color = wordvalence)) +
  geom_hline(yintercept = 0, color = 'grey30') +
  scale_color_manual(breaks = posnegNames, values=posnegWordColours.exp2) +
  geom_quasirandom(dodge.width=.85, alpha = .7) +
  geom_violin(alpha=.7, position = position_dodge(.85)) + theme_pubclean() +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .85),
               geom = 'errorbar',
               size = 1.1, aes(width = .2))  +
  scale_x_discrete(limits = posnegNames) +
  ylab("Self-endorsement\ndifference (T2 - T1)") +
  xlab("Feedback") +
  labs(color = 'Word Valence') + 
  coord_cartesian(ylim = c(-1.7, 1.8)) +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  geom_signif(y_position = c(1.55, 1.55),
              comparisons = list(c("Positive", "Negative")),
              map_signif_level = TRUE,color = 'black', vjust = -.25,
              textsize = pval_size, size = .8,
              annotations = c("p = .014")) +
  geom_signif(y_position = c(1.25, 1.25), color = 'black',
              xmin = c(0.8, 1.8), vjust = -.25,
              textsize = pval_size, size = .8,
              xmax = c(1.2, 2.2),
              annotation = c("n.s.", "p < .0001"), tip_length = 0) +
  geom_signif(y_position = c(.95), color = 'black',
              xmin = c(1.2), vjust = -.2, hjust = 1,
              textsize = pval_size, size = .8,
              xmax = c(2.2),
              annotation = c("p = .03"), tip_length = 0)
f5c

m.reg <- lmer(sretValDiff ~ fbblock*wordvalence + accu + stair +
                (1|sretset) + (1|task) + (1|subj), exp2.df.4)
summary(m.reg)
m.reg <- lmer(sretValDiff ~ fbblock*wordvalence + accu + stair +
                (1|sretset) , exp2.df.4)
summary(m.reg)

plot(m.reg)
plot_model(m.reg, type = ("pred"),
           terms = c("fbblock", "wordvalence"), 
           line.size = 1.5,
           title = '') +
  ylab("Self-endorsement\ndifference (T2 - T1)") +
  xlab("Feedback") +
  labs(colour = 'Word valence') +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  theme_pubclean()

plot_model(m.reg, type = ("pred"),
           terms = c("wordvalence", "fbblock"), 
           line.size = 1.5,
           title = '') +
  scale_color_manual(breaks = posnegNames,
                     values=posnegColours) +
  ylab("Self-endorsement (T2 - T1)") +
  xlab("Word valence") +
  labs(colour = 'Feedback') +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  theme_pubclean()

emm <- emmeans(m.reg, ~wordvalence|fbblock)
contrast(emm)
emm <- emmeans(m.reg, ~fbblock|wordvalence)
contrast(emm)


## preregistered sret pre-post analysis - in Supplementary Results
exp2.df.5 <- exp2.df.4 %>% 
  group_by(subj, wordvalence, task, 
           fbblock, sretset) %>%
  summarise(sretVal = mean(sretValDiff),
            accu = mean(accu),
            stair = mean(stair)) %>%
  pivot_wider(    names_from = wordvalence,
                  values_from = sretVal) %>%
  mutate(sretDiff = Positive - Negative)

m.reg <- lmer(sretDiff ~ fbblock*task + 
                (1|sretset) ,
              exp2.df.5)
summary(m.reg)

m.reg <- lmer(sretDiff ~ fbblock + (1|sretset) ,
              exp2.df.5)
summary(m.reg)

ggplot(exp2.df.5, aes(x = fbblock, y = sretDiff, color = fbblock)) +
  scale_color_manual(breaks = posnegNames, values=posnegColours) +
  geom_quasirandom(dodge.width=.85, alpha = .7) +
  geom_violin(alpha=.7, position = position_dodge(.85)) + theme_pubclean() +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .85),
               geom = 'errorbar',
               size = 1.1, aes(width = .2))  +
  scale_x_discrete(limits = posnegNames) +
  ylab("Self-endorsement double diff.\n(Pos - Neg) & (T2 - T1)") +
  xlab("Feedback") +
  labs(color = 'Word Valence') + 
  coord_cartesian(ylim = c(-2.6, 2)) +
  theme(text = element_text(size=font_size), 
        legend.position = 'none',
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  geom_signif(y_position = c(1.7), color = 'black',
              xmin = c(1), vjust = -.2, #hjust = 1,
              textsize = pval_size, size = .8,
              xmax = c(2),
              annotation = c("p = .04"), tip_length = 0)
