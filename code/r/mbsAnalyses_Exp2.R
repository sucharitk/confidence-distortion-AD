
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
library(BayesFactor)

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


setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/metaBiasShift/")

## read data file
exp2.mbs = read.csv(paste('data/exp2/processed/mbsExp2.csv', sep=''), header = TRUE, sep = ',')
## recode columns as factors
exp2.mbs$subj <- as.factor(exp2.mbs$subj)
exp2.mbs$task <- as.factor(exp2.mbs$task)
exp2.mbs$awarepos[exp2.mbs$awarepos==0] <- 'NA'
exp2.mbs$awarepos[is.nan(exp2.mbs$awarepos)] <- 'NA'
exp2.mbs$awareneg[exp2.mbs$awareneg==0] <- 'NA'
exp2.mbs$awareneg[is.nan(exp2.mbs$awareneg)] <- 'NA'
exp2.mbs$awarepos <- as.factor(exp2.mbs$awarepos)
exp2.mbs$awareneg <- as.factor(exp2.mbs$awareneg)
exp2.mbs$affectpos[is.nan(exp2.mbs$affectpos)] <- 'NA'
exp2.mbs$affectneg[is.nan(exp2.mbs$affectneg)] <- 'NA'
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
         awarepos = recode_factor(awarepos, '1' = 'Yes', '2'='No', 
                                  'NA' = 'No response', 'NaN' = 'No response'),
         awareneg = recode_factor(awareneg, '1' = 'Yes', '2'='No', 
                                  'NA' = 'No response', 'NaN' = 'No response'),
         affectpos = recode_factor(affectpos, '1' = 'Felt better', '2'='Felt worse', 
                                   '3'='No change', '4'='Not sure', 
                                   'NA' = 'No response', 'NaN' = 'No response'),
         affectneg = recode_factor(affectneg, '1' = 'Felt better', '2'='Felt worse', 
                                   '3'='No change', '4'='Not sure', 
                                   'NA' = 'No response', 'NaN' = 'No response'))

# ## check when in block was feedback delivered
# exp2.fb <- filter(exp2.mbs, runnum %in% c(3,5)) %>% 
#   group_by(fbblock, trialnum) %>%
#   summarise(feedback = mean(feedback)) %>%
#   pivot_wider(names_from = fbblock,
#               values_from = feedback) %>%
#   mutate(diff = Positive - Negative) %>% 
#   mutate(expnum = as.factor(2))
# 
# bothexp.fb <- bind_rows(exp1.fb, exp2.fb)
# ggplot(bothexp.fb, aes(x=trialnum, y=diff, colour = expnum)) +
#   geom_smooth(method = 'loess', se=F) +
#   labs(color = 'Experiment number') +
#   ylab('Mean feedback\ndifference (pos-neg)')+
#   xlab('Trial number') +
#   theme_pubclean()

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

exp2.mbs$stair[exp2.mbs$task==taskNames[1]] <- scale(exp2.mbs$stair[exp2.mbs$task==taskNames[1]])
exp2.mbs$stair[exp2.mbs$task==taskNames[2]] <- scale(exp2.mbs$stair[exp2.mbs$task==taskNames[2]])

# exp2.mbs1 <- exp2.mbs

## if to subtract the baseline blocks confidence and spe from the rest of the blocks
  # get the baseline values for the two tasks from runs 1 and 2
  exp2.mbs2 <- filter(exp2.mbs, runnum<3) %>%
    group_by(subj, task, group) %>%
    dplyr::summarise(
      conf.bas = mean(conf),
      spe.bas = mean(spe),
      accu.bas = mean(accu),
      stair.bas = mean(stair)
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
    # mutate(accu_b = accu - accu.bas) %>%
    mutate(conf_b = conf-conf.bas) %>%
    mutate(spe_b = spe-spe.bas) %>%
    mutate(stair_b = stair - stair.bas)
}

## take away the baseline (first two) blocks from the dataset
# exp2.mbs <- exp2.mbs %>%
#   filter(runnum>2) #%>%

exp2.mbs <- exp2.mbs %>%
  mutate(blockType = factor(runnum%%2)) %>%
  mutate(blockType = recode_factor(blockType, "0" = "transf", "1" = "interv"))

levels(exp2.mbs$blockType) <- c(levels(exp2.mbs$blockType), 'base')
exp2.mbs$blockType[exp2.mbs$runnum %in% c(1,2)] <- 'base'

exp2.mbs <- exp2.mbs %>%
  mutate_at(vars(contains("runnum")), as.factor)
exp2.mbs$fbblock <- droplevels(exp2.mbs$fbblock)


## create new columns for trial-by-trial positive and negative feedback
exp2.mbs <- exp2.mbs %>%
  mutate(fbpos=feedback, fbneg=feedback) %>%
  mutate_at("fbpos", ~replace(., feedback==1 & accu==1, 1)) %>%
  mutate_at("fbpos", ~replace(., feedback==1 & accu==0, 0)) %>%
  mutate_at("fbneg", ~replace(., feedback==1 & accu==1, 0)) %>%
  mutate_at("fbneg", ~replace(., feedback==1 & accu==0, 1))  %>%
  mutate_at('feedback', ~replace(., feedback==1 & accu==0, -1)) %>%
  dplyr::select(-c(rt2)) %>%
  mutate(fbblock.notrans = fbblock) %>%
  mutate_at('fbblock.notrans', ~replace(., runnum==4 | runnum==6, 'None'))

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
psych = read.csv(paste('data/exp2/processed/factor_scores.csv', sep=''),
                 header = TRUE, sep = ',')
psych <- psych %>% rename(subj = subjIDs) %>% dplyr::select(-X) %>%
  mutate_at("subj", as.factor)

exp2.mbs <- exp2.mbs %>% left_join(psych)
exp2.mbs.bas <- exp2.mbs.bas %>% left_join(psych)

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
    conf_b = mean(conf_b),
    accu_b = mean(accu)-mean(accu.bas),
    accu = mean(accu),
    scale_rt = mean(scale_rt, na.rm=T),
    rt1 = mean(rt1, na.rm=T),
    fbneg = mean(fbneg),
    fbpos = mean(fbpos),
    spe = mean(spe),
    spe_b = mean(spe_b),
    spe0_perc = mean(spe0_perc),
    spe0_mem = mean(spe0_mem),
    # spe.bas = mean(spe.bas),
    age = mean(age),
    stair = mean(stair),
    stair_b = mean(stair_b),
    AD = mean(AD),
    Compul = mean(Compul),
    SW = mean(SW))


#####################
# exp2.mbs1 <- exp2.mbs1 %>% left_join(psych)


######################################
##### Supp Figure 4. SPE of intervention block

exp2.df.2 <- filter(exp2.df, blockType=="interv")

sf4.left <- ggplot(exp2.df.2, aes(x = fbblock, y = spe_b, color = fbblock)) +
  geom_hline(yintercept = 0, color = 'grey30') +
  facet_wrap(~task) +
  scale_color_manual(breaks = posnegNames, values=posnegColours) +
  theme_pubclean()  +
  # scale_x_discrete(limits = taskNames) +
  geom_quasirandom(dodge.width=.7, alpha = .7) +
  geom_violin(alpha=.6, position = position_dodge(.7)) +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .7),
               geom = 'errorbar',
               size = .9, aes(width = .2)) +
  ylab("Self-performance\nestimate\n(global SPE,\nbaseline-corrected)") +
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
  geom_signif(y_position = c(.55), vjust = -.3, hjust=.4,
              xmin = c(1.1), color = 'black',
              xmax = c(1.9),
              annotation = c("p < .0001"), tip_length = .04,
              textsize = pval_size, size = .8)
sf4.left

sf6.right <- ggplot(exp2.df.2, aes(x = fbblock, y = accu, color = fbblock)) +
  # geom_hline(yintercept = 0, color = 'grey30') +
  facet_wrap(~task) +
  scale_color_manual(breaks = posnegNames, values=posnegColours) +
  theme_pubclean()  +
  # scale_x_discrete(limits = taskNames) +
  geom_quasirandom(dodge.width=.7, alpha = .7) +
  geom_violin(alpha=.6, position = position_dodge(.7)) +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .7),
               geom = 'errorbar',
               size = .9, aes(width = .2)) +
  ylab("\nActual performance\n(accuracy)\n") +
  xlab("Feedback") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.position = 'none') + 
  # coord_cartesian(ylim = c(0,1)) +
  labs(color = 'Feedback') + 
  # labs(tag = "A") +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  geom_signif(y_position = c(.85), vjust = -.3,
              xmin = c(1.1), color = 'black',
              xmax = c(1.9),
              annotation = c("n.s."), tip_length = .04,
              textsize = pval_size, size = .8)
sf6.right

sf6.conf <- ggplot(exp2.df.2, aes(x = fbblock, y = conf_b, color = fbblock)) +
  # geom_hline(yintercept = 0, color = 'grey30') +
  facet_wrap(~task) +
  scale_color_manual(breaks = posnegNames, values=posnegColours) +
  theme_pubclean()  +
  # scale_x_discrete(limits = taskNames) +
  geom_quasirandom(dodge.width=.7, alpha = .7) +
  geom_violin(alpha=.6, position = position_dodge(.7)) +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .7),
               geom = 'errorbar',
               size = .9, aes(width = .2)) +
  ylab("\nLocal confidence\n(baseline-corrected)\n") +
  xlab("Feedback") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        legend.position = 'none') + 
  coord_cartesian(ylim = c(-.45,.45)) +
  labs(color = 'Feedback') + 
  # labs(tag = "A") +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  geom_signif(y_position = c(.35), vjust = -.3, hjust=.4,
              xmin = c(1.1), color = 'black',
              xmax = c(1.9),
              annotation = c("p = .005"), tip_length = .04,
              textsize = pval_size, size = .8)
sf6.conf

##################
### Plot raw data across all blocks
exp2.df11 <- exp2.mbs %>%
  group_by(subj, task, runnum, fbblock.notrans, group) %>%
  dplyr::summarise(
    conf = mean(conf),
    accu = mean(accu, na.rm=T),
    spe = mean(spe),
    stair = mean(stair))

ggplot(exp2.df11, aes(x = runnum, y = spe, color = fbblock.notrans)) +
  facet_wrap(~group, ncol = 2) +
  theme_pubclean()  +
  geom_quasirandom(alpha = .4) +
  geom_violin(alpha=.4) +
  stat_summary(fun.data = mean_cl_boot,
               geom = 'errorbar',
               size = .9, aes(width = .2)) +
  scale_color_manual(breaks = c(posnegNames, 'None'), 
                     values=c(posnegColours, 'grey30')) +
  ylab("Self-performance estimate (global SPE)") +
  xlab("Block number") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  labs(color = 'Feedback') + 
  theme(text = element_text(size=font_size), 
        legend.position = 'none',
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) 

ggplot(exp2.df11, aes(x = runnum, y = accu, color = fbblock.notrans)) +
  facet_wrap(~group, ncol = 2) +
  theme_pubclean()  +
  geom_quasirandom(alpha = .4) +
  geom_violin(alpha=.4) +
  stat_summary(fun.data = mean_cl_boot,
               geom = 'errorbar',
               size = .8, aes(width = .4)) +
  scale_color_manual(breaks = c(posnegNames, 'None'), 
                     values=c(posnegColours, 'grey30')) +
  ylab("Actual performance (accuracy)") +
  xlab("Block number") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  labs(color = 'Feedback') + 
  # coord_cartesian(ylim = c(.5,1)) +
  theme(text = element_text(size=font_size), 
        legend.position = 'none',
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"))

ggplot(exp2.df11, aes(x = runnum, y = conf, color = fbblock.notrans)) +
  facet_wrap(~group, ncol = 2) +
  theme_pubclean()  +
  geom_quasirandom(alpha = .4) +
  geom_violin(alpha=.4) +
  stat_summary(fun.data = mean_cl_boot,
               geom = 'errorbar',
               size = .8, aes(width = .4)) +
  scale_color_manual(breaks = c(posnegNames, 'None'), 
                     values=c(posnegColours, 'grey30')) +
  ylab("Mean confidence") +
  xlab("Block number") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  labs(color = 'Feedback') + 
  theme(text = element_text(size=font_size), 
        legend.position = 'none',
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"))


######################################
##### Supp Figure 5 (lower). Task difficulty during intervention blocks

sf5.bottom <- ggplot(exp2.df.2, aes(x = task, y = stair_b, color = fbblock)) +
  # geom_hline(yintercept = .74, color = 'grey30') +
  scale_color_manual(breaks = posnegNames, values=posnegColours) +
  theme_pubclean()  +
  scale_x_discrete(limits = taskNames) +
  geom_quasirandom(dodge.width=.7, size = .7) +
  geom_violin(alpha=.6, position = position_dodge(.7)) +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .7),
               geom = 'errorbar',
               size = .9, aes(width = .2)) +
  ylab("Mean staircase level\n(z-scored / task,\nbaseline subtracted) ") +
  xlab("Intervention task") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  ggtitle('Exp 2') +
  labs(color = 'Feedback') + 
  theme(text = element_text(size=font_size-4), 
        plot.caption = element_text(size = font_size-4),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  coord_cartesian(ylim = c(-3,2.5))  +
  geom_signif(y_position = c(2, 2), vjust = -.3,
              xmin = c(0.9, 1.9), color = 'black',
              xmax = c(1.1, 2.1),
              annotation = c("n.s.", "n.s."), tip_length = .04,
              textsize = pval_size-1, size = .5)
sf5.bottom

######################################
##### Supp Figure 9. Transfer effect of feedback upon test block confidence

exp2.df.2 <- filter(exp2.mbs, blockType=="transf") %>% 
  mutate(intervTask = relevel(intervTask, "Perception")) %>%
  mutate(testTask = relevel(testTask, "Perception"))

sf9a <- ggplot(exp2.df.2, aes(x = trialnum, y = conf_b, colour = fbblock)) +
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
  summarise(conf = mean(conf_b)) %>%
  pivot_wider(names_from = fbblock,
    values_from = conf) %>% group_by(subj, intervTask, testTask) %>%
  summarise(confPosNeg = mean(Positive) - mean(Negative))%>% 
  mutate(intervTask = relevel(intervTask, "Perception")) %>%
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
  coord_cartesian( ylim = (c(-.3,.5))) +
  theme_pubclean() +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  # geom_signif(map_signif_level = TRUE,
  #             annotations = c("p = .004"),vjust = -.3,
  #             y_position = c(.42),
  #             xmin = 1.03, xmax = 1.97,color= 'black',
  #             textsize = pval_size, size = .8,
  #             tip_length = .03) +
  geom_signif(y_position = c(.4, .32), vjust = -.3,
              xmin = c(.85, 1.85),
              xmax = c(1.15, 2.15),
              annotation = c("p = .002", "p = .999"),  color= 'black',
              textsize = pval_size-1, size = .8,
              tip_length = .0)
sf9b

######################################
##### Supp Figure 14b. transfer of feedback to test block SPE
interv.labs <-  c('Intervention: Perception', 'Intervention: Memory')
exp2.df.2 <- exp2.mbs %>% 
  filter(blockType=="transf") %>% 
  group_by(fbblock, subj, intervTask, testTask) %>%
  summarise(spe_b = mean(spe_b)) %>%
  pivot_wider(names_from = fbblock,
    values_from = spe) %>% group_by(subj, intervTask, testTask) %>%
  summarise(spe = mean(Positive) - mean(Negative)) %>% 
  mutate(intervTask = relevel(intervTask, "Perception")) %>%
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
  coord_cartesian(ylim = c(-.5,.8)) +
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

m.reg <- lmer(spe_b ~ fbblock*task*firstFeedback + accu_b + stair_b +
                (1|subj) + (1|group) + (1|runnum), REML = F, exp2.df.2)
summary(m.reg)
m.reg <- lmer(spe_b ~ fbblock*task*firstFeedback + accu_b + stair_b  +
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
emm <- emmeans(m.reg, spe_bcs = c("fbblock"))
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

m1 <- lmer(spe_b ~ fbblock*task + accu_b + stair_b + (1|subj), 
           REML = F, exp2.df.2)
summary(m1)
m2 <- lmer(spe_b ~ fbblock+task + accu_b + stair_b + (1|subj), 
           REML = F, exp2.df.2)
anova(m1,m2)

library(lmtest)
m1 <- lm(spe_b ~ fbblock + task + accu_b + stair_b, exp2.df.2)
summary(m1)
confint(m1)
m2 <- lm(spe_b ~ task + accu_b + stair_b, exp2.df.2)
lrtest(m2,m1)

cohens_d(spe_b ~ fbblock, data = exp2.df.2)


m1 <- lmer(accu_b ~ fbblock*task + (1|subj), exp2.df.2)
summary(m1)
emm <- emmeans(m1, ~fbblock|task)
contrast(emm)

m1 <- lmer(accu_b ~ fbblock+task + (1|subj), exp2.df.2)
summary(m1)
confint(m1)

m2 <- lmer(accu_b ~ task + (1|subj), exp2.df.2)
anova(m1,m2)
m2 <- lmer(accu_b ~ fbblock*task + (1|subj), exp2.df.2)
anova(m1,m2)


m1.bf <- lmBF(accu_b ~ fbblock+task+subj, whichRandom = c('subj'), 
              iterations = 5000, exp2.df.2)
m2.bf <- lmBF(accu_b ~ task+subj, whichRandom = c('subj'), 
              iterations = 5000, exp2.df.2)

m2.bf/m1.bf

m1 <- lmer(stair_b ~ fbblock*task + (1|subj), exp2.df.2)
summary(m1)
emm <- emmeans(m1, ~fbblock|task)
contrast(emm)

m1 <- lmer(stair_b ~ fbblock+task + (1|subj), exp2.df.2)
summary(m1)
confint(m1)
m2 <- lmer(stair_b ~ task + (1|subj), exp2.df.2)
summary(m2)
anova(m1,m2)
m2 <- lmer(stair_b ~ fbblock*task + (1|subj), exp2.df.2)
anova(m1,m2)


m1.bf <- lmBF(stair_b ~ fbblock+task+subj, whichRandom = c('subj'), 
              iterations = 5000, exp2.df.2)
m2.bf <- lmBF(stair_b ~ task+subj, whichRandom = c('subj'), 
              iterations = 5000, exp2.df.2)

m2.bf/m1.bf

############
### Local confidence during test blocks
## Effect of feedback type (fbblock) on local confidence (test)
exp2.df15 <- filter(exp2.mbs, blockType=="interv")

exp2.mbs.reg = lmer(conf_b ~ fbblock*task + accu + stair_b +
                      (1 + fbblock + accu + stair_b|subj), 
                    exp2.df15)
summary(exp2.mbs.reg)


exp2.mbs.reg = lmer(conf_b ~ fbblock*task + accu + stair_b +
                      (1 + fbblock + accu + stair_b|subj), 
                    exp2.df15)
summary(exp2.mbs.reg)
emm <- emmeans(exp2.mbs.reg, ~fbblock|task)
contrast(emm)

m1 = lmer(conf_b ~ fbblock + task + accu + stair_b +
            (1 + fbblock |subj), REML = F, exp2.df15)
summary(m1)
confint(m1)
m1 = lmer(conf_b ~ fbblock + task + accu + stair_b +
            (1 + fbblock + accu + stair_b|subj), REML = F, exp2.df15)
m2 = lmer(conf_b ~ task + accu + stair_b +
            (1 + fbblock + accu + stair_b|subj), REML = F, exp2.df15)
anova(m1, m2)


######## do local conf and accu predict spe
m1 <- lmer(spe_b ~ conf_b + accu_b + (1|task) + (1|group) + 
             (conf_b +1|subj), REML = F, exp2.df)
summary(m1)
confint(m1)
m2 <- lmer(spe_b ~ accu_b + (1|task) + (1|group) + 
             (conf_b +1|subj), REML = F, exp2.df)
summary(m2)
anova(m1,m2)

m2 <- lmer(spe_b ~ conf_b + (1|task) + (1|group) + 
             (conf_b +1|subj), REML = F, exp2.df)
summary(m2)
anova(m1,m2)

############
############
### Local confidence during feedback blocks

## Effect of feedback type (fbblock) on local confidence (test)
first_n_of_trans <- c(1:20)

exp2.df15 <- filter(exp2.mbs, blockType=="interv",
                    trialnum %in% first_n_of_trans)

####
exp2.mbs.reg = lmer(conf_b ~ fbblock*task + accu + stair_b +
                      (1 + fbblock + accu + stair_b|subj), 
                    exp2.df15)
summary(exp2.mbs.reg)


m1 = lmer(conf_b ~ fbblock+task + accu + stair_b +
                      (1 + fbblock + accu + stair_b|subj), 
          REML = F,exp2.df15)
m2 = lmer(conf_b ~ task + accu + stair_b +
            (1 + fbblock + accu + stair_b|subj), 
          REML = F, exp2.df15)

anova(m1, m2)

############
############
############
### Effect of feedback on local conf (test block)
first_n_of_trans <- c(1:20)

exp2.df3 <- filter(exp2.mbs, blockType=="transf",
                   trialnum %in% first_n_of_trans)

m1 <- lmer(conf_b ~ fbblock*task + accu + stair_b + 
             (1 + fbblock + accu + stair_b|subj), REML = F, exp2.df3)
m2 <- update(m1, ~.-fbblock:task)
anova(m1,m2)

emm <- emmeans(m1, ~fbblock|task)
contrast(emm)
confint(contrast(emm))


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
    scale_rt = mean(scale_rt, na.rm = T))

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
  dplyr::summarise(spe = mean(spe),
                   spe_b = mean(spe_b),
    accu = mean(accu),
    conf = mean(conf),
    scale_rt = mean(scale_rt, na.rm = T),
    stair = mean(stair))

exp2.mbs.reg = lmer(spe_b ~ fbblock*task +
                 accu*task + scale_rt*task + stair*task + 
                 (1|subj), REML = F,
               exp2.df9)
plot(exp2.mbs.reg)
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

m2 <- update(exp2.mbs.reg, ~.-fbblock)
anova(exp2.mbs.reg,m2)

exp2.mbs.reg <- update(exp2.mbs.reg, REML=T)
summary(exp2.mbs.reg)

plot_model(exp2.mbs.reg,
           title="Perception task - interaction",
           show.p=TRUE, show.values = T,
           value.offset = .4, value.size = 3.8) +
  coord_cartesian( ylim = (c(-.65,.65)))

emm <- emmeans(exp2.mbs.reg, pairwise~fbblock)
test(emm)


#######################
####

exp2.df4 <- exp2.mbs %>% group_by(subj, firstFeedback) %>%
  summarise(AD = mean(AD))
reg <- lm(AD ~ firstFeedback, exp2.df4)
summary(reg)
plot(reg)
mean(exp2.df4$AD[exp2.df4$firstFeedback=='Positive'])
mean(exp2.df4$AD[exp2.df4$firstFeedback=='Negative'])

hist(exp2.df4$AD[exp2.df4$firstFeedback=='Positive'])
hist(exp2.df4$AD[exp2.df4$firstFeedback=='Negative'])

t.test(exp2.df4$AD[exp2.df4$firstFeedback=='Positive'], 
       exp2.df4$AD[exp2.df4$firstFeedback=='Negative'])

t.test(AD ~ firstFeedback, exp2.df4)

#######################################
#####3##############################
#### Baseline correlation with mhq

exp2.df4 <- exp2.mbs.bas %>% group_by(task, runnum, group, subj, gender) %>%
  summarise(conf = mean(conf),
            accu = mean(accu, na.rm=T),
            stair = mean(stair),
            scale_rt = mean(scale_rt, na.rm=T),
            spe = mean(spe),
            AD = mean(AD),
            Compul = mean(Compul),
            SW = mean(SW),
            age = mean(age))

# baseline values
exp2.spe0 <- exp2.mbs.bas %>% group_by(group, subj)%>%
  summarise(spe0_perc = mean(spe0_perc),
            spe0_mem = mean(spe0_mem)) %>%
  pivot_longer(cols = starts_with("spe0_"),
    names_to = "task",
    names_prefix = "spe0_",
    values_to = "spe0") %>%  
  mutate(task = recode_factor(task, "perc" = taskNames[1], "mem" = taskNames[2]))

exp2.df4 <- exp2.df4 %>% left_join(exp2.spe0) %>% drop_na()

exp2.df4$confZ <- scale(exp2.df4$conf)
exp2.df4$speZ <- scale(exp2.df4$spe)
exp2.df4$spe0Z <- scale(exp2.df4$spe0)


m.reg = lmer(speZ ~ AD + Compul + SW + age + gender + 
               accu + stair + scale_rt + (AD + Compul + SW +1|task),exp2.df4)
summary(m.reg)
m.reg = lmer(speZ ~ AD + Compul + SW + age + gender + 
               accu + stair + scale_rt + (1|task),exp2.df4)
summary(m.reg)
confint(m.reg)
plot(fitted(m.reg), resid(m.reg))
qqPlot(resid(m.reg))

spe.beta <- summary(m.reg)$coefficients[2:4,1]
spe.ci <- confint(m.reg)
spe.ci <- spe.ci[4:6,]

# exp2.df4$spe.AD.partial <- coef(m.reg)[1] + coef(m.reg)['AD']*exp2.df4$AD + as.numeric(resid(m.reg))
exp2.df4$spe.AD.partial <- keepef(m.reg, fix = 'AD')

ggplot(exp2.df4, aes(x=AD, y=spe.AD.partial, color = task))+
  geom_point()+
  geom_smooth(method = lm) +
  stat_cor(method = "pearson", label.x = -.3, label.y = c(1.8, 2.4), 
           size = 5, aes(color = task))  +
  theme_classic()+
  labs(color = 'Task')+
  ggtitle('Exp 2') +
  ylab('Global SPE\nBaseline blocks\n(partialed, z-scored)') +
  xlab('AD axis') +
  theme(text = element_text(size=font_size-2),
        legend.position = 'top')


m2 = lmer(speZ ~ Compul + SW + age + gender + accu + stair + scale_rt + (1|task) ,exp2.df4)
anova(m.reg, m2)

m2 = lmer(speZ ~ AD + SW + age + gender + accu + stair + scale_rt +(1|task) ,exp2.df4)
anova(m.reg, m2)

m2 = lmer(speZ ~ AD + Compul +age + gender + accu + stair + scale_rt +(1|task) ,exp2.df4)
anova(m.reg, m2)


## local confidence
# m.reg = lmer(confZ ~ AD + Compul + SW + age + gender + 
#                accu + stair + scale_rt + (1+AD+Compul+SW|task), 
#              exp2.df4) # outlier on row 68
# summary(m.reg)
m.reg = lmer(confZ ~ AD + Compul + SW + age + gender + 
               accu + stair + scale_rt + (1+AD|task), 
             exp2.df4) # outlier on row 68
summary(m.reg)
exp2.df4$conf.AD.partial <- keepef(m.reg, fix = 'AD')

ggplot(exp2.df4, aes(x=AD, y=conf.AD.partial, color = task))+
  geom_point()+
  geom_smooth(method = lm) +
  stat_cor(method = "pearson", label.x = -.1, label.y = c(3.4, 4.1), 
           size = 5, aes(color = task))  +
  theme_classic()+
  labs(color = 'Task')+
  ggtitle('Exp 2') +
  ylab('Local confidence\nBaseline blocks\n(partialed, z-scored)') +
  xlab('AD axis') +
  theme(text = element_text(size=font_size-2),
        legend.position = 'top')

m.reg = lm(confZ ~ AD + Compul + SW + age + gender + 
               accu + stair + scale_rt, 
             exp2.df4) # outlier on row 68
summary(m.reg)
confint(m.reg)
plot(fitted(m.reg), resid(m.reg))
qqPlot(resid(m.reg))

lconf.beta <- summary(m.reg)$coefficients[2:4,1]
lconf.ci <- confint(m.reg)
lconf.ci <- lconf.ci[2:4,]


library(lmtest)
m2 = lm(confZ ~ Compul + SW + age + gender + accu + stair+scale_rt ,exp2.df4)
lrtest(m.reg, m2)

m2 = lm(confZ ~ AD + SW + age + gender + accu + stair+scale_rt ,exp2.df4)
lrtest(m.reg, m2)

m2 = lm(confZ ~ AD + Compul +age + gender + accu + stair+ scale_rt,exp2.df4)
lrtest(m.reg, m2)


m.reg = lmer(spe0Z ~ AD + Compul + SW + age + gender + 
               accu + stair + scale_rt +
               (1|task),exp2.df4)
summary(m.reg)

exp2.df4$spe0.AD.partial <- keepef(m.reg, fix = 'AD')

ggplot(exp2.df4, aes(x=AD, y=spe0.AD.partial, color = task))+
  geom_point()+
  geom_smooth(method = lm) +
  stat_cor(method = "pearson", label.x = -.1, label.y = c(3.4, 4.1), 
           size = 5, aes(color = task))  +
  theme_classic()+
  labs(color = 'Task')+
  ggtitle('Exp 2') +
  ylab('Global SPE\nProspective\n(partialed, z-scored)') +
  xlab('AD axis') +
  theme(text = element_text(size=font_size-2),
        legend.position = 'top')

pspe.beta <- summary(m.reg)$coefficients[2:4,1]
pspe.ci <- confint(m.reg)
pspe.ci <- pspe.ci[4:6,]

m2 = lmer(spe0Z ~ Compul + SW + 
               age + gender + 
               accu + stair + scale_rt +
               (1|task),exp2.df4)
anova(m.reg,m2)
m2 = lmer(spe0Z ~ AD + SW + 
            age + gender + 
            accu + stair + scale_rt +
            (1|task),exp2.df4)
anova(m.reg,m2)
m2 = lmer(spe0Z ~ Compul + AD + 
            age + gender + 
            accu + stair + scale_rt +
            (1|task),exp2.df4)
anova(m.reg,m2)


########## correlation of symptom dimension with accuracy
m.reg = lmer(scale(accu) ~ AD + Compul + SW + 
               age + gender +
               (1|task), exp2.df4)
summary(m.reg)
accu.beta <- summary(m.reg)$coefficients[2:4,1]
accu.ci <- confint(m.reg)
accu.ci <- accu.ci[4:6,]


m2 = lmer(scale(accu) ~ AD + Compul + age + gender + (1|task), exp2.df4)
anova(m.reg,m2)
######
#####

betas <- c(lconf.beta, spe.beta, pspe.beta, accu.beta)
cilo <- c(lconf.ci[,1], spe.ci[,1], pspe.ci[,1], accu.ci[,1])
cihi <- c(lconf.ci[,2], spe.ci[,2], pspe.ci[,2], accu.ci[,2])
tdims <- c('AD', 'CIT', 'SW')
tdims <- array(tdims, 12)
# ctype <- c(array('Local',3), 
#            array('Global (retrosp.)',3), 
#            array('Global (prosp.)',3))
ctype <- as.factor(c(array('Local\nconf.',3), array('Global SPE\n(retrosp.)',3), 
                     array('Global SPE\n(prosp.)',3), array('Accuracy\n',3)))

exp2.mhqbeta <- data.frame(betas, tdims, ctype, cilo, cihi)
colnames(exp2.mhqbeta) <- c('betas', 'tdims', 'ctype', 'cilo', 'cihi')
exp2.mhqbeta <- exp2.mhqbeta %>%
  mutate(ctype = factor(ctype, levels = c( 'Local\nconf.',
                                           'Global SPE\n(retrosp.)', 'Global SPE\n(prosp.)',
                                           'Accuracy\n'))) %>%
  mutate(tdims = recode_factor(tdims, 'AD' = 'Anxious-\nDepression',
                               'CIT' = 'Compulsivity &\nIntrusive Thought',
                               'SW' = 'Social\nWithdrawal'))

##################################
##### Figure 2B. Transdiagnostic axes and local-global confidence 
f2b <- ggplot(exp2.mhqbeta, aes(x = tdims, y = betas, fill =ctype)) +
  scale_fill_manual(breaks = c( 'Local\nconf.', 'Global SPE\n(retrosp.)', 
                                'Global SPE\n(prosp.)', 'Accuracy\n'), 
                    values= c('orchid', 'deepskyblue' , 'khaki3', 
                              'darkslategray4')) +
  geom_bar(position=position_dodge(.9), stat='identity',
           width=.8) +
  geom_errorbar(aes(ymin=cilo, ymax=cihi),
                width=.2, size=1,                   # Width of the error bars
                position=position_dodge(.9)) +
  theme_pubclean() +
  labs(fill = '') +
  coord_cartesian(ylim = c(-.55,.35)) +
  ylab('Regression slopes') +
  xlab('') +
  theme(text = element_text(size=font_size-6),
        legend.direction = 'horizontal',
        legend.position = c(.42,1.2),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        plot.margin = margin(t = 70, r = 0, b = 0, l = 5, unit = "pt"))+ 
  theme(legend.text=element_text(size=font_size-10)) +
  geom_signif(y_position = c(-.61, -.56, -.55, .3, .3),
              xmin = c(.66, .89, 1.13, 1.66, 3.34),
              xmax = c(.66, .89, 1.13, 1.66, 3.34),
              annotation = c('****', '****', '****', '*', '*'),
              tip_length = 0, color = 'black', textsize = pval_size-2, size = 0)
f2b



######
## scatter plots of mhq and conf correlations
exp2.df2.3 <- exp2.mbs %>% group_by(subj, task, gender) %>%
  summarise(AD = mean(AD),
            conf = mean(conf),
            spe = mean(spe),
            age = mean(age))
exp2.df2.3$conf[exp2.df2.3$task=='Perception'] <- 
  scale(exp2.df2.3$conf[exp2.df2.3$task=='Perception'])
exp2.df2.3$conf[exp2.df2.3$task=='Memory'] <- 
  scale(exp2.df2.3$conf[exp2.df2.3$task=='Memory'])
exp2.df2.3$spe[exp2.df2.3$task=='Perception'] <- 
  scale(exp2.df2.3$spe[exp2.df2.3$task=='Perception'])
exp2.df2.3$spe[exp2.df2.3$task=='Perception'] <- 
  scale(exp2.df2.3$spe[exp2.df2.3$task=='Perception'])
exp2.df2.3 <- exp2.df2.3 %>% group_by(subj,gender) %>%
  summarise(AD = mean(AD),
            conf = mean(conf),
            spe = mean(spe),
            age = mean(age)) %>%
  drop_na()


m.reg = lm(conf ~ AD + age + gender,exp2.df2.3)
summary(m.reg)
exp2.df2.3$conf.AD.partial <- coef(m.reg)[1] + coef(m.reg)[2]*exp2.df2.3$AD + resid(m.reg)

ggplot(exp2.df2.3, aes(x=AD, y=conf.AD.partial))+
  geom_point()+
  geom_smooth(method = lm) +
  stat_cor(method = "pearson", label.x = 0, label.y = 3, size = 5, color='blue')  +
  theme_classic()+
  ggtitle('Exp 2') +
  ylab('Mean baseline\nconfidence\n(partialed, z-scored)') +
  xlab('AD axis') +
  theme(text = element_text(size=font_size-2))


m.reg = lm(spe ~ AD + age + gender,exp2.df2.3)
summary(m.reg)
exp2.df2.3$spe.AD.partial <- coef(m.reg)[1] + coef(m.reg)[2]*exp2.df2.3$AD + resid(m.reg)

ggplot(exp2.df2.3, aes(x=AD, y=spe.AD.partial))+
  geom_point()+
  geom_smooth(method = lm) +
  stat_cor(method = "pearson", label.x = 0, label.y = 1.4, size = 5, color='blue')  +
  theme_classic()+
  ggtitle('Exp 2') +
  ylab('Mean baseline\nSPE\n(partialed, z-scored)') +
  xlab('AD axis') +
  theme(text = element_text(size=font_size-2))


####################################
####################################
###### model-free predictions
#####################################
# update this folder with the folder location with data
data.folder <- "~/OneDrive - University of Copenhagen/Projects/Experiments/metaBiasShift/data/exp"
expNum <- 2
setwd(paste(data.folder, expNum, sep = ''))
setwd('processed/')
mbsDataExp2 = readMat(paste('mbsDataExp',expNum,'.mat',sep='')) # load data
mbsDataExp2 <- mbsDataExp2[[paste('mbsDataExp',expNum,sep='')]]
mbsDataFit <- mbsDataExp2[,,1]$fitRegZ
# spe.est <- mbsDataFit[,,1]$model.10.q1[,,1]$spe.est
spe.est <- mbsDataFit[,,1]$model.6.q1[,,1]$spe.est
nsubj <- dim(spe.est)[1]
nruns <- dim(spe.est)[2]
subj <- array(1:nsubj, c(nsubj,nruns))
runnum <- t(array(1:nruns, c(nruns,nsubj)))
exp2.fit.data <- data.frame(c(spe.est), c(subj), c(runnum))
names(exp2.fit.data) <- c('spe.fit', 'subj', 'runnum')
exp2.fit.data <- exp2.fit.data %>%
  mutate(subj = as.factor(subj)) %>%
  mutate(runnum = as.factor(runnum))

#####################################
### Figure 4a - plot confidence distortion behavioural effect - relationship between local conf and global spe

exp2.model <- exp2.mbs %>% drop_na() %>% 
  dplyr::select(c(subj, runnum, trialnum, spe, confZ, AD)) %>%
  mutate(confhilo = confZ > 0) %>%
  left_join(exp2.fit.data)

exp2.reg.obs.lo <- lmer(spe ~ AD*confZ + (1|subj) + (1|runnum) + (1|trialnum), 
                        filter(exp2.model, confhilo==F))
summary(exp2.reg.obs.lo)
# plot_model(exp2.reg.obs.lo,type = 'pred', terms= c('confZ', 'AD'))

exp2.reg.fit.lo <- lmer(spe.fit ~ AD*confZ + (1|subj) + (1|runnum) + (1|trialnum), 
                        filter(exp2.model, confhilo==F))
summary(exp2.reg.fit.lo)
# plot_model(exp2.reg.fit.lo,type = 'pred', terms= c('confZ', 'AD'))


exp2.reg.fit <- lmer(spe.fit ~ AD*confZ + (1|subj) + (1|runnum) + (1|trialnum), 
                     filter(exp2.model, confhilo==T))
summary(exp2.reg.fit)
# plot_model(exp2.reg.fit,type = 'pred', terms= c('confZ', 'AD'))

exp2.reg.obs <- lmer(spe ~ AD*confZ + (1|subj) + (1|runnum) + (1|trialnum), 
                     filter(exp2.model, confhilo==T))
summary(exp2.reg.obs)
# plot_model(exp2.reg.obs,type = 'pred', terms= c('confZ', 'AD'))


exp2.model <- exp2.mbs %>% drop_na() %>% 
  dplyr::select(c(subj, runnum, trialnum, spe, confZ, AD, accu, fbblock.notrans)) %>%
  mutate(confhilo = confZ > 0) %>%
  dplyr::group_by(subj) %>%
  mutate(conf.tile = ntile(confZ, 4)) %>%
  group_by(runnum) %>%
  mutate(AD.tile = ntile(AD, 3)) %>%
  mutate_at(c('AD.tile', 'runnum'), as.factor) %>%
  left_join(exp2.fit.data)%>%
  rename(spe.data = spe) %>%
  pivot_longer(cols = c('spe.data', 'spe.fit'),
               values_to = 'spe',
               names_to = 'data.fit') %>%
  mutate_at(c('data.fit'), as.factor)

AD.tile.cents  <- c(mean(exp2.model$AD[exp2.model$AD.tile==1]),
                    mean(exp2.model$AD[exp2.model$AD.tile==3]))

conf.tile.cents <- c(mean(exp2.model$confZ[exp2.model$conf.tile==1]),
                     mean(exp2.model$confZ[exp2.model$conf.tile==2]))

emm.obs.lo <- summary(emmeans(exp2.reg.obs.lo, ~ confZ|AD, 
                              at = list(confZ = conf.tile.cents, AD = AD.tile.cents)))
emm.obs.lo

emm.fit.lo <- summary(emmeans(exp2.reg.fit.lo, ~ confZ|AD, 
                              at = list(confZ = conf.tile.cents, AD = AD.tile.cents)))
emm.fit.lo

conf.tile.cents <- c(mean(exp2.model$confZ[exp2.model$conf.tile==3]),
                     mean(exp2.model$confZ[exp2.model$conf.tile==4]))
emm.obs <- summary(emmeans(exp2.reg.obs, ~ confZ|AD, 
                           at = list(confZ = conf.tile.cents, AD = AD.tile.cents)))
emm.obs

emm.fit <- summary(emmeans(exp2.reg.fit, ~ confZ|AD, 
                           at = list(confZ = conf.tile.cents, AD = AD.tile.cents)))
emm.fit


spe <- c(emm.obs.lo$emmean,   emm.fit.lo$emmean,   emm.obs$emmean,   emm.fit$emmean)
lci <- c(emm.obs.lo$lower.CL, emm.fit.lo$lower.CL, emm.obs$lower.CL, emm.fit$lower.CL)
uci <- c(emm.obs.lo$upper.CL, emm.fit.lo$upper.CL, emm.obs$upper.CL, emm.fit$upper.CL)

conf.tile <- c(rep(c(1:2,1:2),2), rep(c(3:4,3:4),2))
AD.tile <- rep(c(rep('Low',2), rep('High',2)),4)
conftile.hilo <- c(rep(1,8), rep(2,8))
data.fit <- rep(c(rep('Data',4), rep('Fit',4)),2)
exp2.model.df <- data.frame(spe, lci, uci, 
                            conf.tile, AD.tile, conftile.hilo, data.fit)
exp2.model.df$conftile.hilo[exp2.model.df$data.fit=='Fit'] <- 0
exp2.model.df$conftile.hilo <- as.factor(exp2.model.df$conftile.hilo)

ggplot(filter(exp2.model.df),
       aes(x = conf.tile, y = spe, colour = AD.tile, shape = data.fit,
           alpha = data.fit, linetype = conftile.hilo)) +
  scale_alpha_manual(values = c(1,.6)) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', 'royalblue2')) +
  scale_shape_manual(values = c(19,2)) +
  geom_point(position=position_dodge(.5), size=3.2, stroke=.5) +
  geom_errorbar(aes(ymin = lci, ymax = uci), linetype = 1,
                position=position_dodge(.5), linewidth=.6, width=.3) +
  stat_summary(geom="line", position=position_dodge(.5)) +
  scale_linetype_manual(values = c('blank', 'dashed', 'solid')) +
  theme_classic2() +
  labs(color = 'AD score', linetype = 'Confidence level')+
  xlab('Confidence (tiled)') +
  ylab('Global SPE') +
  scale_y_continuous(breaks = c(.5,.6)) +
  theme(text = element_text(size=font_size-6),
        legend.direction = 'vertical',
        # axis.title.y = element_text(hjust=.9),
        legend.text=element_text(size=font_size-6),
        plot.caption = element_text(size = font_size-2),
        legend.position = 'none') +
  geom_segment(aes(x = 3.5, y = .52, xend = 3.5, yend = .58), color = 'grey60') +
  geom_segment(aes(x = 3.4, y = .55, xend = 3.5, yend = .55), color = 'grey60') +
  annotate('text', label = '*', x = 3.38, y = .55, size=6, angle = 90)


### confidence distortion behavioural analysis
# 
# 
exp2.model <- exp2.mbs %>% drop_na() %>%
  dplyr::select(c(subj, runnum, trialnum, spe, confZ, AD, accu)) %>%
  mutate(confhilo = confZ > 0) %>%
  group_by(subj,trialnum) %>%
  mutate(spe.tile = ntile(spe, 6))%>%
  mutate(conf.tile = ntile(confZ, 4)) %>%
  group_by(runnum) %>%
  mutate(AD.tile = ntile(AD, 3)) %>%
  mutate_at(c('AD.tile'), as.factor) %>%
  left_join(exp2.fit.data)
levels(exp2.model$AD.tile) <- c(levels(exp2.model$AD.tile), 'High', 'Low')
exp2.model <- exp2.model %>%
  mutate(AD.tile = recode_factor(AD.tile, '1' = 'Low', '3' = 'High'))
exp2.model$conftile.hilo <- as.factor(ceiling(exp2.model$conf.tile/2))

exp2.reg <- lmer(confZ ~ spe*AD + accu + (1|subj) + (1|runnum)+ (1|trialnum), 
                 data = filter(exp2.model, confhilo==T) %>% drop_na())
summary(exp2.reg)
confint(exp2.reg)
# plot_model(exp2.reg,type = 'pred', terms= c('spe', 'AD'))

m2 <- update(exp2.reg, ~.-spe:AD)
anova(exp2.reg,m2)


##### model-free analysis of feedback distortion

exp2.model2 <- filter(exp2.mbs) %>% drop_na() %>%
  group_by(subj, fbblock.notrans, runnum) %>%
  summarise(spe = mean(spe),
            AD = mean(AD)) %>%
  group_by(runnum) %>%
  mutate(AD.tile = ntile(AD, 3)) %>%
  mutate_at(c('AD.tile'), as.factor) %>%
  mutate(fbblock.notrans = 
           factor(fbblock.notrans, levels = c('Negative','None', 'Positive'))) %>%
  mutate(AD.tile = recode_factor(AD.tile, '1' = 'Low', '3' = 'High'))%>%
  left_join(exp2.fit.data) %>%
  rename(spe.data = spe) %>%
  pivot_longer(cols = c('spe.data', 'spe.fit'),
               values_to = 'spe',
               names_to = 'data.fit') %>%
  mutate_at(c('data.fit'), as.factor)


## Figure 4c - feedback distortion behavioural plot

exp2.fb.obs <- lmer(spe ~ fbblock.notrans*AD + (1|subj) + (1|runnum), filter(exp2.model2, data.fit=='spe.data'))
emm.fb.obs <- summary(emmeans(exp2.fb.obs, ~ fbblock.notrans|AD, at = list(AD = AD.tile.cents)))

exp2.fb.fit <- lmer(spe ~ fbblock.notrans*AD + (1|subj) + (1|runnum), filter(exp2.model2, data.fit=='spe.fit'))
emm.fb.fit <- summary(emmeans(exp2.fb.fit, ~ fbblock.notrans|AD, at = list(AD = AD.tile.cents)))

spe <- c(emm.fb.obs$emmean,   emm.fb.fit$emmean)
lci <- c(emm.fb.obs$lower.CL, emm.fb.fit$lower.CL)
uci <- c(emm.fb.obs$upper.CL, emm.fb.fit$upper.CL)

fb.type <- rep(c('Negative', 'None', 'Positive'),4)
AD.tile <- rep(c(rep('Low',3), rep('High',3)),2)
data.fit <- c(rep('Data',6), rep('Fit',6))
exp2.model.df <- data.frame(spe, lci, uci, 
                            fb.type, AD.tile, data.fit)
exp2.model.df$AD.tile <- as.factor(exp2.model.df$AD.tile)
exp2.model.df <- exp2.model.df %>%
  mutate(AD.tile = factor(AD.tile, levels = c('Low','High'))) 

ggplot(filter(exp2.model.df),
       aes(x = fb.type, y = spe, colour = AD.tile, shape = data.fit,
           alpha = data.fit)) +
  scale_alpha_manual(values = c(1,.6)) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', 'royalblue2')) +
  scale_shape_manual(values = c(19,2)) +
  geom_point(position=position_dodge(.5), size=3.2, stroke=.5) +
  geom_errorbar(aes(ymin = lci, ymax = uci), linetype = 1,
                position=position_dodge(.5), linewidth=.6, width=.3) +
  theme_classic2() +
  labs(color = 'AD score', linetype = 'Confidence level')+
  xlab('Feedback type') +
  ylab('Global SPE') +
  scale_y_continuous(breaks = seq(.4,.7,.1)) +
  theme(text = element_text(size=font_size-6),
        legend.direction = 'vertical',
        legend.text=element_text(size=font_size-6),
        plot.caption = element_text(size = font_size-2),
        legend.position = 'none')  +
  geom_segment(aes(x = 2.1, y = .59, xend = 2.9, yend = .59), color = 'grey60') +
  geom_segment(aes(x = 2.1, y = .55, xend = 2.1, yend = .59), color = 'grey60') +
  geom_segment(aes(x = 2.9, y = .59, xend = 2.9, yend = .63), color = 'grey60') +
  annotate('text', label = 'n.s.', x = 2.5, y = .63, size=5) +
  geom_segment(aes(x = 1.1, y = .48, xend = 1.9, yend = .48), color = 'grey60') +
  geom_segment(aes(x = 1.9, y = .48, xend = 1.9, yend = .52), color = 'grey60') +
  geom_segment(aes(x = 1.1, y = .44, xend = 1.1, yend = .48), color = 'grey60') +
  annotate('text', label = 'n.s.', x = 1.5, y = .51, size=5)


## feedback distortion behavioural analysis

exp2.model2 <- filter(exp2.mbs) %>% drop_na() %>%
  group_by(subj, fbblock.notrans, runnum) %>%
  summarise(spe = mean(spe),
            AD = mean(AD)) %>%
  group_by(runnum) %>%
  mutate(AD.tile = ntile(AD, 3)) %>%
  mutate_at(c('AD.tile'), as.factor) %>%
  mutate(fbblock.notrans = 
           factor(fbblock.notrans, levels = c('Negative','None', 'Positive'))) %>%
  mutate(AD.tile = recode_factor(AD.tile, '1' = 'Low', '3' = 'High'))%>%
  left_join(exp2.fit.data)

exp2.model2$speZ <- c(scale(exp2.model2$spe))


exp2.reg <- lmer(speZ ~ fbblock.notrans*AD + (1+ fbblock.notrans|subj) + (1|runnum), 
                 filter(exp2.model2, fbblock.notrans %in% c('None', 'Positive')))
summary(exp2.reg)
confint(exp2.reg)
plot(exp2.reg)

m2 <- update(exp2.reg, ~.-fbblock.notrans:AD)
anova(exp2.reg, m2)

set.seed(2020)

#### perform feedback analysis 
exp2.pos.AD.m1.bf <- lmBF(speZ ~ fbblock.notrans*AD + subj + runnum + fbblock.notrans:subj, 
                           whichRandom = c('subj', 'runnum'),
                          iterations = 10000,
                           filter(exp2.model2, fbblock.notrans %in% c('None', 'Positive')))
exp2.pos.AD.m2.bf <- lmBF(speZ ~ fbblock.notrans + AD + subj + runnum + fbblock.notrans:subj, 
                           whichRandom = c('subj', 'runnum'),
                          iterations = 10000,
                           filter(exp2.model2, fbblock.notrans %in% c('None', 'Positive')))

exp2.pos.AD.m2.bf/exp2.pos.AD.m1.bf

# exp2.pos.AD.m1.bf/exp2.pos.AD.m2.bf


exp2.reg <- lmer(speZ ~ fbblock.notrans*AD + 
                   (1+fbblock.notrans|subj) + (1|runnum), 
                 filter(exp2.model2, fbblock.notrans %in% c('None', 'Negative')))
summary(exp2.reg)
confint(exp2.reg)
plot(exp2.reg)

m2 <- update(exp2.reg, ~.-fbblock.notrans:AD)
anova(exp2.reg, m2)

exp2.neg.AD.m1.bf <- lmBF(speZ ~ fbblock.notrans*AD + subj + runnum + fbblock.notrans:subj, 
                           whichRandom = c('subj', 'runnum', 'fbblock.notras:subj', 'task'),
                          iterations = 10000,
                           filter(exp2.model2, fbblock.notrans %in% c('None', 'Negative')))
exp2.neg.AD.m2.bf <- lmBF(speZ ~ fbblock.notrans + AD + subj + runnum + fbblock.notrans:subj, 
                           whichRandom = c('subj', 'runnum', 'fbblock.notras:subj', 'task'),
                          iterations = 10000,
                           filter(exp2.model2, fbblock.notrans %in% c('None', 'Negative')))
exp2.neg.AD.m2.bf/exp2.neg.AD.m1.bf

exp2.neg.AD.m1.bf/exp2.neg.AD.m2.bf


####### model free analysis along CIT dimension

exp2.model <- exp2.mbs %>% drop_na() %>% 
  dplyr::select(c(subj, runnum, trialnum, spe, confZ, Compul, AD, accu)) %>%
  mutate(confhilo = confZ > 0) %>%
  group_by(subj,trialnum) %>%
  mutate(conf.tile = ntile(confZ, 4)) %>%
  group_by(runnum) %>%
  mutate(CIT.tile = ntile(Compul, 3)) %>%
  mutate_at(c('CIT.tile'), as.factor)
levels(exp2.model$CIT.tile) <- c(levels(exp2.model$CIT.tile), 'High', 'Low')
exp2.model <- exp2.model %>%
  mutate(CIT.tile = recode_factor(CIT.tile, '1' = 'Low', '3' = 'High'))
exp2.model$conftile.hilo <- as.factor(ceiling(exp2.model$conf.tile/2))

ggplot(filter(exp2.model, CIT.tile!=2),
       aes(x = conf.tile, y = spe, colour = CIT.tile)) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', 'royalblue2')) +
  geom_smooth(method='lm', aes(linetype = conftile.hilo), se = F) +
  scale_linetype_manual(values = c('dotted', 'solid')) +
  stat_summary(fun.data = mean_cl_boot, size=.5, shape = 1) +
  stat_summary(fun.y = mean, size=.5, alpha = .5) +
  labs(color = 'CIT score', linetype = 'Confidence level')+
  theme_classic2() +
  xlab('Confidence (tiled)') +
  ylab('Global SPE') +
  theme(text = element_text(size=font_size-2),
        legend.direction = 'vertical',
        legend.text=element_text(size=font_size-6),
        plot.caption = element_text(size = font_size-2),
        legend.position = 'none')

exp2.reg <- lmer(confZ ~ spe*Compul + accu + (1|subj) + (1|runnum)+ (1|trialnum), 
                 data = filter(exp2.model, confhilo==T) %>% drop_na())
summary(exp2.reg)
m2 <- update(exp2.reg, ~.-spe:Compul)
anova(exp2.reg,m2)


#####Feedback

exp2.temp2 <- filter(exp2.temp) %>% drop_na() %>%
  group_by(subj, fbblock.notrans, task, runnum, AD.tile) %>%
  summarise(spe = mean(spe),
            AD = mean(AD)) %>%
  mutate(fbblock.notrans = factor(fbblock.notrans, levels = c('None', 'Negative','Positive')))
exp2.temp2$speZ <- scale(exp2.temp2$spe)

exp2.reg <- lmer(speZ ~ fbblock.notrans*AD + 
                   (1+fbblock.notrans|subj) + (1|runnum), 
                 filter(exp2.temp2, fbblock.notrans %in% c('None', 'Positive')))
summary(exp2.reg)
plot(exp2.reg)

m2 <- update(exp2.reg, ~.-fbblock.notrans:AD)
anova(exp2.reg, m2)

# BF_BIC = exp((BIC(exp2.reg) - BIC(m2))/2)  # BICs to Bayes factor
set.seed(2020)

exp2.pos.m1.bf <- lmBF(speZ ~ fbblock.notrans*AD + subj + fbblock.notrans:subj+ runnum + task, 
                 whichRandom = c('subj', 'runnum', 'task'),
                 iterations = 50000,
                 filter(exp2.temp2, fbblock.notrans %in% c('None', 'Positive')))
exp2.pos.m2.bf <- lmBF(speZ ~ fbblock.notrans + AD + subj + fbblock.notrans:subj+ runnum + task, 
                 whichRandom = c('subj', 'runnum', 'task'),
                 iterations = 50000,
                 filter(exp2.temp2, fbblock.notrans %in% c('None', 'Positive')))

exp2.pos.m1.bf/exp2.pos.m2.bf


exp2.reg <- lmer(speZ ~ fbblock.notrans*AD + 
                   (1+fbblock.notrans|subj) + (1|runnum), 
                 filter(exp2.temp2, fbblock.notrans %in% c('None', 'Negative')))
summary(exp2.reg)
plot(exp2.reg)

m2 <- update(exp2.reg, ~.-fbblock.notrans:AD)
anova(exp2.reg, m2)

set.seed(2020)

exp2.neg.m1.bf <- lmBF(speZ ~ fbblock.notrans*AD + subj + fbblock.notrans:subj + runnum + task, 
              whichRandom = c('subj', 'runnum', 'task'),
              iterations = 50000,
              filter(exp2.temp2, fbblock.notrans %in% c('None', 'Negative')))
exp2.neg.m2.bf <- lmBF(speZ ~ fbblock.notrans + AD + subj + fbblock.notrans:subj + runnum + task, 
              whichRandom = c('subj', 'runnum', 'task'),
              iterations = 50000,
              filter(exp2.temp2, fbblock.notrans %in% c('None', 'Negative')))
exp2.neg.m1.bf/exp2.neg.m2.bf


exp2.temp2 <- exp2.temp2 %>%
  mutate(fbblock.notrans = factor(fbblock.notrans, 
                                  levels = c('Negative','None', 'Positive')))

ggplot(filter(exp2.temp2, AD.tile!=2), 
       aes(x = fbblock.notrans, y = scale(spe), colour = AD.tile)) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', '#0066cc')) +
  stat_summary(fun.y = mean, size=.8, alpha = .3, 
               position = position_dodge(.5)) + 
  stat_summary(fun.data = mean_cl_boot, size=.8, shape = 1, 
               position = position_dodge(.5)) + 
  xlab('Feedback type') +
  ylab('Global SPE (z-scored)') +
  theme_classic2() +
  theme(text = element_text(size=font_size-2), 
        plot.caption = element_text(size = font_size-2),
        axis.title.y = element_text(hjust=.9),
        legend.position = 'none') +
  geom_segment(aes(x = 2.1, y = .35, xend = 2.9, yend = .35), color = 'grey60') +
  geom_segment(aes(x = 2.1, y = .2, xend = 2.1, yend = .35), color = 'grey60') +
  geom_segment(aes(x = 2.9, y = .35, xend = 2.9, yend = .5), color = 'grey60') +
  annotate('text', label = 'n.s.', x = 2.5, y = .47, size=5) +
  geom_segment(aes(x = 1.1, y = -.35, xend = 1.9, yend = -.35), color = 'grey60') +
  geom_segment(aes(x = 1.9, y = -.2, xend = 1.9, yend = -.35), color = 'grey60') +
  geom_segment(aes(x = 1.1, y = -.35, xend = 1.1, yend = -.5), color = 'grey60') +
  annotate('text', label = 'n.s.', x = 1.5, y = -.23, size=5)




##
####################
####################
#### For sret analysis read in the file with where participants are not excluded in the narrow performance range

library(scales)
setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/metaBiasShift/")

## read data file
exp2.mbs3 = read.csv(paste('data/exp2/processed/mbsExp2_noperfexcl.csv', sep=''),
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
psych = read.csv(paste('data/exp2/processed/factor_scores_noperfexc.csv', sep=''),
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
m2 <- update(m.reg, ~.-wordvalence:AD)
lrtest(m.reg,m2)
m2 <- update(m.reg, ~.-wordvalence:Compul)
lrtest(m.reg,m2)
m2 <- update(m.reg, ~.-wordvalence:SW)
lrtest(m.reg,m2)

##################################
##### Figure 5A. SRET Figure

f5a <- plot_model(m.reg, type = ("pred"),
           terms = c("AD", "wordvalence"), 
           line.size = 1.5, ci.lvl = .95,
           title = '', colors = posnegWordColours.exp2) +
  ylab("Self-endorsement") +
  xlab("Anxious-Depression axis") +
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

m.reg <- glmer(sretVal/4 ~ spe.bas*wordvalence + conf.bas*wordvalence + (1|wordnum) , 
               family = binomial(link = 'logit'), control = glmerControl(optimizer = c('bobyqa')),
               filter(exp2.df.2, sretTimepoint=='1'))
summary(m.reg)
confint(m.reg)
m2 <- update(m.reg, ~.-wordvalence:spe.bas)
lrtest(m.reg,m2)
m2 <- update(m.reg, ~.-wordvalence:conf.bas)
lrtest(m.reg,m2)
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
confint(m.reg)
m2 <- update(m.reg, ~.-fbblock:wordvalence)
anova(m.reg, m2)

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
confint(contrast(emm))
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

m.reg <- lmer(sretDiff ~ fbblock*task + (1|sretset) , exp2.df.5)
summary(m.reg)

m.reg <- lmer(sretDiff ~ fbblock + (1|sretset) , exp2.df.5)
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



############
####### 
############
### A - does meta-sensitivity moderate the effect of feedback
# update this folder with the folder location with data
data.folder <- "/Users/skatyal/OneDrive - University College London/Projects/Experiments/metaBiasShift/data/exp2/processed/"

expNum <- 2
setwd(data.folder)
mbsDataExp2 = readMat(paste('mbsDataExp',expNum,'.mat',sep='')) # load questionnaire data
mbsDataExp2 <- mbsDataExp2[[paste('mbsDataExp',expNum,sep='')]]

mratio <- mbsDataExp2[,,1]$mratio
nsubj <- length(unique(exp2.df$subj))
subj <- array(1:nsubj, c(nsubj, 2))
task <- t(array(c("Perception", "Memory"), c(2, nsubj)))

temp <- data.frame(c(mratio), c(subj), c(task))
colnames(temp) <- c('mratio', 'subj', 'task')
temp <- temp %>% mutate_at('subj', as.factor)

# exp2.df.2 <- exp2.df.2 %>% left_join(temp)
exp2.df <- exp2.df %>% left_join(temp)

exp2.df$fbblock[exp2.df$runnum==4|exp2.df$runnum==6] <- 'None'

m.reg <- lmer(spe_b ~ fbblock*task*mratio + accu_b + stair_b + (1|subj), 
              filter(exp2.df, mratio>0, runnum %in% c(3:6)))
summary(m.reg)
plot(m.reg)
m2 <- update(m.reg, ~.-fbblock:task:mratio)
anova(m.reg,m2)
m.reg <- m2 # accept the reduced model
summary(m.reg)
m2 <- update(m.reg, ~.-fbblock:mratio)
anova(m.reg,m2)



m.reg <- lmer(mratio ~ AD + Compul + SW + age + gender + (1|task) , 
              filter(exp2.df, runnum %in% c(1:2)))
summary(m.reg)
m2 <- update(m.reg, ~.-AD)
anova(m.reg,m2)
m2 <- update(m.reg, ~.-Compul)
anova(m.reg,m2)
m2 <- update(m.reg, ~.-SW)
anova(m.reg,m2)



############
### A8 - plots of debriefing questions
exp2.db <- exp2.df.2 %>% group_by(subj, awarepos, awareneg) %>%
  summarise(AD = mean(AD)) %>%
  pivot_longer(    cols = starts_with("aware"),
                   names_to = "posneg",
                   names_prefix = "aware",
                   values_to = "aware.fb") %>%
  mutate(posneg = recode_factor(posneg, 'pos' = 'Correct', 'neg' = 'Incorrect'))


sum((exp2.df.2 %>% group_by(subj, awareneg))$awareneg=='Yes', na.rm=T)/
  (sum((exp2.df.2 %>% group_by(subj, awareneg))$awareneg=='No', na.rm=T)+
     sum((exp2.df.2 %>% group_by(subj, awareneg))$awareneg=='Yes', na.rm=T))
sum((exp2.df.2 %>% group_by(subj, awarepos))$awarepos=='Yes', na.rm=T)/
  (sum((exp2.df.2 %>% group_by(subj, awarepos))$awarepos=='No', na.rm=T)+
     sum((exp2.df.2 %>% group_by(subj, awarepos))$awarepos=='Yes', na.rm=T))


ggplot(exp2.db, aes(x = aware.fb, fill = posneg)) +
  geom_bar(position = position_dodge()) +
  xlab('Did you notice feedback bias during block?') +
  ylab('Number of participants') +
  theme_classic() +
  ggtitle('Exp 2') +
  labs(fill = 'Feedback block') +
  theme(legend.position = 'top',
        text = element_text(size=16))

t.test(filter(exp2.db, !is.na(aware.fb), posneg=='Correct', aware.fb=='Yes')$AD,
       filter(exp2.db, !is.na(aware.fb), posneg=='Correct', aware.fb=='No')$AD)
t.test(filter(exp2.db, !is.na(aware.fb), posneg=='Incorrect', aware.fb=='Yes')$AD,
       filter(exp2.db, !is.na(aware.fb), posneg=='Incorrect', aware.fb=='No')$AD)


ggplot(filter(exp2.db, !is.na(aware.fb)), aes(x = aware.fb, y = AD, colour = posneg))+
  stat_summary()


exp2.db <- exp2.df %>% group_by(subj, affectpos, affectneg) %>%
  summarise(AD = mean(AD)) %>%
  pivot_longer(    cols = starts_with("affect"),
                   names_to = "posneg",
                   names_prefix = "affect",
                   values_to = "affect.fb") %>%
  mutate(posneg = recode_factor(posneg, 'pos' = 'Correct', 'neg' = 'Incorrect'))

ggplot(exp2.db, aes(x = affect.fb, fill = posneg)) +
  geom_bar(position = position_dodge()) +
  xlab('Did feedback on trial change how you felt?') +
  ylab('Number of participants') +
  theme_classic() +
  ggtitle('Exp 2') +
  labs(fill = 'Feedback on trial')+
  theme(legend.position = 'top',
        text = element_text(size=16))

levels(exp2.db$affect.fb) <- c(levels(exp2.db$affect.fb), 
                               'Not felt better', "Not felt worse")
exp2.db <- exp2.db %>% 
  mutate(affect.better = affect.fb) %>%
  mutate_at('affect.better', ~replace(., .!="Felt better", "Not felt better")) %>% 
  mutate(affect.worse = affect.fb) %>%
  mutate_at('affect.worse', ~replace(., .!="Felt worse", "Not felt worse"))

length(filter(exp2.db, affect.better=="Felt better", posneg == 'Correct')$AD)/
  (length(filter(exp2.db, affect.better=="Felt better", posneg == 'Correct')$AD)+
     length(filter(exp2.db, affect.better=="Not felt better", posneg == 'Correct')$AD))
length(filter(exp2.db, affect.worse=="Felt worse", posneg == 'Incorrect')$AD)/
  (length(filter(exp2.db, affect.worse=="Felt worse", posneg == 'Incorrect')$AD)+
     length(filter(exp2.db, affect.worse=="Not felt worse", posneg == 'Incorrect')$AD))


wilcox.test(filter(exp2.db, affect.better=="Felt better", posneg == 'Correct')$AD,
       filter(exp2.db, affect.better=="Not felt better", posneg == 'Correct')$AD)
wilcox.test(filter(exp2.db, affect.worse=="Felt worse", posneg == 'Incorrect')$AD,
       filter(exp2.db, affect.worse=="Not felt worse", posneg == 'Incorrect')$AD)

ggplot(filter(exp2.db, posneg=='Incorrect'), aes(x = affect.worse, y = AD)) +
  geom_quasirandom(dodge.width=.7, alpha = .7) +
  geom_violin(alpha=.6, position = position_dodge(.7)) +
  stat_summary(fun.data = mean_cl_boot, color = 'black', size = 1, width=.2,
               geom = 'errorbar') +
  xlab('Feeling after incorrect trials') +
  ylab('AD axis score') +
  theme_classic() +
  ggtitle('Exp 2') +
  # coord_cartesian(ylim = c(-.45,.35)) +
  theme(legend.position = 'top',
        text = element_text(size=16)) +
  geom_signif(y_position = c(.4), vjust = -.25, hjust=.5,
              xmin = c(1), color = 'black',
              xmax = c(2),
              annotation = c("p = .0007"), tip_length = .02,
              textsize = pval_size, size = .8)

# 
# levels(exp2.db$affect.fb) <- c(levels(exp2.db$affect.fb), 
#                                'Not felt better', "Not felt worse")
# 
# exp2.db2 <- filter(exp2.db, !is.na(affect.fb), posneg=='Correct') %>%
#   mutate_at('affect.fb', ~replace(., .!="Felt better", "Not felt better"))
# m.reg <- lm(AD ~ affect.fb, exp2.db2)
# summary(m.reg)
# 
# exp2.db2 <- filter(exp2.db, !is.na(affect.fb), posneg=='Incorrect') %>%
#   mutate_at('affect.fb', ~replace(., .!="Felt worse", "Not felt worse"))
# m.reg <- lm(AD ~ affect.fb, exp2.db2)
# summary(m.reg)
# t.test(filter(exp2.db2, affect.fb == 'Felt worse')$AD,
#        filter(exp2.db2, affect.fb == 'Not felt worse')$AD)
# ggplot(exp2.db2, aes(x = affect.fb, y = AD)) +
#   stat_summary(fun.data = mean_cl_boot)

#### Figures for rebuttal
hist(exp2.mbs$conf,
     main = "A.  Local confidence", 
     xlab = "Local confidence")
hist(exp2.mbs$conf_b,
     main = "B.   Local confidence", 
     xlab = "Local confidence (baseline subtracted)")
hist(scale(exp2.df4$conf),
     main = "C.   Mean local confidence", 
     xlab = "Local confidence (z-scored)")
hist(scale(exp2.df4$spe),
     main = "D.   Global SPE", 
     xlab = "Global SPE (z-scored)")



## test if colour of stimulus makes any difference to confidence
exp2.mbs = read.csv(paste('data/exp2/processed/mbsExp2_stim.csv', sep=''), header = TRUE, sep = ',')
exp2.mbs <- exp2.mbs %>%
  mutate_at(c("subj", 'stim'), as.factor)
summary(lmer(conf ~ stim + (1+stim|subj) + (1+stim|runnum) + (1+stim|trialnum),
             filter(exp2.mbs, task==0)))
psych = read.csv(paste('data/exp2/processed/factor_scores.csv', sep=''),
                 header = TRUE, sep = ',')
psych <- psych %>% rename(subj = subjIDs) %>% dplyr::select(-X) %>%
  mutate_at("subj", as.factor)

exp2.mbs <- exp2.mbs %>% left_join(psych)

# summary(lmer(conf ~ AD*stim + (1+AD*stim|subj) + (1+AD*stim|runnum),
#              filter(exp2.mbs, task==0)))

summary(lmer(conf ~ AD*stim + (1+stim|subj) + (1|runnum),
             filter(exp2.mbs, task==0)))

