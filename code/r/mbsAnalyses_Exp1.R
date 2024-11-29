
####################################################################
################# Analysis and Figure generation for Exp 1 #########
####################################################################
################# Katyal, Huys, Dolan, Fleming #####################
### How underconfidence is maintained in anxiety and depression ####
####################################################################
####  This file reproduces all analyses and figures for Exp 1 ######
####################################################################

library(R.matlab)
library(lmerTest) #Linear mixed effect model
library(emmeans) # Least squares means
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm) # for geom_quasirandom()
library(ggpp)
library(sjPlot)
library(see)
library(remef)
library(BayesFactor)
library(lmtest)
# library(mixedpower)

## mediation analysis on mixed models does not work with lmerTest package
## so use only one or the other 
doMediation = F
if (doMediation){
  library(mediation) # mediation analysis
} else {
  library(lmerTest) # anova and summary with p-values
  
}


## setting initial variables and options
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


## set the base directory here which contains a the 'analysis' and 'data' folders
setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/metaBiasShift/")

## read data file and preparing the data
exp1.mbs = read.csv(paste('data/exp1/processed/mbsExp1.csv', sep=''), header = TRUE, sep = ',')
## recode columns as factors
exp1.mbs$subj <- as.factor(exp1.mbs$subj)
exp1.mbs$task <- as.factor(exp1.mbs$task)
exp1.mbs$awarepos[is.nan(exp1.mbs$awarepos)] <- 'NA'
exp1.mbs$awareneg[is.nan(exp1.mbs$awareneg)] <- 'NA'
exp1.mbs$awarepos <- as.factor(exp1.mbs$awarepos)
exp1.mbs$awareneg <- as.factor(exp1.mbs$awareneg)
exp1.mbs$affectpos[is.nan(exp1.mbs$affectpos)] <- 'NA'
exp1.mbs$affectneg[is.nan(exp1.mbs$affectneg)] <- 'NA'
exp1.mbs$affectpos <- as.factor(exp1.mbs$affectpos)
exp1.mbs$affectneg <- as.factor(exp1.mbs$affectneg)
exp1.mbs$gender[is.nan(exp1.mbs$gender)] <- NA
exp1.mbs$gender <- as.factor(exp1.mbs$gender)

exp1.mbs <- exp1.mbs %>%
  mutate_at(vars(contains("fb")), as.factor) %>%
  rename(accu = corr) %>%
  rename(stair = incdec)

exp1.mbs <- exp1.mbs %>%
  mutate(fbblock = recode_factor(fbblock, "1" = posnegNames[1], 
                                 "2" = posnegNames[2], "0" = "None"),
         task = recode_factor(task, "0" = taskNames[1], "1" = taskNames[2]),
         awarepos = recode_factor(awarepos, '1' = 'Yes', '2'='No', 'NA' = 'No response'),
         awareneg = recode_factor(awareneg, '1' = 'Yes', '2'='No', 'NA' = 'No response'),
         affectpos = recode_factor(affectpos, '1' = 'Felt better', '2'='Felt worse', 
                                   '3'='No change', '4'='Not sure',
                                   'NA' = 'No response'),
         affectneg = recode_factor(affectneg, '1' = 'Felt better', '2'='Felt worse', 
                                   '3'='No change', '4'='Not sure',
                                   'NA' = 'No response'))

# ## check when in block was feedback delivered
# exp1.fb <- filter(exp1.mbs, runnum %in% c(3,5)) %>% 
#   group_by(fbblock, trialnum) %>%
#   summarise(feedback = mean(feedback)) %>%
#   pivot_wider(names_from = fbblock,
#               values_from = feedback) %>%
#   mutate(diff = Positive - Negative)%>% 
#   mutate(expnum = as.factor(1))


exp1.mbs <- exp1.mbs %>%
  mutate(endorsediff = endorsepos + (20-endorseneg))

## code highly deviant RTs as NA
max_RT_deviation = 3
rt1_max = c(median(exp1.mbs$rt1[exp1.mbs$task==taskNames[1]], na.rm=T)+ 
              max_RT_deviation*IQR(exp1.mbs$rt1[exp1.mbs$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(exp1.mbs$rt1[exp1.mbs$task==taskNames[2]], na.rm=T)+ 
              max_RT_deviation*IQR(exp1.mbs$rt1[exp1.mbs$task==taskNames[2]], na.rm=T ))
exp1.mbs <- filter(exp1.mbs, trialnum > nLags) %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[1] & .>rt1_max[1],NA)) %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[2] & .>rt1_max[2],NA))

## scale RTs 
exp1.mbs$scale_rt <- exp1.mbs$rt1
exp1.mbs$scale_rt[exp1.mbs$task==taskNames[1]] <- scale(exp1.mbs$rt1[exp1.mbs$task==taskNames[1]])
exp1.mbs$scale_rt[exp1.mbs$task==taskNames[2]] <- scale(exp1.mbs$rt1[exp1.mbs$task==taskNames[2]])

rt2_max = c(median(exp1.mbs$rt2[exp1.mbs$task==taskNames[1]], na.rm=T)+ 
              max_RT_deviation*IQR(exp1.mbs$rt2[exp1.mbs$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(exp1.mbs$rt2[exp1.mbs$task==taskNames[2]], na.rm=T)+ 
              max_RT_deviation*IQR(exp1.mbs$rt2[exp1.mbs$task==taskNames[2]], na.rm=T ))
exp1.mbs <- filter(exp1.mbs, trialnum > nLags) %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[1] & .>rt2_max[1],NA)) %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[2] & .>rt2_max[2],NA))

exp1.mbs$scale_rt2 <- exp1.mbs$rt2
exp1.mbs$scale_rt2[exp1.mbs$task==taskNames[1]] <- scale(exp1.mbs$rt2[exp1.mbs$task==taskNames[1]])
exp1.mbs$scale_rt2[exp1.mbs$task==taskNames[2]] <- scale(exp1.mbs$rt2[exp1.mbs$task==taskNames[2]])

exp1.mbs$stair[exp1.mbs$task==taskNames[1]] <- scale(exp1.mbs$stair[exp1.mbs$task==taskNames[1]])
exp1.mbs$stair[exp1.mbs$task==taskNames[2]] <- scale(exp1.mbs$stair[exp1.mbs$task==taskNames[2]])

# exp1.mbs1 <- exp1.mbs

## to subtract the baseline blocks confidence and spe from the rest of the blocks
if (subtractBaselineConfSpe){
  # get the baseline values for the two tasks from runs 1 and 2
  exp1.mbs2 <- filter(exp1.mbs, runnum<3) %>%
    group_by(subj, task, group) %>%
    dplyr::summarise(
      conf.bas = mean(conf),
      spe.bas = mean(spe),
      accu.bas = mean(accu),
      stair.bas = mean(stair))
  
  # save the baseline data separately
  exp1.mbs.bas <- filter(exp1.mbs, runnum<3) %>%
    group_by(subj, task, runnum, group, trialnum, gender) %>%
    dplyr::summarise(
      conf = mean(conf),
      accu = mean(accu),
      stair = mean(stair),
      scale_rt = mean(scale_rt),
      spe = mean(spe),
      phq = mean(phq),
      gad = mean (gad),
      spin = mean(spin),
      age = mean(age))
  
  # left join the baseline values an subtract them from confidence and spe
  exp1.mbs <- exp1.mbs %>% left_join(exp1.mbs2) %>%
    # mutate(spe = spe) %>% #save untransfromed spe to plot
    mutate(conf_b = conf-conf.bas)%>%
    mutate(spe_b = spe-spe.bas) %>%
    # mutate(accu_b = accu-accu.bas) %>%
    mutate(stair_b = stair - stair.bas)
}

## take away the baseline (first two) blocks from the dataset
# exp1.mbs <- exp1.mbs %>%
#   filter(runnum>2)

exp1.mbs <- exp1.mbs %>%
  mutate(blockType = factor(runnum%%2)) %>%
  mutate(blockType = recode_factor(blockType, "0" = "transf", "1" = "interv"))
levels(exp1.mbs$blockType) <- c(levels(exp1.mbs$blockType), 'base')
exp1.mbs$blockType[exp1.mbs$runnum %in% c(1,2)] <- 'base'

exp1.mbs <- exp1.mbs %>%
  mutate_at(vars(contains("runnum")), as.factor)
exp1.mbs$fbblock <- droplevels(exp1.mbs$fbblock)

## create new columns for trial-by-trial positive and negative feedback
exp1.mbs <- exp1.mbs %>%
  mutate(fbpos=feedback, fbneg=feedback) %>%
  mutate_at("fbpos", ~replace(., feedback==1 & accu==1, 1)) %>%
  mutate_at("fbpos", ~replace(., feedback==1 & accu==0, 0)) %>%
  mutate_at("fbneg", ~replace(., feedback==1 & accu==1, 0)) %>%
  mutate_at("fbneg", ~replace(., feedback==1 & accu==0, 1)) %>%
  mutate_at('feedback', ~replace(., feedback==1 & accu==0, -1)) %>%
  dplyr::select(-c(rt2)) %>%
  mutate(fbblock.notrans = fbblock) %>%
  mutate_at('fbblock.notrans', ~replace(., runnum==4 | runnum==6, 'None'))

## code intervention task by group
exp1.mbs$intervTask = exp1.mbs$group
exp1.mbs$intervTask[exp1.mbs$group %in% c(1,2,3,4)] = taskNames[1]
exp1.mbs$intervTask[exp1.mbs$group %in% c(5,6,7,8)] = taskNames[2]
# 
# transfer task by group (same as task column)
exp1.mbs$testTask = exp1.mbs$group
exp1.mbs$testTask[exp1.mbs$group %in% c(1,2,7,8)] = taskNames[1]
exp1.mbs$testTask[exp1.mbs$group %in% c(3,4,5,6)] = taskNames[2]

## code whether transfer is to same or opposite task
exp1.mbs$transferType = exp1.mbs$group
exp1.mbs$transferType[exp1.mbs$group %in% c(1,2,5,6)] = transfertaskNames[1]
exp1.mbs$transferType[exp1.mbs$group %in% c(3,4,7,8)] = transfertaskNames[2]

## code whether pos or neg is the first block of feedback
exp1.mbs$firstFeedback = exp1.mbs$group
exp1.mbs$firstFeedback[exp1.mbs$group %in% c(1,3,5,7)] = posnegNames[1]
exp1.mbs$firstFeedback[exp1.mbs$group %in% c(2,4,6,8)] = posnegNames[2]

exp1.mbs <- exp1.mbs%>%
  mutate_at(c("group", "transferType", "intervTask", 'testTask'), as.factor) %>%
  mutate(group = recode_factor(group, "1" = 'Group 1', 
                               "2" = 'Group 2', "3" = 'Group 3', "4" = 'Group 4', 
                               "5" = 'Group 5', "6" = 'Group 6', "7" = 'Group 7', 
                               "8" = 'Group 8')) 

# exp1.mbs3 <- exp1.mbs

exp1.mbs <- exp1.mbs %>% drop_na(scale_rt, scale_rt2)

exp1.mbs$scale_rt[as.logical(exp1.mbs$invalid_rt1)] = NA # for trials where participants changed their response, mark rt1 as NA

#####################

exp1.df <- exp1.mbs %>%
  group_by(subj, task, runnum, fbblock, group, 
           blockType, intervTask, transferType, 
           awareneg, awarepos, affectneg, affectpos,
           firstFeedback, testTask, gender) %>%
  dplyr::summarise(
    conf = mean(conf),
    conf_b = mean(conf_b),
    accu = mean(accu),
    accu_b = mean(accu, na.rm=T) - mean(accu.bas),
    scale_rt = mean(scale_rt, na.rm=T),
    rt1 = mean(rt1, na.rm=T),
    fbneg = mean(fbneg),
    fbpos = mean(fbpos),
    spe = mean(spe),
    spe_b = mean(spe_b),
    stair = mean(stair),
    stair_b = mean(stair_b),
    phq = mean(phq),
    gad = mean(gad),
    spin = mean(spin),
    # conf.bas = mean(conf.bas),
    # spe.bas = mean(spe.bas),
    age = mean(age))
  
## check when in block was feedback delivered
exp1.fb <- filter(exp1.mbs, runnum %in% c(3,5)) %>% 
  group_by(subj, fbblock, trialnum) %>%
  summarise(feedback = mean(feedback))
ggplot(exp1.fb, aes(x=trialnum, y=feedback, color=fbblock)) +
  stat_summary()

##### data prep done #####################
##########################################


#################  #################
####### Figure 2
###### effect of feedback on intervention SPEs, accuracy, staircase

exp1.df.2 <- filter(exp1.df, blockType=="interv")

f2a.left <- ggplot(exp1.df.2, aes(x = fbblock, y = spe_b, color = fbblock)) +
  facet_wrap(~task) +
  geom_hline(yintercept = 0, color = 'grey30') +
  scale_color_manual(breaks = posnegNames, values=posnegColours) +
  theme_blank()  +
  # scale_x_discrete(limits = taskNames) +
  geom_quasirandom(dodge.width=.7, size = .5) +
  geom_violin(alpha=.6, position = position_dodge(.7)) +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .7),
               geom = 'errorbar',
               size = .9, aes(width = .2)) +
  ylab("Global Self-performance\nEstimate\n(baseline-corrected)") +
  xlab("Feedback") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  labs(color = 'Feedback') +
  theme(text = element_text(size=font_size-6),
        legend.position = 'none',
        plot.caption = element_text(size = font_size-6),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  coord_cartesian(ylim = c(-.65,.65))  +
  geom_signif(y_position = c(.55), vjust = -.3, hjust=.4,
              xmin = c(1.1), color = 'black',
              xmax = c(1.9),
              annotation = c("p < .0001"), tip_length = .04,
              textsize = pval_size-2, size = .6)
f2a.left

f2a.right <- ggplot(exp1.df.2, aes(x = fbblock, y = accu, color = fbblock)) +
  facet_wrap(~task) +
  scale_color_manual(breaks = posnegNames, values=posnegColours) +
  theme_blank()  +
  # scale_x_discrete(limits = taskNames) +
  geom_quasirandom(dodge.width=.7, size = .5) +
  geom_violin(alpha=.6, position = position_dodge(.7)) +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .7),
               geom = 'errorbar',
               size = .9, aes(width = .2)) +
  ylab("\nActual performance\n(accuracy)") +
  xlab("Feedback") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  labs(color = 'Feedback') + 
  theme(text = element_text(size=font_size-6), 
        legend.position = 'none',
        plot.caption = element_text(size = font_size-6),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  # coord_cartesian(ylim = c(.6,1))  +
  geom_signif(y_position = c(.87), vjust = -.3,
              xmin = c(1.1), color = 'black',
              xmax = c(1.9),
              annotation = c("n.s."), tip_length = .04,
              textsize = pval_size-2, size = .6)
f2a.right

f2a.conf <- ggplot(exp1.df.2, aes(x = fbblock, y = conf_b, color = fbblock)) +
  facet_wrap(~task) +
  geom_hline(yintercept = 0, color = 'grey30') +
  scale_color_manual(breaks = posnegNames, values=posnegColours) +
  theme_blank()  +
  # scale_x_discrete(limits = taskNames) +
  geom_quasirandom(dodge.width=.7, size = .5) +
  geom_violin(alpha=.6, position = position_dodge(.7)) +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .7),
               geom = 'errorbar',
               size = .9, aes(width = .2)) +
  ylab("\nLocal confidence\n(baseline-corrected)") +
  xlab("Feedback") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  labs(color = 'Feedback') + 
  theme(text = element_text(size=font_size-6), 
        legend.position = 'none',
        plot.caption = element_text(size = font_size-6),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  coord_cartesian(ylim = c(-.5,.4))  +
  geom_signif(y_position = c(.32), vjust = -.3, hjust=.4,
              xmin = c(1.1), color = 'black',
              xmax = c(1.9),
              annotation = c("n.s."), tip_length = .04,
              textsize = pval_size-2, size = .6)
f2a.conf

##################
### Plot raw data across all blocks
exp1.df11 <- exp1.mbs %>%
  group_by(subj, task, runnum, fbblock.notrans, group) %>%
  dplyr::summarise(
    conf = mean(conf),
    accu = mean(accu, na.rm=T),
    spe = mean(spe),
    stair = mean(stair))

ggplot(exp1.df11, aes(x = runnum, y = spe, color = fbblock.notrans)) +
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

ggplot(exp1.df11, aes(x = runnum, y = accu, color = fbblock.notrans)) +
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

ggplot(exp1.df11, aes(x = runnum, y = conf, color = fbblock.notrans)) +
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


##################################
####### Supp Figure 5 (upper panel)

################ difficulty level
sf5a.upper <- ggplot(exp1.df.2, aes(x = task, y = stair_b, color = fbblock)) +
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
  ggtitle('Exp 1') +
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
sf5a.upper


##################################
##### Figure 4a. plot transfer of confidence as time series

interv.labs <-  c('Intervention: Perception', 'Intervention: Memory')
names(interv.labs) <- c('Perception', 'Memory')
test.labs <-  c('Test: Perception', 'Test: Memory')
names(test.labs) <- c('Perception', 'Memory')

exp1.df.2 <- filter(exp1.mbs, blockType=="transf") %>% 
  mutate(intervTask = relevel(intervTask, "Perception")) %>%
  mutate(testTask = relevel(testTask, "Perception"))

f4a <- ggplot(exp1.df.2, aes(x = trialnum, y = conf_b, colour = fbblock)) +
  geom_hline(yintercept = 0, color = 'gray') +
  scale_color_manual(breaks = posnegNames, values=posnegColours) +
  scale_linetype_manual(breaks = c('Same', 'Opposite'), 
                        values = c('solid', 'twodash'),
                        labels = c('Same\ntask', 'Opposite\ntask')) +
  facet_grid(testTask~intervTask,
             labeller = labeller(intervTask = interv.labs, testTask = test.labs)) + 
  theme_pubclean() +
  geom_smooth(method='loess', span=.6, 
              aes(linetype=transferType),se=F, size=1.4)+
  ylab("Confidence (test)") +
  xlab("Trial number") +
  labs(linetype = 'Transfer to', colour = 'Interv. block\nfeedback') +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = .5),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15)) +
  guides(color = guide_legend(nrow = 2))
f4a


##################################
##### Figure 4c. transfer of feedback to test block confidence as pos-neg dif

exp1.df.2 <- exp1.mbs %>% 
  filter(blockType=="transf", trialnum %in% c(1:20)) %>% 
  group_by(intervTask, testTask, fbblock, subj) %>%
  summarise(conf = mean(conf)) %>%
  pivot_wider(names_from = fbblock,
    values_from = conf) %>% group_by(subj, intervTask, testTask) %>%
  summarise(confPosNeg = mean(Positive) - mean(Negative)) %>% 
  mutate(intervTask = relevel(intervTask, "Perception")) %>%
  mutate(testTask = relevel(testTask, "Perception"))

f4b <- ggplot(exp1.df.2, aes(x = intervTask, y = confPosNeg, 
                             colour = testTask)) +
  geom_hline(yintercept = 0, color = 'grey30') +
  scale_colour_manual(values=percmemColours, breaks = taskNames) +
  geom_quasirandom(dodge.width=.8, alpha = .7) +
  geom_violin(alpha=.65, position = position_dodge(.8)) +
  stat_summary(fun.data = mean_cl_boot,
               position=position_dodge(width = .8),
               geom = 'errorbar',
               size = .9, aes(width = .2)) +
  coord_cartesian(ylim = c(-.3,.65)) +
  ylab("Confidence (test)\nPositive - Negative") +
  xlab("Intervention task") +
  labs(color = 'Test task') +# ggtitle('Exp 1') +
  theme_pubclean() +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
geom_signif(y_position = c(.58, .32, .3, .3), vjust = -.3,
            xmin = c(.7, 1.1, 1.7,2.1),
            xmax = c(.9, 1.3, 1.9,2.3),
            annotation = c("p = .0004", "p = .030",'p = .78','p = .059'), 
            tip_length = .0, color = 'black',
            textsize = pval_size-1, size = .8)
f4b


##################################
##### Supp Figure 14a. transfer of feedback to test block SPE

exp1.df.2 <- exp1.mbs %>%
  filter(blockType=="transf") %>%
  group_by(fbblock, subj, intervTask, testTask) %>%
  summarise(spe = mean(spe)) %>%
  pivot_wider(
    names_from = fbblock,
    values_from = spe) %>% group_by(subj, intervTask, testTask) %>%
  summarise(spe = mean(Positive) - mean(Negative)) %>%
  mutate(intervTask = factor(intervTask, levels= c('Perception', 'Memory')))

f14a <- ggplot(exp1.df.2, aes(x = intervTask, y = spe, color = testTask)) +
  geom_hline(yintercept = 0, color = 'grey30') +
  scale_color_manual(values=percmemColours, breaks = taskNames) +
  geom_violin(alpha=.8, position = position_dodge(.8)) +
  stat_summary(fun.data = mean_se,
               position=position_dodge(width = .8),
               geom = 'errorbar',
               size = 1, aes(width = .1))  +
  geom_quasirandom(dodge.width=.8,
                   alpha = .15) +
  # scale_x_discrete(limits = taskNames) +
  coord_cartesian(ylim = c(-.6,.5)) +
  ylab("SPE-b (test)\nPositive - Negative") +
  xlab("Intervention task") +
  labs(color = 'Test task') +
  theme_pubclean() +
  theme(text = element_text(size=font_size),
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  geom_signif(y_position = c(.4, .4, .4, .4),
              xmin = c(.7, 1.1, 1.7,2.1),
              xmax = c(.9, 1.3, 1.9, 2.3),
              annotation = c('', '', '', ''),
              tip_length = 0, color = 'black', textsize = pval_size, size = .8)+
  geom_signif(y_position = c(.43, .43), vjust = -.3,
              xmin = c(.8,1.2),
              xmax = c(2.2, 1.8),
              annotation = c('p = .0007', ''),
              tip_length = 0.03, color = 'black', textsize = pval_size, size = .8)

f14a


############
###############################################
### SPE (intervention) analyses
################################################

exp1.df.2 <- filter(exp1.df, blockType=="interv")

## factorial effect of fb on spe

m.reg <- lmer(spe_b ~ fbblock*task*firstFeedback + accu_b*task + stair_b*task + 
                (1|subj) + (1|group) + (1|runnum), exp1.df.2)
summary(m.reg)
m.reg <- lmer(spe_b ~ fbblock*task*firstFeedback + accu_b + stair_b + 
                (1|subj) + (1|group), REML = F, exp1.df.2)
summary(m.reg)
m.reg <- lmer(spe_b ~ fbblock*task*firstFeedback + accu_b + stair_b + 
                (1|subj), REML = F, exp1.df.2)
summary(m.reg)
plot(m.reg)

m2 <- update(m.reg, ~.-fbblock:task:firstFeedback)
anova(m.reg,m2)
m.reg <- m2 # accept the reduced model
summary(m.reg)
m2 <- update(m.reg, ~.-fbblock:task)
anova(m.reg,m2)
m.reg <- m2 # accept the reduced model
summary(m.reg)
m2 <- update(m.reg, ~.-fbblock:firstFeedback)
anova(m.reg,m2)
m.reg <- m2 # accept the reduced model
summary(m.reg)
m2 <- update(m.reg, ~.-task:firstFeedback)
anova(m.reg,m2)
m.reg <- m2 # accept the reduced model
summary(m.reg)
m.reg <- update(m.reg, REML=T)
summary(m.reg)

plot(m.reg)
qqPlot(resid(m.reg))


### as none of the 3- and 2-way interactions are significant evaluate as only main effect

## mixed model with subject as intercept singular for m2 so run without random effects
m1 <- lm(spe_b ~ fbblock + accu_b + stair_b, exp1.df.2)
summary(m1)
confint(m1)

m2 <- lm(spe_b ~ accu_b + stair_b, exp1.df.2)
lrtest(m1,m2)

library(effectsize)
cohens_d(spe_b ~ fbblock, data = exp1.df.2)


m1 <- lmer(accu_b ~ fbblock*task + (1|subj), exp1.df.2)
summary(m1)
emm <- emmeans(m1, ~fbblock|task)
contrast(emm)

m1 <- lmer(accu_b ~ fbblock+task + (1|subj), exp1.df.2)
summary(m1)
confint(m1)
m2 <- lmer(accu_b ~ task + (1|subj), exp1.df.2)
summary(m2)
# m2 <- lmer(accu_b ~ fbblock*task + (1|subj), exp1.df.2)
anova(m1,m2)


m1.bf <- lmBF(accu_b ~ fbblock+task+subj, whichRandom = c('subj'), 
                          iterations = 50000, exp1.df.2)
m2.bf <- lmBF(accu_b ~ task+subj, whichRandom = c('subj'), 
              iterations = 50000, exp1.df.2)

m2.bf/m1.bf


m1 <- lmer(stair_b ~ fbblock*task + (1|subj), exp1.df.2)
summary(m1)
emm <- emmeans(m1, ~fbblock|task)
contrast(emm)

m1 <- lmer(stair_b ~ fbblock+task + (1|subj), exp1.df.2)
summary(m1)
confint(m1)
m2 <- lmer(stair_b ~ task + (1|subj), exp1.df.2)
summary(m2)
anova(m1,m2)
m2 <- lmer(stair_b ~ fbblock*task + (1|subj), exp1.df.2)
anova(m1,m2)


m1.bf <- lmBF(stair_b ~ fbblock+task+subj, whichRandom = c('subj'), 
              iterations = 50000, exp1.df.2)
m2.bf <- lmBF(stair_b ~ task+subj, whichRandom = c('subj'), 
              iterations = 50000, exp1.df.2)

m1.bf/m2.bf
m2.bf/m1.bf
############
### Local confidence during test blocks

## Effect of feedback type (fbblock) on local confidence (test)
exp1.df15 <- filter(exp1.mbs, blockType=="interv")

exp1.mbs.reg = lmer(conf_b ~ fbblock*task + accu + stair_b +
                      (1 + fbblock + accu + stair_b|subj), 
                    exp1.df15)
summary(exp1.mbs.reg)
emm <- emmeans(exp1.mbs.reg, ~fbblock|task)
contrast(emm)

m1 = lmer(conf_b ~ fbblock + task + accu + stair_b +
            (1 + fbblock |subj), 
          REML = F,exp1.df15)
summary(m1)
confint(m1)

m1 = lmer(conf_b ~ fbblock + task + accu + stair_b +
            (1 + fbblock + accu + stair_b|subj), 
          REML = F,exp1.df15)
summary(m1)
confint(m1)
m2 = lmer(conf_b ~ task + accu + stair_b +
            (1 + fbblock + accu + stair_b|subj), REML = F,
          exp1.df15)

anova(m1, m2)



######## do local conf and accu predict spe

m1 <- lmer(spe_b ~ conf_b + accu_b + #(1|task) + (1|group) + 
             (conf_b + 1|subj), REML = F, exp1.df)
summary(m1)
confint(m1)
m2 <- lmer(spe_b ~ accu_b + #(1|task) +(1|group) + 
             (conf_b + 1|subj), REML = F, exp1.df)
summary(m2)
anova(m1,m2)

m2 <- lmer(spe_b ~ conf_b + #(1|task) +(1|group) + 
             (conf_b +1|subj), REML = F, exp1.df)
summary(m2)
anova(m1,m2)


############
############
### Local confidence during test blocks

## Effect of feedback type (fbblock) on local confidence (test)

first_n_of_trans <- c(1:40)
exp1.df3 <- filter(exp1.mbs, blockType=="transf",
             trialnum %in% first_n_of_trans)

####

exp1.mbs.reg = lmer(conf_b ~ fbblock*task*transferType + accu + stair_b +
                      (1 + fbblock + accu + stair_b|subj) + (1|group), exp1.df3)
summary(exp1.mbs.reg)

m2 <- update(exp1.mbs.reg, ~.-fbblock:task:transferType)
anova(exp1.mbs.reg,m2)

summary(exp1.mbs.reg)$coefficients

emm <- emmeans(exp1.mbs.reg, ~fbblock|transferType*task)
contrast(emm)

emm <- emmeans(exp1.mbs.reg, ~fbblock|transferType)
contrast(emm)



####
first_n_of_trans <- c(1:20)
exp1.df3 <- filter(exp1.mbs, blockType=="transf",
                   trialnum %in% first_n_of_trans)

exp1.mbs.reg = lmer(conf_b ~ fbblock*task*transferType + 
                      accu + stair_b +
                      (1 + fbblock + accu + stair_b|subj) + (1|group), exp1.df3)
summary(exp1.mbs.reg)
m2 <- update(exp1.mbs.reg, ~.-fbblock:task:transferType)
anova(exp1.mbs.reg,m2)


summary(exp1.mbs.reg)$coefficients
emm <- emmeans(exp1.mbs.reg, ~fbblock|transferType*task)
contrast(emm)
confint(contrast(emm))

emm <- emmeans(exp1.mbs.reg, ~fbblock|transferType)
contrast(emm)


# ## estimate power for Exp 2
# exp1.df15 = filter(exp1.df3, task=='Memory', transferType=='Opposite')
# exp1.df15$subject = as.numeric(exp1.df15$subj)
# exp1.mbs.reg = lmer(conf ~ fbblock + 
#                  accu + scale_rt + stair +
#                  (1|subj) + (1|group) + (1|runnum) + (1|trialnum), 
#                exp1.df15
#                )
# summary(exp1.mbs.reg)$coefficients
# 
# SESOI <- c(-0.005597832, -0.016, 0.158501287, -0.615647952, 0.021358216) # specify SESOI
# 
# power_SESOI.1 <- mixedpower(model = exp1.mbs.reg, data = exp1.df15,
#                           fixed_effects = c("fbblock", "accu", 'scale_rt', 'stair'),
#                           simvar = "subject", steps = c(150,160,170),
#                           critical_value = 2, n_sim = 100,
#                           SESOI = SESOI, databased = F)



############
## Mediation analysis to test if SPE (intervention block) mediates the effect
## of pos and neg fb on local confidence (test block)

# Mediation analysis cannot be performed on trial-by-trial conf data so average
# over runs -- also because the 

first_n_of_trans <- c(1:40)

exp1.df7 <- filter(exp1.mbs, blockType=="transf",
             trialnum %in% first_n_of_trans) %>%
  group_by(subj, runnum, fbblock, transferType, task, group) %>%
  dplyr::summarise(
    spe = mean(spe),
    accu = mean(accu),
    conf = mean(conf),
    stair = mean(stair),
    scale_rt = mean(scale_rt, na.rm = T)
  )

# num of divisions of intervention block
ndivs_per_run <- 1

exp1.df.temp <- filter(exp1.mbs, blockType=="interv")

# create a new factor by combining run and trialnumber and then divide it into equal divisions
exp1.df.temp <- exp1.df.temp %>%
  mutate(rundiv = factor(ceiling(trialnum/(40/ndivs_per_run))))

exp1.df.temp <- exp1.df.temp %>%
  group_by(rundiv, subj, task, runnum, fbblock, 
           transferType) %>%
  dplyr::summarise(
    conf = mean(conf),
    accu = mean(accu),
    scale_rt = mean(scale_rt, na.rm=T),
    fbneg = mean(fbneg),
    fbpos = mean(fbpos),
    spe = mean(spe),
    stair = mean(stair)
  ) %>% .
  ungroup() %>%
  mutate_at("rundiv", as.factor)

exp1.df.temp <- exp1.df.temp %>% dplyr::select(
  subj, rundiv, conf, spe, accu, stair, scale_rt,
  fbneg, fbpos,fbblock, transferType
) %>%
  pivot_wider(
    names_from = c( rundiv),
    values_from = c(conf, fbneg, fbpos, spe, accu, stair, scale_rt)
  ) 

exp1.df8 <- exp1.df7 %>% left_join(exp1.df.temp)


## mediation combined
med.m1 <- lmer(conf ~ fbpos_1 + fbneg_1 + (1|subj), exp1.df8)
summary(med.m1)

med.m2 <- lmer(spe_1 ~ fbpos_1 + fbneg_1 + (1|subj), exp1.df8)
summary(med.m2)

med.m3 <- lmer(conf ~ fbpos_1 + fbneg_1 + spe_1 + (1|subj), exp1.df8)
summary(med.m3)

med.pos = mediate(med.m2, med.m3, treat='fbpos_1', mediator='spe_1', 
                  boot=F, sims = 5000)
summary(med.pos)

med.pos = mediate(med.m2, med.m3, treat='fbneg_1', mediator='spe_1', 
                  boot=F, sims = 5000)
summary(med.pos)


############
############
############
### SPE (test block) analyses

## Effect of fbblock on spe (test)

exp1.df9 <- filter(exp1.mbs, blockType=="transf") %>%
  group_by(subj, runnum, fbblock, transferType, task, group) %>%
  dplyr::summarise(
    spe = mean(spe),
    spe_b = mean(spe_b),
    accu = mean(accu),
    conf = mean(conf),
    scale_rt = mean(scale_rt, na.rm = T),
    stair = mean(stair),
    phq = mean(phq),
    gad = mean(gad),
    spin = mean(spin))

exp1.mbs.reg = lmer(spe_b ~ fbblock*task*transferType + 
                      accu + scale_rt + stair +
                      (1|subj), REML = F,
                    exp1.df9)
plot(exp1.mbs.reg)
summary(exp1.mbs.reg)

m2 <- update(exp1.mbs.reg, ~.-fbblock:transferType:task)
anova(exp1.mbs.reg,m2)
exp1.mbs.reg <- m2 # accept the reduced model
summary(exp1.mbs.reg)

m2 <- update(exp1.mbs.reg, ~.-transferType:task)
anova(exp1.mbs.reg,m2)
exp1.mbs.reg <- m2 # accept the reduced model
summary(exp1.mbs.reg)

m2 <- update(exp1.mbs.reg, ~.-fbblock:task)
anova(exp1.mbs.reg,m2)
exp1.mbs.reg <- m2 # accept the reduced model
summary(exp1.mbs.reg)

m2 <- update(exp1.mbs.reg, ~.-fbblock:transferType)
anova(exp1.mbs.reg,m2)
exp1.mbs.reg <- m2 # accept the reduced model
summary(exp1.mbs.reg)

exp1.mbs.reg <- update(exp1.mbs.reg, REML=T)
summary(exp1.mbs.reg)

m2 <- update(exp1.mbs.reg, ~.-fbblock)
anova(exp1.mbs.reg,m2)

plot_model(exp1.mbs.reg,
           title="Perception task - interaction",
           show.p=TRUE, show.values = T,
           value.offset = .4, value.size = 3.8) + 
  coord_cartesian(ylim = c(-.65, .65))

emm <- emmeans(exp1.mbs.reg, pairwise~fbblock)
test(emm)

#######################
#### to address question from r4 about any difference in reported AD scores based on the order of feedback

exp1.df4 <- exp1.mbs %>% group_by(subj, firstFeedback) %>%
  summarise(phq = mean(phq),
            gad = mean(gad))
t.test(gad ~ firstFeedback, exp1.df4)
t.test(phq ~ firstFeedback, exp1.df4)

#######################################
###################################
#### baseline correlation with mhq

exp1.df4 <- exp1.mbs.bas %>% group_by(task, runnum, group, subj, gender) %>%
  summarise(conf = mean(conf),
            accu = mean(accu),
            stair = mean(stair),
            scale_rt = mean(scale_rt, na.rm = T),
            spe = mean(spe),
            phq = mean(phq),
            gad = mean(gad),
            spin = mean(spin),
            age = mean(age)) %>%
  drop_na()
exp1.df4$confZ <- scale(exp1.df4$conf)
exp1.df4$speZ <- scale(exp1.df4$spe)

m.reg = lmer(speZ ~ phq + age + gender + accu + stair + scale_rt + 
               ( 1|task),exp1.df4)
summary(m.reg)
confint(m.reg)
qqPlot(resid(m.reg))

# exp1.df4$spe.phq.partial <- coef(m.reg)[1] + coef(m.reg)['phq']*exp1.df4$phq + resid(m.reg)
exp1.df4$spe.phq.partial <- keepef(m.reg, fix = 'phq')

ggplot(exp1.df4, aes(x=phq, y=spe.phq.partial, color = task))+
  geom_point()+
  geom_smooth(method = lm) +
  stat_cor(method = "pearson", label.x = 9, label.y = c(1.8, 2.4), 
           size = 5, aes(color = task))  +
  theme_classic()+
  labs(color = 'Task')+
  ggtitle('Exp 1') +
  ylab('Global SPE\nBaseline blocks\n(partialed, z-scored)') +
  xlab('PHQ-9') +
  theme(text = element_text(size=font_size-2),
        legend.position = 'top')

m2 <- lmer(speZ ~ age + gender + accu + stair + scale_rt + (1|task),exp1.df4)
anova(m.reg,m2)

m.reg = lmer(speZ ~ gad + age + gender + accu + stair + scale_rt + (1|task),exp1.df4)
summary(m.reg)
confint(m.reg)
qqPlot(resid(m.reg))
anova(m.reg,m2)

# exp1.df4$spe.gad.partial <- coef(m.reg)[1] + coef(m.reg)['gad']*exp1.df4$gad + resid(m.reg)
exp1.df4$spe.gad.partial <- keepef(m.reg, fix = 'gad')

ggplot(exp1.df4, aes(x=gad, y=spe.gad.partial, color = task))+
  geom_point()+
  geom_smooth(method = lm) +
  stat_cor(method = "pearson", label.x = 9, label.y = c(1.8, 2.4), 
           size = 5, aes(color = task))  +
  theme_classic()+
  labs(color = 'Task')+
  ggtitle('Exp 1') +
  ylab('Global SPE\nBaseline blocks\n(partialed, z-scored)') +
  xlab('GAD-7') +
  theme(text = element_text(size=font_size-2),
        legend.position = 'top')

m.reg = lmer(speZ ~ spin + accu + stair + age + gender + scale_rt + (1|task),exp1.df4)
anova(m.reg,m2)


## local confidence
m.reg = lmer(confZ ~ phq + age + gender + 
               accu + stair + scale_rt + (1|task), exp1.df4[-c(68),]) # outlier on row 68
summary(m.reg)
confint(m.reg)
qqPlot(resid(m.reg))

m2 = lmer(confZ ~ age + gender + accu + stair + scale_rt + (1|task), exp1.df4[-c(68),]) # outlier on row 68
# lrtest(m.reg,m2)
anova(m.reg,m2)

exp1.df4$conf.phq.partial <- exp1.df4$conf
# exp1.df4$conf.phq.partial[-c(68)] <- coef(m.reg)[1] + 
#   coef(m.reg)['phq']*exp1.df4[-c(68),]$phq + resid(m.reg)
exp1.df4[-c(68),]$conf.phq.partial <- keepef(m.reg, fix = 'phq')

ggplot(exp1.df4[-c(68),], aes(x=phq, y=conf.phq.partial, color = task))+
  geom_point()+
  geom_smooth(method = lm) +
  stat_cor(method = "pearson", label.x = 11, label.y = c(3.4, 4.1), 
           size = 5, aes(color = task))  +
  theme_classic()+
  labs(color = 'Task')+
  ggtitle('Exp 1') +
  ylab('Local confidence\nBaseline blocks\n(partialed, z-scored)') +
  xlab('PHQ-9') +
  theme(text = element_text(size=font_size-2),
        legend.position = 'top')


m.reg = lmer(confZ ~ gad + age + gender + accu + stair+ scale_rt +(1|task), 
             exp1.df4[-c(68),]) # outlier on row 68
summary(m.reg)
confint(m.reg)
qqPlot(resid(m.reg))

m2 = lmer(confZ ~ age + gender + accu + stair+scale_rt + (1|task), exp1.df4[-c(68),]) # outlier on row 68
anova(m.reg,m2)

exp1.df4$conf.gad.partial <- exp1.df4$conf
# exp1.df4$conf.gad.partial[-c(68)] <- coef(m.reg)[1] + 
#   coef(m.reg)['gad']*exp1.df4[-c(68),]$gad + resid(m.reg)
exp1.df4[-c(68),]$conf.gad.partial <- keepef(m.reg, fix = 'gad')

ggplot(exp1.df4[-c(68),], aes(x=gad, y=conf.gad.partial, color = task))+
  geom_point()+
  geom_smooth(method = lm) +
  stat_cor(method = "pearson", label.x = 9, label.y = c(3.4, 4.1), 
           size = 5, aes(color = task))  +
  theme_classic()+
  labs(color = 'Task')+
  ggtitle('Exp 1') +
  ylab('Local confidence\nBaseline blocks\n(partialed, z-scored)') +
  xlab('GAD-7') +
  theme(text = element_text(size=font_size-2),
        legend.position = 'top')


m.reg = lmer(confZ ~ spin + age + gender + accu + stair+ scale_rt +(1|task),exp1.df4[-c(68),])
summary(m.reg)
anova(m.reg,m2)


###

####################################
####################################
###### model-free predictions
#####################################

# update this folder with the folder location with data
data.folder <- "~/OneDrive - University of Copenhagen/Projects/Experiments/metaBiasShift/data/exp"
expNum <- 1
setwd(paste(data.folder, expNum, sep = ''))
setwd('processed/')
mbsDataExp1 = readMat(paste('mbsDataExp',expNum,'.mat',sep='')) # load data
mbsDataExp1 <- mbsDataExp1[[paste('mbsDataExp',expNum,sep='')]]
mbsDataFit <- mbsDataExp1[,,1]$fitRegZ
# spe.est <- mbsDataFit[,,1]$model.10.q2[,,1]$spe.est
spe.est <- mbsDataFit[,,1]$model.6.q2[,,1]$spe.est
nsubj <- dim(spe.est)[1]
nruns <- dim(spe.est)[2]
subj <- array(1:nsubj, c(nsubj,nruns))
runnum <- t(array(1:nruns, c(nruns,nsubj)))
exp1.fit.data <- data.frame(c(spe.est), c(subj), c(runnum))
names(exp1.fit.data) <- c('spe.fit', 'subj', 'runnum')
exp1.fit.data <- exp1.fit.data %>%
  mutate(subj = as.factor(subj)) %>%
  mutate(runnum = as.factor(runnum))


### test

exp1.model <- exp1.mbs %>% drop_na() %>% 
  dplyr::select(c(subj, runnum, trialnum, spe, confZ, gad, phq, accu, fbblock.notrans)) %>%
  mutate(confhilo = confZ > 0) %>%
  left_join(exp1.fit.data)

exp1.reg.obs.lo <- lmer(spe ~ gad*confZ + (1|subj) + (1|runnum), 
                        filter(exp1.model, confhilo==F))
summary(exp1.reg.obs.lo)

exp1.reg.fit.lo <- lmer(spe.fit ~ gad*confZ + (1|subj) + (1|runnum), 
                        filter(exp1.model, confhilo==F))
summary(exp1.reg.fit.lo)


exp1.reg.fit <- lmer(spe.fit ~ gad*confZ + (1|subj) + (1|runnum) + (1|trialnum), 
                 filter(exp1.model, confhilo==T))
summary(exp1.reg.fit)

exp1.reg.obs <- lmer(spe ~ gad*confZ + (1|subj) + (1|runnum), 
                 filter(exp1.model, confhilo==T))
summary(exp1.reg.obs)


exp1.model <- exp1.mbs %>% drop_na() %>% 
  dplyr::select(c(subj, runnum, trialnum, spe, confZ, gad, phq, accu, fbblock.notrans)) %>%
  mutate(confhilo = confZ > 0) %>%
  dplyr::group_by(subj) %>%
  mutate(conf.tile = ntile(confZ, 4)) %>%
  group_by(runnum) %>%
  mutate(AD.tile = ntile(gad, 3)) %>%
  mutate_at(c('AD.tile', 'runnum'), as.factor) %>%
  left_join(exp1.fit.data)%>%
  rename(spe.data = spe) %>%
  pivot_longer(cols = c('spe.data', 'spe.fit'),
               values_to = 'spe',
               names_to = 'data.fit') %>%
  mutate_at(c('data.fit'), as.factor)

AD.tile.cents  <- c(mean(exp1.model$gad[exp1.model$AD.tile==1]),
                    mean(exp1.model$gad[exp1.model$AD.tile==3]))

conf.tile.cents <- c(mean(exp1.model$confZ[exp1.model$conf.tile==1]),
                     mean(exp1.model$confZ[exp1.model$conf.tile==2]))

emm.obs.lo <- summary(emmeans(exp1.reg.obs.lo, ~ confZ|gad, 
                              at = list(confZ = conf.tile.cents, gad = AD.tile.cents)))
emm.obs.lo

emm.fit.lo <- summary(emmeans(exp1.reg.fit.lo, ~ confZ|gad, 
                              at = list(confZ = conf.tile.cents, gad = AD.tile.cents)))
emm.fit.lo

conf.tile.cents <- c(mean(exp1.model$confZ[exp1.model$conf.tile==3]),
                     mean(exp1.model$confZ[exp1.model$conf.tile==4]))
# conf.tile.cents <- c(1,2.5)
emm.obs <- summary(emmeans(exp1.reg.obs, ~ confZ|gad, 
                       at = list(confZ = conf.tile.cents, gad = AD.tile.cents)))
emm.obs

emm.fit <- summary(emmeans(exp1.reg.fit, ~ confZ|gad, 
                           at = list(confZ = conf.tile.cents, gad = AD.tile.cents)))
emm.fit


spe <- c(emm.obs.lo$emmean,   emm.fit.lo$emmean,   emm.obs$emmean,   emm.fit$emmean)
lci <- c(emm.obs.lo$lower.CL, emm.fit.lo$lower.CL, emm.obs$lower.CL, emm.fit$lower.CL)
uci <- c(emm.obs.lo$upper.CL, emm.fit.lo$upper.CL, emm.obs$upper.CL, emm.fit$upper.CL)

conf.tile <- c(rep(c(1:2,1:2),2), rep(c(3:4,3:4),2))
AD.tile <- rep(c(rep('Low',2), rep('High',2)),4)
conftile.hilo <- c(rep(1,8), rep(2,8))
data.fit <- rep(c(rep('Data',4), rep('Fit',4)),2)
exp1.model.df <- data.frame(spe, lci, uci, 
                            conf.tile, AD.tile, conftile.hilo, data.fit)
exp1.model.df$conftile.hilo[exp1.model.df$data.fit=='Fit'] <- 0
exp1.model.df$conftile.hilo <- as.factor(exp1.model.df$conftile.hilo)

ggplot(filter(exp1.model.df),
       aes(x = conf.tile, y = spe, colour = AD.tile, shape = data.fit, 
           alpha = data.fit, linetype = conftile.hilo)) +
  scale_alpha_manual(values = c(1,.6)) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', 'royalblue2')) +
  scale_shape_manual(values = c(19,2)) +
  geom_point(position=position_dodge(.5), size=3.2, stroke = .5) +
  geom_errorbar(aes(ymin = lci, ymax = uci), linetype = 1,
                position=position_dodge(.5), linewidth=.6, width=.3) +
  stat_summary(geom="line", position=position_dodge(.5), size=.5) +
  scale_linetype_manual(values = c('blank', 'dashed', 'solid')) +
  theme_classic2() +
  labs(color = 'AD score', linetype = 'Confidence level')+
  xlab('Confidence (tiled)') +
  ylab('Global SPE') +
  scale_y_continuous(breaks = c(.5,.6)) +
  theme(text = element_text(size=font_size-6),
        legend.direction = 'vertical',
        legend.text=element_text(size=font_size-6),
        plot.caption = element_text(size = font_size-2),
        legend.position = 'none') +
  geom_segment(aes(x = 3.6, y = .555, xend = 3.6, yend = .605), color = 'grey60') +
  geom_segment(aes(x = 3.45, y = .58, xend = 3.6, yend = .58), color = 'grey60') +
  annotate('text', label = '****', x = 3.45, y = .58, size=6, angle = 90)



### original 
exp1.model <- exp1.mbs %>% drop_na() %>%
  dplyr::select(c(subj, runnum, trialnum, spe, confZ, gad, phq, accu, fbblock.notrans)) %>%
  mutate(confhilo = confZ > 0) %>%
  group_by(subj,trialnum) %>%
  mutate(spe.tile = ntile(spe, 6))%>%
  mutate(conf.tile = ntile(confZ, 4)) %>%
  group_by(runnum) %>%
  mutate(AD.tile = ntile(gad, 3)) %>%
  mutate_at(c('AD.tile', 'runnum'), as.factor) %>%
  left_join(exp1.fit.data)

# 

# exp1.reg <- lmer(confZ ~ spe*gad + accu + (1|subj) + (1|runnum)+ (1+spe*gad|trialnum), 
#                  filter(exp1.model, confhilo==T))
exp1.reg <- lmer(confZ ~ spe*gad + accu + (1|subj) + (1|runnum)+ (1|trialnum), 
                 filter(exp1.model, confhilo==T))
summary(exp1.reg)
confint(exp1.reg)
plot(exp1.reg)
# plot_model(exp1.reg,type = 'pred', terms= c('spe', 'gad'))

m2 <- update(exp1.reg, ~.-spe:gad)
anova(exp1.reg,m2)


exp1.reg <- lmer(confZ ~ spe*phq + accu + (1|subj) + (1|runnum)+ (1|trialnum), 
                 filter(exp1.model, confhilo==T))
summary(exp1.reg)
confint(exp1.reg)
# plot_model(exp1.reg,type = 'pred', terms= c('spe', 'phq'))

m2 <- update(exp1.reg, ~.-spe:phq)
anova(exp1.reg,m2)


##### model-free analysis of feedback distortion

exp1.model2 <- filter(exp1.mbs) %>% drop_na() %>%
  group_by(subj, fbblock.notrans, runnum) %>%
  summarise(spe = mean(spe),
            gad = mean(gad),
            phq = mean(phq)) %>%
  group_by(runnum) %>%
  mutate(AD.tile = ntile(gad, 3)) %>%
  mutate_at(c('AD.tile'), as.factor) %>%
  mutate(fbblock.notrans = 
           factor(fbblock.notrans, levels = c('Negative','None', 'Positive'))) %>%
  mutate(AD.tile = recode_factor(AD.tile, '1' = 'Low', '3' = 'High')) %>%
  left_join(exp1.fit.data) %>%
  rename(spe.data = spe) %>%
  pivot_longer(cols = c('spe.data', 'spe.fit'),
               values_to = 'spe',
               names_to = 'data.fit') %>%
  mutate_at(c('data.fit'), as.factor)

AD.tile.cents  <- c(mean(exp1.model2$gad[exp1.model2$AD.tile=='Low']),
                    mean(exp1.model2$gad[exp1.model2$AD.tile=='High']))


exp1.fb.obs <- lmer(spe ~ fbblock.notrans*gad + (1|subj) + (1|runnum), filter(exp1.model2, data.fit=='spe.data'))
emm.fb.obs <- summary(emmeans(exp1.fb.obs, ~ fbblock.notrans|gad, at = list(gad = AD.tile.cents)))

exp1.fb.fit <- lmer(spe ~ fbblock.notrans*gad + (1|subj) + (1|runnum), filter(exp1.model2, data.fit=='spe.fit'))
emm.fb.fit <- summary(emmeans(exp1.fb.fit, ~ fbblock.notrans|gad, at = list(gad = AD.tile.cents)))

spe <- c(emm.fb.obs$emmean,   emm.fb.fit$emmean)
lci <- c(emm.fb.obs$lower.CL, emm.fb.fit$lower.CL)
uci <- c(emm.fb.obs$upper.CL, emm.fb.fit$upper.CL)

fb.type <- rep(c('Negative', 'None', 'Positive'),4)
AD.tile <- rep(c(rep('Low',3), rep('High',3)),2)
data.fit <- c(rep('Data',6), rep('Fit',6))
exp1.model.df <- data.frame(spe, lci, uci, 
                            fb.type, AD.tile, data.fit)
exp1.model.df$AD.tile <- as.factor(exp1.model.df$AD.tile)
exp1.model.df <- exp1.model.df %>%
  mutate(AD.tile = factor(AD.tile, levels = c('Low','High'))) 

ggplot(filter(exp1.model.df),
       aes(x = fb.type, y = spe, colour = AD.tile, shape = data.fit, alpha = data.fit)) +
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
  geom_segment(aes(x = 2.1, y = .61, xend = 2.9, yend = .61), color = 'grey60') +
  geom_segment(aes(x = 2.1, y = .58, xend = 2.1, yend = .61), color = 'grey60') +
  geom_segment(aes(x = 2.9, y = .61, xend = 2.9, yend = .64), color = 'grey60') +
  annotate('text', label = 'n.s.', x = 2.5, y = .64, size=5) +
  geom_segment(aes(x = 1.1, y = .50, xend = 1.9, yend = .50), color = 'grey60') +
  geom_segment(aes(x = 1.9, y = .50, xend = 1.9, yend = .54), color = 'grey60') +
  geom_segment(aes(x = 1.1, y = .46, xend = 1.1, yend = .50), color = 'grey60') +
  annotate('text', label = 'n.s.', x = 1.5, y = .53, size=5)



exp1.model2 <- filter(exp1.mbs) %>% drop_na() %>%
  group_by(subj, fbblock.notrans, runnum) %>%
  summarise(spe = mean(spe),
            gad = mean(gad),
            phq = mean(phq)) %>%
  group_by(runnum) %>%
  mutate(AD.tile = ntile(gad, 3)) %>%
  mutate_at(c('AD.tile'), as.factor) %>%
  mutate(fbblock.notrans = 
           factor(fbblock.notrans, levels = c('Negative','None', 'Positive'))) %>%
  mutate(AD.tile = recode_factor(AD.tile, '1' = 'Low', '3' = 'High')) %>%
  left_join(exp1.fit.data) 
exp1.model2$speZ <- c(scale(exp1.model2$spe))

# exp1.reg <- lmer(spe.fit ~ fbblock.notrans*gad + (1|subj) + (1|runnum), 
#                  filter(exp1.model2, fbblock.notrans %in% c('None', 'Positive')))
# summary(exp1.reg)

exp1.reg <- lmer(speZ ~ fbblock.notrans*gad + (1|subj) + (1|runnum), 
                 filter(exp1.model2, fbblock.notrans %in% c('None', 'Positive')))
summary(exp1.reg)
confint(exp1.reg)
plot(exp1.reg)

m2 <- update(exp1.reg, ~.-fbblock.notrans:gad)
anova(exp1.reg, m2)

set.seed(2020)

#### perform feedback analysis 
exp1.pos.gad.m1.bf <- lmBF(speZ ~ fbblock.notrans*gad + subj + runnum, 
                           whichRandom = c('subj', 'runnum'),
                           iterations = 10000,
                           filter(exp1.model2, fbblock.notrans %in% c('None', 'Positive')))
exp1.pos.gad.m2.bf <- lmBF(speZ ~ fbblock.notrans + gad + subj + runnum, 
                           whichRandom = c('subj', 'runnum'),
                           iterations = 10000,
                           filter(exp1.model2, fbblock.notrans %in% c('None', 'Positive')))

exp1.pos.gad.m1.bf/exp1.pos.gad.m2.bf
exp1.pos.gad.m2.bf/exp1.pos.gad.m1.bf

exp1.reg <- lmer(speZ ~ fbblock.notrans*gad + (1+fbblock.notrans|subj) + (1|runnum), 
                 filter(exp1.model2, fbblock.notrans %in% c('None', 'Negative')))
summary(exp1.reg)
confint(exp1.reg)
plot(exp1.reg)

m2 <- update(exp1.reg, ~.-fbblock.notrans:gad)
anova(exp1.reg, m2)

exp1.neg.gad.m1.bf <- lmBF(speZ ~ fbblock.notrans*gad + subj + runnum + fbblock.notrans:subj, 
                           whichRandom = c('subj', 'runnum'),
                           iterations = 10000,
                           filter(exp1.model2, fbblock.notrans %in% c('None', 'Negative')))
exp1.neg.gad.m2.bf <- lmBF(speZ ~ fbblock.notrans + gad + subj + runnum + fbblock.notrans:subj, 
                           whichRandom = c('subj', 'runnum'),
                           iterations = 10000,
                           filter(exp1.model2, fbblock.notrans %in% c('None', 'Negative')))

exp1.neg.gad.m2.bf/exp1.neg.gad.m1.bf
exp1.neg.gad.m1.bf/exp1.neg.gad.m2.bf


### now do for phq
exp1.reg <- lmer(speZ ~ fbblock.notrans*phq + (1+fbblock.notrans|subj) + (1|runnum), 
                 filter(exp1.model2, fbblock.notrans %in% c('None', 'Positive')))
summary(exp1.reg)
confint(exp1.reg)

plot(exp1.reg)

m2 <- update(exp1.reg, ~.-fbblock.notrans:phq)
anova(exp1.reg, m2)

# BF_BIC = exp((BIC(exp1.reg) - BIC(m2))/2)  # BICs to Bayes factor

exp1.pos.phq.m1.bf <- lmBF(speZ ~ fbblock.notrans*phq + subj + runnum + fbblock.notrans:subj, 
                           whichRandom = c('subj', 'runnum'),
                           iterations = 10000,
                           filter(exp1.model2, fbblock.notrans %in% c('None', 'Positive')))
exp1.pos.phq.m2.bf <- lmBF(speZ ~ fbblock.notrans + phq + subj + runnum + fbblock.notrans:subj, 
                           whichRandom = c('subj', 'runnum'), 
                           iterations = 10000,
                           filter(exp1.model2, fbblock.notrans %in% c('None', 'Positive')))
exp1.pos.phq.m2.bf/exp1.pos.phq.m1.bf

exp1.pos.phq.m1.bf/exp1.pos.phq.m2.bf


exp1.reg <- lmer(speZ ~ fbblock.notrans*phq + (1+fbblock.notrans|subj) + (1|runnum), 
                 filter(exp1.model2, fbblock.notrans %in% c('None', 'Negative')))
summary(exp1.reg)
confint(exp1.reg)
plot(exp1.reg)

m2 <- update(exp1.reg, ~.-fbblock.notrans:phq)
anova(exp1.reg, m2)

set.seed(2020)

exp1.neg.phq.m1.bf <- lmBF(speZ ~ fbblock.notrans*phq + subj + runnum + fbblock.notrans:subj, 
                           whichRandom = c('subj', 'runnum', 'fbblock.notrans:subj'),
                           iterations = 10000,
                           filter(exp1.model2, fbblock.notrans %in% c('None', 'Negative')))
exp1.neg.phq.m2.bf <- lmBF(speZ ~ fbblock.notrans + phq + subj + runnum + fbblock.notrans:subj, 
                           whichRandom = c('subj', 'runnum', 'fbblock.notrans:subj'), 
                           iterations = 10000,
                           filter(exp1.model2, fbblock.notrans %in% c('None', 'Negative')))
exp1.neg.phq.m2.bf/exp1.neg.phq.m1.bf
exp1.neg.phq.m1.bf/exp1.neg.phq.m2.bf


##### Supp Figure 9
## do the same plot for high and low values of gad based on clinical cutoff
exp1.model <- exp1.mbs %>% drop_na() %>% 
  dplyr::select(c(subj, runnum, trialnum, spe, confZ, gad, phq, accu)) %>%
  mutate(confhilo = confZ > 0) %>%
  group_by(subj) %>%
  mutate(conf.tile = ntile(confZ, 4)) %>%
  group_by(runnum) %>%
  mutate(gad.cc = cut(gad, breaks = c(0,5,10,20))) %>%
  mutate(phq.cc = cut(phq, breaks = c(0,5,10,25))) %>%
  mutate_at(c('runnum', 'gad.cc', 'phq.cc'), as.factor)
exp1.model$conftile.hilo <- as.factor(ceiling(exp1.model$conf.tile/2))
exp1.model <- exp1.model %>%
  mutate(conftile.hilo = recode_factor(conftile.hilo, '1' = 'Low', '2' = 'High')) %>%
  mutate(gad.cc = recode_factor(gad.cc, '(0,5]' = 'Minimal', '(10,20]' = 'Moderate-to-Severe'))  %>%
  mutate(phq.cc = recode_factor(phq.cc, '(0,5]' = 'Minimal', '(10,25]' = 'Moderate-to-Severe')) 


exp1.reg <- lmer(spe ~ confZ*gad.cc + accu + (1|subj) + (1|runnum), 
                 filter(exp1.model, gad.cc %in% c('Minimal', 'Moderate-to-Severe'), 
                        confhilo==F))
summary(exp1.reg)

conf.tile.cents <- c(mean(exp1.model$confZ[exp1.model$conf.tile==1]),
                     mean(exp1.model$confZ[exp1.model$conf.tile==2]))

emm.lo <- summary(emmeans(exp1.reg, ~ confZ|gad.cc, 
                              at = list(confZ = conf.tile.cents)))
emm.lo

exp1.reg <- lmer(spe ~ confZ*gad.cc + accu + (1|subj) + (1|runnum), 
                 filter(exp1.model, gad.cc %in% c('Minimal', 'Moderate-to-Severe'), 
                        confhilo==T))
summary(exp1.reg)

conf.tile.cents <- c(mean(exp1.model$confZ[exp1.model$conf.tile==3]),
                     mean(exp1.model$confZ[exp1.model$conf.tile==4]))
emm.hi <- summary(emmeans(exp1.reg, ~ confZ|gad.cc, 
                           at = list(confZ = conf.tile.cents)))
emm.hi

spe <- c(emm.lo$emmean,   emm.hi$emmean)
lci <- c(emm.lo$lower.CL, emm.hi$lower.CL)
uci <- c(emm.lo$upper.CL, emm.hi$upper.CL)

conf.tile <- c(c(1:2,1:2), c(3:4,3:4))
gad.cc <- rep(c(rep('Minimal',2), rep('Moderate-to-Severe',2)),2)
conftile.hilo <- c(rep('Low',4), rep('High',4))
exp1.model.df <- data.frame(spe, lci, uci, 
                            conf.tile, gad.cc, conftile.hilo)
exp1.model.df$conftile.hilo <- as.factor(exp1.model.df$conftile.hilo)

ggplot(filter(exp1.model.df, gad.cc %in% c('Minimal', 'Moderate-to-Severe')), 
       aes(x = conf.tile, y = spe, colour = gad.cc, linetype = conftile.hilo)) +
  scale_color_manual(breaks = c('Minimal', 'Moderate-to-Severe'),
                     values = c('tan3', 'royalblue2')) +
  geom_point(position=position_dodge(.2), size=3.2, stroke=.5) +
  geom_errorbar(aes(ymin = lci, ymax = uci), linetype = 1,
                position=position_dodge(.2), linewidth=.6, width=.3) +
  stat_summary(geom="line", position=position_dodge(.2)) +
  scale_linetype_manual(values = c('solid', 'dashed')) +
  theme_classic2() +
  labs(color = 'GAD-7 score', linetype = 'Confidence level')+
  xlab('Confidence (tiled)') +
  ylab('Global SPE') +
  scale_y_continuous(breaks = c(.5,.6)) +
  theme(text = element_text(size=font_size-6),
        legend.direction = 'vertical',
        # axis.title.y = element_text(hjust=.9),
        legend.text=element_text(size=font_size-10),
        plot.caption = element_text(size = font_size-4),
        legend.position = 'top') +
  geom_segment(aes(x = 3.7, y = .52, xend = 3.7, yend = .58), color = 'grey60') +
  geom_segment(aes(x = 3.6, y = .55, xend = 3.7, yend = .55), color = 'grey60') +
  annotate('text', label = 'p = .016', x = 3.3, y = .55, size=4, angle = 0)



exp1.reg <- lmer(spe ~ confZ*phq.cc + accu + (1|subj) + (1|runnum), 
                 filter(exp1.model, phq.cc %in% c('Minimal', 'Moderate-to-Severe'), 
                        confhilo==F))
summary(exp1.reg)

conf.tile.cents <- c(mean(exp1.model$confZ[exp1.model$conf.tile==1]),
                     mean(exp1.model$confZ[exp1.model$conf.tile==2]))

emm.lo <- summary(emmeans(exp1.reg, ~ confZ|phq.cc, 
                          at = list(confZ = conf.tile.cents)))
emm.lo

exp1.reg <- lmer(spe ~ confZ*phq.cc + accu + (1|subj) + (1|runnum), 
                 filter(exp1.model, phq.cc %in% c('Minimal', 'Moderate-to-Severe'), 
                        confhilo==T))
summary(exp1.reg)

conf.tile.cents <- c(mean(exp1.model$confZ[exp1.model$conf.tile==3]),
                     mean(exp1.model$confZ[exp1.model$conf.tile==4]))
emm.hi <- summary(emmeans(exp1.reg, ~ confZ|phq.cc, 
                          at = list(confZ = conf.tile.cents)))
emm.hi

spe <- c(emm.lo$emmean,   emm.hi$emmean)
lci <- c(emm.lo$lower.CL, emm.hi$lower.CL)
uci <- c(emm.lo$upper.CL, emm.hi$upper.CL)

conf.tile <- c(c(1:2,1:2), c(3:4,3:4))
phq.cc <- rep(c(rep('Minimal',2), rep('Moderate-to-Severe',2)),2)
conftile.hilo <- c(rep('Low',4), rep('High',4))
exp1.model.df <- data.frame(spe, lci, uci, 
                            conf.tile, phq.cc, conftile.hilo)
exp1.model.df$conftile.hilo <- as.factor(exp1.model.df$conftile.hilo)

ggplot(filter(exp1.model.df, phq.cc %in% c('Minimal', 'Moderate-to-Severe')), 
       aes(x = conf.tile, y = spe, colour = phq.cc, linetype = conftile.hilo)) +
  scale_color_manual(breaks = c('Minimal', 'Moderate-to-Severe'),
                     values = c('tan3', 'royalblue2')) +
  geom_point(position=position_dodge(.2), size=3.2, stroke=.5) +
  geom_errorbar(aes(ymin = lci, ymax = uci), linetype = 1,
                position=position_dodge(.2), linewidth=.6, width=.3) +
  stat_summary(geom="line", position=position_dodge(.0)) +
  scale_linetype_manual(values = c('solid', 'dashed')) +
  theme_classic2() +
  labs(color = 'PHQ-9 score', linetype = 'Confidence level')+
  xlab('Confidence (tiled)') +
  ylab('Global SPE') +
  scale_y_continuous(breaks = c(.5,.6)) +
  theme(text = element_text(size=font_size-6),
        legend.direction = 'vertical',
        # axis.title.y = element_text(hjust=.9),
        legend.text=element_text(size=font_size-10),
        plot.caption = element_text(size = font_size-4),
        legend.position = 'top') +
  geom_segment(aes(x = 3.75, y = .55, xend = 3.75, yend = .60), color = 'grey60') +
  geom_segment(aes(x = 3.65, y = .57, xend = 3.75, yend = .57), color = 'grey60') +
  annotate('text', label = 'p = .067', x = 3.35, y = .57, size=4, angle = 0)



exp1.reg <- lmer(confZ ~ spe*gad.cc + accu + (1|subj) + (1|runnum), 
                 filter(exp1.model, gad.cc %in% c('Minimal', 'Moderate-to-Severe'), 
                        confhilo==T))
summary(exp1.reg)
plot(exp1.reg)
plot_model(exp1.reg,type = 'pred', terms= c('spe', 'gad.cc'))
m2 <- update(exp1.reg, ~.-spe:gad.cc)
anova(exp1.reg,m2)


exp1.reg <- lmer(confZ ~ spe*phq.cc + accu + (1|subj) + (1|runnum), 
                 filter(exp1.model, phq.cc %in% c('Minimal', 'Moderate-to-Severe'), 
                        confhilo==T))
summary(exp1.reg)
plot(exp1.reg)
plot_model(exp1.reg,type = 'pred', terms= c('spe', 'phq.cc'))
m2 <- update(exp1.reg, ~.-spe:phq.cc)
anova(exp1.reg,m2)


exp1.reg <- lmer(confZ ~ spe*confhilo*gad + (1+spe*confhilo*gad|subj) +
                   (1|task) + (1|runnum), 
                 control = lmerControl(optimizer = c('bobyqa'), calc.derivs = F),
                 exp1.model)
summary(exp1.reg)

# # reduce the random effects till model converges
# exp1.reg <- lmer(confZ ~ spe*gad + (1+spe*gad|subj), filter(exp1.model, confhilo=='High'))
# summary(exp1.reg)
# plot_model(exp1.reg,type = 'pred', terms= c('spe', 'gad'))

exp1.reg <- lmer(confZ ~ spe*confhilo*gad + (1+spe+confhilo|subj), exp1.model)
summary(exp1.reg)
# plot_model(exp1.reg,type = 'pred', terms= c('spe', 'gad'))
plot_model(exp1.reg,type = 'pred', terms= c('spe', 'gad', 'confhilo'))

m2 <- update(exp1.reg, ~.-spe:confhilo:gad)
anova(exp1.reg, m2)

# # run the reverse model (regression upon SPE) to get partialled values of SPE for plotting 
# exp1.reg <- lmer(speZ ~ confZ*gad + feedback*gad + (1+feedback|subj), exp1.model)
# summary(exp1.reg)
# plot_model(exp1.reg, type = c('pred'), terms = c('confZ', 'gad')) +
#   theme_classic()

# run the reverse model (regression upon SPE) to get partialled values of SPE for plotting 
exp1.reg <- lmer(speZ ~ confZ*confhilo*gad + feedback*gad + (1+feedback|subj), exp1.model)
summary(exp1.reg)
plot(exp1.reg)
plot_model(exp1.reg, type = c('pred'), terms = c('confZ', 'gad', 'confhilo')) +
  theme_classic()

exp1.model$spe.partial.3way <- remef(exp1.reg, keep.intercept = T, ran = 'all',
                                    fix = c('gad:feedbackNegative',
                                            'gad:feedbackPositive'))

##### Figure 4A
ggplot(filter(exp1.model, AD.tile!=2), 
       aes(x = conf.tile, y = spe.partial.3way, colour = AD.tile)) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', 'royalblue2')) +
  stat_summary(geom='line', aes(linetype = conftile.hilo), size=1) +
  stat_summary(fun.data = mean_cl_boot, size=.5, shape = 1) + 
  stat_summary(fun.y = mean, size=.5, alpha = .5) + 
  labs(color = 'AD score', linetype = 'Confidence level')+
  theme_classic2() +
  xlab('Confidence (tiled)') +
  ylab('Global SPE (z-scored)') +
  theme(text = element_text(size=font_size-2), 
        legend.direction = 'vertical',
        legend.text=element_text(size=font_size-6),
        plot.caption = element_text(size = font_size-2),
        axis.title.y = element_text(hjust = .9),
        legend.position = 'none') +
  geom_segment(aes(x = 1.5, y = 0, xend = 3.5, yend = 0), color = 'black') +
  geom_segment(aes(x = 1.5, y = -.05, xend = 1.5, yend = .05), color = 'black') +
  geom_segment(aes(x = 3.5, y = -.05, xend = 3.5, yend = .05), color = 'black') +
  annotate('text', label = 'p < .0001', x = 2.5, y = .06, size=5)

##### Supp Figure 9
## do the same plot for high and low values of gad based on clinical cutoff
ggplot(filter(exp1.model, gad.cc %in% c('(0,5]', '(10,20]')), 
       aes(x = conf.tile, y = spe.partial.3way, colour = gad.cc)) +
  scale_color_manual(breaks = c('(0,5]', '(10,20]'),
                     values = c('tan3', 'royalblue2')) +
  stat_summary(geom='line', aes(linetype = conftile.hilo), size=1) +
  stat_summary(fun.data = mean_cl_boot, size=.5, shape = 1) + 
  stat_summary(fun.y = mean, size=.5, alpha = .5) + 
  labs(color = 'GAD-7 score', linetype = 'Confidence level')+
  theme_classic2() +
  xlab('Confidence (tiled)') +
  ylab('Global SPE (z-scored)') +
  theme(text = element_text(size=font_size-2), 
        legend.direction = 'vertical',
        legend.text=element_text(size=font_size-6),
        plot.caption = element_text(size = font_size-2),
        legend.position = 'top') 

ggplot(filter(exp1.model, phq.cc %in% c('(0,5]', '(10,25]')), 
       aes(x = conf.tile, y = spe.partial.3way, colour = phq.cc)) +
  scale_color_manual(breaks = c('(0,5]', '(10,25]'),
                     values = c('tan3', 'royalblue2')) +
  stat_summary(geom='line', aes(linetype = conftile.hilo), size=1) +
  stat_summary(fun.data = mean_cl_boot, size=.5, shape = 1) + 
  stat_summary(fun.y = mean, size=.5, alpha = .5) + 
  labs(color = 'PHQ-9 score', linetype = 'Confidence level')+
  theme_classic2() +
  xlab('Confidence (tiled)') +
  ylab('Global SPE (z-scored)') +
  theme(text = element_text(size=font_size-2), 
        legend.direction = 'vertical',
        legend.text=element_text(size=font_size-6),
        plot.caption = element_text(size = font_size-2),
        legend.position = 'top')


### also do the same analysis for depression score PHQ (above was for anxiety GAD)
exp1.reg <- lmer(confZ ~ spe*confhilo*phq + (1+spe+confhilo|subj), exp1.model)
summary(exp1.reg)
m2 <- update(exp1.reg, ~.-spe:confhilo:phq)
anova(exp1.reg, m2)


exp1.reg <- lmer(speZ ~ confZ*confhilo*phq + feedback*phq +
                   (1+ feedback|subj), exp1.model)
exp1.model$spe.partial.3way <- remef(exp1.reg, keep.intercept = T, ran = 'all',
                                    fix = c('phq:feedbackNegative',
                                            'phq:feedbackPositive'))

ggplot(filter(exp1.model, phq.tile!=2), 
       aes(x = conf.tile, y = spe.partial.3way, colour = phq.tile)) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', 'royalblue2')) +
  stat_summary(geom='line', aes(linetype = conftile.hilo), size=1) +
  stat_summary(fun.data = mean_cl_boot, size=.5, shape = 1) + 
  stat_summary(fun.y = mean, size=.5, alpha = .5) + 
  labs(color = 'AD score', linetype = 'Confidence level')+
  theme_classic2() +
  xlab('Confidence (tiled)') +
  ylab('Global SPE (z-scored)') +
  theme(text = element_text(size=font_size-2), 
        legend.direction = 'vertical',
        axis.title.y = element_text(hjust=.9),
        legend.text=element_text(size=font_size-6),
        plot.caption = element_text(size = font_size-2),
        legend.position = 'none') +
  geom_segment(aes(x = 1.5, y = .01, xend = 3.5, yend = .01), color = 'grey60') +
  geom_segment(aes(x = 1.5, y = -.03, xend = 1.5, yend = .05), color = 'grey60') +
  geom_segment(aes(x = 3.5, y = -.03, xend = 3.5, yend = .05), color = 'grey60') +
  annotate('text', label = 'p = .0006', x = 2.5, y = .05, size=5)



plot(exp1.reg)

# 
# ### how does phq relate feedback in forming spe(intervention)
# ### old analysis for preregistration
# m.reg = lmer(spe ~ accu*task +
#                conf*task +
#                fbpos*task*phq +
#                fbneg*task*phq +
#                scale_rt*task +
#                stair*task +
#                (1|subj),
#              REML=F,exp1.df)
# summary(m.reg)
# # plot(m.reg)
# 
# m2 <- update(m.reg, ~.-task:fbpos:phq)
# anova(m.reg,m2)
# m.reg <- m2 # accept the reduced model
# summary(m.reg)
# 
# m2 <- update(m.reg, ~.-task:fbneg:phq)
# anova(m.reg,m2)
# m.reg <- m2 # accept the reduced model
# summary(m.reg)
# 
# m2 <- update(m.reg, ~.-task:fbpos)
# anova(m.reg,m2)
# m.reg <- m2 # accept the reduced model
# summary(m.reg)
# 
# m2 <- update(m.reg, ~.-task:stair)
# anova(m.reg,m2)
# m.reg <- m2 # accept the reduced model
# summary(m.reg)
# 
# m2 <- update(m.reg, ~.-task:scale_rt)
# anova(m.reg,m2)
# m.reg <- m2 # accept the reduced model
# summary(m.reg)
# 
# m2 <- update(m.reg, ~.-task:fbneg)
# anova(m.reg,m2)
# m.reg <- m2 # accept the reduced model
# summary(m.reg)
# 
# m2 <- update(m.reg, ~.-task:conf)
# anova(m.reg,m2)
# m.reg <- m2 # accept the reduced model
# summary(m.reg)
# 
# m2 <- update(m.reg, ~.-task:accu)
# anova(m.reg,m2)
# m.reg <- m2 # accept the reduced model
# summary(m.reg)
# 
# m2 <- update(m.reg, ~.-phq:fbpos)
# anova(m.reg,m2)
# m.reg <- m2 # accept the reduced model
# summary(m.reg)
# 
# m.reg <- update(m.reg, REML=T)
# summary(m.reg)
# # plot(m.reg)
# 
# plot_model(m.reg,
#            title="spe (intervention)",
#            show.p=TRUE, show.values = T,
#            value.offset = .4, value.size = 3.8,
#            vline.color = 'grey') + 
#   coord_cartesian(ylim = c(-.9, .9)) +
#   labs(tag = 'A') +
#   theme(text = element_text(size=font_size), 
#         plot.caption = element_text(size = font_size),
#         plot.tag = element_text(size = tag_size, face = "bold")) +
#   theme_pubclean()
# 
# f4a <- plot_model(m.reg, type = ("pred"),
#            terms = c("fbneg", "phq"), 
#            line.size = 1.5, ci.lvl=.95,
#            title = '') +
#   ylab("SPE (intervention)") +
#   xlab("Proportion negative feedback") +
#   labs(colour = 'PHQ-9 score', tag = 'A') +
#   theme_pubclean() +
#   annotate('text', x=.3, y=-.1, label="p = .0053") +
#   theme(text = element_text(size=font_size), 
#         plot.caption = element_text(size = font_size),
#         plot.tag = element_text(size = tag_size, face = "bold"))
# f4a
# 
# emm <- emtrends(m.reg, var = 'fbneg', c('phq'), at=list(phq=c(1,6,12)))
# test(emm)
# contrast(emm)

# 
# # Power estimation for Exp 2
# 
# exp1.df$subject = as.numeric(exp1.df$subj)
# FLPmodel = lmer(spe ~
#                fbneg*phq + 
#                (1|subject),
#              REML=F,exp1.df)
# summary(FLPmodel)
# 
# SESOI <- c(0.011485, -0.632253, 0.001377, -0.023) # specify SESOI
# 
# power_SESOI <- mixedpower(model = FLPmodel, data = exp1.df,
#                           fixed_effects = c("fbneg", "phq"),
#                           simvar = "subject", steps = c(320,330,340),
#                           critical_value = 2, n_sim = 1500,
#                           SESOI = SESOI, databased = F)
# #             310   320   330  mode    effect
# # fbneg     1.000 1.000 1.000 SESOI     fbneg
# # phq       0.242 0.263 0.271 SESOI       phq
# # fbneg:phq 0.886 0.900 0.883 SESOI fbneg:phq
# # 
# # 2000 iterations
# #           300   310   320   330   340  mode    effect
# # fbneg     1.000 1.000 1.000 1.000 1.000 SESOI     fbneg
# # phq       0.230 0.248 0.254 0.309 0.283 SESOI       phq
# # fbneg:phq 0.855 0.897 0.885 0.895 0.916 SESOI fbneg:phq
# 
# # 1500 iterations
#             # 320       330       340  mode    effect
# # fbneg     1.0000000 1.0000000 1.0000000 SESOI     fbneg
# # phq       0.2573333 0.2753333 0.2773333 SESOI       phq
# # fbneg:phq 0.9033333 0.8866667 0.8973333 SESOI fbneg:phq


##########
#####
# For SRET analyses load the file with the less strict inclusion criterion

setwd("~/OneDrive - University of Copenhagen/Projects/Experiments/metaBiasShift/")

## read data file
exp1.mbs3 = read.csv(paste('data/exp1/processed/mbsExp1_noperfexcl.csv', sep=''),
                     header = TRUE, sep = ',')
## recode columns as factors
exp1.mbs3$subj <- as.factor(exp1.mbs3$subj)
exp1.mbs3$task <- as.factor(exp1.mbs3$task)
exp1.mbs3$gender[is.nan(exp1.mbs3$gender)] <- NA
exp1.mbs3$gender <- as.factor(exp1.mbs3$gender)

exp1.mbs3 <- exp1.mbs3 %>%
  mutate_at(vars(contains("fb")), as.factor) %>%
  rename(accu = corr) %>%
  rename(stair = incdec)

exp1.mbs3 <- exp1.mbs3 %>%
  mutate(fbblock = recode_factor(fbblock, "1" = posnegNames[1], 
                                 "2" = posnegNames[2], "0" = "None"),
         task = recode_factor(task, "0" = taskNames[1], "1" = taskNames[2]))

## code highly deviant RTs as NA
max_RT_deviation = 3
rt1_max = c(median(exp1.mbs3$rt1[exp1.mbs3$task==taskNames[1]], na.rm=T)+ 
              max_RT_deviation*IQR(exp1.mbs3$rt1[exp1.mbs3$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(exp1.mbs3$rt1[exp1.mbs3$task==taskNames[2]], na.rm=T)+ 
              max_RT_deviation*IQR(exp1.mbs3$rt1[exp1.mbs3$task==taskNames[2]], na.rm=T ))
exp1.mbs3 <- filter(exp1.mbs3, trialnum > nLags) %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[1] & .>rt1_max[1],NA)) %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[2] & .>rt1_max[2],NA))

## scale RTs 
exp1.mbs3$scale_rt <- exp1.mbs3$rt1
exp1.mbs3$scale_rt[exp1.mbs3$task==taskNames[1]] <- scale(exp1.mbs3$rt1[exp1.mbs3$task==taskNames[1]])
exp1.mbs3$scale_rt[exp1.mbs3$task==taskNames[2]] <- scale(exp1.mbs3$rt1[exp1.mbs3$task==taskNames[2]])

rt2_max = c(median(exp1.mbs3$rt2[exp1.mbs3$task==taskNames[1]], na.rm=T)+ 
              max_RT_deviation*IQR(exp1.mbs3$rt2[exp1.mbs3$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(exp1.mbs3$rt2[exp1.mbs3$task==taskNames[2]], na.rm=T)+ 
              max_RT_deviation*IQR(exp1.mbs3$rt2[exp1.mbs3$task==taskNames[2]], na.rm=T ))
exp1.mbs3 <- filter(exp1.mbs3, trialnum > nLags) %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[1] & .>rt2_max[1],NA)) %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[2] & .>rt2_max[2],NA))

exp1.mbs3$scale_rt2 <- exp1.mbs3$rt2
exp1.mbs3$scale_rt2[exp1.mbs3$task==taskNames[1]] <- scale(exp1.mbs3$rt2[exp1.mbs3$task==taskNames[1]])
exp1.mbs3$scale_rt2[exp1.mbs3$task==taskNames[2]] <- scale(exp1.mbs3$rt2[exp1.mbs3$task==taskNames[2]])

## if to subtract the baseline blocks confidence and spe from the rest of the blocks
# get the baseline values for the two tasks from runs 1 and 2
exp1.mbs3.bas <- filter(exp1.mbs3, runnum<3) %>%
  group_by(subj, task, group) %>%
  dplyr::summarise(
    conf.bas = mean(conf),
    spe.bas = mean(spe),
    accu.bas = mean(accu)
  )

exp1.mbs3 <- exp1.mbs3 %>% left_join(exp1.mbs3.bas)

if (subtractBaselineConfSpe){
  # left join the baseline values an subtract them from confidence and spe
  exp1.mbs3 <- exp1.mbs3 %>%
    mutate(accu_b = accu - accu.bas) %>%
    mutate(conf = conf-conf.bas)%>%
    mutate(spe = spe-spe.bas) 
}

## take away the baseline (first two) blocks from the dataset
exp1.mbs3 <- exp1.mbs3 %>%
  filter(runnum>2) #%>%
# mutate(runnum = runnum-2)
exp1.mbs3 <- exp1.mbs3 %>%
  mutate(blockType = factor(runnum%%2)) %>%
  mutate(blockType = recode_factor(blockType, "0" = "transf", "1" = "interv"))
exp1.mbs3 <- exp1.mbs3 %>%
  mutate_at(vars(contains("runnum")), as.factor)
exp1.mbs3$fbblock <- droplevels(exp1.mbs3$fbblock)

## create new columns for trial-by-trial positive and negative feedback
exp1.mbs3 <- exp1.mbs3 %>%
  mutate(fbpos=feedback, fbneg=feedback) %>%
  mutate_at("fbpos", ~replace(., feedback==1 & accu==1, 1)) %>%
  mutate_at("fbpos", ~replace(., feedback==1 & accu==0, 0)) %>%
  mutate_at("fbneg", ~replace(., feedback==1 & accu==1, 0)) %>%
  mutate_at("fbneg", ~replace(., feedback==1 & accu==0, 1)) 


# give both interv and transf run the same value of fbblock (feedback block type)
exp1.mbs3$fbblock[exp1.mbs3$runnum==2] = exp1.mbs3$fbblock[exp1.mbs3$runnum==1]
exp1.mbs3$fbblock[exp1.mbs3$runnum==4] = exp1.mbs3$fbblock[exp1.mbs3$runnum==3]

## code intervention task by group
exp1.mbs3$intervTask = exp1.mbs3$group
exp1.mbs3$intervTask[exp1.mbs3$group %in% c(1,2,3,4)] = taskNames[1]
exp1.mbs3$intervTask[exp1.mbs3$group %in% c(5,6,7,8)] = taskNames[2]
# 
# transfer task by group (same as task column)
exp1.mbs3$testTask = exp1.mbs3$group
exp1.mbs3$testTask[exp1.mbs3$group %in% c(1,2,7,8)] = taskNames[1]
exp1.mbs3$testTask[exp1.mbs3$group %in% c(3,4,5,6)] = taskNames[2]

## code whether transfer is to same or opposite task
exp1.mbs3$transferType = exp1.mbs3$group
exp1.mbs3$transferType[exp1.mbs3$group %in% c(1,2,5,6)] = transfertaskNames[1]
exp1.mbs3$transferType[exp1.mbs3$group %in% c(3,4,7,8)] = transfertaskNames[2]

## code whether pos or neg is the first block of feedback
exp1.mbs3$firstFeedback = exp1.mbs3$group
exp1.mbs3$firstFeedback[exp1.mbs3$group %in% c(1,3,5,7)] = posnegNames[1]
exp1.mbs3$firstFeedback[exp1.mbs3$group %in% c(2,4,6,8)] = posnegNames[2]

# exp1.mbs3$group <- as.factor(exp1.mbs3$group)
exp1.mbs3 <- exp1.mbs3 %>%
  mutate_at(c("group", "intervTask", "firstFeedback", 'testTask'), as.factor) %>%
  mutate(group = recode_factor(group, "1" = 'Group 1', 
                               "2" = 'Group 2', "3" = 'Group 3', "4" = 'Group 4', 
                               "5" = 'Group 5', "6" = 'Group 6', "7" = 'Group 7', 
                               "8" = 'Group 8')) 

############
############
############
### A6 - baseline correlation with sret
posnegWordColours =  c('springgreen3', 'maroon2')

exp1.df6 <- filter(exp1.mbs3, trialnum==1, runnum==3) %>%
  pivot_longer(cols = starts_with("posword_"),
               names_to = "wordnum",
               names_prefix = "posword_",
               values_to = "word_pos") %>%
  dplyr::select(-starts_with("posword_")) 

exp1.temp <- filter(exp1.mbs3, trialnum==1, runnum==3) %>%
  pivot_longer(cols = starts_with("negword_"),
               names_to = "wordnum",
               names_prefix = "negword_",
               values_to = "word_neg") %>%
  dplyr::select(-starts_with("negword_"))

exp1.df6 <- exp1.df6 %>% left_join(exp1.temp) %>%
  dplyr::select(-starts_with("posword_"))  %>%
  dplyr::select(-starts_with("negword_")) %>%
  pivot_longer(cols = starts_with("word_"),
               names_to = "wordvalence",
               names_prefix = "word_",
               values_to = "word_endorse")

exp1.df6$wordnum = as.numeric(exp1.df6$wordnum)
exp1.df6$wordnum[exp1.df6$wordvalence=='neg'] = exp1.df6$wordnum[exp1.df6$wordvalence=='neg'] + 20
exp1.df6$wordnum = as.factor(exp1.df6$wordnum)
exp1.df6$word_endorse = as.factor(exp1.df6$word_endorse)

exp1.df6 <- exp1.df6 %>% 
  mutate(word_endorse = recode_factor(word_endorse, "0" = "No", "1" = "Yes")) %>%
  mutate(wordvalence = recode_factor(wordvalence, "pos" = "Positive", "neg" = "Negative")) %>%
  group_by(word_endorse, wordvalence, gender, group, wordnum, subj) %>%
  summarise(
    phq = mean(phq),
    gad = mean(gad),
    spin = mean(spin),
    scale_rt = mean(scale_rt, na.rm = T),
    accu = mean(accu),
    conf = mean(conf),
    spe.bas = mean(spe.bas),
    conf.bas = mean(conf.bas),
    age = mean(age)
  )


m.reg <- glmer(word_endorse ~ phq*wordvalence + 
               ( 1|wordnum), exp1.df6,
             family = binomial)
summary(m.reg)
m2 <- update(m.reg, ~.-wordvalence:phq)
lrtest(m.reg,m2)

##################################
##### Supp Figure 10 - interaction of SRET word valence with MHQ and Confidence

sf10a <- plot_model(m.reg, type = ("pred"),
           terms = c("phq [all]", "wordvalence"), 
           line.size = 1.5, ci.lvl = .95,
           title = '', colors = posnegWordColours) +
  theme_pubclean() +
  ylab("Self-endorsement") +
  xlab("PHQ-9") +
  labs(colour = 'Word valence') +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        legend.position="none") +
  guides(color = guide_legend(nrow = 1))
sf10a

m.reg <- glmer(word_endorse ~ 
                 gad*wordvalence +
                 (1|wordnum), exp1.df6,
               family = binomial)
summary(m.reg)

m2 <- update(m.reg, ~.-wordvalence:gad)
lrtest(m.reg,m2)
sf10b <- plot_model(m.reg, type = ("pred"),
           terms = c("gad [all]", "wordvalence"), 
           line.size = 1.5,ci.lvl = .95,
           title = '', colors = posnegWordColours) +
  theme_pubclean()+
  ylab("Self-endorsement") +
  xlab("GAD-7") +
  labs(colour = 'Word valence') +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        legend.position="none") 
sf10b

m.reg <- glmer(word_endorse ~ 
                 spin*wordvalence +
                 (1|wordnum), exp1.df6,
               family = binomial)
summary(m.reg)

m2 <- update(m.reg, ~.-wordvalence:spin)
lrtest(m.reg,m2)

sf10c <- plot_model(m.reg, type = ("pred"),
           terms = c("spin [all]", "wordvalence"), 
           line.size = 1.5,ci.lvl = .95,
           title = '', colors = posnegWordColours)  +
  theme_pubclean()+
  coord_cartesian(xlim = c(0,12))+
  ylab("Self-endorsement") +
  xlab("mini-SPIN") +
  labs(colour = 'Word valence') +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        legend.position="none") 
sf10c


m.reg <- glmer(word_endorse ~ spe.bas*wordvalence + conf.bas*wordvalence + 
                 (1|wordnum), 
               exp1.df6, control=glmerControl(optimizer="bobyqa"),
               family = binomial)
summary(m.reg)
confint(m.reg)
m2 <- update(m.reg, ~.-wordvalence:spe.bas)
lrtest(m.reg,m2)

m2 <- update(m.reg, ~.-wordvalence:conf.bas)
lrtest(m.reg,m2)

sf10d<-plot_model(m.reg, type = ("pred"),
           terms = c("spe.bas [all]", "wordvalence"), 
           line.size = 1.5, ci.lvl = .95,
           title = '', colors = posnegWordColours)  +
  theme_pubclean()+
  ylab("Self-endorsement") +
  xlab("SPE (baseline)") +
  labs(colour = 'Word valence') +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        legend.position = 'none')
sf10d
emm <- emtrends(m.reg, ~wordvalence, "spe.bas")
test(emm)

sf10e<-plot_model(m.reg, type = ("pred"),
           terms = c("conf.bas [all]", "wordvalence"), 
           line.size = 1.5, ci.lvl = .95,
           title = '', colors = posnegWordColours) +
  theme_pubclean()+
  ylab("Self-endorsement") +
  xlab("Confidence (baseline)") +
  labs(colour = 'Word valence') +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        legend.position = 'none')
sf10e
emm <- emtrends(m.reg, ~wordvalence, "conf.bas")
test(emm)


############
####### 
############
### A7 - does meta-sensitivity moderate the effect of feedback
# update this folder with the folder location with data
data.folder <- "/Users/skatyal/OneDrive - University College London/Projects/Experiments/metaBiasShift/data/exp1/processed/"

expNum <- 1
setwd(data.folder)
mbsDataExp1 = readMat(paste('mbsDataExp',expNum,'.mat',sep='')) # load questionnaire data
mbsDataExp1 <- mbsDataExp1[[paste('mbsDataExp',expNum,sep='')]]

mratio <- mbsDataExp1[,,1]$mratio
nsubj <- length(unique(exp1.df.2$subj))
subj <- array(1:nsubj, c(nsubj, 2))
task <- t(array(c("Perception", "Memory"), c(2, nsubj)))

temp <- data.frame(c(mratio), c(subj), c(task))
colnames(temp) <- c('mratio', 'subj', 'task')
temp <- temp %>% mutate_at('subj', as.factor)

# exp1.df.2 <- exp1.df.2 %>% left_join(temp)
exp1.df <- exp1.df %>% left_join(temp)

exp1.df$fbblock[exp1.df$runnum==4|exp1.df$runnum==6] <- 'None'

m.reg <- lmer(spe_b ~ fbblock*task*mratio + accu_b + stair_b + (1|subj), 
              filter(exp1.df, mratio>0, runnum %in% c(3:6)))
summary(m.reg)
plot(m.reg)
m2 <- update(m.reg, ~.-fbblock:task:mratio)
anova(m.reg,m2)
m.reg <- m2 # accept the reduced model
summary(m.reg)
m2 <- update(m.reg, ~.-fbblock:mratio)
anova(m.reg,m2)


m.reg <- lmer(mratio ~ phq + age + gender + (1|task) , 
              filter(exp1.df, runnum %in% c(1:2)))
summary(m.reg)
m2 <- update(m.reg, ~.-phq)
anova(m.reg,m2)

m.reg <- lmer(mratio ~ gad + age + gender + (1|task) , 
              filter(exp1.df, runnum %in% c(1:2)))
summary(m.reg)
m2 <- update(m.reg, ~.-gad)
anova(m.reg,m2)




############
### A8 - histogram of mhq
exp1.mhq <- exp1.df.2 %>% group_by(subj) %>%
  summarise(phq = mean(phq),
            gad = mean(gad))

exp1.mhq$phq.cutoff <- exp1.mhq$phq
exp1.mhq$phq.cutoff <- 'Minimal'
exp1.mhq$phq.cutoff[exp1.mhq$phq>4] <- 'Mild'
exp1.mhq$phq.cutoff[exp1.mhq$phq>9] <- 'Moderate'
exp1.mhq$phq.cutoff[exp1.mhq$phq>14] <- 'Moderately\nsevere'
exp1.mhq$phq.cutoff[exp1.mhq$phq>19] <- 'Severe'

exp1.mhq$gad.cutoff <- exp1.mhq$phq
exp1.mhq$gad.cutoff <- 'Minimal'
exp1.mhq$gad.cutoff[exp1.mhq$gad>4] <- 'Mild'
exp1.mhq$gad.cutoff[exp1.mhq$gad>9] <- 'Moderate'
exp1.mhq$gad.cutoff[exp1.mhq$gad>14] <- 'Severe'

exp1.mhq <- exp1.mhq %>% 
  mutate_at(c('phq.cutoff', 'gad.cutoff'), as.factor) %>%
  mutate(phq.cutoff = relevel(phq.cutoff, 'Minimal')) %>%
  mutate(gad.cutoff = relevel(gad.cutoff, 'Minimal'))

ggplot(exp1.mhq, aes(x = phq, fill = phq.cutoff)) +
  geom_histogram(bins = 20, colour = 'grey40') +
  scale_fill_brewer(palette = "RdYlGn", direction = -1)+
  geom_vline(xintercept = median(exp1.mhq$phq), color = 'black', size=1) +
  xlab('PHQ-9 score') +
  ylab('Number of participants') +
  labs(fill = 'Depression severity') +
  theme_classic() +
  theme(text = element_text(size=16), 
        plot.caption = element_text(size = 16),
        plot.tag = element_text(size = tag_size, face = "bold"), 
        legend.position = 'top',
        legend.direction = 'horizontal') + 
  guides(fill = guide_legend(title.position = "top", 
                             hjust = 0.5,
                             title.hjust = 0.5)) 


ggplot(exp1.mhq, aes(x = gad, fill = gad.cutoff)) + 
  geom_histogram(bins = 21, colour = 'grey40') +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  geom_vline(xintercept = median(exp1.mhq$gad), color = 'black', size=1) +
  xlab('GAD-7 score') +
  ylab('Number of participants') +
  labs(fill = 'Anxiety severity') +
  theme_classic() +
  theme(text = element_text(size=16), 
        plot.caption = element_text(size = 16),
        plot.tag = element_text(size = tag_size, face = "bold"), 
        legend.position = 'top') + 
  guides(fill = guide_legend(title.position = "top", 
                              hjust = 0.5,
                              title.hjust = 0.5)) 


sum(exp1.mhq$phq>9) /length(exp1.mhq$phq)
sum(exp1.mhq$gad>9) /length(exp1.mhq$gad)
############
### A8 - plots of debriefing questions
exp1.db <- exp1.df.2 %>% group_by(subj, awarepos, awareneg) %>%
  summarise(phq = mean(phq),
            gad = mean(gad)) %>%
  pivot_longer(    cols = starts_with("aware"),
                   names_to = "posneg",
                   names_prefix = "aware",
                   values_to = "aware.fb") %>%
  mutate(posneg = recode_factor(posneg, 'pos' = 'Correct', 'neg' = 'Incorrect'))

sum((exp1.df.2 %>% group_by(subj, awareneg))$awareneg=='Yes', na.rm=T)/
  (sum((exp1.df.2 %>% group_by(subj, awareneg))$awareneg=='No', na.rm=T)+
     sum((exp1.df.2 %>% group_by(subj, awareneg))$awareneg=='Yes', na.rm=T))
sum((exp1.df.2 %>% group_by(subj, awarepos))$awarepos=='Yes', na.rm=T)/
  (sum((exp1.df.2 %>% group_by(subj, awarepos))$awarepos=='No', na.rm=T)+
     sum((exp1.df.2 %>% group_by(subj, awarepos))$awarepos=='Yes', na.rm=T))

ggplot(exp1.db, aes(x = aware.fb, fill = posneg)) +
  geom_bar(position = position_dodge()) +
  xlab('Did you notice feedback bias during block?') +
  ylab('Number of participants') +
  theme_classic() +
  ggtitle('Exp 1') +
  labs(fill = 'Feedback block') +
  theme(legend.position = 'top',
        text = element_text(size=16))

t.test(filter(exp1.db, !is.na(aware.fb), posneg=='Correct', aware.fb=='Yes')$phq,
       filter(exp1.db, !is.na(aware.fb), posneg=='Correct', aware.fb=='No')$phq)
t.test(filter(exp1.db, !is.na(aware.fb), posneg=='Incorrect', aware.fb=='Yes')$phq,
       filter(exp1.db, !is.na(aware.fb), posneg=='Incorrect', aware.fb=='No')$phq)



ggplot(filter(exp1.db, !is.na(aware.fb)), aes(x = aware.fb, y = phq, colour = posneg))+
  stat_summary()


exp1.db <- exp1.df.2 %>% group_by(subj, affectpos, affectneg) %>%
  summarise(phq = mean(phq),
            gad = mean(gad)) %>%
  pivot_longer(    cols = starts_with("affect"),
                   names_to = "posneg",
                   names_prefix = "affect",
                   values_to = "affect.fb") %>%
  mutate(posneg = recode_factor(posneg, 'pos' = 'Correct', 'neg' = 'Incorrect'))

ggplot(exp1.db, aes(x = affect.fb, fill = posneg)) +
  geom_bar(position = position_dodge()) +
  xlab('Did feedback on trial change how you felt?') +
  ylab('Number of participants') +
  theme_classic() +
  ggtitle('Exp 1') +
  labs(fill = 'Feedback on trial')+
  theme(legend.position = 'top',
        text = element_text(size=16))

levels(exp1.db$affect.fb) <- c(levels(exp1.db$affect.fb), 
                               'Not felt better', "Not felt worse")
exp1.db <- exp1.db %>% 
  mutate(affect.better = affect.fb) %>%
  mutate_at('affect.better', ~replace(., .!="Felt better", "Not felt better")) %>% 
  mutate(affect.worse = affect.fb) %>%
  mutate_at('affect.worse', ~replace(., .!="Felt worse", "Not felt worse"))

length(filter(exp1.db, affect.better=="Felt better", posneg == 'Correct')$phq)/
  (length(filter(exp1.db, affect.better=="Felt better", posneg == 'Correct')$phq)+
     length(filter(exp1.db, affect.better=="Not felt better", posneg == 'Correct')$phq))
length(filter(exp1.db, affect.worse=="Felt worse", posneg == 'Incorrect')$phq)/
  (length(filter(exp1.db, affect.worse=="Felt worse", posneg == 'Incorrect')$phq)+
     length(filter(exp1.db, affect.worse=="Not felt worse", posneg == 'Incorrect')$phq))


wilcox.test(filter(exp1.db, affect.better=="Felt better", posneg == 'Correct')$phq,
       filter(exp1.db, affect.better=="Not felt better", posneg == 'Correct')$phq)
wilcox.test(filter(exp1.db, affect.worse=="Felt worse", posneg == 'Incorrect')$phq,
       filter(exp1.db, affect.worse=="Not felt worse", posneg == 'Incorrect')$phq)

wilcox.test(filter(exp1.db, affect.better=="Felt better", posneg == 'Correct')$gad,
       filter(exp1.db, affect.better=="Not felt better", posneg == 'Correct')$gad)
wilcox.test(filter(exp1.db, affect.worse=="Felt worse", posneg == 'Incorrect')$gad,
       filter(exp1.db, affect.worse=="Not felt worse", posneg == 'Incorrect')$gad)


ggplot(filter(exp1.db, posneg=='Incorrect'), aes(x = affect.worse, y = gad)) +
  geom_quasirandom(dodge.width=.7, alpha = .7) +
  geom_violin(alpha=.6, position = position_dodge(.7)) +
  stat_summary(fun.data = mean_cl_boot, color = 'black', size = 1, width=.2,
               geom = 'errorbar') +
  xlab('Feeling after incorrect trials') +
  ylab('GAD-7 score') +
  theme_classic() +
  ggtitle('Exp 1') +
  # coord_cartesian(ylim = c(3,8)) +
  theme(legend.position = 'top',
        text = element_text(size=16)) +
  geom_signif(y_position = c(7.5), vjust = -.25, hjust=.5,
              xmin = c(1), color = 'black',
              xmax = c(2),
              annotation = c("p = .015"), tip_length = .01,
              textsize = pval_size, size = .8)


###########################
##########################
# baseline correlation of SPE and SRET to split the words into 2 sets for Exp 2

library(lme4) #Linear mixed effect model
library(emmeans) # Least squares means
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm) # for geom_quasirandom()
library(ggpp)
library(sjPlot)
library(see)
library(lmerTest)

font_size <- 20
tag_size <- 22
nLags <- 0
taskNames = c('Perception', 'Memory')
transfertaskNames = c('Same', 'Opposite')
posnegNames = c('Positive', 'Negative')
posnegColours = c("chartreuse3", "coral3")

if (Sys.info()["sysname"]=='Darwin'){
  setwd("~/OneDrive - University College London/Projects/Experiments/metaBiasShift/")
} else{
  setwd("C:/Users/skatyal/OneDrive - University College London/Projects/Experiments/metaBiasShift/")
}

## read data file
exp1.mbs.bas = read.csv(paste('data/exp1/mbsExp1_allsubj.csv', sep=''),
                   header = TRUE, sep = ',')
## recode columns as factors
exp1.mbs.bas$subj <- as.factor(exp1.mbs.bas$subj)
exp1.mbs.bas$task <- as.factor(exp1.mbs.bas$task)
exp1.mbs.bas$gender[is.na(exp1.mbs.bas$gender)] <- NA
exp1.mbs.bas$gender <- as.factor(exp1.mbs.bas$gender)
exp1.mbs.bas$endorsedif <- exp1.mbs.bas$endorsepos - exp1.mbs.bas$endorseneg

exp1.mbs.bas <- exp1.mbs.bas %>%
  mutate_at(vars(contains("fb")), as.factor) %>%
  dplyr::rename(accu = corr) %>%
  dplyr::rename(stair = incdec)

exp1.mbs.bas <- exp1.mbs.bas %>%
  mutate(fbblock = recode_factor(fbblock, "1" = "pos", "2" = "neg", "0" = "neut"),
         task = recode_factor(task, "0" = taskNames[1], "1" = taskNames[2]),
         gender = recode_factor(gender, '1' = 'female', '2' = 'male')
  )

exp1.mbs.bas$age <- exp1.mbs.bas$age/10

max_RT_deviation = 3
rt1_max = c(median(exp1.mbs.bas$rt1[exp1.mbs.bas$task==taskNames[1]], na.rm=T)+ 
              max_RT_deviation*IQR(exp1.mbs.bas$rt1[exp1.mbs.bas$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(exp1.mbs.bas$rt1[exp1.mbs.bas$task==taskNames[2]], na.rm=T)+ 
              max_RT_deviation*IQR(exp1.mbs.bas$rt1[exp1.mbs.bas$task==taskNames[2]], na.rm=T ))
exp1.mbs.bas <- filter(exp1.mbs.bas, trialnum > nLags) %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[1] & .>rt1_max[1],NA)) %>%
  mutate_at(vars(contains("rt1")), ~replace(.,task==taskNames[2] & .>rt1_max[2],NA))

## scale RTs 
exp1.mbs.bas$scale_rt <- exp1.mbs.bas$rt1
exp1.mbs.bas$scale_rt[exp1.mbs.bas$task==taskNames[1]] <- scale(exp1.mbs.bas$rt1[exp1.mbs.bas$task==taskNames[1]])
exp1.mbs.bas$scale_rt[exp1.mbs.bas$task==taskNames[2]] <- scale(exp1.mbs.bas$rt1[exp1.mbs.bas$task==taskNames[2]])

rt2_max = c(median(exp1.mbs.bas$rt2[exp1.mbs.bas$task==taskNames[1]], na.rm=T)+ 
              max_RT_deviation*IQR(exp1.mbs.bas$rt2[exp1.mbs.bas$task==taskNames[1]], na.rm=T ), # max rt for perception task
            median(exp1.mbs.bas$rt2[exp1.mbs.bas$task==taskNames[2]], na.rm=T)+ 
              max_RT_deviation*IQR(exp1.mbs.bas$rt2[exp1.mbs.bas$task==taskNames[2]], na.rm=T ))
exp1.mbs.bas <- filter(exp1.mbs.bas, trialnum > nLags) %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[1] & .>rt2_max[1],NA)) %>%
  mutate_at(vars(contains("rt2")), ~replace(.,task==taskNames[2] & .>rt2_max[2],NA))

# exp1.mbs.bas <- exp1.mbs.bas %>% drop_na(scale_rt, scale_rt2)
exp1.mbs.bas <- exp1.mbs.bas %>% drop_na(scale_rt)

######################################
#### regress each word on baseline spe
###
##
#

exp1.mbs.bas.df <- filter(exp1.mbs.bas, runnum<3) %>%
  group_by(subj, group, task, gender, runnum) %>%
  dplyr::summarise_at(
    vars(conf:scale_rt),
    .funs = c('mean'), na.rm = T
  ) 

regpos = array(dim = 20)
for (i in c(1:20)){
  chari = as.character(i)
  m.reg = lmer(paste(
    "spe ~ posword_",chari, "+accu + stair +gender + age + (1|task)", 
    sep = ""
  ), exp1.mbs.bas.df)
  regpos[i] <- summary(m.reg)$coefficient[2,4]
}
mean(regpos)

regneg = array(dim = 20)
for (i in c(1:20)){
  chari = as.character(i)
  m.reg = lmer(paste(
    "spe ~ negword_",chari, "+accu + stair +gender + age + (1|task)", 
    sep = ""
  ), exp1.mbs.bas.df)
  regneg[i] <- summary(m.reg)$coefficient[2,4]
}
mean(regneg)

regposconf = array(dim = 20)
for (i in c(1:20)){
  chari = as.character(i)
  m.reg = lmer(paste(
    "conf ~ posword_",chari, "+accu + stair +gender + age + (1|task)",
    sep = ""
  ), exp1.mbs.bas.df)
  regposconf[i] <- summary(m.reg)$coefficient[2,4]
}
mean(regposconf)

regnegconf = array(dim = 20)
for (i in c(1:20)){
  chari = as.character(i)
  m.reg = lmer(paste(
    "conf ~ negword_",chari, "+accu + stair +gender + age + (1|task)",
    sep = ""
  ), exp1.mbs.bas.df)
  regnegconf[i] <- summary(m.reg)$coefficient[2,4]
}
mean(regnegconf)

regposphq = array(dim = 20)
for (i in c(1:20)){
  chari = as.character(i)
  m.reg = lm(paste(
    "phq ~ posword_",chari, " +gender + age",
    sep = ""
  ), exp1.mbs.bas.df)
  regposphq[i] <- summary(m.reg)$coefficient[2,3]
}
mean(regposphq)

regnegphq = array(dim = 20)
for (i in c(1:20)){
  chari = as.character(i)
  m.reg = lm(paste(
    "phq ~ negword_",chari, " +gender + age",
    sep = ""
  ), exp1.mbs.bas.df)
  regnegphq[i] <- summary(m.reg)$coefficient[2,3]
}
mean(regnegphq)

regposgad = array(dim = 20)
for (i in c(1:20)){
  chari = as.character(i)
  m.reg = lm(paste(
    "gad ~ posword_",chari, " +gender + age",
    sep = ""
  ), exp1.mbs.bas.df)
  regposgad[i] <- summary(m.reg)$coefficient[2,3]
}
mean(regposgad)

regneggad = array(dim = 20)
for (i in c(1:20)){
  chari = as.character(i)
  m.reg = lm(paste(
    "gad ~ negword_",chari, " +gender + age",
    sep = ""
  ), exp1.mbs.bas.df)
  regneggad[i] <- summary(m.reg)$coefficient[2,3]
}
mean(regneggad)


library(wordcloud)
library(treemap)
library(hrbrthemes) # theme ispum

## read data file
sretWords = read.csv(paste('data/sretWords.txt', sep=''),
                   header = F, sep = ',')
sretWords$tValSPE = NA
sretWords$tValConf = NA
sretWords$tValPhq = NA
sretWords$tValGad = NA
sretWords$valence = NA
sretWords$tValSPE[1:20] = regpos
sretWords$tValSPE[21:40] = regneg
sretWords$tValConf[1:20] = regposconf
sretWords$tValConf[21:40] = regnegconf
sretWords$tValPhq[1:20] = regposphq
sretWords$tValPhq[21:40] = regnegphq
sretWords$tValGad[1:20] = regposgad
sretWords$tValGad[21:40] = regneggad
sretWords$valence[1:20] = 'pos'
sretWords$valence[21:40] = 'neg'
sretWords$valence[41:42] = 'neut'
sretWords <- sretWords %>% rename(words = V1)

write.csv(sretWords, "data/wordTValues.csv")

##################################
##### Supp Figure 2. Regression of individual SRET words on baseline SPE

filter(sretWords, valence!='neut') %>%
  arrange(tValSPE) %>%
  mutate(words=factor(words, words)) %>%
  ggplot( aes(x=words, y=tValSPE) ) +
  theme_minimal()+
  geom_segment( aes(x=words ,xend=words, y=0, yend=tValSPE), colour="grey") +
  geom_point(size=3, aes(colour = valence))+#, color=c("#69b3a2", "#e42143")) +
  coord_flip()  +
  ylab("t value on SPE") +
  xlab('Words') +
  labs(colour = 'valence')+
  theme(text = element_text(size=16), 
        plot.caption = element_text(size = 18),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(r=15),
    # axis.text=element_text(size=16),
    axis.title=element_text(size=18),
    legend.position="none"
  ) 


## test if colour of stimulus makes any difference to confidence
exp1.mbs = read.csv(paste('data/exp1/processed/mbsExp1_stim.csv', sep=''), header = TRUE, sep = ',')
summary(lmer(conf ~ stim + (1+stim|subj) + (1+stim|runnum) + (1+stim|trialnum),
             filter(exp1.mbs, task==0)))

summary(lmer(phq ~ stim + (1+stim|subj) + (1+stim|runnum) + (1+stim|trialnum),
             filter(exp1.mbs, task==0)))
