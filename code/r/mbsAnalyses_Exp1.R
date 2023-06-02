####################################################################
####################################################################
################# Analysis and Figure generation for Exp 1 #################
################# Katyal et al ##################################
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
# library(mixedpower)

## mediation analysis on mixed models does not work with lmerTest package
## so use only one or the other 
doMediation = F
if (doMediation){
  library(mediation) # mediation analysis
} else {
  library(lmerTest) # anova and summary with p-values
  
}

## set the base directory here which contains a the 'analysis' and 'data' folders
if (Sys.info()["sysname"]=='Darwin'){
  setwd("~/OneDrive - University College London/Projects/Experiments/metaBiasShift/")
} else{
  setwd("C:/Users/skatyal/OneDrive - University College London/Projects/Experiments/metaBiasShift/")
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


## read data file and preparing the data
exp1.mbs = read.csv(paste('data/exp1/mbsExp1.csv', sep=''),
               header = TRUE, sep = ',')
## recode columns as factors
exp1.mbs$subj <- as.factor(exp1.mbs$subj)
exp1.mbs$task <- as.factor(exp1.mbs$task)
exp1.mbs$awarepos[is.nan(exp1.mbs$awarepos)] <- NA
exp1.mbs$awareneg[is.nan(exp1.mbs$awareneg)] <- NA
exp1.mbs$awarepos <- as.factor(exp1.mbs$awarepos)
exp1.mbs$awareneg <- as.factor(exp1.mbs$awareneg)
exp1.mbs$affectpos[is.nan(exp1.mbs$affectpos)] <- NA
exp1.mbs$affectneg[is.nan(exp1.mbs$affectneg)] <- NA
exp1.mbs$affectpos[exp1.mbs$affectpos!=1] <- 0
exp1.mbs$affectneg[exp1.mbs$affectneg!=2] <- 0
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
         awarepos = recode_factor(awarepos, '1' = 'Yes', '2'='No'),
         awareneg = recode_factor(awareneg, '1' = 'Yes', '2'='No'),
         affectpos = recode_factor(affectpos, '1' = 'Yes', '0'='No'),
         affectneg = recode_factor(affectneg, '2' = 'Yes', '0'='No')
  )

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

exp1.mbs1 <- exp1.mbs

## to subtract the baseline blocks confidence and spe from the rest of the blocks
if (subtractBaselineConfSpe){
  # get the baseline values for the two tasks from runs 1 and 2
  exp1.mbs2 <- filter(exp1.mbs, runnum<3) %>%
    group_by(subj, task, group) %>%
    dplyr::summarise(
      conf.bas = mean(conf),
      spe.bas = mean(spe),
      accu.bas = mean(accu)
    )
  
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
      age = mean(age)
    )
  
  # left join the baseline values an subtract them from confidence and spe
  exp1.mbs <- exp1.mbs %>% left_join(exp1.mbs2) %>%
    mutate(spe_ut = spe) %>% #save untransfromed spe to plot
    mutate(conf = conf-conf.bas)%>%
    mutate(spe = spe-spe.bas) %>%
    mutate(accu_b = accu-accu.bas)
}

## take away the baseline (first two) blocks from the dataset
exp1.mbs <- exp1.mbs %>%
  filter(runnum>2) #%>%
  # mutate(runnum = runnum-2)
exp1.mbs <- exp1.mbs %>%
  mutate(blockType = factor(runnum%%2)) %>%
  mutate(blockType = recode_factor(blockType, "0" = "transf", "1" = "interv"))
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
  dplyr::select(-c(feedback, rt2))

# give both interv and transf run the same value of fbblock (feedback block type)
exp1.mbs$fbblock[exp1.mbs$runnum==2] = exp1.mbs$fbblock[exp1.mbs$runnum==1]
exp1.mbs$fbblock[exp1.mbs$runnum==4] = exp1.mbs$fbblock[exp1.mbs$runnum==3]

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

exp1.mbs3 <- exp1.mbs

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
    accu = mean(accu, na.rm=T),
    accu_b = mean(accu_b),
    scale_rt = mean(scale_rt, na.rm=T),
    rt1 = mean(rt1, na.rm=T),
    fbneg = mean(fbneg),
    fbpos = mean(fbpos),
    spe = mean(spe),
    spe_ut = mean(spe_ut),
    stair = mean(stair),
    phq = mean(phq),
    gad = mean(gad),
    spin = mean(spin),
    conf.bas = mean(conf.bas),
    spe.bas = mean(spe.bas),
    age = mean(age)
  )

##### data prep done #####################
##########################################

#################  #################
####### Figure 2
###### effect of feedback on intervention SPEs, accuracy, staircase

exp1.df.2 <- filter(exp1.df, blockType=="interv")

f2a.left <- ggplot(exp1.df.2, aes(x = fbblock, y = spe, color = fbblock)) +
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
  ylab("Self-performance estimate\n(global SPE)") +
  xlab("Feedback") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  labs(color = 'Feedback') + 
  theme(text = element_text(size=font_size), 
        legend.position = 'none',
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  coord_cartesian(ylim = c(-.65,.65))  +
  geom_signif(y_position = c(.56), vjust = -.3, hjust=.4,
              xmin = c(1.1), color = 'black',
              xmax = c(1.9),
              annotation = c("p < .0001"), tip_length = .04,
              textsize = pval_size, size = .8)
f2a.left


f2a.right <- ggplot(exp1.df.2, aes(x = fbblock, y = accu_b, color = fbblock)) +
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
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  labs(color = 'Feedback') + 
  theme(text = element_text(size=font_size), 
        legend.position = 'none',
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold")) +
  # coord_cartesian(ylim = c(-.65,.65))  +
  geom_signif(y_position = c(.2), vjust = -.3,
              xmin = c(1.1), color = 'black',
              xmax = c(1.9),
              annotation = c("n.s."), tip_length = .04,
              textsize = pval_size, size = .8)
f2a.right


m.reg <- lmer(accu ~ fbblock*task + (1|subj),
              exp1.df.2)
summary(m.reg)
m.reg <- lmer(accu ~ fbblock + (1|subj),
              exp1.df.2)
summary(m.reg)
m.reg <- lmer(stair ~ fbblock*task + (1|subj),
              exp1.df.2)
summary(m.reg)
m.reg <- lmer(stair ~ fbblock + (1|subj),
              exp1.df.2)
summary(m.reg)


##################################
####### Supp Figure 5 (upper panel)

sf5a.upper <- ggplot(exp1.df.2, aes(x = task, y = stair, color = fbblock)) +
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
  ggtitle('Exp 1') +
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

f4a <- ggplot(exp1.df.2, aes(x = trialnum, y = conf, colour = fbblock)) +
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
  filter(blockType=="transf") %>% 
  group_by(intervTask, testTask, fbblock, subj) %>%
  summarise(conf = mean(conf)) %>%
  pivot_wider(
    names_from = fbblock,
    values_from = conf
  ) %>% group_by(subj, intervTask, testTask) %>%
  summarise(
    confPosNeg = mean(Positive) - mean(Negative)
  ) %>% mutate(intervTask = relevel(intervTask, "Perception")) %>%
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
  coord_cartesian(ylim = c(-.25,.48)) +
  ylab("Confidence (test)\nPositive - Negative") +
  xlab("Intervention task") +
  labs(color = 'Test task') +# ggtitle('Exp 1') +
  theme_pubclean() +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
geom_signif(y_position = c(.4, .26, .26, .26), vjust = -.3,
            xmin = c(.7, 1.1, 1.7,2.1),
            xmax = c(.9, 1.3, 1.9,2.3),
            annotation = c("p < .0001", "p = .01",'p = .56','p = .004'), 
            tip_length = .0, color = 'black',
            textsize = pval_size, size = .8)
f4b


##################################
##### Supp Figure 14a. transfer of feedback to test block SPE

exp1.df.2 <- exp1.mbs %>%
  filter(blockType=="transf") %>%
  group_by(fbblock, subj, intervTask, testTask) %>%
  summarise(spe = mean(spe)) %>%
  pivot_wider(
    names_from = fbblock,
    values_from = spe
  ) %>% group_by(subj, intervTask, testTask) %>%
  summarise(
    spe = mean(Positive) - mean(Negative)
  ) %>%
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


###########

############
###############################################
### SPE (intervention) analyses
################################################

exp1.df.2 <- filter(exp1.df, blockType=="interv")

## factorial effect of fb on spe

m.reg <- lmer(spe ~ fbblock*task*firstFeedback + 
                accu + conf + stair + scale_rt +
                (1|subj) + (1|group) + (1|runnum), exp1.df.2)
summary(m.reg)
m.reg <- lmer(spe ~ fbblock*task*firstFeedback + 
                accu + conf + stair + scale_rt +
                (1|subj) + (1|group), REML = F, exp1.df.2)
summary(m.reg)
m.reg <- lmer(spe ~ fbblock*task*firstFeedback + 
                accu + conf + stair + scale_rt +
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

emmeans(m.reg, pairwise ~ fbblock)$contrasts

emm <- emmeans(m.reg, specs = c("fbblock"))
test(emm)

m.reg <- lmer(spe ~ fbblock + 
                accu + conf + stair + scale_rt +
                (1|subj), REML = F, exp1.df.2)
summary(m.reg)
plot(m.reg)

m2 <- update(m.reg, ~.-fbblock:task:firstFeedback)
anova(m.reg,m2)
m.reg <- m2 # accept the reduced model
summary(m.reg)


############
############
### Local confidence during test blocks

## Effect of feedback type (fbblock) on local confidence (test)

first_n_of_trans <- c(1:40)
exp1.df3 <- filter(exp1.mbs, blockType=="transf",
             trialnum %in% first_n_of_trans)

####
#

exp1.mbs.reg = lmer(conf ~ fbblock*task*transferType + 
                 accu + scale_rt + stair +
                 (1|subj) + (1|group) + (1|runnum) + (1|trialnum), exp1.df3)

summary(exp1.mbs.reg)
exp1.mbs.reg = lmer(conf ~ fbblock*task*transferType +
                 accu + scale_rt + stair +
                 (1|subj) + (1|group) + (1|trialnum), exp1.df3)

summary(exp1.mbs.reg)

summary(exp1.mbs.reg)$coefficients
emm <- emmeans(exp1.mbs.reg, ~fbblock|transferType*task)
contrast(emm)

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
  ) %>% 
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


## mediation analysis on only transfer to the same domain
med.m1 <- lmer(conf ~ fbpos_1 + fbneg_1 + (1|subj), filter(exp1.df8, transferType=="Same"))
summary(med.m1)

med.m2 <- lmer(spe_1 ~ fbpos_1 + fbneg_1 + (1|subj), filter(exp1.df8, transferType=="Same"))
summary(med.m2)

med.m3 <- lmer(conf ~ fbpos_1 + fbneg_1 + spe_1 + (1|subj), filter(exp1.df8, transferType=="Same"))
summary(med.m3)

med.pos = mediate(med.m2, med.m3, treat='fbpos_1', mediator='spe_1', 
                  boot=F, sims = 5000)
summary(med.pos)

med.pos = mediate(med.m2, med.m3, treat='fbneg_1', mediator='spe_1', 
                  boot=F, sims = 5000)
summary(med.pos)


## mediation analysis on only transfer to the opposite domain
med.m1 <- lmer(conf ~ fbpos_1 + fbneg_1 + (1|subj), 
               filter(exp1.df8, transferType=="Opposite"))
summary(med.m1)

med.m2 <- lmer(spe_1 ~ fbpos_1 + fbneg_1 + (1|subj), 
               filter(exp1.df8, transferType=="Opposite"))
summary(med.m2)

med.m3 <- lmer(conf ~ fbpos_1 + fbneg_1 + spe_1 + (1|subj), 
               filter(exp1.df8, transferType=="Opposite"))
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
    accu = mean(accu),
    conf = mean(conf),
    scale_rt = mean(scale_rt, na.rm = T),
    stair = mean(stair),
    phq = mean(phq),
    gad = mean(gad),
    spin = mean(spin)
  )

exp1.mbs.reg = lmer(spe ~ fbblock*task*transferType + 
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

plot_model(exp1.mbs.reg,
           title="Perception task - interaction",
           show.p=TRUE, show.values = T,
           value.offset = .4, value.size = 3.8) + 
  coord_cartesian(ylim = c(-.65, .65))

emm <- emmeans(exp1.mbs.reg, pairwise~fbblock)
test(emm)

############
## Mediation analysis - if intervention block SPE mediates transfer of feedback
## to test block SPE

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
  ) %>% 
  ungroup() %>%
  mutate_at("rundiv", as.factor)


exp1.df.temp <- exp1.df.temp %>% dplyr::select(
  subj, rundiv, conf, accu, scale_rt, spe,
  fbneg, fbpos,stair, fbblock, transferType
) %>%
  pivot_wider(
    names_from = c( rundiv),
    values_from = c(conf, accu, scale_rt, fbneg, fbpos, stair, spe)
  ) 

exp1.df10 <- exp1.df9 %>% left_join(exp1.df.temp)


## mediation analyses

med.m1 <- lmer(spe ~ fbpos_1 + fbneg_1 + (1|subj), exp1.df10)
summary(med.m1)

med.m2 <- lmer(spe_1 ~ fbpos_1 + fbneg_1 + (1|subj), exp1.df10)
summary(med.m2)

med.m3 <- lmer(spe ~ fbpos_1 + fbneg_1 + spe_1 + (1|subj), exp1.df10)
summary(med.m3)

med.pos = mediate(med.m2, med.m3, treat='fbpos_1', mediator='spe_1', 
                  boot=F, sims = 5000)
summary(med.pos)

med.pos = mediate(med.m2, med.m3, treat='fbneg_1', mediator='spe_1', 
                  boot=F, sims = 5000)
summary(med.pos)



#######################################
###################################
#### baseline correlation with mhq

exp1.df4 <- exp1.mbs.bas %>% group_by(task, runnum, group, subj, gender) %>%
  summarise(conf = mean(conf),
            accu = mean(accu),
            stair = mean(stair),
            scale_rt = mean(scale_rt),
            spe = mean(spe),
            phq = mean(phq),
            gad = mean(gad),
            spin = mean(spin),
            age = mean(age))

m.reg = lmer(conf ~ phq + 
               age + gender +
               accu + stair + 
             (1|task),exp1.df4)
summary(m.reg)

m.reg = lmer(conf ~ gad + 
               accu + stair + 
               age + gender +
               (1|task),exp1.df4)
summary(m.reg)

m.reg = lmer(conf ~ spin + 
               age + gender +
               accu + stair + 
               (1|task),exp1.df4)
summary(m.reg)


m.reg = lmer(spe ~ phq + 
               age + gender +
               accu + stair + 
               (1|task),exp1.df4)
summary(m.reg)

m.reg = lmer(spe ~ gad + 
               age + gender +
               accu + stair + 
               (1|task),exp1.df4)
summary(m.reg)

m.reg = lmer(spe ~ spin + 
               accu + stair + 
               age + gender +
               (1|task),exp1.df4)
summary(m.reg)


####################################
####################################
###### model predictions with mhq

# setwd("/Users/skatyal/OneDrive - University College London/Projects/Experiments/metaBiasShift/data/exp1")
# mbsDataExp1 = readMat('mbsDataExp1.mat') # load questionnaire data
# mbsDataExp1 <- mbsDataExp1$mbsDataExp1
# exp1.speFit <- mbsDataExp1[,,1]$fitZ[,,1]$model.4[,,1]$spe.est
# nsubjrun <- dim(exp1.speFit)
# subj <- kronecker(matrix(1,1,nsubjrun[2]), c(1:nsubjrun[1]))
# runnum <- t(kronecker(matrix(1,1,nsubjrun[1]), c(1:nsubjrun[2])))
# exp1.speFit <- data.frame(c(exp1.speFit), c(subj), c(runnum))
# colnames(exp1.speFit) <- c('speFit', 'subj', 'runnum')
# exp1.speFit <- exp1.speFit %>%
#   mutate_at(c('subj'), as.factor)

exp1.df2 <- exp1.mbs1 %>%
  select(c(subj, runnum, phq, gad, fbblock, task, spe, feedback, accu, 
           conf, confZ, group)) %>%
  mutate(firstFeedback = ifelse(group %in% c(1,3,5,7), 'Positive', 'Negative')) %>%
  group_by(subj, runnum, fbblock, task, firstFeedback) %>%
  summarise(spe = mean(spe),
            conf = mean(conf),
            confZ = mean(confZ),
            fbpos = sum(feedback==1 & accu==1),
            fbneg = sum(feedback==1 & accu==0),
            fbboth = sum(feedback!=0),
            phq = mean(phq),
            gad = mean(gad)) %>%
  ungroup() %>%
  mutate(phqTile = ntile(phq,2)) %>%
  mutate(gadTile = ntile(gad, 2)) %>%
  mutate_at(c('phqTile', 'gadTile'), as.factor) %>%
  mutate(gadTile = recode_factor(gadTile, '1' = 'Low', '2' = 'High')) %>%
  mutate(phqTile = recode_factor(phqTile, '1' = 'Low', '2' = 'High')) %>%
  # mutate(gadTile = recode_factor(gadTile, '1' = 'Low', '2' = 'Medium', '3' = 'High')) %>%
  # mutate(phqTile = recode_factor(phqTile, '1' = 'Low', '2' = 'Medium', '3' = 'High')) %>%
  mutate_at(c('runnum'), as.numeric) #%>%left_join(exp1.speFit)

exp1.df2$fbblock[exp1.df2$runnum==4|exp1.df2$runnum==6] <- 'None'

exp1.df2.bas <- filter(exp1.df2, runnum %in% c(1,2)) %>%
  group_by(subj) %>%
  summarise(spe.bas = mean(spe),
            # speFit.bas = mean(speFit)
            )

exp1.df3 <- exp1.df2 %>% left_join(exp1.df2.bas) %>%
  group_by(subj, runnum, fbblock, phqTile, gadTile, firstFeedback) %>%
  summarise(
    spe_b = mean(spe) - mean(spe.bas),
    # speFit_b = mean(speFit) - mean(speFit.bas),
    spe = mean(spe),
    conf = mean(conf),
    confZ = mean(confZ),
    # speFit = mean(speFit),
    fbpos = mean(fbpos),
    fbneg = mean(fbneg),
    fbboth = mean(fbboth),
    fbdif = mean(fbpos)-mean(fbneg),
    phq = mean(phq),
    gad = mean(gad)
  ) 


##################################
##### Figure 4A and 4B. model predictions

f4a <- ggplot(exp1.df3 %>% filter(runnum %in% c(3,5)), 
       aes(x = fbblock, y = spe_b, color = phqTile)) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', 'royalblue3')) +
  stat_summary(fun.data = mean_se, size = .6, 
               position=position_dodge(width = .1), alpha=.8) +
  stat_summary(fun = mean, size = 1.5, geom = 'line', aes(group=phqTile),
               position=position_dodge(width = .1), alpha=.8) +
  theme_pubclean() +
  theme(text = element_text(size=font_size), 
        plot.caption = element_text(size = font_size),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none') +
  coord_cartesian(ylim=c(-.2,.15)) +
  ylab("SPE-b") +
  xlab("Feedback") + labs(color = 'PHQ score')
f4a


exp1.df2 <- exp1.mbs1 %>%
  mutate(phqTile = ntile(phq, 2)) %>%
  mutate(gadTile = ntile(gad, 2)) %>%
  drop_na(scale_rt) %>%
  group_by(subj, runnum) %>%
  summarise(spe = mean(spe),
            confZ = mean(confZ),
            phqTile = median(phqTile),
            gadTile = median(gadTile)
  ) %>%
  mutate(confTile = ntile(confZ, 6)) %>%
  mutate(hiloConf = ifelse(confTile<=2, 0, 1)) %>%
  mutate_at(c('gadTile', 'hiloConf', 'phqTile'), as.factor) %>%
  mutate(hiloConf = recode_factor(hiloConf, '0' = 'lo', '1' = 'hi')) %>%
  mutate(phqTile = recode_factor(phqTile, '1' = 'Low', '2' ='High')) %>%
  mutate(gadTile = recode_factor(gadTile, '1' = 'Low', '2'  = 'High'))# %>%# left_join(exp1.speFit) 

f4b <- ggplot(exp1.df2 %>% filter(runnum %in% c(1,2,4,6)), 
       aes(x = runnum, y = spe, color = phqTile)) +
  stat_summary(fun.data = mean_se, size = .8,
               position=position_dodge(width = .1), alpha=.8) +
  scale_color_manual(breaks = c('Low', 'High'),
                     values = c('tan3', 'royalblue3')) +
  stat_summary(fun = mean, geom = 'line', size=1.5,
               position=position_dodge(width = .1), alpha=.8) +
  theme_pubclean() +
  theme(text = element_text(size=font_size),
        plot.caption = element_text(size = font_size),
        # legend.title = element_text(color = 'blue'),
        plot.tag = element_text(size = tag_size, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none') +
  ylab("SPE") +
  coord_cartesian(ylim=c(.47,.65)) +
  scale_x_continuous(breaks = c(1,2,4,6))+
  xlab("Block #") +
  labs(colour = 'PHQ score', linetype = 'Conf level')
f4b



# 
# ### how does phq relate feedback in forming spe(intervention)
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

if (Sys.info()["sysname"]=='Darwin'){
  setwd("~/OneDrive - University College London/Projects/Experiments/metaBiasShift/")
} else{
  setwd("C:/Users/skatyal/OneDrive - University College London/Projects/Experiments/metaBiasShift/")
}
## read data file
exp1.mbs3 = read.csv(paste('data/exp1/mbsExp1_noperfexcl.csv', sep=''),
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
         task = recode_factor(task, "0" = taskNames[1], "1" = taskNames[2])
  )

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
               (1|wordnum), exp1.df6,
             family = binomial)
summary(m.reg)

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

