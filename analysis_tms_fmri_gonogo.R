setwd('/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project')

df=read.csv('data/tms-plip-out_061523_edit.csv')

img_vars=c('act_gonogo_nogovsgo_013136', 'act_gonogo_nogovsgo_533294', 'act_gonogo_nogovsgo_291082', 'ppi_gonogo_thrvneu_013136_291082', 'ppi_gonogo_thrvneu_291082_013136','ppi_gonogo_thrvneu_533294_291082', 'ppi_gonogo_thrvneu_291082_533294')
img_vars_rename=c('act_gonogo_LdlPFC', 'act_gonogo_RdlPFC', 'act_gonogo_MdACC', 'ppi_gonogo_LdlPFC_MdACC', 'ppi_gonogo_MdACC_LdlPFC', 'ppi_gonogo_RdlPFC_MdACC', 'ppi_gonogo_MdACC_RdlPFC')

library(data.table)
setnames(df, img_vars, img_vars_rename)

#### COMBAT correction ####

# Import controls
hc=read.csv('data/data_hc.csv')

# Make hc and clinical have similar variables
img_vars_hc=c('id', 'scanner', 'act.gng.nogo_vs_go.013136_Left_dlPFC', 'act.gng.nogo_vs_go.533294_Right_dlPFC', 'act.gng.nogo_vs_go.291082_Medial_dACC', 'ppi.gng.nogo_vs_go.013136_Left_dlPFC.291082_Medial_dACC', 'ppi.gng.nogo_vs_go.291082_Medial_dACC.013136_Left_dlPFC', 'ppi.gng.nogo_vs_go.533294_Right_dlPFC.291082_Medial_dACC', 'ppi.gng.nogo_vs_go.291082_Medial_dACC.533294_Right_dlPFC')
img_vars_rename_hc=c('Subjects', 'Site', 'act_gonogo_LdlPFC', 'act_gonogo_RdlPFC', 'act_gonogo_MdACC', 'ppi_gonogo_LdlPFC_MdACC', 'ppi_gonogo_MdACC_LdlPFC', 'ppi_gonogo_RdlPFC_MdACC', 'ppi_gonogo_MdACC_RdlPFC')
setnames(hc, img_vars_hc, img_vars_rename_hc)
hc$Site=factor(hc$Site, labels = c('Sidney', 'Palo Alto Discovery'))
hc$Visit=0
df$group=1

# Merge the hc and clinical
img_only=rbind(df[, c('Visit', 'group', img_vars_rename_hc)], hc[, c('Visit', 'group', img_vars_rename_hc)])

# COMBAT
library(ez.combat)
img_only_adjusted_combat  =  img_only[complete.cases(img_only[c('Visit', 'Subjects', 'Site', 'group', img_vars_rename)]), ]
cb=ez.combat(img_only_adjusted_combat,
             'Site',
             adjust.var = img_vars_rename,
             output = c("overwrite"),
             use.eb = TRUE,
             verbose = TRUE)
img_only_adjusted_combat=cb$df

# Scale clinical imaging variables to controls
for (var in img_vars_rename){
  mn=mean(img_only_adjusted_combat[img_only_adjusted_combat$group==0, var])
  sdev=sd(img_only_adjusted_combat[img_only_adjusted_combat$group==0, var])
  img_only_adjusted_combat[, var]=(img_only_adjusted_combat[, var]-mn)/sdev
}

# Keep only clinical
img_only_adjusted_combat=img_only_adjusted_combat[img_only_adjusted_combat$group==1, ]

#### MERGING OF DATA 

# Add questionnaire data
ques=read.csv('data/TMS Data 2.9.23.csv')
data_adjusted_combat=merge(img_only_adjusted_combat, ques, by = c('Subjects', 'Visit', 'Site'), all = TRUE)
qids=read.csv('data/clinical_edit.csv')
data_adjusted_combat=merge(data_adjusted_combat, qids, by = c('Subjects', 'Visit'), all = TRUE)

# Add behavioral data
beh=read.csv('data/PsychoStats_mod.csv')
data_adjusted_combat=merge(data_adjusted_combat, beh, by = c('Subjects', 'Visit'), all = TRUE)

# Add motion
fd=read.csv('data/fd.csv')
data_adjusted_combat=merge(data_adjusted_combat, fd, by = c('Subjects', 'Visit'), all = TRUE)

# Merge Site variable
library(dplyr)
library(purrr)
data_adjusted_combat  =  data_adjusted_combat %>%
  mutate(Site = coalesce(Site.x, Site.y))

# Import dates as assessments
dates=read.csv('data/all_dates.csv')

# Calculate time from baseline
s0=as.Date(dates[, 'Symptom_date_0'], format = '%m/%d/%y')
for (var in c('Symptom_date_0', 'Symptom_date_1', 'Symptom_date_2')){
  dates[, var] =  difftime(as.Date(dates[, var], format = '%m/%d/%y'), s0 , units = c("days"))
}
s0=as.Date(dates[, 'WebNeuro_date_0'], format = '%m/%d/%y')
for (var in c('WebNeuro_date_0', 'WebNeuro_date_1', 'WebNeuro_date_2')){
  dates[, var] =  difftime(as.Date(dates[, var], format = '%m/%d/%y'), s0 , units = c("days"))
}
s0=as.Date(dates[, 'MRI_date_0'], format = '%m/%d/%y')
for (var in c('MRI_date_0', 'MRI_date_1', 'MRI_date_2')){
  dates[, var] =  difftime(as.Date(dates[, var], format = '%m/%d/%y'), s0 , units = c("days"))
}

# Add the times for each session
for (r in 1:nrow(data_adjusted_combat)){
  sub=data_adjusted_combat[r, 'Subjects']
  vis=data_adjusted_combat[r, 'Visit']
  data_adjusted_combat[r, 'Symptom_timediff']=as.numeric(dates[dates$Subjects==sub, paste('Symptom_date_', vis, sep='')])
  data_adjusted_combat[r, 'WebNeuro_timediff']=as.numeric(dates[dates$Subjects==sub, paste('WebNeuro_date_', vis, sep='')])
  data_adjusted_combat[r, 'MRI_timediff']=as.numeric(dates[dates$Subjects==sub, paste('MRI_date_', vis, sep='')])
}

# Keep only people in the first half of the sample
subs_orig=c('1976-320', '1976-322', '1976-318', '2000-012', '2000-018', '2000-030', '1976-312', '2000-023', '1976-311', '1991-64', '2000-027', '1984-5', '2000-016', '1976-323', '1976-315', '1976-319', '1976-313', '2000-029', '1976-325', '1991-72', '1991-92', '1991-43', '1976-317', '2000-022', '1991-107', '2000-010', '2000-031', '1991-70', '1984-11', '1984-10', '1991-103', '1991-114', '1984-7', '1991-90', '1984-26', '1976-309', '2000-011', '1991-58', '1976-308', '1991-93', '1991-74', '2000-013', '1976-314')
data_adjusted_combat=data_adjusted_combat[data_adjusted_combat$Subjects %in% subs_orig, ]

# Change Visit to factor
data_adjusted_combat$Visit=factor(data_adjusted_combat$Visit, labels = c('Baseline', 'Early Treatment', 'Post-Treatment'))

# Write data
write.csv(data_adjusted_combat, 'data/df_merged.csv', row.names = F)

#### GENERATE NEW VARIABLES ####

# Define change, responders, remitters and 0 split
for (sub in unique(data_adjusted_combat$Subjects)){
  
  # Store ppi_gonogo_LdlPFC_MdACC change from baseline
  ppi_bl=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Baseline', 'ppi_gonogo_LdlPFC_MdACC']
  ppi_v1=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'ppi_gonogo_LdlPFC_MdACC']
  ppi_v2=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'ppi_gonogo_LdlPFC_MdACC']
  data_adjusted_combat[data_adjusted_combat$Subjects==sub  && data_adjusted_combat$Visit=='Baseline', 'ppi_diff_bl']= 0
  if (length(ppi_bl)>0 & length(ppi_v1)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'ppi_diff_bl']= ppi_v1-ppi_bl}
  if (length(ppi_bl)>0 & length(ppi_v2)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'ppi_diff_bl']= ppi_v2-ppi_bl}
  if (length(ppi_bl)>0 ){data_adjusted_combat[data_adjusted_combat$Subjects==sub, 'ppi_bl']= ppi_bl}
  
  # Store QIDS change from baseline
  q_total_bl=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Baseline', 'q_total']
  q_total_v1=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'q_total']
  q_total_v2=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'q_total']
  data_adjusted_combat[data_adjusted_combat$Subjects==sub  && data_adjusted_combat$Visit=='Baseline', 'q_total_diff_bl']= 0
  if (length(q_total_bl)>0 & length(q_total_v1)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'q_total_diff_bl']= q_total_v1-q_total_bl}
  if (length(q_total_bl)>0 & length(q_total_v2)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'q_total_diff_bl']= q_total_v2-q_total_bl}
  if (length(q_total_bl)>0 ){data_adjusted_combat[data_adjusted_combat$Subjects==sub, 'q_total_bl']=q_total_bl}
  
  # Store errors change from baseline
  g2errk_norm_bl=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Baseline', 'g2errk_norm']
  g2errk_norm_v1=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'g2errk_norm']
  g2errk_norm_v2=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'g2errk_norm']
  data_adjusted_combat[data_adjusted_combat$Subjects==sub  && data_adjusted_combat$Visit=='Baseline', 'g2errk_norm_diff_bl']= 0
  if (length(g2errk_norm_bl)>0 & length(g2errk_norm_v1)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'g2errk_norm_diff_bl']= g2errk_norm_v1-g2errk_norm_bl}
  if (length(g2errk_norm_bl)>0 & length(g2errk_norm_v2)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'g2errk_norm_diff_bl']= g2errk_norm_v2-g2errk_norm_bl}
  if (length(g2errk_norm_bl)>0 ){data_adjusted_combat[data_adjusted_combat$Subjects==sub, 'g2errk_norm_bl']=g2errk_norm_bl}
  
  # Store QIDS change from baseline
  qids_bl=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Baseline', 'q_total']
  qids_fu=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'q_total']
  if (length(qids_bl)>0 & length(qids_fu)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub, 'qids_diff']=qids_fu-qids_bl}
  if (length(qids_bl)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub, 'qids_bl']=qids_bl}
  if (length(qids_fu)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub, 'qids_fu']=qids_fu}
  
  # Zero split of Biotype values
  bl_conn=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Baseline', 'ppi_gonogo_LdlPFC_MdACC']
  if (length(bl_conn)>0 && bl_conn>=0 && !is.na(bl_conn)){data_adjusted_combat[data_adjusted_combat$Subjects==sub, 'split_0']=1}
  if (length(bl_conn)>0 && bl_conn<0 && !is.na(bl_conn)){data_adjusted_combat[data_adjusted_combat$Subjects==sub, 'split_0']=0}
  
  # Store QIDS10 change from baseline
  q_10_bl=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Baseline', 'q10']
  q_10_v1=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'q10']
  q_10_v2=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'q10']
  data_adjusted_combat[data_adjusted_combat$Subjects==sub  && data_adjusted_combat$Visit=='Baseline', 'q_10_diff_bl']= 0
  if (length(q_total_bl)>0 & length(q_total_v1)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'q_10_diff_bl']= q_10_v1-q_10_bl}
  if (length(q_total_bl)>0 & length(q_total_v2)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'q_10_diff_bl']= q_10_v2-q_10_bl}
  if (length(q_total_bl)>0 ){data_adjusted_combat[data_adjusted_combat$Subjects==sub, 'q_10_bl']=q_10_bl}
  
}  

data_adjusted_combat$remit=factor((data_adjusted_combat$qids_fu<=5)*1)
data_adjusted_combat$qids_perc=(data_adjusted_combat$qids_bl-data_adjusted_combat$qids_fu)/data_adjusted_combat$qids_bl
data_adjusted_combat$respond=factor((data_adjusted_combat$qids_perc>0.5)*1)
data_adjusted_combat$split_0  =  factor(data_adjusted_combat$split_0, levels = c(0, 1), labels = c("C -", "C +"))

#### Connectivity and time ####

# Mixed model
library(lme4)
library(lmerTest)
library(emmeans)

# Compare at baseline
bl=data_adjusted_combat[data_adjusted_combat$Visit=='Baseline', ]
wilcox.test(data=bl, ppi_gonogo_LdlPFC_MdACC~split_0)

# Comparison between groups for demographics
t.test(data=bl, sessage ~split_0)
chisq.test(bl$Gender, bl$split_0)

# Fit the linear mixed effects model
model  =  lmer(data=data_adjusted_combat, ppi_diff_bl ~ split_0*Visit + MRI_timediff + (1 | Subjects))

# Display the model summary
summary(model)
anova(model)

# Extract the estimated marginal means and their standard errors
emmeans_results  =  emmeans(model, ~Visit|split_0)
emmeans_summary  =  summary(emmeans_results)
emmeans_contrast  =  pairs(emmeans_results, adjust = "none")

# Create a dataframe from the summary object
emmeans_df  =  data.frame(emmeans_summary)

#### Behavior and time ####

# Compare at baseline
bl=data_adjusted_combat[data_adjusted_combat$Visit=='Baseline', ]
t.test(data=bl, g2errk_norm_bl~split_0)

# Fit the linear mixed effects model
model  =  lmer(data=data_adjusted_combat, g2errk_norm_diff_bl ~ split_0*Visit + WebNeuro_timediff + (1 | Subjects))

# Display the model summary
summary(model)
anova(model)

# Extract the estimated marginal means and their standard errors
emmeans_results  =  emmeans(model, ~Visit|split_0)
emmeans_summary  =  summary(emmeans_results)

# Create a dataframe from the summary object
emmeans_df  =  data.frame(emmeans_summary)
emmeans_results  =  emmeans(model, ~Visit|split_0)
emmeans_contrast  =  pairs(emmeans_results, adjust = "none")

#### QIDS and time ####

# Compare at baseline
bl=data_adjusted_combat[data_adjusted_combat$Visit=='Baseline', ]
wilcox.test(data=bl, q_total~split_0)

# Fit the linear mixed effects model
model  =  lmer(data=data_adjusted_combat, q_total_diff_bl ~ split_0*Visit + Symptom_timediff + (1 | Subjects))

# Display the model summary
summary(model)
anova(model)

# Extract the estimated marginal means and their standard errors
emmeans_results  =  emmeans(model, ~Visit|split_0)
emmeans_summary  =  summary(emmeans_results)

# Create a dataframe from the summary object
emmeans_df  =  data.frame(emmeans_summary)
emmeans_results  =  emmeans(model, ~Visit)
emmeans_contrast  =  pairs(emmeans_results, adjust = "none")

#### Correlations

# Correlation between PPI change and GNG error change
cor.test(data_adjusted_combat$ppi_diff_bl, data_adjusted_combat$g2errk_norm_diff_bl, use = 'pairwise.complete.obs')

# Correlation between in scanner GNG errors change and GNG errors
data_adjusted_combat$scanner_errors=data_adjusted_combat$Commission.Errors+data_adjusted_combat$Omission.Errors
temp=data_adjusted_combat[data_adjusted_combat$scanner_errors<50 & data_adjusted_combat$g2errk<50, ]
cor.test(temp$scanner_errors, temp$g2errk, use = 'pairwise.complete.obs', method = 'spearman')

#### Responders and remitters
data_adjusted_combat_bl=data_adjusted_combat[data_adjusted_combat$Visit=='Baseline', ]

# Responders
table(data_adjusted_combat_bl$respond, data_adjusted_combat_bl$split_0)
chisq.test(data_adjusted_combat_bl$respond, data_adjusted_combat_bl$split_0)

# Remitters
table(data_adjusted_combat_bl$remit, data_adjusted_combat_bl$split_0)
chisq.test(data_adjusted_combat_bl$remit, data_adjusted_combat_bl$split_0)

# Save the data frame
write.csv(data_adjusted_combat, 'data/data_adjusted_combat.csv', row.names = F)

#### Models with covariates ####

# Activation
model  =  lmer(data=data_adjusted_combat, ppi_diff_bl ~ split_0*Visit + act_gonogo_LdlPFC + MRI_timediff + (1 | Subjects))
anova(model)
emmeans_results  =  emmeans(model, ~Visit|split_0)
emmeans_contrast  =  pairs(emmeans_results, adjust = "none")

model  =  lmer(data=data_adjusted_combat, g2errk_norm_diff_bl ~ split_0*Visit + act_gonogo_LdlPFC + WebNeuro_timediff + (1 | Subjects))
anova(model)
emmeans_results  =  emmeans(model, ~Visit|split_0)
emmeans_contrast  =  pairs(emmeans_results, adjust = "none")

model  =  lmer(data=data_adjusted_combat, q_total_diff_bl ~ split_0*Visit + act_gonogo_LdlPFC + Symptom_timediff + (1 | Subjects))
anova(model)
emmeans_results  =  emmeans(model, ~Visit)
emmeans_contrast  =  pairs(emmeans_results, adjust = "none")


# Motion
model  =  lmer(data=data_adjusted_combat, ppi_diff_bl ~ split_0*Visit + FD_03 + MRI_timediff + (1 | Subjects))
anova(model)
emmeans_results  =  emmeans(model, ~Visit|split_0)
emmeans_contrast  =  pairs(emmeans_results, adjust = "none")

model  =  lmer(data=data_adjusted_combat, g2errk_norm_diff_bl ~ split_0*Visit + FD_03 + WebNeuro_timediff + (1 | Subjects))
anova(model)
emmeans_results  =  emmeans(model, ~Visit|split_0)
emmeans_contrast  =  pairs(emmeans_results, adjust = "none")

model  =  lmer(data=data_adjusted_combat, q_total_diff_bl ~ split_0*Visit + FD_03 + Symptom_timediff + (1 | Subjects))
anova(model)
emmeans_results  =  emmeans(model, ~Visit)
emmeans_contrast  =  pairs(emmeans_results, adjust = "none")


#### Behavior

# Fit the linear mixed effects model
model  =  lmer(data=data_adjusted_combat, g2errk_norm_diff_bl ~ split_0*Visit + WebNeuro_timediff + (1 | Subjects))

# Extract the estimated marginal means and their standard errors
emmeans_results  =  emmeans(model, ~Visit|split_0)
emmeans_summary  =  summary(emmeans_results)

# Create a dataframe from the summary object
emmeans_df  =  data.frame(emmeans_summary)
emmeans_results  =  emmeans(model, ~Visit|split_0)

#### QIDS and time ####

# Fit the linear mixed effects model
model  =  lmer(data=data_adjusted_combat, q_total_diff_bl ~ split_0*Visit + Symptom_timediff + (1 | Subjects))

# Extract the estimated marginal means and their standard errors
emmeans_results  =  emmeans(model, ~Visit|split_0)
emmeans_summary  =  summary(emmeans_results)

# Create a dataframe from the summary object
emmeans_df  =  data.frame(emmeans_summary)
emmeans_results  =  emmeans(model, ~Visit)

#### MEDIATION MODEL

library("mediation")

set.seed(123124)

temp=temp  =  data_adjusted_combat %>%
  filter(!is.na(split_0) & !is.na(ppi_diff_bl) & !is.na(g2errk_norm_diff_bl))

fit.mediator=lm(ppi_diff_bl~split_0,temp)
summary(fit.mediator)
fit.dv=lm(g2errk_norm_diff_bl~ppi_diff_bl+split_0,temp)
summary(fit.dv)

results = mediate(fit.mediator, fit.dv, treat='split_0', mediator='ppi_diff_bl', boot=F, control.value = '>=0', treat.value = '<0')
summary(results)

#### COMPARE GONOGO PERFORMANCE TO OTHER BEHAVIORS

# Create new webneuro variables
bl$cog_con=rowMeans(bl[, c('g2fpk_norm', 'g2avrtk_norm', 'g2fnk_norm')], na.rm = T)
bl$verbal_mem=rowMeans(bl[, c('ctmlearn_norm', 'ctmrec4_norm', 'ctmsco13_norm')], na.rm = T)
bl$working_mem=rowMeans(bl[, c('digitot_norm', 'digitsp_norm')], na.rm = T)
bl$sustained_att=rowMeans(bl[, c('wmfpk_norm', 'wmfnk_norm', 'wmrtk_norm')], na.rm = T)

# Correlation between gonogo and other variables
cor.test(bl$g2errk_norm_bl, bl$verbal_mem, method='spearman')
cor.test(bl$g2errk_norm_bl, bl$working_mem, method='spearman')
cor.test(bl$g2errk_norm_bl, bl$sustained_att, method='spearman')

# Comparison between groups
t.test(data=bl, verbal_mem ~split_0)
t.test(data=bl, working_mem ~split_0)
t.test(data=bl, sustained_att ~split_0)

# Correlation between measures that are not gonogo
data_adjusted_combat$cog_con=rowMeans(data_adjusted_combat[, c('g2fpk_norm', 'g2avrtk_norm', 'g2fnk_norm')], na.rm = T)
data_adjusted_combat$verbal_mem=rowMeans(data_adjusted_combat[, c('ctmlearn_norm', 'ctmrec4_norm', 'ctmsco13_norm')], na.rm = T)
data_adjusted_combat$working_mem=rowMeans(data_adjusted_combat[, c('digitot_norm', 'digitsp_norm')], na.rm = T)
data_adjusted_combat$sustained_att=rowMeans(data_adjusted_combat[, c('wmfpk_norm', 'wmfnk_norm', 'wmrtk_norm')], na.rm = T)
for (sub in unique(data_adjusted_combat$Subjects)){
  
  # Store change from baseline
  verbal_mem_bl=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Baseline', 'verbal_mem']
  verbal_mem_v1=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'verbal_mem']
  verbal_mem_v2=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'verbal_mem']
  data_adjusted_combat[data_adjusted_combat$Subjects==sub  && data_adjusted_combat$Visit=='Baseline', 'verbal_mem_diff_bl']= 0
  if (length(verbal_mem_bl)>0 & length(verbal_mem_v1)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'verbal_mem_diff_bl']= verbal_mem_v1-verbal_mem_bl}
  if (length(verbal_mem_bl)>0 & length(verbal_mem_v2)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'verbal_mem_diff_bl']= verbal_mem_v2-verbal_mem_bl}
  if (length(verbal_mem_bl)>0 ){data_adjusted_combat[data_adjusted_combat$Subjects==sub, 'verbal_mem_bl']=verbal_mem_bl}
  
  # Store change from baseline
  working_mem_bl=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Baseline', 'working_mem']
  working_mem_v1=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'working_mem']
  working_mem_v2=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'working_mem']
  data_adjusted_combat[data_adjusted_combat$Subjects==sub  && data_adjusted_combat$Visit=='Baseline', 'working_mem_diff_bl']= 0
  if (length(working_mem_bl)>0 & length(working_mem_v1)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'working_mem_diff_bl']= working_mem_v1-working_mem_bl}
  if (length(working_mem_bl)>0 & length(working_mem_v2)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'working_mem_diff_bl']= working_mem_v2-working_mem_bl}
  if (length(working_mem_bl)>0 ){data_adjusted_combat[data_adjusted_combat$Subjects==sub, 'working_mem_bl']=working_mem_bl}
  
  # Store change from baseline
  sustained_att_bl=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Baseline', 'sustained_att']
  sustained_att_v1=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'sustained_att']
  sustained_att_v2=data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'sustained_att']
  data_adjusted_combat[data_adjusted_combat$Subjects==sub  && data_adjusted_combat$Visit=='Baseline', 'sustained_att_diff_bl']= 0
  if (length(sustained_att_bl)>0 & length(sustained_att_v1)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Early Treatment', 'sustained_att_diff_bl']= sustained_att_v1-sustained_att_bl}
  if (length(sustained_att_bl)>0 & length(sustained_att_v2)>0){data_adjusted_combat[data_adjusted_combat$Subjects==sub & data_adjusted_combat$Visit=='Post-Treatment', 'sustained_att_diff_bl']= sustained_att_v2-sustained_att_bl}
  if (length(sustained_att_bl)>0 ){data_adjusted_combat[data_adjusted_combat$Subjects==sub, 'sustained_att_bl']=sustained_att_bl}
  
}  

cor.test(data_adjusted_combat$ppi_diff_bl, data_adjusted_combat$verbal_mem_diff_bl, use = 'pairwise.complete.obs')
cor.test(data_adjusted_combat$ppi_diff_bl, data_adjusted_combat$working_mem_diff_bl, use = 'pairwise.complete.obs')
cor.test(data_adjusted_combat$ppi_diff_bl, data_adjusted_combat$sustained_att_diff_bl, use = 'pairwise.complete.obs')

#### PLOTS ####

library(ggplot2) 

#### Behavior at baseline

# Compute means and standard errors
bl_summary  =  bl %>%
  group_by(split_0) %>%
  summarise(
    mean_value = mean(g2errk_norm_bl, na.rm = TRUE),
    se_value = sd(g2errk_norm_bl, na.rm = TRUE) / sqrt(n())
  )

# Create the bar plot
png(filename="/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/gng_bl_bar.png", width = 500, height = 600)
ggplot(bl_summary, aes(x = split_0, y = mean_value, fill = split_0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.2) +
  theme_minimal() +
  xlab('Biotype') +
  ylab('Go-NoGo Performance at Baseline') +
  theme(legend.position = "none") + 
  theme(text = element_text(size = 24)) +
  scale_fill_manual(values = c("#B75A65", "gray40")) +
  ylim(c(-1, 1))
dev.off()

### Connectivity raw

# Calculate means and standard errors
data_means  =  data_adjusted_combat_temp %>%
  group_by(Visit, split_0) %>%
  summarise(mean = mean(ppi_gonogo_LdlPFC_MdACC, na.rm = TRUE),
            se = sd(ppi_gonogo_LdlPFC_MdACC, na.rm = TRUE) / sqrt(n()), 
            .groups = "drop")

# Line plot
png(filename="/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/conn_raw_mean.png", width = 600, height = 600)
ggplot(data_means, aes(x = Visit, y = mean, color = split_0)) +
  geom_line(aes(group = split_0), size = 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, size=1) +
  labs(y = "Connectivity", x = "") +
  theme_minimal() +
  theme(text = element_text(size = 24))+
  theme(legend.position = "top") +
  scale_color_manual(name = "Biotype", values = c("#B75A65", "gray40"))+ 
  ylim(c(-1, 1))
dev.off()

# Dots and lines plot
png(filename="/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/conn_raw_mean_dots.png", width = 1000, height = 600)
ggplot() +
  geom_jitter(data = data_adjusted_combat_temp, aes(x = Visit, y = ppi_gonogo_LdlPFC_MdACC, color = split_0), 
              size = 2, alpha = 0.5, position = position_dodge(0.8)) +
  geom_line(data = data_adjusted_combat_temp, aes(x = Visit, y = ppi_gonogo_LdlPFC_MdACC, group = interaction(Subjects, split_0), color = split_0), 
            alpha = 0.5, size = 0.5) +
  geom_line(data = data_means, aes(x = Visit, y = mean, group = split_0, color = split_0), size = 2) +
  geom_errorbar(data = data_means, aes(x = Visit, y = mean, ymin = mean - se, ymax = mean + se, color = split_0), 
                width = 0.2, size = 1) +
  facet_wrap(~ split_0) +
  labs(y = 'Connectivity', color = 'Biotype', x = '') +
  theme_minimal() +
  theme(legend.position = "top", text = element_text(size = 24)) +
  scale_color_manual(values = c("#B75A65", "gray40")) +
  scale_fill_manual(values = c("#B75A65", "gray40")) 
dev.off()

#### Behavior raw

# Calculate means and standard errors
data_means  =  data_adjusted_combat_temp %>%
  group_by(Visit, split_0) %>%
  summarise(mean = mean(g2errk_norm, na.rm = TRUE),
            se = sd(g2errk_norm, na.rm = TRUE) / sqrt(n()), 
            .groups = "drop")

# Line plot
png(filename="/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/beh_raw_mean.png", width = 600, height = 600)
ggplot(data_means, aes(x = Visit, y = mean, color = split_0)) +
  geom_line(aes(group = split_0), size = 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, size=1) +
  labs(y = "Go-NoGo Performance", x = "") +
  theme_minimal() +
  theme(text = element_text(size = 24))+
  theme(legend.position = "top") +
  scale_color_manual(name = "Biotype", values = c("#B75A65", "gray40"))+ 
  ylim(c(-1, 1))
dev.off()

# Dots and lines plot
png(filename="/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/beh_raw_mean_dots.png", width = 1000, height = 600)
ggplot() +
  geom_jitter(data = data_adjusted_combat_temp, aes(x = Visit, y = g2errk_norm, color = split_0), 
              size = 2, alpha = 0.5, position = position_dodge(0.8)) +
  geom_line(data = data_adjusted_combat_temp, aes(x = Visit, y = g2errk_norm, group = interaction(Subjects, split_0), color = split_0), 
            alpha = 0.5, size = 0.5) +
  geom_line(data = data_means, aes(x = Visit, y = mean, group = split_0, color = split_0), size = 2) +
  geom_errorbar(data = data_means, aes(x = Visit, y = mean, ymin = mean - se, ymax = mean + se, color = split_0), 
                width = 0.2, size = 1) +
  facet_wrap(~ split_0) +
  labs(y = 'Go-NoGo Performance', color = 'Biotype', x = '') +
  theme_minimal() +
  theme(legend.position = "top", text = element_text(size = 24)) +
  scale_color_manual(values = c("#B75A65", "gray40")) +
  scale_fill_manual(values = c("#B75A65", "gray40")) 
dev.off()


#### QIDS raw

# Calculate means and standard errors
data_means  =  data_adjusted_combat_temp %>%
  group_by(Visit, split_0) %>%
  summarise(mean = mean(q_total, na.rm = TRUE),
            se = sd(q_total, na.rm = TRUE) / sqrt(n()), 
            .groups = "drop")

# Line plot
png(filename="/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/qids_raw_mean.png", width = 600, height = 600)
ggplot(data_means, aes(x = Visit, y = mean, color = split_0)) +
  geom_line(aes(group = split_0), size = 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, size=1) +
  labs(y = "QIDS total", x = "") +
  theme_minimal() +
  theme(text = element_text(size = 24))+
  theme(legend.position = "top") +
  scale_color_manual(name = "Biotype", values = c("#B75A65", "gray40"))+ 
  ylim(c(0, 20))
dev.off()

# Dots and lines plot
png(filename="/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/qids_raw_mean_dots.png", width = 1000, height = 600)
ggplot() +
  geom_jitter(data = data_adjusted_combat_temp, aes(x = Visit, y = q_total, color = split_0), 
              size = 2, alpha = 0.5, position = position_dodge(0.8)) +
  geom_line(data = data_adjusted_combat_temp, aes(x = Visit, y = q_total, group = interaction(Subjects, split_0), color = split_0), 
            alpha = 0.5, size = 0.5) +
  geom_line(data = data_means, aes(x = Visit, y = mean, group = split_0, color = split_0), size = 2) +
  geom_errorbar(data = data_means, aes(x = Visit, y = mean, ymin = mean - se, ymax = mean + se, color = split_0), 
                width = 0.2, size = 1) +
  facet_wrap(~ split_0) +
  labs(y = 'QIDS total', color = 'Biotype', x = '') +
  theme_minimal() +
  theme(legend.position = "top", text = element_text(size = 24)) +
  scale_color_manual(values = c("#B75A65", "gray40")) +
  scale_fill_manual(values = c("#B75A65", "gray40")) 
dev.off()


#### Change in connectivity

# Line plot
data_means  =  temp %>%
  group_by(Visit, split_0) %>%
  summarise(mean = mean(ppi_diff_bl, na.rm = TRUE),
            se = sd(ppi_diff_bl, na.rm = TRUE) / sqrt(n()), 
            .groups = "drop")

# Create the raw line plot
png(filename="/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/conn_raw_diff_mean.png", width = 800, height = 600)
ggplot(data_means, aes(x = Visit, y = mean, color = split_0)) +
  geom_line(aes(group = split_0), size = 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, size=1) +
  labs(y = "Connectivity Change", x = "") +
  theme_minimal() +
  theme(text = element_text(size = 24))+
  theme(legend.position = "top") +
  scale_color_manual(name = "Biotype", values = c("#B75A65", "gray40"))
dev.off()

# Dots and lines plot
png(filename="/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/conn_raw_diff_mean_dots.png", width = 1000, height = 600)
ggplot() +
  geom_jitter(data = data_adjusted_combat_temp, aes(x = Visit, y = ppi_diff_bl, color = split_0), 
              size = 2, alpha = 0.5, position = position_dodge(0.8)) +
  geom_line(data = data_adjusted_combat_temp, aes(x = Visit, y = ppi_diff_bl, group = interaction(Subjects, split_0), color = split_0), 
            alpha = 0.5, size = 0.5) +
  geom_line(data = data_means, aes(x = Visit, y = mean, group = split_0, color = split_0), size = 2) +
  geom_errorbar(data = data_means, aes(x = Visit, y = mean, ymin = mean - se, ymax = mean + se, color = split_0), 
                width = 0.2, size = 1) +
  facet_wrap(~ split_0) +
  labs(y = 'Connectivity Change', color = 'Biotype', x = '') +
  theme_minimal() +
  theme(legend.position = "top", text = element_text(size = 24)) +
  scale_color_manual(values = c("#B75A65", "gray40")) +
  scale_fill_manual(values = c("#B75A65", "gray40")) 
dev.off()

#### Change in Behavior

# Calculate means and standard errors
data_means  =  temp %>%
  group_by(Visit, split_0) %>%
  summarise(mean = mean(g2errk_norm_diff_bl, na.rm = TRUE),
            se = sd(g2errk_norm_diff_bl, na.rm = TRUE) / sqrt(n()), 
            .groups = "drop")

# Line plot
png(filename="/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/beh_raw_diff.png", width = 800, height = 600)
ggplot(data_means, aes(x = Visit, y = mean, color = split_0)) +
  geom_line(aes(group = split_0), size = 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, size=1) +
  labs(y = "Go-NoGo Performance Change", x = "") +
  theme_minimal() +
  theme(text = element_text(size = 24))+
  theme(legend.position = "top") +
  scale_color_manual(name = "Biotype", values = c("#B75A65", "gray40"))
dev.off()

# Dots and lines plot
png(filename="/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/beh_raw_diff_mean_dots.png", width = 1000, height = 600)
ggplot() +
  geom_jitter(data = data_adjusted_combat_temp, aes(x = Visit, y = g2errk_norm_diff_bl, color = split_0), 
              size = 2, alpha = 0.5, position = position_dodge(0.8)) +
  geom_line(data = data_adjusted_combat_temp, aes(x = Visit, y = g2errk_norm_diff_bl, group = interaction(Subjects, split_0), color = split_0), 
            alpha = 0.5, size = 0.5) +
  geom_line(data = data_means, aes(x = Visit, y = mean, group = split_0, color = split_0), size = 2) +
  geom_errorbar(data = data_means, aes(x = Visit, y = mean, ymin = mean - se, ymax = mean + se, color = split_0), 
                width = 0.2, size = 1) +
  facet_wrap(~ split_0) +
  labs(y = 'Go-NoGo Performance Change', color = 'Biotype', x = '') +
  theme_minimal() +
  theme(legend.position = "top", text = element_text(size = 24)) +
  scale_color_manual(values = c("#B75A65", "gray40")) +
  scale_fill_manual(values = c("#B75A65", "gray40")) 
dev.off()


#### Change in QIDS

# Calculate means and standard errors
data_means  =  temp %>%
  group_by(Visit, split_0) %>%
  summarise(mean = mean(q_total_diff_bl, na.rm = TRUE),
            se = sd(q_total_diff_bl, na.rm = TRUE) / sqrt(n()), 
            .groups = "drop")

# Line plot
png(filename="/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/qids_raw_diff.png", width = 800, height = 600)
ggplot(data_means, aes(x = Visit, y = mean, color = split_0)) +
  geom_line(aes(group = split_0), size = 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, size=1) +
  labs(y = "QIDS Total Change", x = "") +
  theme_minimal() +
  theme(text = element_text(size = 24))+
  theme(legend.position = "top") +
  scale_color_manual(name = "Biotype", values = c("#B75A65", "gray40"))
dev.off()

# Dots and lines plot
png(filename="/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/qids_raw_diff_mean_dots.png", width = 1000, height = 600)
ggplot() +
  geom_jitter(data = data_adjusted_combat_temp, aes(x = Visit, y = q_total_diff_bl, color = split_0), 
              size = 2, alpha = 0.5, position = position_dodge(0.8)) +
  geom_line(data = data_adjusted_combat_temp, aes(x = Visit, y = q_total_diff_bl, group = interaction(Subjects, split_0), color = split_0), 
            alpha = 0.5, size = 0.5) +
  geom_line(data = data_means, aes(x = Visit, y = mean, group = split_0, color = split_0), size = 2) +
  geom_errorbar(data = data_means, aes(x = Visit, y = mean, ymin = mean - se, ymax = mean + se, color = split_0), 
                width = 0.2, size = 1) +
  facet_wrap(~ split_0) +
  labs(y = 'QIDS Total Change', color = 'Biotype', x = '') +
  theme_minimal() +
  theme(legend.position = "top", text = element_text(size = 24)) +
  scale_color_manual(values = c("#B75A65", "gray40")) +
  scale_fill_manual(values = c("#B75A65", "gray40")) 
dev.off()


#### Plots for each QIDS symptom ####

library(tidyverse)
library(scales)

# Rename the variables
data_adjusted_combat  =  rename(data_adjusted_combat, QIDS01 = q1, QIDS02 = q2, QIDS03 = q3, QIDS04 = q4,
                                QIDS05 = q5, QIDS06 = q6, QIDS07 = q7, QIDS08 = q8, QIDS09 = q9, QIDS10 = q10,
                                QIDS11 = q11, QIDS12 = q12, QIDS13 = q13, QIDS14 = q14, QIDS15 = q15, QIDS16 = q16)

# Reshape the data into a long format
long_data  =  data_adjusted_combat %>%
  pivot_longer(
    cols = c("QIDS01", "QIDS02", "QIDS03", "QIDS04", "QIDS05", "QIDS06", "QIDS07", "QIDS08", 
             "QIDS09", "QIDS10", "QIDS11", "QIDS12", "QIDS13", "QIDS14", "QIDS15", "QIDS16"), 
    names_to = "q_var", 
    values_to = "q_value"
  ) %>%
  group_by(q_var, Visit, split_0)

# Calculate means and standard errors
data_means  =  long_data %>%
  group_by(q_var, Visit, split_0) %>%
  summarise(mean = mean(q_value, na.rm = TRUE),
            se = sd(q_value, na.rm = TRUE) / sqrt(n()), 
            .groups = "drop")

# Create the plot
png(file=paste("/Users/ltozzi/My Drive (ltozzi@stanford.edu)/Projects/tms_gonogo_project/plots/qids_indiv.png", sep=''),width=800, height=800)
ggplot(long_data, aes(x = Visit, y = mean, color = split_0, group = split_0)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, size = 1) +
  facet_wrap(~ q_var, scales = 'free_y') +
  ylim(c(0, 3)) +
  labs(y = 'QIDS Item Response', color = 'Biotype', x = '') +
  theme_minimal() +
  theme(legend.position = "top",
        text = element_text(size = 24),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + # Adjust x-axis labels
  scale_color_manual(values = c("#B75A65", "gray40"))
dev.off()
