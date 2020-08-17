setwd('/home/chandan/Downloads/')

library("survival")
library("survminer")



#loading the dataset
clinical <- read.csv('/home/chandan/Downloads/outp.csv')



#categorising the methylation columns

clinical$cg23527067 <- ifelse(clinical$cg23527067 > 0.8, 1, 0)
clinical$cg16977035 <- ifelse(clinical$cg16977035 > 0.8, 1, 0)
clinical$cg23412777 <- ifelse(clinical$cg23412777 > 0.8, 1, 0)
#changing the dead and alive status
clinical$vital_stat <- ifelse(clinical$vital_stat %in% c('Dead'), 1, 0)


ind_keep <- grep('days_to_death',colnames(clinical))
death <- as.matrix(clinical[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}

ind_keep <- grep('days_to_last_followup',colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if ( sum (is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}

all_clin <- data.frame(death_collapsed,fl_collapsed)
colnames(all_clin) <- c('death_days', 'followUp_days')

all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}

clinical$new_death <- all_clin$new_death

fit <- survfit(Surv(clinical$new_death, clinical$vital_stat) ~ clinical$cg23527067 + clinical$cg16977035 + clinical$cg23412777, data = clinical)

ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 500,     # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Non Methylated", "Methylated"),    # change legend labels.
  palette = 
    c("#E7B800", "#2E9FDF"), # custom color palettes.
  title = "cg23527067",
)

ggsurvplot(fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           title = "cg23527067",
           legend.labs = 
             c("Non Methylated", "Methylated"),  
           fun = "cumhaz")

ggsurv <- ggsurvplot(fit, conf.int = TRUE,
                     ggtheme = theme_bw(), fun = "event")

ggsurv$plot +theme_bw() + 
  theme (legend.position = "right")+
  facet_grid(cg16977035 ~ cg23412777)



surv_diff <- survdiff(Surv(clinical$new_death, clinical$vital_stat) ~ clinical$cg23527067, data = clinical)

