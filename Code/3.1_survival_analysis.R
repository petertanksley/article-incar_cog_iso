source("0_packages.R")

#color palettes






if(!file.exists("hrs_surv_ind.rds") | !file.exists("hrs_surv_dep.rds")){
  hrs_full <- import("hrs_full_analytic.rds") 
  
  #=Data set up for survival analysis==========================================
  
  #make survival model recodes (age as time-scale approach)
  hrs_surv <- hrs_full %>% 
    group_by(hhidpn) %>% 
    fill(dod_yr, .direction="updown") %>% 
    ungroup() %>% 
    mutate(study_age = year-birthyr,
           firstiw_age = firstiw-birthyr,
           dod_age = as_numeric(dod_yr)-birthyr) %>% 
    group_by(hhidpn) %>% 
    mutate(cog_first = ifelse(cog_2cat_num==1, study_age, NA),
           cog_first = pmin(cog_first)) %>%
    fill(cog_first, .direction="updown") %>% 
    mutate(cog_ever = max(cog_2cat_num),
           cog_surv_age = ifelse(cog_ever>0, cog_first, max(study_age))) %>% 
    ungroup() 
    # filter(age>=50) #removed 2,032 rows (4%), 53,313 rows remaining
  
  
  #time-independent 
  hrs_surv_ind <- hrs_surv %>% 
    distinct(hhidpn, .keep_all=TRUE) %>% #removed 48,404 rows (87%), 6,951 rows remaining
    select(-c(year, stroke_ever, age, cog_2cat, cog_2cat_num, cogfunction, edu_yrs)) %>% 
    filter(cog_surv_age>firstiw_age) #removed 290 rows (4%), 6,661 rows remaining
  
  #format time-independent variables
  hrs_surv_ind_fmt <- tmerge(data1 = hrs_surv_ind,
                             data2 = hrs_surv_ind,
                             id=hhidpn,
                             event=event(cog_surv_age, cog_ever))
  
  #time-dependent
  hrs_surv_dep <- hrs_surv %>% 
    distinct(hhidpn, study_age, .keep_all=TRUE) %>% 
    select(hhidpn, study_age, stroke_ever)
  
  #merge
  hrs_surv_final <- tmerge(data1 = hrs_surv_ind_fmt,
                           data2 = hrs_surv_dep,
                           id=hhidpn,
                           stroke=tdc(study_age, stroke_ever)) %>% 
    filter(tstart>0) #removed 6,661 rows (14%), 42,015 rows remaining
  
  export(hrs_surv_final, "hrs_surv_dep.rds")
  export(hrs_surv_ind_fmt, "hrs_surv_ind.rds")
} else {
  
  hrs_surv_ind<- import("hrs_surv_ind.rds")
  hrs_surv_dep<- import("hrs_surv_dep.rds")
}

#=Fit Cox model================================================================

#INCARCERATION
#time-dependent covariates
cox1 <- coxph(Surv(tstart, tstop, event) ~ factor(incar_ever) +
                factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) + 
                factor(smoke_ever) + factor(stroke) +
                strata(study), 
              data=hrs_surv_dep, id=hhidpn)
cox1_res <- tidy(cox1, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "incar_ever")
#extract estimates for visuals
cox1_res_incar <- cox1_res %>% 
  filter(term=="incar_everIncarcerated") %>% 
  pull(estimate) 

#APOE4
#time-dependent covariates
cox2 <- coxph(Surv(tstart, tstop, event) ~ factor(apoe_info99_4ct) +
                factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) + 
                factor(smoke_ever) + factor(stroke) +
                strata(study), 
              data=hrs_surv_dep, id=hhidpn)
cox2_res <- tidy(cox2, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "apoe_4")
#extract estimates for visuals
cox2_res_apoe1 <- cox2_res %>% 
  filter(term=="factor(apoe_info99_4ct)one copy") %>% 
  pull(estimate) 
cox2_res_apoe2 <- cox2_res %>% 
  filter(term=="factor(apoe_info99_4ct)two copies") %>% 
  pull(estimate) 

#INCARCERATION + APOE4
#time-dependent covariates
cox3 <- coxph(Surv(tstart, tstop, event) ~ factor(incar_ever) + factor(apoe_info99_4ct) +
                factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) + 
                factor(smoke_ever) + factor(stroke) +
                strata(study), 
              data=hrs_surv_dep, id=hhidpn)
cox3_res <- tidy(cox3, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "incar_apoe_4")


#INCARCERATION x APOE4
#time-dependent covariates
cox4 <- coxph(Surv(tstart, tstop, event) ~ factor(incar_ever)*factor(apoe_info99_4ct) +
                factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) + 
                factor(smoke_ever) + factor(stroke) +
                strata(study), 
              data=hrs_surv_dep, id=hhidpn)
cox4_res <- tidy(cox4, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "incar_x_apoe_4") 

#bind results
cox_all_res <- bind_rows(cox1_res,
                         cox2_res,
                         cox3_res,
                         cox4_res)

export(c("cox1", "cox2", "cox3", "cox4"), "../output/results/main_results_surv_models.rdata")
export(cox_all_res, "../output/results/main_results_surv.csv")

#=Fit and visualize basic survival curves======================================

#INCARCERATION
surv_ind_incar <- survfit(Surv(study_age, event) ~ incar_ever, data=hrs_surv_ind)
logrank_incar <- survdiff(Surv(study_age, event) ~ incar_ever, data=hrs_surv_ind) %>% 
  glance() %>% 
  pull(statistic)

survplot1 <- ggsurvplot(surv_ind_incar,
                        conf.int = TRUE,
                        pval = TRUE, 
                        pval.method = TRUE,
                        risk.table = "nrisk_cumevents",
                        fontsize=6,
                        cumevents = TRUE,
                        xlim=c(50,100),
                        break.x.by=10,
                        # surv.median.line = "hv",
                        censor.shape=NA,
                        palette = c("darkblue", "darkred"),
                        cumcensor = TRUE) 

#median survival age================================#
closest<-function(var,val){
  var[which(abs(var-val)==min(abs(var-val)))] 
}

median_surv_nojail <- survplot1$data.survplot %>% 
  filter(incar_ever=="Not Incarcerated") %>% 
  filter(surv==closest(surv, .5)) %>% 
  pull(time)

median_surv_jail <- survplot1$data.survplot %>% 
  filter(incar_ever=="Incarcerated") %>% 
  filter(surv==closest(surv, .5)) %>% 
  pull(time)
#==================================================#


logrank_incar <- expression(paste("Log-rank: ", chi^2, "(1)=189.1; ", italic("P"), "<0.001"))
cox_incar <- expression(paste("Cox: HR=1.33; ", italic("P"), "<0.001"))

plot1 <- survplot1$plot +
  labs(y="Survival probability\n(no cognitive impairment)",
       x="Age",
       tag = "B.") +
  theme(legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 22),
        legend.direction = "vertical",
        legend.position = c(.85, .85),
        legend.key.size = unit(1, "cm"),
        axis.title.x.bottom = element_text(size=22, face = "bold"),
        axis.title.y.left = element_text(size=22, face = "bold", vjust = -15),
        axis.text.y.left = element_text(size = 18),
        axis.text.x.bottom = element_text(size=18)) +
  geom_segment(aes(x=45, xend=79, y=.5, yend=.5), linewidth=1, linetype="dashed") +
  geom_segment(aes(x=median_surv_nojail, xend=median_surv_nojail, y=0, yend=.5),  linewidth=1, linetype="dashed") +
  geom_segment(aes(x=median_surv_jail, xend=median_surv_jail, y=0, yend=.5),  linewidth=1, linetype="dashed") +
  scale_y_continuous(expand = c(0,0)) +
  # annotate("text", x=60, y=.25, size=7, fontface="bold", label=expression(atop(textstyle("Log-rank"), paste(italic("p"), "<0.0001")))) +
  annotate("text", x=48, y=.1, size=7, fontface="bold", label=logrank_incar, hjust=0) +
  annotate("text", x=48, y=.05, size=7, fontface="bold", label=cox_incar, hjust=0) +
  scale_fill_manual(name="Lifetime incarceration",
                    values = c("darkblue", "darkred"),
                    labels = c("Never-incarcerated", "Incarcerated")) +
  scale_color_manual(name="Lifetime incarceration",
                     values = c("darkblue", "darkred"),
                    labels = c("Never-incarcerated", "Incarcerated")) 
plot1

tab1 <- survplot1$table +
  labs(y="",
       x="") +
  scale_y_discrete(labels=c("Incarcerated", "Never-\nincarcerated")) +
  theme(axis.text.y.left = element_text(size = 22, face = "bold", hjust = .5),
        plot.title       = element_text(size = 22, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y.left = element_line(colour = NA),
        panel.background = element_rect(color = "black", linewidth = 1))


  
survplot_incar <- plot1 / tab1 + plot_layout(heights = c(3,1)) 
survplot_incar
ggsave("../output/figures/survplot_incar.png", survplot_incar, width = 15, height = 8)


#APOE-4
surv_ind_apoe4 <- survfit(Surv(study_age, event) ~ apoe_info99_4ct, data=hrs_surv_ind)
logrank_apoe <- survdiff(Surv(study_age, event) ~ apoe_info99_4ct, data=hrs_surv_ind) %>% 
  glance() %>% 
  pull(statistic)

survplot2 <- ggsurvplot(surv_ind_apoe4,
                        conf.int = TRUE,
                        pval = TRUE, 
                        pval.method = TRUE,
                        risk.table = "nrisk_cumevents",
                        fontsize=6,
                        cumevents = TRUE,
                        xlim=c(50,100),
                        break.x.by=10,
                        # surv.median.line = "hv",
                        censor.shape=NA,
                        palette = c("darkblue", "darkslateblue", "darkred"),
                        cumcensor = TRUE) 

#median survival age================================#
closest<-function(var,val){
  var[which(abs(var-val)==min(abs(var-val)))] 
}

median_surv_zero <- survplot2$data.survplot %>% 
  filter(apoe_info99_4ct=="zero copies") %>% 
  filter(surv==closest(surv, .5)) %>% 
  pull(time)

median_surv_one <- survplot2$data.survplot %>% 
  filter(apoe_info99_4ct=="one copy") %>% 
  filter(surv==closest(surv, .5)) %>% 
  pull(time)

median_surv_two <- survplot2$data.survplot %>% 
  filter(apoe_info99_4ct=="two copies") %>% 
  filter(surv==closest(surv, .5)) %>% 
  pull(time)
#==================================================#

logrank_apoe <- expression(paste("Log-rank: ", chi^2, "(2)=45.3; ", italic("P"), "<0.001"))
cox_apoe1 <- expression(paste("Cox: One copy HR=1.19; ", italic("P"), "<0.001"))
cox_apoe2 <- expression(paste("Cox: Two copies HR=1.56; ", italic("P"), "<0.001"))

plot2 <- survplot2$plot +
  labs(y="Survival probability\n(no cognitive impairment)",
       x="Age",
       tag = "A.") +
  theme(legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 22),
        legend.direction = "vertical",
        legend.position = c(.85, .85),
        legend.key.size = unit(1, "cm"),
        axis.title.x.bottom = element_text(size=22, face = "bold"),
        axis.title.y.left = element_text(size=22, face = "bold", vjust = -15),
        axis.text.y.left   = element_text(size = 18),
        axis.text.x.bottom = element_text(size = 18)) +
  geom_segment(aes(x=45, xend=79, y=.5, yend=.5), linewidth=1, linetype="dashed") +
  geom_segment(aes(x=median_surv_zero, xend=median_surv_zero, y=0, yend=.5),  linewidth=1, linetype="dashed") +
  geom_segment(aes(x=median_surv_one, xend=median_surv_one, y=0, yend=.5),  linewidth=1, linetype="dashed") +
  geom_segment(aes(x=median_surv_two, xend=median_surv_two, y=0, yend=.5),  linewidth=1, linetype="dashed") +
  scale_y_continuous(expand = c(0,0)) +
  annotate("text", x=48, y=.15, size=7, fontface="bold", hjust=0, label=logrank_apoe) +
  annotate("text", x=48, y=.1, size=7, fontface="bold", hjust=0, label=cox_apoe1) +
  annotate("text", x=48, y=.05, size=7, fontface="bold", hjust=0, label=cox_apoe2) +
  scale_fill_manual(name=expression(paste(bolditalic("APOE-\u03b54"), bold(" genotype"))),
                    values = c("darkblue", "darkslateblue", "darkred"),
                    labels = c("Zero copies", "One copy", "Two copies")) +
  scale_color_manual(name=expression(paste(bolditalic("APOE-\u03b54"), bold(" genotype"))),
                     values = c("darkblue",  "darkslateblue", "darkred"),
                     labels = c("Zero copies", "One copy", "Two copies"))
plot2

tab2 <- survplot2$table +
  labs(y="",
       x="") +
  scale_y_discrete(labels=c("Two copies", "One copy", "Zero copies")) +
  theme(axis.title.y.left = element_blank(),
        axis.text.y.left  = element_text(size = 22, face = "bold", hjust = .5),
        plot.title        = element_text(size = 22, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y.left = element_line(colour = NA),
        panel.background = element_rect(color = "black", linewidth = 1))

survplot_apoe4 <- plot2 / tab2 + plot_layout(heights = c(3,1)) 
# survplot_apoe4
ggsave("../output/figures/survplot_apoe4.png", survplot_apoe4, width = 15, height = 8)

#=Combine plots========================================#

# design1="
# AAACCC
# AAACCC
# AAACCC
# AAACCC
# AAACCC
# BBBDDD
# "
# 
# combined1 <- plot1 + tab1 + plot2 + tab2 + 
#   plot_layout(design = design1, tag_level = "new") + 
#   plot_annotation(tag_levels = "A", tag_suffix=".") &
#   theme(plot.tag = element_text(size = 24, face = "bold"))
# combined1
# ggsave("../output/figures/survplot_combined1.png", combined1, width = 28, height =10)

design2="
AAA
AAA
AAA
AAA
AAA
BBB
CCC
CCC
CCC
CCC
CCC
DDD"

combined2 <- plot2 + tab2 + plot1 + tab1 + 
  plot_layout(design = design2, tag_level = "new") & 
  # plot_annotation(tag_levels = "A", tag_suffix=".") 
  theme(plot.tag = element_text(size = 24, face = "bold"))
# combined2
ggsave("../output/figures/survplot_combined2.tiff", combined2, width = 16, height =20, dpi = 700)

#=test of the proportionality assumption=====================================

#examine Schoenfeld residuals (looking for no time trends)
#visual examination will be pursued due to size of data (too easy to reject null)
cox.zph(cox1) %>% ggcoxzph()
cox.zph(cox2) %>% ggcoxzph()
cox.zph(cox3) %>% ggcoxzph()
cox.zph(cox4) %>% ggcoxzph()




