source("0_packages.R")

hrs <- import("hrs_full_analytic.rds")

#get case/observation numbers
cases <- n_distinct(hrs$hhidpn)
obs <- nrow(hrs)

#set table theme
set_gtsummary_theme(theme_gtsummary_compact())

#=Table 1 - Descriptives=======================================================

#generate descriptives table
hrs_labs <- hrs %>% 
  select(hhidpn,
         study, 
         sex, 
         age,
         year,
         race_ethn,
         edu,
         cog_2cat,
         incar_ever, 
         lone_scale_pro,
         soc_iso_index_pro,
         stroke_ever,
         smoke_ever,
         apoe_info99_4ct,
         social_origins,
  ) %>% 
  mutate(incar_ever = fct_recode(incar_ever, 
                                 "Yes" = "Incarcerated", 
                                 "No" = "Not Incarcerated"),
         smoke_ever = fct_recode(as.factor(smoke_ever),
                                 "Yes" = "1",
                                 "No" = "0"),
         smoke_ever = fct_relevel(smoke_ever, "Yes"),
         edu = fct_recode(edu, 
                          "High school or more" = "hs or more",
                          "Less than high school" = "less than hs"),
         apoe_info99_4ct = fct_recode(apoe_info99_4ct,
                                      "Zero copies" = "zero copies",
                                      "One copy" = "one copy", 
                                      "Two copies" = "two copies"),
         stroke_ever=fct_recode(as_factor(stroke_ever), "No"="0", "Yes"="1"),
         stroke_ever = fct_relevel(stroke_ever, "Yes")) %>% 
  mutate(study = fct_drop(study)) %>%
  mutate(social_origins = as_numeric(social_origins)) %>% 
  var_labels(study           = "HRS cohort", 
             sex             = "Self-reported sex", 
             race_ethn       = "Race/Ethnicity",
             edu             = "HS completion",
             age             = "Age",
             year            = "Study year",
             cog_2cat        = "Cognitive function",
             incar_ever      = "Lifetime incarceration",
             lone_scale_pro  = "Loneliness",
             soc_iso_index_pro = "Social isolation",
             stroke_ever     = "History of stroke",
             smoke_ever      = "Ever smoker",
             apoe_info99_4ct = "APOE-4 count",
             social_origins  = "Social origins index"
  ) 

#time-independent variables
tab1_top <- hrs_labs %>% 
  select(-c(age, year, cog_2cat, stroke_ever, lone_scale_pro, soc_iso_index_pro)) %>% 
  # mutate(ad_pgs_blk = ifelse(race_ethn=="Black", ad_pgs, NA),
  #        ad_pgs_wht = ifelse(race_ethn=="White", ad_pgs, NA)) %>% 
  # var_labels(ad_pgs_blk = "Polygenic Index for AD (Black)",
  #            ad_pgs_wht = "Polygenic Index for AD (White)") %>% 
  distinct(hhidpn, .keep_all=TRUE) %>% 
  select(-c(hhidpn, 
            # ad_pgs
            )) %>% 
  tbl_summary(by=incar_ever,
              type = list(c(#ad_pgs_blk, ad_pgs_wht, 
                            social_origins) ~ "continuous",
                          c(study, sex, race_ethn, edu, apoe_info99_4ct, smoke_ever) ~ "categorical"),
              statistic = all_continuous() ~ "{mean} ({sd})") %>% 
  add_p() %>% 
  bold_p() %>% 
  add_overall() %>% 
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Ever incarcerated?**") 
  # remove_row_type(variables = c(ad_pgs_blk, ad_pgs_wht), type = "missing")

#time-dependent variables
tab1_bottom <- hrs_labs %>%
  select(hhidpn, age, cog_2cat, stroke_ever, lone_scale_pro, soc_iso_index_pro, year, incar_ever) %>%
  mutate(year=as_factor(year)) %>%
  group_by(hhidpn) %>% 
  mutate(age=mean(as.double(age))) %>% 
  var_labels(age = "Age (mean)") %>% 
  ungroup() %>% 
  tbl_summary(by=incar_ever,
              type = list(age ~ "continuous",
                          stroke_ever ~ "categorical",
                          lone_scale_pro ~ "continuous", 
                          soc_iso_index_pro ~ "continuous"),
              statistic = all_continuous() ~ "{mean} ({sd})",
              include = -hhidpn,
              digits = all_continuous() ~ 2) %>%
  add_p(include=-c(stroke_ever, cog_2cat)) %>%
  bold_p() %>% 
  add_overall() %>%
  modify_header(label ~ "**Variable**") %>% 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Ever incarcerated?**") 
  

#stack tables
tab1_combined <- tbl_stack(tbls = list(tab1_top, tab1_bottom),
          group_header = c(glue("Individuals (N = {style_number(cases,big.mark=',')})"), 
                           glue("Observations (N = {style_number(obs,big.mark=',')})"))) %>% 
  modify_caption('**Table 1**. Individual and observation-level descriptive statistics of the Health and Retirement Study.')
tab1_combined


#get pvalues for cognitive status and history of stroke
cog_pval <- glmer(cog_2cat_num ~ factor(incar_ever) + (1|hhidpn), 
      data=hrs, 
      family=binomial, 
      control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) %>% 
  tidy() %>% 
  filter(term=="factor(incar_ever)Incarcerated") %>% 
  pull(p.value)
# cog_pval
# [1] 1.274484e-05

strok_pval <- glmer(stroke_ever ~ factor(incar_ever) + (1|hhidpn), 
                  data=hrs, 
                  family=binomial, 
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) %>% 
  tidy() %>% 
  filter(term=="factor(incar_ever)Incarcerated") %>% 
  pull(p.value)
# strok_pval
# [1] 0.3086792



lone_pval <- lmerTest::lmer(lone_scale_pro ~ factor(incar_ever) + (1|hhidpn), 
                    data=hrs) %>% 
  tidy() %>%
  filter(term=="factor(incar_ever)Incarcerated") %>%
  pull(p.value)
# lone_pval
# [1] 1.572525e-05

soc_iso_pval <- lmerTest::lmer(soc_iso_index_pro ~ factor(incar_ever) + (1|hhidpn), 
                               data=hrs) %>% 
  tidy() %>%
  filter(term=="factor(incar_ever)Incarcerated") %>%
  pull(p.value)
# lone_pval
# [1] 1.572525e-05

#export .html table
gtsave(as_gt(tab1_combined), "../output/results/tab1_descriptives.html")

# #export .docx table
# tab1_combined %>%
#   as_flex_table() %>%
#   flextable::save_as_docx(path="../output/results/tab1_descriptives.docx")


#=Table 2 - Main results=======================================================

main_results <- import_list("../output/results/tab3_main_results.rdata")
list2env(main_results, .GlobalEnv)
rm(main_results)

#function to calculate delta r2
delta_r2 <- function(model, base=m0){
  r2_0 <- r2(m0) %>% unlist()
  r2_plus1 <- r2(model) %>% unlist()
  delta_r2 <- r2_plus1[2] - r2_0[2]
  print(delta_r2)
}


tab21 <- tbl_regression(m11, exponentiate=TRUE,
                        include=c("factor(apoe_info99_4ct)"),
                        label=c("factor(apoe_info99_4ct)"="APOE-4 allele count"))
tab22 <- tbl_regression(m12, exponentiate=TRUE,
                        include=c("factor(incar_ever)"),
                        label=c("factor(incar_ever)"="Lifetime incarceration"))
tab23 <- tbl_regression(m13, exponentiate=TRUE,
                        include=c("factor(incar_ever)","factor(apoe_info99_4ct)"),
                        label=list("factor(apoe_info99_4ct)"="APOE-4 allele count",
                                   "factor(incar_ever)"="Lifetime incarceration"))
tab24 <- tbl_regression(m14, exponentiate=TRUE,
                        include=c("factor(incar_ever)","factor(apoe_info99_4ct)",
                                  "factor(incar_ever):factor(apoe_info99_4ct)"),
                        label=list("factor(apoe_info99_4ct)"="APOE-4 allele count",
                                   "factor(incar_ever)"="Lifetime incarceration",
                                   "factor(incar_ever):factor(apoe_info99_4ct)"="Lifetime Incarceration*APOE-4 allele count"))

tab2_mods <- list(tab21,
                  tab22,
                  tab23,
                  tab24)
#write function to update gtsumamry tables
tab_updates <- function(x){
  x %>%
    add_significance_stars(hide_ci = FALSE, hide_p = TRUE, hide_se = TRUE) %>%
    remove_row_type(type = "reference") %>%
    modify_table_body(~ .x %>%
                        dplyr::mutate(ci = ifelse(!is.na(ci), paste0("[", ci, "]"), "")))
}

tab1_mods_update <- lapply(tab2_mods, tab_updates)

tab2_all <- tbl_merge(tab1_mods_update, tab_spanner = paste("Model 2.", 1:4, sep = "")) %>%
  modify_header(label = "**Variable**") %>%
  modify_caption(glue('**Table 2**. Mixed effect Poisson regression of cognitive status on APOE-4 genotype and lifetime incarceration')) %>%
  modify_footnote(label ~ "All models also adjusted for age, sex, race/ethnicity, high school completion, smoking history, stroke history, social origins index, and HRS cohort.")
tab2_all

tab2_all %>%
  as_gt() %>%
  gtsave("../output/results/tab2_main_res.html")

# tab2_all %>%
#   as_flex_table() %>%
#   flextable::save_as_docx(path="../output/results/tab2_main_res.docx")

#=Table S1 - polygenic risk (EUR)==============================================

hrs_pgs_eur <- import("hrs_pgs_analytic.rds") %>% filter(race_ethn=="White" & !is.na(ad_pgs))

#get case/observation numbers
cases_eur <- n_distinct(hrs_pgs_eur$hhidpn)
obs_eur <- nrow(hrs_pgs_eur)


poly_eur_results <- import_list("../output/results/tab_s1_eur.rdata")
list2env(poly_eur_results, .GlobalEnv)

#function to calculate delta r2
delta_r2 <- function(model, base=m0){
  r2_0 <- r2(m0) %>% unlist()
  r2_plus1 <- r2(model) %>% unlist()
  delta_r2 <- r2_plus1[2] - r2_0[2]
  print(delta_r2)
}


tab_s1 <- tbl_regression(m21, exponentiate=TRUE,
                         include=c("scale(ad_pgs_resid)"),
                         label=c("scale(ad_pgs_resid)"="PGS AD"))
tab_s2 <- tbl_regression(m22, exponentiate=TRUE,
                         include=c("factor(apoe_info99_4ct)"),
                         label=c("factor(apoe_info99_4ct)"="APOE-4 allele count"))
tab_s3 <- tbl_regression(m23, exponentiate=TRUE,
                         include=c("scale(ad_pgs_resid)","factor(apoe_info99_4ct)"),
                         label=c("scale(ad_pgs_resid)"="PGS AD",
                                 "factor(apoe_info99_4ct)"="APOE-4 allele count"))
tab_s4 <- tbl_regression(m24, exponentiate=TRUE,
                         include=c("factor(incar_ever)"),
                         label=c("factor(incar_ever)"="Lifetime incarceration"))
tab_s5 <- tbl_regression(m25, exponentiate=TRUE,
                         include=c("factor(incar_ever)","scale(ad_pgs_resid)","factor(apoe_info99_4ct)"),
                         label=list("scale(ad_pgs_resid)"="PGS AD",
                                    "factor(apoe_info99_4ct)"="APOE-4 allele count",
                                    "factor(incar_ever)"="Lifetime incarceration"))
tab_s6 <- tbl_regression(m26, exponentiate=TRUE,
                         include=c("factor(incar_ever)","scale(ad_pgs_resid)","factor(apoe_info99_4ct)",
                                   "factor(incar_ever):scale(ad_pgs_resid)",
                                   "factor(incar_ever):factor(apoe_info99_4ct)"),
                         label=list("scale(ad_pgs_resid)"="PGS AD",
                                    "factor(apoe_info99_4ct)"="APOE-4 allele count",
                                    "factor(incar_ever)"="Lifetime incarceration",
                                    "factor(incar_ever):factor(apoe_info99_4ct)"="PGS AD",
                                    "factor(incar_ever):factor(apoe_info99_4ct)"="Lifetime Incarceration*APOE-4 allele count"
                         ))

tab_s1_mods <- list(tab_s1,
                    tab_s2,
                    tab_s3,
                    tab_s4,
                    tab_s5,
                    tab_s6)

#write function to update gtsumamry tables
tab_updates <- function(x){
  x %>%
    add_significance_stars(hide_ci = FALSE, hide_p = TRUE, hide_se = TRUE) %>%
    remove_row_type(type = "reference") %>%
    modify_table_body(~ .x %>%
                        dplyr::mutate(ci = ifelse(!is.na(ci), paste0("[", ci, "]"), "")))
}

tab_s1_mods_update <- lapply(tab_s1_mods, tab_updates)

tab_s1_all <- tbl_merge(tab_s1_mods_update, tab_spanner = paste("Model S", 1:6, sep = "")) %>%
  modify_header(label = "**Variable**") %>%
  modify_caption(glue('**Table S1**. Mixed effect Poisson regression of cognitive status on mono-/polygenic risk and lifetime incarceration among White and Black HRS participants')) %>%
  modify_footnote(label ~ "All models also adjusted for age, sex, high school completion, smoking history, stroke history, social origins index, and HRS cohort.")
tab_s1_all

rm(list = setdiff(ls(), c("tab_s1_all", "obs_eur", "cases_eur")))

tab_s1_all %>%
  as_gt() %>%
  gtsave("../output/results/tab_s1_eur.html")

# tab_s1_all %>%
#   as_flex_table() %>%
#   flextable::save_as_docx(path="../output/results/tab_s1_eur.docx")



#=Table S1 - polygenic risk (AFR)==============================================

hrs_pgs_afr <- import("hrs_pgs_analytic.rds") %>% filter(race_ethn=="Black" & !is.na(ad_pgs))

#get case/observation numbers
cases_afr <- n_distinct(hrs_pgs_afr$hhidpn)
obs_afr <- nrow(hrs_pgs_afr)


poly_afr_results <- import_list("../output/results/tab_s1_afr.rdata")
list2env(poly_afr_results, .GlobalEnv)

#function to calculate delta r2
delta_r2 <- function(model, base=m0){
  r2_0 <- r2(m0) %>% unlist()
  r2_plus1 <- r2(model) %>% unlist()
  delta_r2 <- r2_plus1[2] - r2_0[2]
  print(delta_r2)
}


tab_s1 <- tbl_regression(m31, exponentiate=TRUE,
                         include=c("scale(ad_pgs_resid)"),
                         label=c("scale(ad_pgs_resid)"="PGS AD"))
tab_s2 <- tbl_regression(m32, exponentiate=TRUE,
                         include=c("factor(apoe_info99_4ct)"),
                         label=c("factor(apoe_info99_4ct)"="APOE-4 allele count"))
tab_s3 <- tbl_regression(m33, exponentiate=TRUE,
                         include=c("scale(ad_pgs_resid)","factor(apoe_info99_4ct)"),
                         label=c("scale(ad_pgs_resid)"="PGS AD",
                                 "factor(apoe_info99_4ct)"="APOE-4 allele count"))
tab_s4 <- tbl_regression(m34, exponentiate=TRUE,
                         include=c("factor(incar_ever)"),
                         label=c("factor(incar_ever)"="Lifetime incarceration"))
tab_s5 <- tbl_regression(m35, exponentiate=TRUE,
                         include=c("factor(incar_ever)","scale(ad_pgs_resid)","factor(apoe_info99_4ct)"),
                         label=list("scale(ad_pgs_resid)"="PGS AD",
                                    "factor(apoe_info99_4ct)"="APOE-4 allele count",
                                    "factor(incar_ever)"="Lifetime incarceration"))
tab_s6 <- tbl_regression(m36, exponentiate=TRUE,
                         include=c("factor(incar_ever)","scale(ad_pgs_resid)","factor(apoe_info99_4ct)",
                                   "factor(incar_ever):scale(ad_pgs_resid)",
                                   "factor(incar_ever):factor(apoe_info99_4ct)"),
                         label=list("scale(ad_pgs_resid)"="PGS AD",
                                    "factor(apoe_info99_4ct)"="APOE-4 allele count",
                                    "factor(incar_ever)"="Lifetime incarceration",
                                    "factor(incar_ever):factor(apoe_info99_4ct)"="PGS AD",
                                    "factor(incar_ever):factor(apoe_info99_4ct)"="Lifetime Incarceration*APOE-4 allele count"
                         ))

tab_s2_mods <- list(tab_s1,
                    tab_s2,
                    tab_s3,
                    tab_s4,
                    tab_s5,
                    tab_s6)

#write function to update gtsumamry tables
tab_updates <- function(x){
  x %>%
    add_significance_stars(hide_ci = FALSE, hide_p = TRUE, hide_se = TRUE) %>%
    remove_row_type(type = "reference") %>%
    modify_table_body(~ .x %>%
                        dplyr::mutate(ci = ifelse(!is.na(ci), paste0("[", ci, "]"), "")))
}

tab_s2_mods_update <- lapply(tab_s2_mods, tab_updates)

tab_s2_all <- tbl_merge(tab_s2_mods_update, tab_spanner = paste("Model S1.", 1:6, sep = ""))
tab_s2_all

tab_s2_all %>%
  as_gt() %>%
  gtsave("../output/results/tab_s2_afr.html")

# tab_s2_all %>%
#   as_flex_table() %>%
#   flextable::save_as_docx(path="../output/results/tab_s2_afr.docx")


#stack EUR/AFR tables
tab_s1_combined <- tbl_stack(list(tab_s1_all, tab_s2_all),
          group_header = c(glue('White genotyped sample (Observations={format(obs_eur, big.mark=",")}; Cases={format(cases_eur, big.mark=",")})'),
                           glue('Black genotyped sample (Observations={format(obs_afr, big.mark=",")}; Cases={format(cases_afr, big.mark=",")})')))

tab_s1_combined %>%
  as_gt() %>%
  gtsave("../output/results/tab_s1_combined.html")

# tab_s1_combined %>%
#   as_flex_table() %>%
#   flextable::save_as_docx(path="../output/results/tab_s1_combined.docx")

#=Table S2 - Survival results==================================================

# hrs <- import("hrs_full_analytic.rds")


surv_results <- import_list("../output/results/main_results_surv_models.rdata")
list2env(surv_results, .GlobalEnv)
rm(surv_results)

tab_s2_cox1 <- tbl_regression(cox2, exponentiate = TRUE,
                              include=c("factor(apoe_info99_4ct)one copy",
                                        "factor(apoe_info99_4ct)two copies"),
                              label = list("factor(apoe_info99_4ct)one copy" = "One APOE-4 allele",
                                           "factor(apoe_info99_4ct)two copies" = "Two APOE-4 alleles"))
tab_s2_cox2 <- tbl_regression(cox1, exponentiate = TRUE,
               include = "factor(incar_ever)Incarcerated",
               label = list("factor(incar_ever)Incarcerated" = "Lifetime Incarceration"))
tab_s2_cox3 <- tbl_regression(cox3, exponentiate = TRUE,
                              include=c("factor(apoe_info99_4ct)one copy",
                                        "factor(apoe_info99_4ct)two copies",
                                        "factor(incar_ever)Incarcerated"),
                              label = list("factor(apoe_info99_4ct)one copy" = "One APOE-4 allele",
                                           "factor(apoe_info99_4ct)two copies" = "Two APOE-4 alleles",
                                           "factor(incar_ever)Incarcerated" = "Lifetime Incarceration"))
tab_s2_cox4 <- tbl_regression(cox4, exponentiate = TRUE,
                              include=c("factor(apoe_info99_4ct)one copy",
                                        "factor(apoe_info99_4ct)two copies",
                                        "factor(incar_ever)Incarcerated",
                                        "factor(incar_ever)Incarcerated:factor(apoe_info99_4ct)one copy",
                                        "factor(incar_ever)Incarcerated:factor(apoe_info99_4ct)two copies"),
                              label=list("factor(apoe_info99_4ct)one copy" = "One APOE-4 allele",
                                         "factor(apoe_info99_4ct)two copies" = "Two APOE-4 alleles",
                                         "factor(incar_ever)Incarcerated" = "Lifetime Incarceration",
                                         "factor(incar_ever)Incarcerated:factor(apoe_info99_4ct)one copy" = 
                                           "Lifetime Incarceration*One APOE-4 allele",
                                         "factor(incar_ever)Incarcerated:factor(apoe_info99_4ct)two copies" = 
                                           "Lifetime Incarceration*Two APOE-4 alleles"))
#group summary tables
tab_s2_mods <- list(tab_s2_cox1,
                    tab_s2_cox2,
                    tab_s2_cox3,
                    tab_s2_cox4)

#function to update tables
tab_updates <- function(x){
  x %>%
    add_significance_stars(hide_ci = FALSE, hide_p = TRUE, hide_se = TRUE) %>%
    remove_row_type(type = "reference") %>%
    modify_table_body(~ .x %>%
                        dplyr::mutate(ci = ifelse(!is.na(ci), paste0("[", ci, "]"), "")))
}

#apply updates
tab_s2_mods_update <- lapply(tab_s2_mods, tab_updates)

#merge models
tab_s2_all <- tbl_merge(tab_s2_mods_update, tab_spanner = paste("Model S1.", 1:4, sep = "")) %>% 
  modify_caption(glue("**Table S1.** Cox proportional hazard model of lifetime incarceration and APOE-4 genotype on first cognitive impairment.")) %>% 
  modify_footnote(label ~ "All models adjusted for sex, race/ethnicity, high school completion, smoking history, stroke history, and social origins index. All models were stratified by HRS cohort.")
tab_s2_all


#save out table
tab_s2_all %>%
  as_gt() %>%
  gtsave("../output/results/tab_s2_all.html")

# tab_s1_combined %>%
#   as_flex_table() %>%
#   flextable::save_as_docx(path="../output/results/tab_s1_combined.docx")

