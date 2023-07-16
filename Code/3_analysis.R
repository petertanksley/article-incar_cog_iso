source("0_packages.R")

hrs_full <- import("hrs_full_analytic.rds") 
hrs_pgs_eur <- import("hrs_pgs_analytic.rds") %>% filter(race_ethn=="White")
hrs_pgs_afr <- import("hrs_pgs_analytic.rds") %>% filter(race_ethn=="Black")

#=Main analysis================================================================

#set model formulas 
m10_base <- formula("cog_2cat_num ~ 
                   factor(sex) + factor(race_ethn) + scale(age) + factor(edu) + scale(social_origins) + factor(smoke_ever) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m11_apoe <- formula("cog_2cat_num ~ factor(apoe_info99_4ct) +
                   factor(sex) + factor(race_ethn) + scale(age) + factor(edu) + scale(social_origins) + factor(smoke_ever) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m12_incar <- formula("cog_2cat_num ~ factor(incar_ever) + 
                   factor(sex) + factor(race_ethn) + scale(age) + factor(edu) + scale(social_origins) + factor(smoke_ever) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m13_main_eff <- formula("cog_2cat_num ~ factor(incar_ever) + factor(apoe_info99_4ct) +
                   factor(sex) + factor(race_ethn) + scale(age) + factor(edu) + scale(social_origins) + factor(smoke_ever) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m14_int_eff <- formula("cog_2cat_num ~ factor(incar_ever)*factor(apoe_info99_4ct) +
                   factor(sex) + factor(race_ethn) + scale(age) + factor(edu) + scale(social_origins) + factor(smoke_ever) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")

#run models 
m10 <- glmer(m10_base,      data=hrs_full, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m11 <- glmer(m11_apoe,      data=hrs_full, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m12 <- glmer(m12_incar,     data=hrs_full, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m13 <- glmer(m13_main_eff,  data=hrs_full, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m14 <- glmer(m14_int_eff,   data=hrs_full, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

#check linear probabily model
# lmer(m14_int_eff,   data=hrs_full) %>% broom.mixed::tidy(conf.int=TRUE) %>% view()

rio::export(c("m10","m11","m12","m13","m14"),
       "../output/results/tab3_main_results.rdata")

#table models
res_10 <- tidy(m10, exponentiate=TRUE) %>% mutate(model="model base") %>% mutate(glance(m10)) 
res_11 <- tidy(m11, exponentiate=TRUE) %>% mutate(model="model 1")    %>% mutate(glance(m11)) 
res_12 <- tidy(m12, exponentiate=TRUE) %>% mutate(model="model 2")    %>% mutate(glance(m12))
res_13 <- tidy(m13, exponentiate=TRUE) %>% mutate(model="model 3")    %>% mutate(glance(m13))
res_14 <- tidy(m14, exponentiate=TRUE) %>% mutate(model="model 4")    %>% mutate(glance(m14))


main_results <- bind_rows(res_10,
                          res_11,
                          res_12,
                          res_13,
                          res_14)

#export results
rio::export(main_results, "../output/results/tab3_main_results.csv")

#=PGS analysis (EUR)================================================================

#set model formulas
m20_base <- formula("cog_2cat_num ~
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m21_pgs <- formula("cog_2cat_num ~ scale(ad_pgs_resid) +
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m22_apoe <- formula("cog_2cat_num ~ factor(apoe_info99_4ct) +
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m23_pgs_apoe <- formula("cog_2cat_num ~ scale(ad_pgs_resid) + factor(apoe_info99_4ct) +
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m24_incar <- formula("cog_2cat_num ~ factor(incar_ever) +
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m25_main_eff <- formula("cog_2cat_num ~ factor(incar_ever) + scale(ad_pgs_resid) + factor(apoe_info99_4ct) +
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m26_int_eff <- formula("cog_2cat_num ~ factor(incar_ever)*scale(ad_pgs_resid) + factor(incar_ever)*factor(apoe_info99_4ct) +
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")

#run models
m20 <- glmer(m20_base,     data=hrs_pgs_eur, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m21 <- glmer(m21_pgs,      data=hrs_pgs_eur, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m22 <- glmer(m22_apoe,     data=hrs_pgs_eur, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m23 <- glmer(m23_pgs_apoe, data=hrs_pgs_eur, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m24 <- glmer(m24_incar,    data=hrs_pgs_eur, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m25 <- glmer(m25_main_eff, data=hrs_pgs_eur, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m26 <- glmer(m26_int_eff,  data=hrs_pgs_eur, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

export(c("m20","m21","m22","m23","m24","m25","m26"),
            "../output/results/tab_s1_eur.rdata")


#table models
res_20 <- tidy(m20, exponentiate=TRUE) %>% mutate(model="model base") %>% mutate(glance(m20))
res_21 <- tidy(m21, exponentiate=TRUE) %>% mutate(model="model 1")    %>% mutate(glance(m21))
res_22 <- tidy(m22, exponentiate=TRUE) %>% mutate(model="model 2")    %>% mutate(glance(m22))
res_23 <- tidy(m23, exponentiate=TRUE) %>% mutate(model="model 3")    %>% mutate(glance(m23))
res_24 <- tidy(m24, exponentiate=TRUE) %>% mutate(model="model 4")    %>% mutate(glance(m24))
res_25 <- tidy(m25, exponentiate=TRUE) %>% mutate(model="model 5")    %>% mutate(glance(m25))
res_26 <- tidy(m26, exponentiate=TRUE) %>% mutate(model="model 6")    %>% mutate(glance(m26))

main_results2 <- bind_rows(res_20,
                          res_21,
                          res_22,
                          res_23,
                          res_24,
                          res_25,
                          res_26)

#export results
export(main_results2, "../output/results/tab_s1_eur.csv")

#=PGS analysis (AFR)================================================================

#set model formulas
m30_base <- formula("cog_2cat_num ~
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m31_pgs <- formula("cog_2cat_num ~ scale(ad_pgs_resid) +
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m32_apoe <- formula("cog_2cat_num ~ factor(apoe_info99_4ct) +
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m33_pgs_apoe <- formula("cog_2cat_num ~ scale(ad_pgs_resid) + factor(apoe_info99_4ct) +
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m34_incar <- formula("cog_2cat_num ~ factor(incar_ever) +
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m35_main_eff <- formula("cog_2cat_num ~ factor(incar_ever) + scale(ad_pgs_resid) + factor(apoe_info99_4ct) +
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")
m36_int_eff <- formula("cog_2cat_num ~ factor(incar_ever)*scale(ad_pgs_resid) + factor(incar_ever)*factor(apoe_info99_4ct) +
                   factor(sex) + scale(age) + scale(social_origins) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")

#run models
m30 <- glmer(m30_base,     data=hrs_pgs_afr, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m31 <- glmer(m31_pgs,      data=hrs_pgs_afr, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m32 <- glmer(m32_apoe,     data=hrs_pgs_afr, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m33 <- glmer(m33_pgs_apoe, data=hrs_pgs_afr, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m34 <- glmer(m34_incar,    data=hrs_pgs_afr, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m35 <- glmer(m35_main_eff, data=hrs_pgs_afr, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
m36 <- glmer(m36_int_eff,  data=hrs_pgs_afr, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

export(c("m30","m31","m32","m33","m34","m35","m36"),
       "../output/results/tab_s1_afr.rdata")


#table models
res_30 <- tidy(m30, exponentiate=TRUE) %>% mutate(model="model base") %>% mutate(glance(m30))
res_31 <- tidy(m31, exponentiate=TRUE) %>% mutate(model="model 1")    %>% mutate(glance(m31))
res_32 <- tidy(m32, exponentiate=TRUE) %>% mutate(model="model 2")    %>% mutate(glance(m32))
res_33 <- tidy(m33, exponentiate=TRUE) %>% mutate(model="model 3")    %>% mutate(glance(m33))
res_34 <- tidy(m34, exponentiate=TRUE) %>% mutate(model="model 4")    %>% mutate(glance(m34))
res_35 <- tidy(m35, exponentiate=TRUE) %>% mutate(model="model 5")    %>% mutate(glance(m35))
res_36 <- tidy(m36, exponentiate=TRUE) %>% mutate(model="model 6")    %>% mutate(glance(m36))

main_results3 <- bind_rows(res_30,
                           res_31,
                           res_32,
                           res_33,
                           res_34,
                           res_35,
                           res_36)

#export results
export(main_results3, "../output/results/tab_s1_afr.csv")

#=sub-group analysis===========================================================

#baseline model
m40_main_eff <- formula("cog_2cat_num ~ factor(incar_ever) + factor(apoe_info99_4ct) +
                   factor(sex) + factor(race_ethn) + scale(age) + factor(edu) + factor(social_origins) + factor(smoke_ever) + factor(stroke_ever) + factor(study) +
                   (1|hhidpn)")

# #Sex
# m41_incar_sex <- formula("cog_2cat_num ~ factor(incar_ever)*factor(sex) + factor(apoe_info99_4ct)+
#                     factor(race_ethn) + scale(age) + factor(edu) + scale(social_origins) + factor(smoke_ever) + factor(stroke_ever) + factor(study) +
#                    (1|hhidpn)")
# m42_apoe_sex <- formula("cog_2cat_num ~ factor(incar_ever) + factor(apoe_info99_4ct)*factor(sex) +
#                     factor(race_ethn) + scale(age) + factor(edu) + scale(social_origins) + factor(smoke_ever) + factor(stroke_ever) + factor(study) +
#                    (1|hhidpn)")
# 
# #Race/ethnicity
# m43_incar_race <- formula("cog_2cat_num ~ factor(incar_ever)*factor(race_ethn) + factor(apoe_info99_4ct) +
#                     factor(sex) + scale(age) + factor(edu) + scale(social_origins) + factor(smoke_ever) + factor(stroke_ever) + factor(study) +
#                    (1|hhidpn)")
# m44_apoe_race <- formula("cog_2cat_num ~ factor(incar_ever) + factor(apoe_info99_4ct)*factor(race_ethn) +
#                     factor(sex) + scale(age) + factor(edu) + scale(social_origins) + factor(smoke_ever) + factor(stroke_ever) + factor(study) +
#                    (1|hhidpn)")
# 
# #Edu
# m45_incar_edu <- formula("cog_2cat_num ~ factor(incar_ever)*factor(edu) + factor(apoe_info99_4ct) +
#                     factor(sex) + factor(race_ethn) + scale(age) + scale(social_origins) + factor(smoke_ever) + factor(stroke_ever) + factor(study) +
#                    (1|hhidpn)")
# m46_apoe_edu <- formula("cog_2cat_num ~ factor(incar_ever) + factor(apoe_info99_4ct)*factor(edu) +
#                     factor(sex) + factor(race_ethn) + scale(age) + scale(social_origins) + factor(smoke_ever) + factor(stroke_ever) + factor(study) +
#                    (1|hhidpn)")

#run models 
m40_main_eff   <- glmer(m40_main_eff,   data=hrs_full, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# m41_incar_sex  <- glmer(m41_incar_sex,  data=hrs_full, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# m42_apoe_sex   <- glmer(m42_apoe_sex,   data=hrs_full, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# m43_incar_race <- glmer(m43_incar_race, data=hrs_full, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# m44_apoe_race  <- glmer(m44_apoe_race,  data=hrs_full, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# m45_incar_edu  <- glmer(m45_incar_edu,  data=hrs_full, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# m46_apoe_edu   <- glmer(m46_apoe_edu,   data=hrs_full, family=poisson(link="log"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# 
# 
# #export raw results
# export(c("m40_main_eff","m41_incar_sex","m42_apoe_sex",
#          "m43_incar_race","m44_apoe_race","m45_incar_edu",
#          "m46_apoe_edu"), "../output/results/res_s2_raw.rdata")

# #table models
# res_40 <- tidy(m40_main_eff,   exponentiate=TRUE) %>% mutate(model="base")       %>% mutate(glance(m40_main_eff)) 
# res_41 <- tidy(m41_incar_sex,  exponentiate=TRUE) %>% mutate(model="incar_sex")  %>% mutate(glance(m41_incar_sex)) 
# res_42 <- tidy(m42_apoe_sex,   exponentiate=TRUE) %>% mutate(model="apoe_sex")   %>% mutate(glance(m42_apoe_sex)) 
# res_43 <- tidy(m43_incar_race, exponentiate=TRUE) %>% mutate(model="incar_race") %>% mutate(glance(m43_incar_race)) 
# res_44 <- tidy(m44_apoe_race,  exponentiate=TRUE) %>% mutate(model="apoe_race")  %>% mutate(glance(m44_apoe_race)) 
# res_45 <- tidy(m45_incar_edu,  exponentiate=TRUE) %>% mutate(model="incar_edu")  %>% mutate(glance(m45_incar_edu)) 
# res_46 <- tidy(m46_apoe_edu,   exponentiate=TRUE) %>% mutate(model="apoe_edu")   %>% mutate(glance(m46_apoe_edu)) 
# 
# #bind results
# subgroup_res <- bind_rows(res_40,
#                           res_41,
#                           res_42,
#                           res_43,
#                           res_44,
#                           res_45,
#                           res_46
# )
# 
# #export results
# export(subgroup_res, "../output/results/tab_s2_subgroups.csv")

#import
res <- import_list("../output/results/res_s2_raw.rdata")
list2env(res, .GlobalEnv)

m40_emmeans_i    <- emmeans(m40_main_eff, specs = pairwise ~ incar_ever,                 #type = "response", 
                            adjust="fdr")
m40_emmeans_a    <- emmeans(m40_main_eff, specs = pairwise ~ apoe_info99_4ct,            type = "response", adjust="fdr")
m40_emmeans_isex <- emmeans(m40_main_eff, specs = pairwise ~ incar_ever      |sex,       type = "response", adjust="fdr")
m40_emmeans_asex <- emmeans(m40_main_eff, specs = pairwise ~ apoe_info99_4ct |sex,       type = "response", adjust="fdr")
m40_emmeans_irac <- emmeans(m40_main_eff, specs = pairwise ~ incar_ever      |race_ethn, type = "response", adjust="fdr")
m40_emmeans_arac <- emmeans(m40_main_eff, specs = pairwise ~ apoe_info99_4ct |race_ethn, type = "response", adjust="fdr")
m40_emmeans_iedu <- emmeans(m40_main_eff, specs = pairwise ~ incar_ever      |edu,       type = "response", adjust="fdr")
m40_emmeans_aedu <- emmeans(m40_main_eff, specs = pairwise ~ apoe_info99_4ct |edu,       type = "response", adjust="fdr")

# m41_emmeans <- emmeans(m41_incar_sex,  specs = pairwise ~ incar_ever     |sex,       type = "response", adjust="fdr")
# m42_emmeans <- emmeans(m42_apoe_sex,   specs = pairwise ~ apoe_info99_4ct|sex,       type = "response", adjust="fdr")
# m43_emmeans <- emmeans(m43_incar_race, specs = pairwise ~ incar_ever     |race_ethn, type = "response", adjust="fdr")
# m44_emmeans <- emmeans(m44_apoe_race,  specs = pairwise ~ apoe_info99_4ct|race_ethn, type = "response", adjust="fdr")
# m45_emmeans <- emmeans(m45_incar_edu,  specs = pairwise ~ incar_ever     |edu,       type = "response", adjust="fdr")
# m46_emmeans <- emmeans(m46_apoe_edu,   specs = pairwise ~ apoe_info99_4ct|edu,       type = "response", adjust="fdr")

# m41_emmeans
# m42_emmeans
# m43_emmeans
# m44_emmeans
# m45_emmeans
# m46_emmeans

# rio::export(c("m41_emmeans",
#               "m42_emmeans",
#               "m43_emmeans",
#               "m44_emmeans",
#               "m45_emmeans",
#               "m46_emmeans"),
#             "../output/results/tab_s1_emmeans.rdata")

rio::export(c("m40_emmeans_i",
              "m40_emmeans_a",
              "m40_emmeans_isex",
              "m40_emmeans_asex",
              "m40_emmeans_irac",
              "m40_emmeans_arac",
              "m40_emmeans_iedu",
              "m40_emmeans_aedu"),
            "../output/results/tab_s1_emmeans.rdata")

#=END==========================================================================











