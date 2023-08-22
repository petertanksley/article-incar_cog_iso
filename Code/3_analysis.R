source("0_packages.R")
set.seed(793101)

hrs <- import("hrs_full_analytic_surv.rds")

#=STEPS========================================================================
# 1. Estimate total effect model (coxph without mediators)
# 2. Estimate mediator model (exposure predicting individual mediators)
# 3. Test independence (exposure and mediators predicting additional mediator)
# 4. Estimate weighted mediation model
# 5. Compute bootstrap confidence intervals 
# Sensitivity tests
# 1. Combine mediators
# 2. Non-linearities, interactions, and misclassification
#=============================================================================#

#=Step 1 - total effect model==================================================
cox1 <- coxph(Surv(event=cog_ever, time = cog_surv_age) ~ factor(incar_ever) +
                scale(age_base) + factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) +
                factor(smoke_ever) + factor(stroke_ever) +
                strata(study),
              data=hrs, id=hhidpn)
cox_toteff <- tidy(cox1, exponentiate = TRUE)
export(cox_toteff, "../output/results/1_coxtoteff.csv")


#=Step 2 - mediator models=====================================================

#lone <- incar
med_lone <- lm(scale(lonely_pro_m) ~ factor(incar_ever) +
                 scale(age_base) + factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) +
                 factor(smoke_ever) + factor(stroke_ever) + factor(study),
               data=hrs)
med_lone_res <- tidy(med_lone)
export(med_lone_res, "../output/results/2.1_med_lone.csv")

#soc_iso <- incar
med_soc_iso <- lm(scale(soc_iso_pro_m) ~ factor(incar_ever) +
                    scale(age_base) + factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) +
                    factor(smoke_ever) + factor(stroke_ever) + factor(study),
                  data=hrs)
med_soc_iso_res <- tidy(med_soc_iso)
export(med_soc_iso_res, "../output/results/2.2_med_soc_iso.csv")

# #=Step 3 - test independence of mediators======================================
# #(this is only needed for testing multiple paths simultaneously)
# 
# #lone <- soc_iso + incar
# med_loneXsoc_iso <- lm(scale(lonely_pro_m) ~ scale(soc_iso_pro_m) + factor(incar_ever) +
#                          scale(age_base) + factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) +
#                          factor(smoke_ever) + factor(stroke_ever) + factor(study),
#                        data=hrs)
# med_loneXsoc_iso_res <- tidy(med_loneXsoc_iso)
# 
# #highly sig. effect for soc_iso on loneliness
# #this indicates that the mediators are "intertwined". Estimation of indirect effects
# #should account for this. I will explore options for combining mediators.
# 
# 
#=Step 4 - full weighted mediation model=======================================

# NOTE: function neWeight creates expanded datasets with which to estimate NDE 
# and NIE. The function creates a new ID variable, and two new exposure 
# variables: 
# exp0 (indexing the NDE)
# exp1 (indexing the NIE) 


#LONELINESS
if (!file.exists("../output/results/3.1_full_med_lone.csv")){
  #write function to estimate mediation models for one mediator
  doEffectDecomp_lonely = function(data, index) {
    dat <- data[index,]
    
    #mediator model
    med <- paste0("scale(lonely_pro_m)")
    covars <- c("factor(incar_ever)", "scale(age_base)", "factor(sex)","factor(race_ethn)","factor(edu)",
                "scale(social_origins)","factor(smoke_ever)","factor(stroke_ever)","factor(study)")
    formula <- paste(med, "~", paste(covars, collapse = "+"))
    
    #expand data
    m_mod <- glm(as.formula(formula),
                 family = gaussian(),
                 data=dat)
    exp_data <- neWeight(m_mod)
    
    #extract weights
    m_weights <- weights(exp_data)
    
    #fit weighted model
    cox <- coxph(Surv(event = cog_ever, time = cog_surv_age) ~ factor(incar_ever1) + factor(incar_ever0) + 
                   scale(age_base) + factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) + 
                   factor(smoke_ever) + factor(stroke_ever) + strata(study),
                 data = exp_data, weights = m_weights)
    
    # Return values: Estimates for total, direct, indirect effect 
    TE = exp(sum(coef(cox)[c('factor(incar_ever1)Incarcerated', 'factor(incar_ever0)Incarcerated')]))
    DE = exp(unname(coef(cox)['factor(incar_ever0)Incarcerated']))
    IE = exp(sum(coef(cox)['factor(incar_ever1)Incarcerated']))
    PM = log(IE) / log(TE)
    
    med_res <- tibble(term = c("total", "direct", "indirect", "pro_med"),
                      estimate = c(TE, DE, IE, PM))
    
    return(c(exp(coef(cox)), TE=TE, DE=DE, IE=IE, PM=PM))
  }
  
  #res_lonely <- doEffectDecomp_lonely(data = hrs)
  
  #bootstrap estimation function (takes a long time)
  boot_lonely <- boot(data=hrs, doEffectDecomp_lonely, R=1000, parallel = "multicore", ncpus = 4)
  #tidy dataframe
  lonely_res  <- tidy(boot_lonely,  conf.int = TRUE, conf.method = "perc")
  #export
  export(lonely_res, "../output/results/3.1_full_med_lone.csv")
  
  #clear environment
  rm(boot_lonely)
  
} else {
  lonely_res <- import("../output/results/3.1_full_med_lone.csv")
}


#=SOCIAL ISOLATION
if(!file.exists("../output/results/3.2_full_med_soc_iso.csv")) {
  #write function for social isolation
  doEffectDecomp_soc_iso = function(data, index) {
    
    dat <- data[index,]
    #mediator model
    med <- paste0("scale(soc_iso_pro_m)")
    covars <- c("factor(incar_ever)", "scale(age_base)", "factor(sex)","factor(race_ethn)","factor(edu)",
                "scale(social_origins)","factor(smoke_ever)","factor(stroke_ever)","factor(study)")
    formula <- paste(med, "~", paste(covars, collapse = "+"))
    
    #expand data
    m_mod <- glm(as.formula(formula),
                 family = gaussian(),
                 data=dat)
    exp_data <- neWeight(m_mod)
    
    #extract weights
    m_weights <- weights(exp_data)
    
    #fit weighted model
    cox <- coxph(Surv(event = cog_ever, time = cog_surv_age) ~ factor(incar_ever1) + factor(incar_ever0) + 
                   scale(age_base) + factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) + 
                   factor(smoke_ever) + factor(stroke_ever) + strata(study),
                 data = exp_data, weights = m_weights)
    
    # Return values: Estimates for total, direct, indirect effect 
    TE = exp(sum(coef(cox)[c('factor(incar_ever1)Incarcerated', 'factor(incar_ever0)Incarcerated')]))
    DE = exp(unname(coef(cox)['factor(incar_ever0)Incarcerated']))
    IE = exp(sum(coef(cox)['factor(incar_ever1)Incarcerated']))
    PM = log(IE) / log(TE)
    
    return(c(exp(coef(cox)), TE=TE, DE=DE, IE=IE, PM=PM))
  }
  
  # res_soc_iso <- doEffectDecomp_soc_iso(data = hrs)
  
  #bootstrap estimation function (takes a long time)
  boot_soc_iso <- boot(data=hrs, doEffectDecomp_soc_iso, R=1000, parallel = "multicore", ncpus = 4)
  #tidy dataframe
  soc_iso_res <- tidy(boot_soc_iso, conf.int = TRUE, conf.method = "perc")
  #export
  export(soc_iso_res, "../output/results/3.2_full_med_soc_iso.csv")
  
  #clean environment
  rm(boot_soc_iso)
  
} else {
  soc_iso_res <- import("../output/results/3.2_full_med_soc_iso.csv")
}


#=Sensitivity analyses=========================================================

#=Test 1: non-linearity, interactions in mediator models
#lone <- incar
med_lone <- lm(scale(lonely_pro_m) ~ factor(incar_ever) * 
                 (scale(age_base) + factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) +
                    factor(smoke_ever) + factor(stroke_ever) + factor(study)) + I(age_base^2),
               data=hrs)
med_lone_res <- tidy(med_lone)
export(med_lone_res, "../output/results/s1.1_med_lone.csv")

#soc_iso <- incar
med_soc_iso <- lm(scale(soc_iso_pro_m) ~ factor(incar_ever) * 
                    (scale(age_base) + factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) +
                       factor(smoke_ever) + factor(stroke_ever) + factor(study)) + I(age_base^2),
                  data=hrs)
med_soc_iso_res <- tidy(med_soc_iso)
export(med_soc_iso_res, "../output/results/s1.2_med_soc_iso.csv")
#Incarceration*Black is nominally significant at P=0.0019

#SOCIAL ISOLATION + interaction
med <- paste0("scale(soc_iso_pro_m)")
covars <- c("factor(incar_ever)", "scale(age_base)", "factor(sex)","factor(race_ethn)","factor(edu)",
            "scale(social_origins)","factor(smoke_ever)","factor(stroke_ever)","factor(study)",
            "factor(incar_ever):factor(race_ethn)")
formula <- paste(med, "~", paste(covars, collapse = "+"))

#expand data
m_mod <- glm(as.formula(formula),
             family = gaussian(),
             data=hrs)
exp_data <- neWeight(m_mod)

#extract weights
m_weights <- weights(exp_data)

#fit weighted model
cox <- coxph(Surv(event = cog_ever, time = cog_surv_age) ~ factor(incar_ever1) + factor(incar_ever0) + 
               scale(age_base) + factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) + 
               factor(smoke_ever) + factor(stroke_ever) + strata(study),
             data = exp_data, weights = m_weights)
s1.2_med_soc_iso_int <- tidy(cox, exponentiate = TRUE)
export(s1.2_med_soc_iso_int, "../output/results/s1.2_med_soc_iso_int.csv")
#results don't really change with adjustment

#=Test 2: visually inspect Schoenfeld residuals

#PRIMARY OUTCOME MODEL
cox <- coxph(Surv(event = cog_ever, time = cog_surv_age) ~ factor(incar_ever) +  
               scale(age_base) + factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) + 
               factor(smoke_ever) + factor(stroke_ever) + strata(study),
             data = hrs)
cox.zph(cox) %>% ggcoxzph()

#LONELINESS
med <- paste0("scale(lonely_pro_m)")
covars <- c("factor(incar_ever)", "scale(age_base)", "factor(sex)","factor(race_ethn)","factor(edu)",
            "scale(social_origins)","factor(smoke_ever)","factor(stroke_ever)","factor(study)")
formula <- paste(med, "~", paste(covars, collapse = "+"))

#expand data
m_mod <- glm(as.formula(formula),
             family = gaussian(),
             data=hrs)
exp_data <- neWeight(m_mod)

#extract weights
m_weights <- weights(exp_data)

#fit weighted model
cox <- coxph(Surv(event = cog_ever, time = cog_surv_age) ~ factor(incar_ever1) + factor(incar_ever0) + 
               scale(age_base) + factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) + 
               factor(smoke_ever) + factor(stroke_ever) + strata(study),
             data = exp_data, weights = m_weights)
cox.zph(cox) %>% ggcoxzph()


#SOCIAL ISOLATION
med <- paste0("scale(soc_iso_pro_m)")
covars <- c("factor(incar_ever)", "scale(age_base)", "factor(sex)","factor(race_ethn)","factor(edu)",
            "scale(social_origins)","factor(smoke_ever)","factor(stroke_ever)","factor(study)")
formula <- paste(med, "~", paste(covars, collapse = "+"))

#expand data
m_mod <- glm(as.formula(formula),
             family = gaussian(),
             data=hrs)
exp_data <- neWeight(m_mod)

#extract weights
m_weights <- weights(exp_data)

#fit weighted model
cox <- coxph(Surv(event = cog_ever, time = cog_surv_age) ~ factor(incar_ever1) + factor(incar_ever0) + 
               scale(age_base) + factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) + 
               factor(smoke_ever) + factor(stroke_ever) + strata(study),
             data = exp_data, weights = m_weights)
cox.zph(cox) %>% ggcoxzph()


#=END==========================================================================






