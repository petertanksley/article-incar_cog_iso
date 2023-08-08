
hrs <- import("hrs_full_analytic.rds")
hrs_surv <- import("hrs_full_analytic_surv.rds")


test <- glmer(factor(lone_class) ~ factor(incar_ever) +
                factor(sex) + factor(race_ethn) + scale(age) + factor(edu) + scale(social_origins) +
                factor(smoke_ever) + factor(study) + (1|hhidpn),
              family = binomial(),
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
              data=hrs)

# #total effect
# tot_eff <- coxph(Surv(cog_surv_age, cog_ever) ~ incar_ever_lgl +
#                    factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) +
#                    factor(smoke_ever) + strata(study),
#                  data = hrs,
#                  id=hhidpn)
# 
# d=hrs

#effect decomposition function from: https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-14-9
#=single mediator - social isolation


doEffectDecomp = function(d) {
  #step #1: replicate exposure variable, predict mediators
  d=hrs_surv
  
  d$incar_ever_lgl_TEMP <- d$incar_ever_lgl
  m_soc_iso <- glm(soc_iso_class_lgl ~ incar_ever_lgl_TEMP +
                     factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) +
                     factor(smoke_ever) + factor(study),
                   family = binomial(),
                   data = d)
  
  # Step 2: Replicate data with different exposures for loneliness
  d1 = d2 = d
  d1$med_soc_iso = d1$incar_ever_lgl
  d2$med_soc_iso = !d2$incar_ever_lgl
  newd = rbind(d1,d2)
  
  # Step 3: Compute weights for social isolation
  newd$incar_ever_lgl_TEMP = newd$incar_ever_lgl
  w = predict(m_soc_iso, newdata = newd, type = "response")
  direct = ifelse(newd$soc_iso_class_lgl, w, 1-w)
  newd$incar_ever_lgl_TEMP = newd$med_soc_iso
  w = predict(m_soc_iso, newdata = newd, type = "response")
  indirect = ifelse(newd$soc_iso_class_lgl, w, 1-w)
  newd$w_soc_iso_class_lgl = indirect/direct
  
  # Step 4: Weighted Cox Model
  newd$weight = newd$w_soc_iso_class_lgl
  cox = coxph(Surv(cog_surv_age, cog_ever) ~ incar_ever_lgl + med_soc_iso +
                factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) +
                factor(smoke_ever) + strata(study),
              weights = weight,
              data = newd)
  
  # Return value: Estimates for total, direct, indirect effects
  TE = exp(sum(coef(cox)[c('incar_ever_lglTRUE','med_soc_isoTRUE')])) 
  DE = exp(unname(coef(cox)['incar_ever_lglTRUE']))
  IE = exp(sum(coef(cox)[c('med_soc_isoTRUE')]))
  PM = log(IE) / log(TE)
  
  return(c(exp(coef(cox)), TE=TE, DE=DE, IE=IE, PM=PM))
}

doEffectDecomp(d=d)
HRs = replicate(10, doEffectDecomp(d)) 
apply(HRs, 1, quantile, c(0.025, 0.975))

#=Both mediators===============================================================

doEffectDecomp = function(d) {
  #step #1: replicate exposure variable, predict mediators
  d=hrs
  
  d$incar_ever_lgl_TEMP <- d$incar_ever_lgl
  m_lone <- glm(lone_class_lgl ~ incar_ever_lgl_TEMP +
                  factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) +
                  factor(smoke_ever) + factor(study),
                family = binomial(),
                data = d)
  m_soc_iso <- glm(soc_iso_class_lgl ~ incar_ever_lgl_TEMP +
                     factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) +
                     factor(smoke_ever) + factor(study),
                   family = binomial(),
                   data = d)
  
  # Step 2a: Replicate data with different exposures for loneliness
  d1 = d2 = d
  d1$med_lone = d1$incar_ever_lgl
  d2$med_lone = !d2$incar_ever_lgl
  newd = rbind(d1,d2)
  
  # Step 2b: Replicate data with different exposures for social isolation
  d1 = d2 = newd
  d1$med_soc_iso = d1$incar_ever_lgl
  d2$med_soc_iso = !d2$incar_ever_lgl
  newd = rbind(d1,d2)
  
  # Step 3a: Compute weights for loneliness
  newd$incar_ever_lgl_TEMP = newd$incar_ever_lgl
  w = predict(m_lone, newdata = newd, type = "response")
  direct = ifelse(newd$lone_class_lgl, w, 1-w)
  newd$incar_ever_lgl_TEMP = newd$med_lone
  w = predict(m_lone, newdata = newd, type = "response")
  indirect = ifelse(newd$lone_class_lgl, w, 1-w)
  newd$w_lone_class_lgl = indirect/direct
  
  # Step 3b: Compute weights for social isolation
  newd$incar_ever_lgl_TEMP = newd$incar_ever_lgl
  w = predict(m_soc_iso, newdata = newd, type = "response")
  direct = ifelse(newd$soc_iso_class_lgl, w, 1-w)
  newd$incar_ever_lgl_TEMP = newd$med_soc_iso
  w = predict(m_soc_iso, newdata = newd, type = "response")
  indirect = ifelse(newd$soc_iso_class_lgl, w, 1-w)
  newd$w_soc_iso_class_lgl = indirect/direct
  
  # Step 4: Weighted Cox Model
  newd$weight = newd$w_lone_class_lgl * newd$w_soc_iso_class_lgl
  cox = coxph(Surv(cog_surv_age, cog_ever) ~ incar_ever_lgl + med_lone + med_soc_iso +
                factor(sex) + factor(race_ethn) + factor(edu) + scale(social_origins) +
                factor(smoke_ever) + strata(study),
              weights = weight,
              data = newd)
  
  # Return value: Estimates for total, direct, indirect effects
  TE = exp(sum(coef(cox)[c('incar_ever_lglTRUE', 'med_loneTRUE', 'med_soc_isoTRUE')])) 
  DE = exp(unname(coef(cox)['incar_ever_lglTRUE']))
  IE = exp(sum(coef(cox)[c('med_loneTRUE', 'med_soc_isoTRUE')]))
  PM = log(IE) / log(TE)
  
  return(c(exp(coef(cox)), TE=TE, DE=DE, IE=IE, PM=PM))
}

doEffectDecomp(d=d)
HRs = replicate(100, doEffectDecomp(d)) 
apply(HRs, 1, quantile, c(0.025, 0.975))

#=confidence intervals







# #check correlations
# 
# hrs %>% 
#   group_by(year) %>% 
#   summarize(cor=cor(as.numeric(lone_scale_pro_cat), as.numeric(soc_iso_index_pro_cat)))
# # year      cor
# # <dbl>    <dbl>
# #  2004    0.163
# #  2006    0.182
# #  2008    0.160
# #  2010    0.130
# #  2012    0.158
# #  2014    0.147
# #  2016    0.151
# #  2018    0.159
# #  2020    0.135

