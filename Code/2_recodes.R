source("0_packages.R")

hrs_merged <- import("hrs_merged.rds")

conflicted::conflict_prefer_all("dplyr", quiet = TRUE)
conflicted::conflict_prefer_all("tidylog", quiet = TRUE)

#=Recodes======================================================================
# age, two category version of cogfunction
#=============================================================================#


hrs_recodes <- hrs_merged %>% 
  #two-category cognitive function variable
  mutate(cog_2cat = case_when((cogfunction=="normal") ~ "Normal",
                              (cogfunction %in% c("cind", "demented")) ~ "Impaired",
                              TRUE ~ NA_character_),
         cog_2cat_num = case_when((cogfunction=="normal") ~ 0,
                                  (cogfunction %in% c("cind", "demented")) ~ 1,
                                  TRUE ~ NA_real_)) %>% 
  #binary edu
  mutate(edu = ifelse(edu_yrs>=12, "hs or more", "less than hs"),
         edu = fct_relevel(edu, "hs or more", "less than hs")) %>% 
  #numeric loneliness/social isolation
  mutate(lone_scale_pro_num = case_when((lone_scale_pro_cat=="lonely") ~ 1, 
                                        (lone_scale_pro_cat=="not lonely") ~ 0,
                                        TRUE ~ NA_real_)) %>% 
  mutate(soc_iso_index_pro_num = case_when((soc_iso_index_pro_cat=="isolated") ~ 1, 
                                        (lone_scale_pro_cat=="not isolated") ~ 0,
                                        TRUE ~ NA_real_)) %>% 
  #summarize loneliness/social isolation within individuals
  group_by(hhidpn) %>% 
  mutate(lonely_ct = sum(lone_scale_pro_num, na.rm = TRUE),
         lonely_m = ifelse(lonely_ct==0, 0, mean(lone_scale_pro_num, na.rm = TRUE)),
         lonely_pro_m = mean(lone_scale_pro, na.rm = TRUE),
         lonely_per = case_when((lonely_ct>=2) ~ "persistent", 
                                (lonely_ct<2) ~ "not",
                                TRUE ~ NA_character_)) %>% 
  mutate(soc_iso_ct = sum(soc_iso_index_pro_num, na.rm = TRUE),
         soc_iso_m = ifelse(soc_iso_ct==0, 0, mean(soc_iso_index_pro_num, na.rm = TRUE)),
         soc_iso_pro_m = mean(soc_iso_index_pro, na.rm = TRUE),
         soc_iso_per = case_when((soc_iso_ct>=2) ~ "persistent",
                                 (soc_iso_ct<2) ~ "not",
                                 TRUE ~ NA_character_)) %>% 
  ungroup() %>% 
  #fix non-valid values
  mutate(across(c(lonely_pro_m, soc_iso_pro_m), ~ifelse(.=='NaN', NA, .)))


#create analytic sample
hrs_full <- hrs_recodes %>% 
  select(hhidpn,
         study, race_ethn, sex, birthyr, year, firstiw, dod_yr,
         starts_with("cog_2cat"), cogfunction,
         incar_ever, incar_time_3cat,
         lone_scale_pro, lone_scale_pro_cat, lonely_ct, lonely_per, lonely_pro_m, #starts_with("lone_"),
         soc_iso_index_pro, soc_iso_index_pro_cat, soc_iso_ct, soc_iso_per, soc_iso_pro_m, #ends_with("_bin"),
         stroke_ever,
         apoe_info99_4ct,
         social_origins,
         edu_yrs, edu,
         smoke_ever) %>%
  drop_na(-c(dod_yr)) #removed 290,847 rows (90%), 33,766 rows remaining


# #check case count
# hrs_full %>% count(cases = n_distinct(hhidpn))
# # cases     n
# # 11355 33766


#export analytic dataframes
rio::export(hrs_full, "hrs_full_analytic.rds")


#pivot dataframe wider to "fold" observations into obs_1, obs_2,...,obs_n format

  
#make survival model recodes (age as time-scale approach)
hrs_surv <- hrs_full %>%
  mutate(study_age = year-birthyr,
         firstiw_age = firstiw-birthyr,
         dod_age = as_numeric(dod_yr)-birthyr) %>%
  group_by(hhidpn) %>% 
  mutate(age_base = min(study_age)) %>% #age at baseline
  fill(dod_yr, .direction="updown") %>% #fill DOD
  mutate(cog_first = ifelse(cog_2cat_num==1, study_age, NA),
         cog_first = pmin(cog_first)) %>%
  fill(cog_first, .direction="updown") %>%
  mutate(cog_ever = max(cog_2cat_num),
         cog_surv_age = ifelse(cog_ever>0, cog_first, max(study_age))) %>%
  ungroup()
  



#time-independent HRS data
hrs_ind <- hrs_surv %>%
  distinct(hhidpn, .keep_all=TRUE) %>% #removed 22,411 rows (66%), 11,355 rows remaining
  filter(cog_surv_age>firstiw_age) %>%  #removed 114 rows (1%), 11,241 rows remaining
  select(hhidpn, study, age_base, incar_ever, cog_ever, cog_surv_age, 
         smoke_ever, stroke_ever, sex, race_ethn, edu, social_origins,
         lonely_pro_m, soc_iso_pro_m)

#export analytic dataframes
rio::export(hrs_ind, "hrs_full_analytic_surv.rds")

#clean up environment
rm(list=ls())

