source("0_packages.R")

hrs_merged <- import("hrs_merged.rds")



#=Recodes======================================================================
# age, two category version of cogfunction, drop HRS/LBB cohorts (too few) 
#=============================================================================#


hrs_recodes <- hrs_merged %>% 
  mutate(age = year-birthyr) %>% 
  mutate(cog_2cat = case_when((cogfunction=="normal") ~ "Normal",
                              (cogfunction %in% c("cind", "demented")) ~ "Impaired",
                              TRUE ~ NA_character_),
         cog_2cat_num = case_when((cogfunction=="normal") ~ 0,
                                  (cogfunction %in% c("cind", "demented")) ~ 1,
                                  TRUE ~ NA_real_)) %>% 
  mutate(edu = ifelse(edu_yrs>=12, "hs or more", "less than hs"),
         edu = fct_relevel(edu, "hs or more", "less than hs")) %>% 
  filter(!study %in% c("HRS", "LBB")) #removed 137,126 rows (53%), 122,510 rows remaining

#create analytic sample
hrs_study_vars <- hrs_recodes %>% 
  select(hhidpn,
         study, race_ethn, sex, birthyr, year, age, firstiw, 
         # dod_yr,
         cogfunction,
         # ad_pgs, starts_with("pc"),
         incar_ever, 
         # incar_time_3cat,
         stroke_ever,
         apoe_info99_4ct,
         social_origins,
         # edu_yrs, 
         edu,
         smoke_ever)

hrs_study_vars %>% count(cases = n_distinct(hhidpn))


tiff("../output/figures/missing_upset.tiff", units="in", width=10, height=8, res=300)
naniar::gg_miss_upset(hrs_study_vars, 
                      nsets = n_var_miss(hrs_study_vars)) 
dev.off()

