source("0_packages.R")

#=import cleaned HRS data======================================================
#path to cleaned data directory
hrs_clean_dir <- "../../../../Desktop/Data/HRS/cleaned_variables"

#DV: Cognitive status
cog_stat <- import(glue("{hrs_clean_dir}/dementia_status/dementia_status_clean.rds"))

#IVs: incarceration history
incar   <- import(glue("{hrs_clean_dir}/incarceration/incarceration_cleaned.rds")) 

#Meds: loneliness & social isolation
lone <- import(glue("{hrs_clean_dir}/loneliness/loneliness_clean.rds")) %>% select(hhidpn, year, lone_scale_pro)
soc_iso <- import(glue("{hrs_clean_dir}/social_isolation/social_iso_clean.rds")) %>% select(hhidpn, year, soc_iso_index_pro)

#CVs: demographics; APOE status; history of stroke; social origins; educational attainment; smoking status
demo     <- import(glue("{hrs_clean_dir}/demographics/demo_cleaned.rds"))           
apoe     <- import(glue("{hrs_clean_dir}/serotonin_apoe/serotonin_apoe_clean.rds")) 
stroke   <- import(glue("{hrs_clean_dir}/stroke/stroke_clean.rds"))   
smoke    <- import(glue("{hrs_clean_dir}/smoke/smoke_clean.rds")) %>% distinct(hhidpn, .keep_all=TRUE) %>% select(hhidpn, smoke_ever)
soc_orig <- import(glue("{hrs_clean_dir}/social_origins/social_origin_clean.rds"))  
# ses      <- import(glue("{hrs_clean_dir}ses/ses_clean.rds"))
edu      <- import(glue("{hrs_clean_dir}/education/edu_tracker_clean.rds"))
# income   <- import(glue("{hrs_clean_dir}income/income_clean.rds"))    
death    <- import(glue("{hrs_clean_dir}/death/death_clean.rds"))

#=merge cleaned HRS data=======================================================

# convert all year columns to factor
year_num <- function(x) {
  x %>%
    clean_names() %>% 
    mutate(year = as_numeric(year))
}

# convert all hhidpn columns to factor
id_fct <- function(x) {
  x %>% 
    clean_names() %>% 
    mutate(hhidpn = as_factor(hhidpn))
}

#longitudinal dataframes: convert year to numeric
dfs_years <- list(cog_stat, stroke, death, lone, soc_iso
                  ) %>% 
  lapply(., year_num)

#all dataframes: convert id for factor
dfs_ids <- c(dfs_years, list(incar, demo, apoe, soc_orig, edu, smoke
                             )) %>% 
  lapply(., id_fct)

#clear up memory
rm(list = setdiff(ls(), "dfs_ids"))

#merge all dataframes
hrs_merged <- reduce(dfs_ids, full_join)
rm(dfs_ids)

#=Select study variables; drop observations who died before interview

#select study variables
hrs_merged_studyvars <- hrs_merged %>% 
  select(hhidpn, year, firstiw,
         race_ethn, sex, study, birthyr,
         dod_yr, alive,
         cogfunction,
         incar_ever, incar_time_3cat,
         lone_scale_pro, 
         soc_iso_index_pro,
         stroke_ever,
         apoe_info99_4ct,
         social_origins,
         # ses,
         smoke_ever,
         edu_yrs,
         # income_hh,
         ) %>% 
  # filter(race_ethn %in% c("White", "Black")) %>% 
  filter(firstiw<=year) %>% #removed 255,512 rows (27%), 704,836 rows remaining
  filter(as_numeric(dod_yr)>=year | is.na(dod_yr)) %>% #removed 118,886 rows (17%), 585,950 rows remaining
  filter(as_numeric(dod_yr)>2012 | is.na(dod_yr))  #removed 88,183 rows (15%), 497,767 rows remaining
  # select(-firstiw)
#drop rows with missing on all columns (ignore time-stable variables)
var_list <- c("cogfunction",
              "incar_ever", "incar_time_3cat",
              "lone_scale_pro",
              "soc_iso_index_pro",
              "stroke_ever", "smoke_ever",
              "apoe_info99_4ct")

hrs_merged_studyvars_sparse <- hrs_merged_studyvars %>% 
  filter(!if_all(all_of(var_list), is.na)) #removed 720 rows (<1%), 497,047 rows remaining

#export final merged dataframe
export(hrs_merged_studyvars_sparse, "hrs_merged.rds")
rm(list=ls())

#=END==========================================================================

