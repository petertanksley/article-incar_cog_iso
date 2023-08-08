if(system.file(package="pacman")=="") utils::install.packages("pacman")
library(pacman)

# p_load_gh(c("zabore/condsurv","zabore/ezfun"))
p_load(rio,
       tidyverse,
       # tidylog,
       ###
       boot,
       broom,
       broom.mixed,
       corrr,
       emmeans,
       flextable,
       geepack,
       ggforce,
       ggridges,
       glue,
       gt,
       gtExtras,
       gtsummary,
       janitor,
       lme4,
       naniar,
       patchwork,
       purrr,
       renv,
       sjlabelled,
       webshot2,
       #survival packages
       # condsurv,
       # ezfun,
       medflex,
       survival,
       survminer,
       tidycmprsk)

conflicted::conflict_prefer_all("dplyr", quiet = TRUE)
# conflicted::conflict_prefer_all("tidylog", quiet = TRUE)

# # renv::init()
# renv::snapshot()




