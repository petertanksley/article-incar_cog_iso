if(system.file(package="pacman")=="") utils::install.packages("pacman")
library(pacman)

# p_load_gh(c("zabore/condsurv","zabore/ezfun"))
p_load(rio,
       tidyverse,
       tidylog,
       ###
       broom,
       broom.mixed,
       corrr,
       emmeans,
       flextable,
       ggforce,
       ggridges,
       glue,
       gt,
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
       survival,
       survminer,
       tidycmprsk)

# # renv::init()
# renv::snapshot()




