source("0_packages.R")

hrs <- import("hrs_full_analytic.rds")

theme_set(theme_classic())

#=venn: race, sex, incar, apoe=================================================

p_load(ggalluvial)

hrs %>% 
  distinct(hhidpn, .keep_all=TRUE) %>% 
  # arrange(incar_ever) %>% 
  ggplot(aes(axis4=race_ethn, axis3=fct_rev(sex), axis1=fct_rev(incar_ever), axis2=fct_rev(apoe_info99_4ct))) +
  geom_alluvium(aes(fill=incar_ever, alpha=incar_ever)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme(legend.position = "bottom") +
  scale_fill_manual(name = "",
                    values = c("grey90", "red")) +
  scale_alpha_manual(values = c(1, .8)) +
  scale_x_discrete(limits = c("Incarcerated", "APOE-4 Count", "Sex", "Race/Ethnicity")) +
  guides(alpha="none")
  # coord_polar()

#=age x incarceration==========================================================

age_mean <- mean(hrs$age)

incar_age <- ggplot(hrs, aes(incar_ever, y=age, col=incar_ever)) +
  geom_point(position = position_jitternormal(sd_x = .05, seed = 973101), alpha=.25, size=1.5) +
  geom_hline(yintercept = age_mean, linetype="dashed") +
  geom_violin(alpha=.15, fill="grey", color="black", linewidth=.75) +
  geom_boxplot(width=.25, color="black", size=.5, outlier.colour = NA) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.ticks.x = element_blank()) +
  labs(x="",
       y="Years of Age") +
  guides(color="none") +
  scale_color_manual(values = c("grey", "red"))

ggsave("../output/figures/1_incar_age.png", height = 5, width = 7)


#=cognitive function by year, by cohort========================================
cog_apoe4 <- hrs %>% 
  mutate(study = fct_recode(study,
                            "AHEAD\n<1924" = "AHEAD",
                           "CODA\n1924-30" = "CODA",
                            "HRS\n1931-41" = "HRS",
                           "WB\n1942-47"   = "WB",
                           "EBB\n1948-53"  = "EBB",
                           "MBB\n1954-59"  = "MBB")) %>% 
  group_by(year, study, apoe_info99_4ct
  ) %>% 
  summarise(tot=n(),
            count=sum(cog_2cat_num)) %>%
  mutate(prop = count/tot) %>% 
  mutate(apoe_info99_4ct = fct_recode(apoe_info99_4ct,
                                      "Zero copies" = "zero copies",
                                      "One copy" = "one copy",
                                      "Two copies" = "two copies")) %>% 
  ggplot(aes(year, prop, fill=fct_rev(apoe_info99_4ct))) +
  stat_smooth(
    geom = 'area', method = 'loess', span = 1/3,
    alpha = .75) +
  theme(legend.position = c(.15,.2),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text.y.right = element_text(size = 14, angle = 0),
        legend.background = element_rect(fill = "grey95", color = "grey20"),
        panel.spacing.y = unit(2, "lines")) +
  labs(title = 'Percent of APEO-4 allele carriers with "impaired" congnitive status by study year and HRS cohort',
       x="Study year",
       y='Percent cognitively impaired') +
  scale_x_continuous(breaks = seq(1998, 2018, 4)) +
  scale_y_continuous(breaks = c(0, .5, 1),
                     labels = scales::percent_format(),
                     limits = c(0,1)) +
  scale_fill_viridis_d(name="APOE-4 allele count",
                    option = "A") +
  coord_cartesian(clip = "off", expand = 0) +
  facet_grid(study~.) 
  
cog_apoe4

ggsave("../output/figures/2_cogstat_apoe4.png", height = 7, width = 10)

#=Subgroup analysis===========================================================

subgroup_res <- import_list("../output/results/tab_s1_emmeans.rdata")
list2env(subgroup_res, .GlobalEnv)

#group all EMMs
emmeans4i <- m40_emmeans_i$emmeans    %>% as_tibble() %>% mutate(across(c(rate, SE, asymp.LCL, asymp.UCL), ~exp(.)), var = incar_ever,      strata = "ref")     %>% select(-c(incar_ever))
emmeans4a <- m40_emmeans_a$emmeans    %>% as_tibble() %>% mutate(across(c(rate, SE, asymp.LCL, asymp.UCL), ~exp(.)), var = apoe_info99_4ct, strata = "ref")     %>% select(-c(apoe_info99_4ct))
emmeans41 <- m40_emmeans_isex$emmeans %>% as_tibble() %>% mutate(across(c(rate, SE, asymp.LCL, asymp.UCL), ~exp(.)), var = incar_ever,      strata = sex)       %>% select(-c(sex, incar_ever))
emmeans42 <- m40_emmeans_asex$emmeans %>% as_tibble() %>% mutate(across(c(rate, SE, asymp.LCL, asymp.UCL), ~exp(.)), var = apoe_info99_4ct, strata = sex)       %>% select(-c(sex,apoe_info99_4ct))
emmeans43 <- m40_emmeans_irac$emmeans %>% as_tibble() %>% mutate(across(c(rate, SE, asymp.LCL, asymp.UCL), ~exp(.)), var = incar_ever,      strata = race_ethn) %>% select(-c(race_ethn, incar_ever))
emmeans44 <- m40_emmeans_arac$emmeans %>% as_tibble() %>% mutate(across(c(rate, SE, asymp.LCL, asymp.UCL), ~exp(.)), var = apoe_info99_4ct, strata = race_ethn) %>% select(-c(race_ethn,apoe_info99_4ct))
emmeans45 <- m40_emmeans_iedu$emmeans %>% as_tibble() %>% mutate(across(c(rate, SE, asymp.LCL, asymp.UCL), ~exp(.)), var = incar_ever,      strata = edu)       %>% select(-c(edu, incar_ever))
emmeans46 <- m40_emmeans_aedu$emmeans %>% as_tibble() %>% mutate(across(c(rate, SE, asymp.LCL, asymp.UCL), ~exp(.)), var = apoe_info99_4ct, strata = edu)       %>% select(-c(edu,apoe_info99_4ct))

emmeans_all4 <- bind_rows(emmeans4i,
                          emmeans4a,
                          emmeans41,
                          emmeans42,
                          emmeans43,
                          emmeans44,
                          emmeans45,
                          emmeans46) 
rm(emmeans4i,
   emmeans4a,
   emmeans41,
   emmeans42,
   emmeans43,
   emmeans44,
   emmeans45,
   emmeans46)


# #group all contrasts
# contast41 <- m41_emmeans$contrasts %>% as_tibble() %>% mutate(strata = sex)       %>% select(-sex)
# contast42 <- m42_emmeans$contrasts %>% as_tibble() %>% mutate(strata = sex)       %>% select(-sex)
# contast43 <- m43_emmeans$contrasts %>% as_tibble() %>% mutate(strata = race_ethn) %>% select(-race_ethn)
# contast44 <- m44_emmeans$contrasts %>% as_tibble() %>% mutate(strata = race_ethn) %>% select(-race_ethn)
# contast45 <- m45_emmeans$contrasts %>% as_tibble() %>% mutate(strata = edu)       %>% select(-edu)
# contast46 <- m46_emmeans$contrasts %>% as_tibble() %>% mutate(strata = edu)       %>% select(-edu)
# 
# constrats_all4 <- bind_rows(contast41,
#                             contast42,
#                             contast43,
#                             contast44,
#                             contast45,
#                             contast46) %>% 
#   mutate(strata_group = case_when((strata %in% c("Male", "Female")) ~ "Sex",
#                            (strata %in% c("Black", "White", "Hispanic", "Other")) ~ "Race/Ethnicity",
#                            (strata %in% c("hs or more", "less than hs")) ~ "Education"))
# rm(contast41,
#    contast42,
#    contast43,
#    contast44,
#    contast45,
#    contast46)

#visualize
theme_set(theme_classic2())

refs_i <- emmeans_all4 %>% filter(strata_group=="Ref." & var_group=='Lifetime Incarceration') %>% pull(rate)
refs_a <- emmeans_all4 %>% filter(strata_group=="Ref." & !var_group=='Lifetime Incarceration') %>% pull(rate)

subgroup_plot <- emmeans_all4 %>% 
  mutate(strata_group = case_when((strata %in% c("Male", "Female")) ~ "Sex",
                                  (strata %in% c("Black", "White", "Hispanic", "Other")) ~ "Race",
                                  (strata %in% c("hs or more", "less than hs")) ~ "Education",
                                  TRUE ~ "Ref."),
         strata_group = fct_relevel(strata_group, "Ref.", "Sex", "Race", "Education"),
         var_group = case_when((var %in% c("Incarcerated", "Not Incarcerated")) ~ "Lifetime Incarceration",
                               (var %in% c("zero copies", "one copy", "two copies")) ~ "APOE-4 allele count"),
         var_group = fct_relevel(var_group, "Lifetime Incarceration")) %>% 
  
  filter(!strata_group=="Ref.") %>% 
ggplot(aes(rate, strata, col=var)) +
  geom_vline(xintercept = 1, linetype="dashed", color="grey") +
  geom_vline(data=. %>% filter(var_group=='Lifetime Incarceration'),  aes(xintercept=refs_i[1], col=var), linetype="dotted", show.legend = F) +
  geom_vline(data=. %>% filter(var_group=='Lifetime Incarceration'),  aes(xintercept=refs_i[2], col=var), linetype="dotted", show.legend = F) +
  geom_vline(data=. %>% filter(!var_group=='Lifetime Incarceration'), aes(xintercept=refs_a[1], col=var), linetype="dotted", show.legend = F) +
  geom_vline(data=. %>% filter(!var_group=='Lifetime Incarceration'), aes(xintercept=refs_a[2], col=var), linetype="dotted", show.legend = F) +
  geom_vline(data=. %>% filter(!var_group=='Lifetime Incarceration'), aes(xintercept=refs_a[3], col=var), linetype="dotted", show.legend = F) +
  geom_point(position = position_dodge(width = .5), size=2) +
  geom_linerange(aes(xmin=asymp.LCL, xmax=asymp.UCL), position = position_dodge(width = .5), linewidth=1) +
  theme(strip.placement = "outside",
        strip.background = element_rect(color = NA, fill = NA),
        panel.background = element_rect(color = "black"),
        axis.ticks.y = element_blank()) +
  labs(title = "Estimated Marginal Means of IRR within Strata",
       y="",
       x="IRR") +
  scale_y_discrete(limits=rev) +
  facet_grid(strata_group~var_group, scales = "free", switch = "y", space = "free_y") +
  coord_cartesian(xlim = c(1,2))


ggsave("../output/figures/fig_subgroup.png", height = 7, width = 10)

#=theoretical model===========================================================



# p_load_gh("ricardo-bion/ggradar")
p_load(ggradar)

data <- tibble(group = c("Genetic risk model",
                          "Environmental risk model",
                          "G+E model",
                          "G\u00D7E model"),
                Excess = c( 0, 0,0,.5),
                Environ. = c(  .25,1,1,  1),
                Genetic = c(1,  .25,1,  1)) %>% 
  mutate(group = fct_relevel(group, "Genetic risk model",
                             "Environmental risk model",
                             "G+E model",
                             "G\u00D7E model"))

group_cols <- c("darkred", "darkblue", "darkslateblue", "darkorchid4")

models <- ggradar(data,
        background.circle.colour = "white",
        gridline.min.linetype = 1,
        gridline.mid.linetype = 1,
        gridline.max.linetype = 1,
        gridline.min.colour = "grey",
        gridline.mid.colour = "grey",
        gridline.max.colour = "grey",
        plot.legend = FALSE,
        fill = TRUE,
        fill.alpha = .5,
        values.radar = c("", "", ""),
        group.point.size = 3,
        group.colours = group_cols) +
  facet_wrap(~group, ncol=2) +
  theme(strip.background = element_rect(fill = "black", color = "black"),
        strip.text = element_text(colour = "white", face = "bold")) +
  coord_fixed(clip = "off")

models
ggsave("../output/figures/theoretical_mods.tiff",
       models, height = 7, width = 9, dpi = 300)
