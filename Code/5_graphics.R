source("0_packages.R")

p_load(ggrepel,
       ggridges,
       ggtext)

hrs <- import("hrs_full_analytic.rds")
hrs_surv <- import("hrs_full_analytic_surv.rds")


theme_set(theme_classic())


#=Loneliness & social isolation overtime======================================

hrs %>% 
  group_by(year, study) %>% 
  summarise(lone = mean(lone_scale_pro, na.rm=TRUE),
            lone_sd = sd(lone_scale_pro, na.rm=TRUE),
            soc_iso = mean(soc_iso_index_pro, na.rm=TRUE),
            soc_iso_sd = sd(soc_iso_index_pro, na.rm=TRUE)) %>% 
  ggplot(aes(year, lone, col=study)) +
  geom_line()

hrs %>% 
  # group_by(year) %>%
  mutate(across(c(lone_scale_pro, soc_iso_index_pro), ~scale(., center=TRUE, scale=TRUE)[,1])) %>%
  group_by(year, study) %>% 
  summarise(lone = mean(lone_scale_pro, na.rm=TRUE),
            lone_sd = sd(lone_scale_pro, na.rm=TRUE),
            soc_iso = mean(soc_iso_index_pro, na.rm=TRUE),
            soc_iso_sd = sd(soc_iso_index_pro, na.rm=TRUE)) %>% 
  ggplot(aes(year, log(lone+1), col=study)) +
  geom_smooth(method = "lm")



hrs %>% 
  filter(!study=="MBB") %>% 
  # group_by(year) %>%
  # mutate(across(c(lone_scale_pro, soc_iso_index_pro), ~scale(., center=TRUE, scale=TRUE)[,1])) %>%
  group_by(year, study) %>% 
  summarise(lone = mean(lone_scale_pro, na.rm=TRUE),
            lone_sd = sd(lone_scale_pro, na.rm=TRUE),
            soc_iso = mean(soc_iso_index_pro, na.rm=TRUE),
            soc_iso_sd = sd(soc_iso_index_pro, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(study) %>% 
  mutate(soc_iso_index08 = soc_iso[year==2008],
         soc_iso_rel = soc_iso-soc_iso_index08) %>% 
  mutate(name_lab = ifelse(year==2020, as.character(study), NA)) %>% 
  #plot
  ggplot(aes(year, soc_iso_rel)) +
  geom_vline(xintercept = 2008, linetype="dashed", color="grey30") +
  geom_segment(aes(x=2004, xend=2020, y=0, yend=0), linetype="dashed", color="grey30") +
  geom_line(aes(col=study), linewidth=2, show.legend = F) +
  annotate("label", x=2008, y=.5, label="2008", size=5) +
  geom_label_repel(aes(label=name_lab),
                   size = 5,
                   direction = "y",
                   xlim = c(2020.8, NA),
                   hjust = 0,
                   segment.size = .7,
                   segment.alpha = .5,
                   segment.linetype = "dotted",
                   box.padding = .4,
                   segment.curvature = -0.1,
                   segment.ncp = 3,
                   segment.angle = 20) +
  theme(plot.title = element_text(face = "bold")) +
  labs(title = "Social isolation scores (0-6) in the HRS participants from 2004-2020, indexed to 2008",
       y="Change in social isolation score",
       x="Year") +
  scale_x_continuous(limits = c(2004, 2023),
                     breaks = seq(2005, 2020, by = 5)) +
  scale_y_continuous(limits = c(-1,1),
                     breaks = seq(-1, 1, .5),
                     labels=c("-1", "-.5", "0", "+.5", "+1")) +
  scale_color_viridis_d() +
  coord_cartesian(expand = FALSE)




hrs_surv %>% 
  ggplot(aes(y=study, x=lonely_pro_m, fill=factor(after_stat(quantile))))   +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = c(.025, .95), quantile_lines = TRUE
  ) +
  scale_fill_viridis_d(name = "Quartiles")



lexis_plot <- hrs %>% 
  mutate(age = year-birthyr) %>% 
  mutate(year_2 = case_when((year %in% c(2004, 2005)) ~ "2004-05",
                            (year %in% c(2006, 200))  ~ "2006-07",
                            (year %in% c(2008, 200))  ~ "2008-09",
                            (year %in% c(2010, 2011)) ~ "2010-11",
                            (year %in% c(2012, 2013)) ~ "2012-13",
                            (year %in% c(2014, 2015)) ~ "2014-15",
                            (year %in% c(2016, 2017)) ~ "2016-17",
                            (year %in% c(2018, 2019)) ~ "2018-19",
                            (year == 2020) ~ "2020")) %>% 
  group_by(age, year_2) %>% 
  summarise(soc_iso_ageyear = mean(soc_iso_index_pro)) %>% 
  group_by(age) %>% 
  mutate(soc_agestd = scale(soc_iso_ageyear)[,1]) %>% 
  ggplot(aes(year_2, age)) +
  geom_tile(aes(fill=soc_agestd)) +
  theme(plot.title = element_markdown(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 35, hjust=1),
        axis.title.y = element_text(angle = 0, vjust = .5),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        legend.position = "bottom",
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 14)) +
  scale_y_continuous(limits = c(50, 100)) +
  labs(title = str_wrap("Lexis plot of social isolation scores (age-standardized) 
                        for detecting possible <br><span style = 'color:blue;'>period</span> and 
       <span style = 'color:red3;'>cohort</span> effects among HRS participants."),
       x="",
       y="Age") +
  coord_fixed(ratio = .175, expand = F) +
  # scale_fill_viridis_c(name="Change in age-standardized\nsocial isolation score (SD)",
  #                      labels=c("+2", "+1", "0", "-1", "-2"),
  #                      direction = 1,
  #                      option = "B") +
  scale_fill_gradient2(name="Change in age-standardized\nsocial isolation score (SD)",
                       high = "darkred", mid = "white", low = "navy",
                       labels=c("-2", "-1", "0", "+1", "+2"))+
  geom_vline(xintercept = "2008-09", color="blue", size=1, linetype="dashed") +
  geom_segment(aes(x="2010-11", xend="2020", y=50, yend=58), color="red3", size=1, linetype="dashed")
  
lexis_plot
ggsave("../output/lexis_plot.png", lexis_plot, height = 10, width = 10)


