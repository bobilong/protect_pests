library(mgcv)
library(ggplot2)
library(terra)
library(data.table)
library(dplyr)
library(MatchIt)

country_shp <- rnaturalearthdata::countries50 %>% vect() %>% .[,c('iso_a3','formal_en','continent')]
  as.data.frame()
countryRaster <- rasterize(country_shp, globalRaster, field = 'iso_a3')

#=================match result=================
match_table <- fread("./data/match_table.csv")

match_table <- match_table %>%
  mutate(
    allPest = Animals_records + Dieases_records + Plants_records,
    any_pest = as.integer(allPest > 0),
    
    gbif_effort_nonpest = pmax(gbif_effort - allPest, 0),
    log_gbif_effort = log1p(gbif_effort_nonpest),
    
    carbon_decile = cut(
      carbon_percentile,
      breaks = seq(0, 1, by = 0.1),
      include.lowest = TRUE,
      labels = c(
        "0–10", "10–20", "20–30", "30–40", "40–50",
        "50–60", "60–70", "70–80", "80–90", "90–100"
      )
    )
  )

match_df <- match_table %>%
  mutate(
    treat = case_when(
      pa_fraction >= 0.5 ~ 1L,
      pa_fraction == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),
    ecoregion = factor(ecoregion),
    iso_a3 = factor(iso_a3),
    
    log_human_foot = log1p(human_foot),
    log_ele = log1p(ele + abs(min(ele, na.rm = TRUE)) + 1)
  ) %>%
  filter(
    !is.na(treat),
    treeCover > 10,
    !is.na(carbon_percentile),
    !is.na(gbif_effort_nonpest),
    !is.na(iso_a3)
  )

df_match_country <- match_df %>%
  group_by(ecoregion, iso_a3) %>%
  filter(
    any(treat == 1),
    any(treat == 0)
  ) %>%
  ungroup()

m_country <- matchit(
  treat ~
    carbon_percentile +
    treeCover +
    human_foot +
    bio01 +
    bio02 +
    bio04 +
    bio15 +
    ele +
    FFI,
  exact = ~ ecoregion + iso_a3,
  data = df_match_country,
  method = "nearest",
  distance = "glm",
  caliper = 0.2,
  ratio = 1,
  replace = FALSE,
  discard = "both"
)
# saveRDS(m_country, file = "./data/matchit_model.rds")
matchit_data <- match.data(m_country)
fwrite(matchit_data, "./data/matchit_data.csv")

m_country_effort <- matchit(
  treat ~
    carbon_percentile +
    treeCover +
    human_foot +
    bio01 +
    bio02 +
    bio04 +
    bio15 +
    ele +
    FFI +
    log_gbif_effort,
  exact = ~ ecoregion + iso_a3,
  data = df_match_country,
  method = "nearest",
  distance = "glm",
  caliper = 0.2,
  ratio = 1,
  replace = FALSE,
  discard = "both"
)

matchit_data_effort <- match.data(m_country_effort)
# fwrite(matchit_data_effort, "./data/matchit_data_effort.csv")

#==============modeling=======================
library(mgcv)

matchit_data2 <- fread("./data/matchit_data.csv") %>%
  as.data.frame() %>%
  mutate(
    any_pest = as.integer(allPest > 0),
    
    treat = factor(
      treat,
      levels = c(0, 1),
      labels = c("Outside PA", "Inside PA")
    ),
    
    carbon_decile = cut(
      carbon_percentile,
      breaks = seq(0, 1, by = 0.1),
      include.lowest = TRUE,
      labels = c(
        "0–10", "10–20", "20–30", "30–40", "40–50",
        "50–60", "60–70", "70–80", "80–90", "90–100"
      )
    ),
    
    gbif_effort_nonpest = pmax(gbif_effort - allPest, 0),
    log_gbif_effort = log1p(gbif_effort_nonpest),
    
    ecoregion = factor(ecoregion),
    iso_a3 = factor(iso_a3)
  ) %>%
  filter(
    !is.na(any_pest),
    !is.na(treat),
    !is.na(carbon_percentile),
    !is.na(carbon_decile),
    !is.na(weights),
    !is.na(log_gbif_effort)
  )


decile_labs <- c(
  "0–10", "10–20", "20–30", "30–40", "40–50",
  "50–60", "60–70", "70–80", "80–90", "90–100"
)


gam_fig2 <- bam(
  any_pest ~
    treat +
    s(carbon_percentile, by = treat, bs = "cr", k = 6) +
    s(log_gbif_effort, bs = "cr", k = 10) +
    s(bio01, bs = "cr", k = 10) +
    s(bio02, bs = "cr", k = 10) +
    s(bio04, bs = "cr", k = 10) +
    s(bio15, bs = "cr", k = 10) +
    s(human_foot, bs = "cr", k = 20) +
    s(ele, bs = "cr", k = 20) +
    s(treecover2000, bs = "cr", k = 10) +
    s(FFI, bs = "cr", k = 10) +
    s(x, y, k = 80),
  data = matchit_data2,
  weights = weights,
  family = binomial(link = "logit"),
  method = "fREML",
  discrete = TRUE,
  select = TRUE
)

saveRDS(gam_fig2, file = "./data/gam_fig2_model.rds")

summary(gam_fig2)

set.seed(123)

gam_fig2 <- readRDS("./data/gam_fig2_model.rds")
n_sim <- 1000

beta <- coef(gam_fig2)
Vb <- vcov(gam_fig2, unconditional = TRUE)

sim_beta <- MASS::mvrnorm(n_sim, beta, Vb)

new_in <- matchit_data2
new_out <- matchit_data2

new_in$treat <- factor(
  "Inside PA",
  levels = levels(matchit_data2$treat)
)

new_out$treat <- factor(
  "Outside PA",
  levels = levels(matchit_data2$treat)
)

X_in <- predict(gam_fig2, newdata = new_in, type = "lpmatrix")
X_out <- predict(gam_fig2, newdata = new_out, type = "lpmatrix")

pred_in_sims <- plogis(X_in %*% t(sim_beta))
pred_out_sims <- plogis(X_out %*% t(sim_beta))

diff_sims <- pred_in_sims - pred_out_sims


decile_keys <- matchit_data2 %>%
  distinct(carbon_decile) %>%
  arrange(carbon_decile)

fig2_list <- lapply(seq_len(nrow(decile_keys)), function(i) {

  decile_i <- decile_keys$carbon_decile[i]

  idx <- which(
    matchit_data2$carbon_decile == decile_i &
      matchit_data2$treat == "Inside PA"
  )

  if (length(idx) < 20) return(NULL)

  w <- matchit_data2$weights[idx]

  in_sim <- apply(
    pred_in_sims[idx, , drop = FALSE],
    2,
    weighted.mean,
    w = w,
    na.rm = TRUE
  )

  out_sim <- apply(
    pred_out_sims[idx, , drop = FALSE],
    2,
    weighted.mean,
    w = w,
    na.rm = TRUE
  )

  diff_sim <- in_sim - out_sim

  # simulation relative difference
  rel_sim <- diff_sim / pmax(out_sim, 1e-6) * 100

  data.frame(
    carbon_decile = decile_i,
    n = length(idx),

    pred_inside = mean(in_sim),
    inside_lwr = quantile(in_sim, 0.025, na.rm = TRUE),
    inside_upr = quantile(in_sim, 0.975, na.rm = TRUE),

    pred_outside = mean(out_sim),
    outside_lwr = quantile(out_sim, 0.025, na.rm = TRUE),
    outside_upr = quantile(out_sim, 0.975, na.rm = TRUE),

    diff = mean(diff_sim),
    diff_lwr = quantile(diff_sim, 0.025, na.rm = TRUE),
    diff_upr = quantile(diff_sim, 0.975, na.rm = TRUE),

    rel_diff = mean(rel_sim, na.rm = TRUE),
    rel_diff_lwr = quantile(rel_sim, 0.025, na.rm = TRUE),
    rel_diff_upr = quantile(rel_sim, 0.975, na.rm = TRUE),

    rel_diff_plugin = (mean(in_sim) - mean(out_sim)) / pmax(mean(out_sim), 1e-6) * 100
  )
})

fig2_df <- bind_rows(fig2_list) %>%
  mutate(
    carbon_decile = factor(
      carbon_decile,
      levels = rev(decile_labs)
    ),
    significant_abs = ifelse(diff_lwr > 0 | diff_upr < 0,
                             "Significant", "Not significant"),
    significant_rel = ifelse(rel_diff_lwr > 0 | rel_diff_upr < 0,
                             "Significant", "Not significant")
  )

fwrite(fig2_df, "./data/fig2_df.csv", row.names = FALSE)

fig2_pred_long <- fig2_df %>%
  select(
    carbon_decile,
    pred_inside, inside_lwr, inside_upr,
    pred_outside, outside_lwr, outside_upr
  ) %>%
  tidyr::pivot_longer(
    cols = -carbon_decile,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    pa_status = case_when(
      grepl("inside", metric) ~ "Inside PA",
      grepl("outside", metric) ~ "Outside PA"
    ),
    stat = case_when(
      grepl("^pred", metric) ~ "estimate",
      grepl("lwr", metric) ~ "lwr",
      grepl("upr", metric) ~ "upr"
    )
  ) %>%
  select(carbon_decile, pa_status, stat, value) %>%
  tidyr::pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  mutate(
    pa_status = factor(
      pa_status,
      levels = c("Outside PA", "Inside PA")
    )
  )
fwrite(fig2_pred_long, "./data/fig2_pred_long.csv", row.names = FALSE)


fig2_pred_long <- fread("./data/fig2_pred_long.csv")

pd <- position_dodge(width = 0.55)
fig2_pred_long$carbon_decile <- factor(
  fig2_pred_long$carbon_decile,
  levels = rev(decile_labs)
)
fig2a <- ggplot(
  fig2_pred_long,
  aes(
    x = estimate,
    y = carbon_decile,
    color = pa_status
  )
) +
  geom_point(
    size = 2.5,
    position = pd
  ) +
  geom_errorbarh(
    aes(xmin = lwr, xmax = upr),
    height = 0.15,
    linewidth = 0.65,
    position = pd
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_color_manual(
    values = c(
      "Outside PA" = "#4f81bd",
      "Inside PA" = "#c43c2f"
    ),
    labels = c(
      "Outside PA" = "Control",
      "Inside PA" = "Inside PA"
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    # legend.position = "top"
  ) +
  labs(
    x = "Predicted probability of forest pest records (%)",
    y = "Carbon percentile within ecoregion"
  )

fig2a


fig2_df <- fread("./data/fig2_df.csv")
fig2_df$carbon_decile <- factor(
  fig2_df$carbon_decile,
  levels = rev(decile_labs)
)

fig2_df <- fig2_df %>%
  mutate(
    rel_diff = (pred_inside - pred_outside) / pred_outside * 100,
    
    # 下面这个 CI 是近似的，不是严格 simulation CI
    rel_diff_lwr = (inside_lwr - outside_upr) / outside_upr * 100,
    rel_diff_upr = (inside_upr - outside_lwr) / outside_lwr * 100,
    
    significant_rel = ifelse(rel_diff_lwr > 0 | rel_diff_upr < 0,
                             "Significant", "Not significant")
  )


fig2b <- ggplot(
  fig2_df,
  aes(
    x = rel_diff,
    y = carbon_decile
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "black",
    linewidth = 0.5
  ) +
  geom_point(
    aes(color = significant_rel),
    size = 2.5
  ) +
  geom_errorbarh(
    aes(xmin = rel_diff_lwr, xmax = rel_diff_upr, color = significant_rel),
    height = 0.15,
    linewidth = 0.65
  ) +
  scale_x_continuous(
    labels = function(x) paste0(round(x, 0), "%"),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_color_manual(
    values = c(
      "Significant" = "black",
      "Not significant" = "grey60"
    ),
    guide = "none"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 11)
  ) +
  labs(
    x = "Difference in pest record probability (%)",
    # y = 'Carbon percentile within biome'
    y=''
  )

fig2b


fig2a + fig2b +
  plot_layout(
    ncol = 2,
    widths = c(1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title = element_blank()
  )


#==sup matching quality=============================
bal_plot_df <- fread("./data/balance_plot_data.csv")

covariate_labels <- c(
  "carbon_percentile" = "Forest carbon",
  "treeCover"         = "Tree canopy cover",
  "human_foot"        = "Human Footprint Index",
  "bio01"             = "Mean temperature (BIO1)",
  "bio02"             = "Mean diurnal temperature range (BIO2)",
  "bio04"             = "Temperature seasonality (BIO4)",
  "bio15"             = "Precipitation seasonality (BIO15)",
  "ele"               = "Elevation",
  "FFI"               = "Forest Fragmentation Index"
)

bal_plot_df <- bal_plot_df %>%
  mutate(
    covariate_full = recode(covariate, !!!covariate_labels),
    covariate_full = factor(
      covariate_full,
      levels = rev(unname(covariate_labels))
    )
  )


figs4 <- ggplot(
  bal_plot_df,
  aes(x = abs(smd), y = covariate_full, color = stage)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "grey50"
  ) +
  geom_vline(
    xintercept = 0.25,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_point(size = 2.6) +
  scale_color_manual(
    values = c(
      "Before matching" = "#c43c2f",
      "After matching" = "#4f81bd"
    )
  ) +
  theme_bw() +
  labs(
    x = "Absolute standardized mean difference",
    y = NULL,
    color = NULL
  ) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# p_balance

balance_summary <- bal_plot_df %>%
  group_by(stage) %>%
  summarise(
    mean_abs_smd = mean(abs(smd), na.rm = TRUE),
    median_abs_smd = median(abs(smd), na.rm = TRUE),
    max_abs_smd = max(abs(smd), na.rm = TRUE),
    n_over_025 = sum(abs(smd) > 0.25, na.rm = TRUE),
    n_over_01 = sum(abs(smd) > 0.1, na.rm = TRUE),
    .groups = "drop"
  )

balance_summary


#============sup fig5===================
area_effect_decile <- fread('./data/area_effect_by_carbon_decile.csv')

p_area_decile <- ggplot(
  area_effect_decile,
  aes(x = carbon_decile, y = protected_forest_area_mha/100)
) +
  geom_col(fill = "grey55", width = 0.7) +
  theme_bw() +
  labs(
    x = "Carbon percentile within biome",
    y = "Matched protected forest area (Mha)"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_area_decile


p_absdiff_decile <- ggplot(
  area_effect_decile,
  aes(x = carbon_decile, y = mean_abs_diff)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_point(size = 2.6, color = "#756BB1") +
  geom_line(aes(group = 1), linewidth = 0.8, color = "#2C7FB8") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  theme_bw() +
  labs(
    # x = "Carbon percentile within biome",
    y = "Absolute probability difference"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
  )


p_absdiff_decile


biome_effect_summary <- fread('/root/autodl-tmp/pests/protect_pests/data/biome_effect_summary.csv')
forest_biomes <- c(
  "Tropical & Subtropical Moist Broadleaf Forests",
  "Tropical & Subtropical Dry Broadleaf Forests",
  "Tropical & Subtropical Coniferous Forests",
  "Temperate Broadleaf & Mixed Forests",
  "Temperate Conifer Forests",
  "Boreal Forests/Taiga",
  "Mediterranean Forests, Woodlands & Scrub",
  "Mangroves"
)

biome_lookup <- tibble(
  ecoregion_code = 1:14,
  biome_name = c(
    "Tropical & Subtropical Moist Broadleaf Forests",
    "Tropical & Subtropical Dry Broadleaf Forests",
    "Tropical & Subtropical Coniferous Forests",
    "Temperate Broadleaf & Mixed Forests",
    "Temperate Conifer Forests",
    "Boreal Forests/Taiga",
    "Tropical & Subtropical Grasslands, Savannas & Shrublands",
    "Temperate Grasslands, Savannas & Shrublands",
    "Flooded Grasslands & Savannas",
    "Montane Grasslands & Shrublands",
    "Tundra",
    "Mediterranean Forests, Woodlands & Scrub",
    "Mangroves",
    "Deserts & Xeric Shrublands"
    
  )
)

biome_effect_summary_forest <- biome_effect_summary %>%
  mutate(
    ecoregion_code = as.integer(as.character(ecoregion))
  ) %>%
  left_join(biome_lookup, by = "ecoregion_code") %>%
  filter(biome_name %in% forest_biomes) %>%
  mutate(
    biome_name = factor(biome_name, levels = forest_biomes)
  ) %>%
  arrange(biome_name)



p_biome_effect <- ggplot(
  biome_effect_summary_forest,
  aes(x = reorder(biome_name, abs_diff), y = abs_diff)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey45"
  ) +
  geom_point(
    aes(size = n),
    color = "#2C7FB8",
    alpha = 0.75,
    # show.legend = F
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_size_continuous(range = c(2, 6)) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Absolute probability difference",
    size = "Matched grids"
  ) +
  theme(
    axis.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    # axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "bottom"
  )

p_biome_effect


figS5 <- (
  p_biome_effect |
    plot_spacer() |
    (p_absdiff_decile / p_area_decile)
) +
  plot_layout(
    widths = c(1, 0.035, 1)
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.01, 0.99)
  )
figS5
