library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(forcats)
library(countrycode)
library(stringr)
library(patchwork)

# =========================
# Load data
# =========================
attr_summary <- fread("/root/autodl-tmp/pests/protect_pests/data/fig3_attribute_pest_visibility_summary.csv")
size_limits <- range(attr_summary$protected_forest_area_mha, na.rm = TRUE)
size_breaks <- pretty(size_limits, n = 3)

common_size_scale <- scale_size_continuous(
  range = c(2.5, 7),
  limits = size_limits,
  breaks = size_breaks,
  labels = scales::comma,
  name = "Forest area (Mha)"
)

# =========================
# Fig. 3a Governance
# =========================

p3a_gov <- attr_summary %>%
  filter(attribute == "Governance") %>%
  mutate(
    level = forcats::fct_reorder(level, area_weighted_pest_prob)
  ) %>%
  ggplot(
    aes(
      x = area_weighted_pest_prob,
      y = level
    )
  ) +
  geom_point(
    aes(size = protected_forest_area_mha),
    color = "#2C7FB8",
    alpha = 0.85
  ) +
  common_size_scale +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Governance"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

# =========================
# Fig. 3b IUCN category
# =========================

iucn_order <- c(
  "Ia", "Ib", "II", "III",
  "IV", "V", "VI", "Not reported"
)

p3b_iucn <- attr_summary %>%
  filter(attribute == "IUCN category") %>%
  mutate(
    level = factor(level, levels = rev(iucn_order))
  ) %>% 
  ggplot(
    aes(
      x = area_weighted_pest_prob,
      y = level
    )
  ) +
  geom_point(
    aes(size = protected_forest_area_mha),
    color = "#2C7FB8",
    alpha = 0.85
  ) +
  common_size_scale +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Pest-record probability",
    y = NULL,
    title = "IUCN category"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )
p3b_iucn

# =========================
# fig3c  Country-level summary
# =========================
country_visibility <- fread("/root/autodl-tmp/pests/protect_pests/data/fig3_country_pest_visibility_summary.csv")
country_top15 <- country_visibility %>%
  slice_head(n = 15) %>%
  mutate(
    country = forcats::fct_reorder(
      country,
      elevated_visibility_area_mha
    )
  )

p3c_country <- ggplot(
  country_top15,
  aes(
    x = elevated_visibility_area_mha,
    y = country
  )
) +
  geom_col(
    fill = "#D7301F",
    alpha = 0.85,
    width = 0.72
  ) +
  scale_x_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = "High-probability forest area (Mha)",
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p3c_country

p3a_gov <- p3a_gov +
  scale_x_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  theme(
    plot.margin = margin(t = 5.5, r = 12, b = 5.5, l = 5.5)
  )



(p3a_gov / p3b_iucn) +
  plot_layout(
    # widths = c(1.2, 1),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.01, 0.99),
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  ) | p3c_country


#============sup fig6===================
pa_level_fli <- fread(
  "./data/pa_level_FLII.csv"
)
p_iucn_fli <- pa_level_fli %>%
  mutate(
    IUCN_group = factor(
      IUCN_group,
      levels = rev(iucn_order)
    )
  ) %>%
  ggplot(
    aes(
      x = mean_FLII,
      y = IUCN_group
    )
  ) +
  geom_boxplot(
    fill = "#BFD3E6",
    color = "#2C7FB8",
    width = 0.65,
    outlier.shape = NA,
    linewidth = 0.35
  ) +
  geom_point(
    data = iucn_fli_summary %>%
      mutate(
        IUCN_group = factor(
          IUCN_group,
          levels = rev(iucn_order)
        )
      ),
    aes(
      x = area_weighted_FLII,
      y = IUCN_group,
      size = protected_forest_area_mha
    ),
    inherit.aes = FALSE,
    color = "#D7301F",
    alpha = 0.9
  ) +
  scale_x_continuous(
    limits = c(0, 10),
    breaks = seq(0, 10, by = 2),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_size_continuous(
    range = c(2.5, 7),
    name = "Forest area (Mha)"
  ) +
  labs(
    x = "Forest Landscape Integrity Index",
    y = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

p_iucn_fli

pred_pest_fli <- fread(
  "./data/predicted_pest_probability_by_FLII.csv"
)

p_fli_pest <- ggplot(
  pred_pest_fli,
  aes(
    x = FLII,
    y = estimate
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    fill = "#9ECAE1",
    alpha = 0.45
  ) +
  geom_line(
    color = "#2C7FB8",
    linewidth = 0.9
  ) +
  scale_x_continuous(
    breaks = seq(0, 10, by = 2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  labs(
    x = "Forest Landscape Integrity Index",
    y = "Pest-record probability"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

p_fli_pest


pred_effort_fli <- fread(
  "./data/predicted_sampling_effort_by_FLII.csv"
)
p_fli_effort <- ggplot(
  pred_effort_fli,
  aes(
    x = FLII,
    y = estimate
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    fill = "#CBC9E2",
    alpha = 0.5
  ) +
  geom_line(
    color = "#756BB1",
    linewidth = 0.9
  ) +
  scale_x_continuous(
    breaks = seq(0, 10, by = 2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Forest Landscape Integrity Index",
    y = "Background sampling effort"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

p_fli_effort


fli_supp_figure <-
  p_iucn_fli |
  (p_fli_pest / p_fli_effort)

fli_supp_figure <-
  fli_supp_figure +
  plot_layout(
    widths = c(1.05, 1)
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(
      # face = "bold",
      size = 13
    ),
    plot.tag.position = c(0.01, 0.99)
  )

fli_supp_figure





