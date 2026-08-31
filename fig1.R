.libPaths('/home/qyx/R/x86_64-pc-linux-gnu-library/4.3')
#=================load library=================
library(terra)
library(dplyr)
library(data.table)
library(ggplot2)
library(sf)
library(tidyr)
library(stringr)
library(purrr)
library(scales)
library(tidyterra)
library(patchwork)


#=================load data=================
continentPolygon <- rnaturalearthdata::coastline50 %>% vect()
crs <- '+proj=longlat +datum=WGS84'
globalRaster <- rast(vals=NA,nrows=1800, ncols=3600,xmin=-180, xmax=180,ymin=-90, ymax=90,crs=crs)

#plot
treeCover   <- rast('./data/hasenForestCover_resample.tif')
pa_fraction <- rast('./data/pa_fraction_combined.tif')

pa_mask <- ifel(pa_fraction > 0, 1, NA)
# protected forest
pa_forest_mask <- ifel(treeCover > 10 & pa_fraction > 0, 1, NA)

pa_forest_rast <- c(pa_mask,pa_forest_mask) |> sum(na.rm=T)

greenland <- rnaturalearth::ne_countries(country = "Greenland", scale = "medium", returnclass = "sf")
pa_forest_rast_filter <- mask(pa_forest_rast, greenland,inverse=T)

pa_forest_only <- terra::ifel(
  pa_forest_rast_filter == 2,
  2,
  NA
)


fig1a <- ggplot() +
  geom_spatvector(
    data = continentPolygon,
    fill = "grey92",
    color = "grey75",
    linewidth = 0.12
  ) +
  geom_spatraster(
    data = pa_forest_only
  ) +
  scale_fill_gradientn(
    name = NULL,
    colours = c("#5E7D4A", "#5E7D4A"),
    limits = c(1.5, 2.5),
    breaks = 2,
    labels = "Protected forests",
    na.value = "transparent",
    guide = guide_legend(
      order = 1,
      direction = "horizontal"
    )
  ) +
  geom_point(
    data = pest_df,
    aes(x = x, y = y, color = pest_class_num),
    shape = 15,
    size = 0.1,
    alpha = 0.95
  ) +
  scale_color_gradientn(
    colors = pest_cols,
    name = " Forest pest records",
    limits = c(1, 4),
    breaks = c(1, 2, 3, 4),
    labels = c("1", "10", "20", ">20"),
    guide = guide_colorbar(
      title.position = "right",
      label.position = "bottom",
      direction = "horizontal",
      barwidth = unit(2.2, "cm"),
      barheight = unit(0.25, "cm")
    )
  ) +
  coord_sf(
    crs = crs,
    xlim = c(-160, 163.5),
    ylim = c(-56, 83),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_blank(),

    legend.position = "bottom",
    legend.box = "horizontal",

    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

#fig1b
compotion_table <- fread('/root/autodl-tmp/pests/protect_pests/data/compotion_table.csv')

decile_labs <- c(
  "0–10", "10–20", "20–30", "30–40", "40–50",
  "50–60", "60–70", "70–80", "80–90", "90–100"
)

pa_forest <- compotion_table %>%
  filter(
    pa_fraction > 0,
    treeCover > 10
  ) %>%
  mutate(
    allPest = Animals_records + Dieases_records + Plants_records,
    any_pest = as.integer(allPest > 0),
    carbon_decile = cut(
      carbon_percentile,
      breaks = seq(0, 1, by = 0.1),
      include.lowest = TRUE,
      labels = decile_labs
    )
  )

fig1b_df <- pa_forest %>%
  group_by(carbon_decile) %>%
  summarise(
    n = n(),
    prop_with_pest = mean(any_pest, na.rm = TRUE),
    se = sqrt(prop_with_pest * (1 - prop_with_pest) / n),
    ci_low = pmax(0, prop_with_pest - 1.96 * se),
    ci_high = pmin(1, prop_with_pest + 1.96 * se),
    .groups = "drop"
  )


fig1c_df <- pa_forest %>%
  filter(allPest > 0) %>%
  mutate(
    log_records = log10(allPest + 1)
  ) %>%
  group_by(carbon_decile) %>%
  summarise(
    n_positive = n(),
    median_log_records = median(log_records, na.rm = TRUE),
    q25 = quantile(log_records, 0.25, na.rm = TRUE),
    q75 = quantile(log_records, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

fig1b_df$carbon_decile <- factor(fig1b_df$carbon_decile, levels = rev(decile_labs))
pa_forest$carbon_decile <- factor(pa_forest$carbon_decile, levels = rev(decile_labs))

fig1b_h <- ggplot(fig1b_df, aes(y = carbon_decile, x = prop_with_pest)) +
  geom_col(
    width = 0.75,
    fill = "#810F7C",
    color = "white",
    linewidth = 0.25
  ) +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    height = 0.18,
    linewidth = 0.35,
    color = "grey25"
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11)
  ) +
  labs(
    x = "Protected forest grids with pest records (%)",
    y = "Carbon percentile within biome"
  )

fig1c_h <- pa_forest %>%
  filter(allPest > 0) %>%
  mutate(
    log_records = log10(allPest + 1)
  ) %>%
  ggplot(aes(y = carbon_decile, x = log_records)) +
  geom_boxplot(
    fill = "#8C96C6",
    color = "#4D004B",
    width = 0.65,
    outlier.shape = NA,
    linewidth = 0.35
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11)
  ) +
  labs(
    x = "Observation intensity",
    y = NULL
  )+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

(fig1b_h | fig1c_h) +
  plot_layout(widths = c(1, 1.1))

#==========sup1=================
pest_var <- "allPest"
pa_forest <- compotion_table %>%
  filter(pa_fraction > 0.5, treeCover > 10) %>%
  mutate(
    pest_value = .data[[pest_var]],
    any_pest = as.integer(pest_value > 0),
    gbif_effort_nonpest = pmax(gbif_effort - allPest, 0),

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

figS1a <- ggplot() +
  geom_spatvector(
    data = continentPolygon,
    fill = "grey92",
    color = "grey75",
    linewidth = 0.1
  ) +
  geom_raster(
    data = pa_forest,
    aes(
      x = x,
      y = y,
      fill = log10(gbif_effort_nonpest + 1)
    )
  ) +
  scale_fill_viridis_c(
    name = "log10(GBIF sampling effort + 1)",
    na.value = "transparent"
  ) +
  coord_sf(
    crs = crs,
    xlim = c(-160, 163.5),
    ylim = c(-56, 83),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) +
  labs(x = NULL, y = NULL)

figS1a


sampling_presence_summary <- pa_forest %>%
  mutate(
    any_sampling = as.integer(gbif_effort_nonpest > 0)
  ) %>%
  group_by(carbon_decile) %>%
  summarise(
    n = n(),
    prop_with_sampling = mean(any_sampling, na.rm = TRUE),
    se = sqrt(prop_with_sampling * (1 - prop_with_sampling) / n),
    ci_low = pmax(0, prop_with_sampling - 1.96 * se),
    ci_high = pmin(1, prop_with_sampling + 1.96 * se),
    .groups = "drop"
  )

figS1b <- ggplot(
  sampling_presence_summary,
  aes(x = carbon_decile, y = prop_with_sampling)
) +
  geom_col(
    width = 0.75,
    fill = "#3182bd",
    color = "white",
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.25,
    linewidth = 0.35
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Carbon percentile within ecoregion",
    y = "Protected forest grids with GBIF records (%)"
  )

figS1b


sampling_intensity_summary <- pa_forest %>%
  group_by(carbon_decile) %>%
  summarise(
    n = n(),
    median_effort = median(log10(gbif_effort_nonpest + 1), na.rm = TRUE),
    q25 = quantile(log10(gbif_effort_nonpest + 1), 0.25, na.rm = TRUE),
    q75 = quantile(log10(gbif_effort_nonpest + 1), 0.75, na.rm = TRUE),
    .groups = "drop"
  )

figS1c <- ggplot(
  sampling_intensity_summary,
  aes(x = carbon_decile, y = median_effort, group = 1)
) +
  geom_ribbon(
    aes(ymin = q25, ymax = q75),
    alpha = 0.25
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Carbon percentile within ecoregion",
    y = "Observation effort"
  )

figS1c


#==========sup2=================
compotion_table <- fread('./data/compotion_table.csv')
compotion_table <- compotion_table %>%
  mutate(
    allPest = Animals_records + Dieases_records + Plants_records
  )

  pa_biome <- compotion_table %>%
  filter(
    pa_fraction > 0,
    treeCover > 10,
    ecoregion!='N/A',
    !is.na(carbon)
  ) %>%
  group_by(ecoregion) %>%
  ungroup() %>%
  mutate(
    carbon_class = case_when(
      carbon_percentile <= 0.25 ~ 1,
      carbon_percentile <= 0.50 ~ 2,
      carbon_percentile <= 0.75 ~ 3,
      carbon_percentile <= 1.00 ~ 4,
      TRUE ~ NA_real_
    ),
    pest_class = case_when(
      allPest == 0 ~ 1,
      allPest == 1 ~ 2,
      allPest <= 20 ~ 3,
      allPest > 20 ~ 4,
      TRUE ~ NA_real_
    ),
    group = paste(carbon_class, "-", pest_class),
    biome_name = as.character(ecoregion)
  )


bivariate_color_scale_4 <- tibble(
  "1 - 1" = "#f7f7f7",
  "2 - 1" = "#d9f0d3",
  "3 - 1" = "#a6dba0",
  "4 - 1" = "#5aae61",

  "1 - 2" = "#e7d4e8",
  "2 - 2" = "#c2e5d3",
  "3 - 2" = "#80cdc1",
  "4 - 2" = "#35978f",

  "1 - 3" = "#c2a5cf",
  "2 - 3" = "#a6bddb",
  "3 - 3" = "#74add1",
  "4 - 3" = "#2b8cbe",

  "1 - 4" = "#9970ab",
  "2 - 4" = "#8073ac",
  "3 - 4" = "#6a51a3",
  "4 - 4" = "#3f007d"
) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "group",
    values_to = "fill"
  )


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


pa_biome_plot <- pa_biome %>%
  left_join(bivariate_color_scale_4, by = "group") %>%
  filter(
    !is.na(fill),
    biome_name %in% forest_biomes,
    biome_name != "N/A"
  ) %>%
  mutate(
    biome_name = factor(biome_name, levels = forest_biomes)
  )

biome_labels <- setNames(
  paste0(letters[seq_along(forest_biomes)], ")"),
  forest_biomes
)



ggplot() +
  geom_spatvector(
    data = continentPolygon,
    fill = "grey94",
    color = "grey80",
    linewidth = 0.08
  ) +
  geom_point(
    data = pa_biome_plot,
    aes(x = x, y = y, color = fill),
    shape = 15,
    size = 0.08,
    alpha = 0.95
  ) +
  scale_color_identity() +
  facet_wrap(
    ~ biome_name,
    ncol = 2,
    labeller = labeller(biome_name = biome_labels)
  ) +
  coord_sf(
    crs = crs,
    xlim = c(-160, 163.5),
    ylim = c(-56, 83),
    expand = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.placement = "inside",
    strip.clip = "off",
    strip.text = element_text(
      size = 8,        
      face = "bold",
      hjust = 0,         
      vjust = 1,       
      margin = margin(t = 1, r = 0, b = 0, l = 2)
    ),

    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),

    panel.spacing = unit(0.1, "lines"),

    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )


ggplot() +
  geom_spatvector(
    data = continentPolygon,
    fill = "grey94",
    color = "grey80",
    linewidth = 0.08
  ) +
  geom_point(
    data = pa_biome_plot,
    aes(x = x, y = y, color = fill),
    shape = 15,
    size = 0.08,
    alpha = 0.95
  ) +
  scale_color_identity() +
  facet_wrap(
    ~ biome_name,
    ncol = 3
  ) +
  coord_sf(
    crs = crs,
    xlim = c(-160, 163.5),
    ylim = c(-56, 83),
    expand = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.placement = "inside",
    strip.clip = "off",
    strip.text = element_text(
      size = 8,          
      face = "bold",
      hjust = 0,     
      vjust = 1,      
      margin = margin(t = 1, r = 0, b = 0, l = 2)
    ),

    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),

    panel.spacing = unit(0.08, "lines"),

    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )


legend_df <- bivariate_color_scale_4 %>%
  tidyr::separate(
    group,
    into = c("carbon_class", "pest_class"),
    sep = " - ",
    convert = TRUE
  )

bivar_legend <- ggplot(legend_df) +
  geom_tile(
    aes(x = carbon_class, y = pest_class, fill = fill),
    color = "white",
    linewidth = 0.25
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("0–25", "25–50", "50–75", "75–100"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = 1:4,
    labels = c("0", "1", "2–20", ">20"),
    expand = c(0, 0)
  ) +
  coord_fixed() +
    labs(
    x = "Carbon percentile⟶",
    y = "Pest records⟶"
  )  +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_blank(),
    plot.margin = margin(2, 2, 2, 2)
  )

bivar_legend


#==============sup3==================
biome_fig1b_df <- fread("./data/biome_carbon_decile_pest_coverage.csv")
# biome 排序：按 protected forest grid 数从大到小
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

label_data <- tibble(
  biome_name = forest_biomes,
  label_name = paste0(letters[seq_along(forest_biomes)], ")")
)

biome_fig1b_df_forest <- biome_fig1b_df %>%
  mutate(
    biome_name = as.character(biome_name)
  ) %>%
  filter(
    biome_name %in% forest_biomes,
    biome_name != "N/A"
  ) %>%
  left_join(label_data, by = "biome_name") %>%
  mutate(
    biome_name = factor(biome_name, levels = forest_biomes),
    carbon_decile = factor(carbon_decile, levels = decile_labs),
    label_name = factor(label_name, levels = label_data$label_name)
  )


figS3 <- ggplot(
  biome_fig1b_df_forest,
  aes(x = carbon_decile, y = prop_with_pest, group = 1)
) +
  geom_ribbon(
    aes(ymin = ci_low, ymax = ci_high),
    fill = "#8C96C6",
    alpha = 0.25,
    na.rm = TRUE
  ) +
  geom_line(
    color = "#4D004B",
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  geom_point(
    color = "#4D004B",
    size = 0.9,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ biome_name,
    ncol = 2,
    scales = "free_y"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 8.5, face = "bold", hjust = 0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.25, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    # axis.title = element_text(size = 9)
  ) +
  labs(
    x = "Carbon percentile within biome",
    y = "Protected forest grids with pest records (%)"
  )





