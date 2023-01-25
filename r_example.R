
# Required Packages -------------------------------------------------------
library(tidyverse)
library(here)
library(broom)
library(vegan)

# Import Data -------------------------------------------------------------
## Sites
sites_df <- readr::read_csv(here::here("data",
                          "zms_thesis-sites_2017-06-18.csv"),
                          col_types = cols(
                            station_id = col_character(),
                            lat = col_double(),
                            long = col_double()
                          ))
## Macroinvertebrates
macros_df <- readr::read_csv(here::here("data",
                                       "zms_thesis-macro_2017-06-18.csv"),
                            col_types = cols(
                              unique_id = col_character(),
                              lake = col_character(),
                              station_id = col_character(),
                              sample_number = col_double(),
                              lat = col_double(),
                              long = col_double(),
                              agency_code = col_character(),
                              method = col_character(),
                              date = col_character(),
                              count = col_double(),
                              life_stage = col_character(),
                              final_id = col_character(),
                              taxon_level = col_character(),
                              phylum = col_character(),
                              subphylum = col_character(),
                              class = col_character(),
                              subclass = col_character(),
                              order = col_character(),
                              suborder = col_character(),
                              family = col_character(),
                              subfamily = col_character(),
                              tribe = col_character(),
                              genus = col_character(),
                              species = col_character(),
                              picked = col_double(),
                              squares = col_double(),
                              proportion = col_double()
                            )
                          )

## In Situ Data
insitu_df <- readr::read_csv(here::here("data",
                                       "zms_thesis-env_2017-06-18.csv"),
                             col_types = cols(
                               unique_id = col_character(),
                               wind = col_character(),
                               oncolite = col_character(),
                               disturbance = col_character(),
                               substrate_size_d50 = col_double(),
                               substrate = col_character(),
                               distance_fram_shore = col_double(),
                               temperature = col_double(),
                               conductivity = col_double(),
                               ph = col_double(),
                               macrophyte = col_double()
                             )
                            )

# Joins -------------------------------------------------------------------

all_df <- purrr::reduce(
  .x = list(sites_df, macros_df, insitu_df),
  left_join
) |>
  mutate(lake = case_when(
    lake %in% "caz" ~ "Cazenovia",
    lake %in% "onon" ~ "Onondaga",
    lake %in% "ot" ~ "Otisco",
    TRUE ~ "ERROR"
  ),
  lake = factor(lake, levels = c("Onondaga", "Otisco", "Cazenovia"))
  )


# Macro Matrix ------------------------------------------------------------

macros_wide <- macros_df |>
  select(unique_id,
         final_id,
         count = picked) |>
  distinct() |>
  pivot_wider(
    names_from = final_id,
    values_from = count
  )

macros_wide[is.na(macros_wide)] <- 0

macros_matrix <- as.matrix(macros_wide[!names(macros_wide) %in% "unique_id"])
rownames(macros_matrix) <- macros_wide$unique_id
# Plots -------------------------------------------------------------------

(
  spcond_plot <- all_df |>
    select(station_id,
           lake,
           conductivity) |>
    distinct() |>
    tidyr::drop_na() |>
    ggplot(aes(x = lake, y = conductivity)) +
    geom_boxplot() +
    xlab("Lake") +
    ylab("Specific Conductance (µS/cm at 25°C)") +
    theme_bw()
)



# Heatmap -----------------------------------------------------------------

plotree <- hclust(vegdist(macros_matrix), "average")
sptree <- hclust(vegdist(t(macros_matrix), "raup"), "average")
heatmap_plot <- tabasco(macros_matrix, plotree, sptree)

# ANOVA -------------------------------------------------------------------

order_df <- all_df |>
  select(
    unique_id,
    station_id,
    lake,
    order,
    count = picked
  ) |>
  group_by(station_id, unique_id, lake, order) |>
  summarize(count = sum(count),
            .groups = "drop") |>
  group_by(station_id, lake) |>
  mutate(total = sum(count),
         relative_abundance = count / total * 100) |>
  ungroup() |>
  complete(order, nesting(lake, unique_id, station_id, total))

order_df[is.na(order_df)] <- 0

mean_rel_abund <- order_df |>
  group_by(station_id, lake, order) |>
  summarize(relative_abundance = mean(relative_abundance),
            .groups = "drop")

cote_df <- mean_rel_abund |>
  filter(order %in% c("coleoptera", "odonata", "ephemeroptera", "trichoptera")) |>
  group_by(lake, station_id) |>
  summarize(relative_abundance = sum(relative_abundance),
          .groups = "drop")


anova <- aov(formula = relative_abundance ~ lake, data = cote_df)
summary(anova)
anova_tidy <- tidy(anova)
(tukey <- TukeyHSD(anova))
tukey_tidy <- tidy(tukey)

(cote_plot <- ggplot(cote_df, aes(lake, relative_abundance)) +
  geom_boxplot() +
  xlab("Lake") +
  ylab("Relative Abundance (%)") +
  theme_bw())
