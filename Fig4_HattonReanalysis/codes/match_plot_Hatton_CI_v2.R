library(brms)
library(posterior)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(data.table)
library(zoo)
library(scales)

# Load models
rds_path <- "rds/Hatton_v1.2.rds"
results <- readRDS(rds_path)

# Prepare simulation range
logN_grid <- seq(-2, 4, length.out = 300)
N_grid <- 10^logN_grid

# Function to compute predicted growth from posterior draws
get_curve_draws <- function(model, dataset_name, x_range) {
  draws <- as_draws_df(model, variables = c("b_r_Intercept", "b_a_Intercept", "b_K_Intercept"))
  if (nrow(draws) == 0) return(NULL)
  
  draws <- draws %>%
    mutate(draw_id = row_number()) %>%
    filter(is.finite(b_r_Intercept), is.finite(b_a_Intercept), is.finite(b_K_Intercept))
  
  # Generate grid within empirical range of x_gm
  logN_grid <- seq(log10(x_range[1]), log10(x_range[2]), length.out = 300)
  N_grid <- 10^logN_grid
  
  expand_grid(draw_id = draws$draw_id, N = N_grid) %>%
    left_join(draws, by = "draw_id") %>%
    mutate(
      sign_a = ifelse(b_a_Intercept >= 0, 1, -1),
      dNdt = b_r_Intercept * sign_a * N * (1 - (N / b_K_Intercept)^b_a_Intercept),
      log10_N = log10(N),
      log10_dNdt = log10(pmax(dNdt, 1e-8))  # avoid -Inf
    ) %>%
    filter(is.finite(log10_dNdt), log10_dNdt > -5, log10_dNdt < 5) %>%
    group_by(log10_N) %>%
    summarise(
      mean_raw = mean(log10_dNdt, na.rm = TRUE),
      lower_raw = quantile(log10_dNdt, 0.025, na.rm = TRUE),
      upper_raw = quantile(log10_dNdt, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mean = zoo::rollapply(mean_raw, width = 7, FUN = mean, fill = NA, align = "center"),
      lower = zoo::rollapply(lower_raw, width = 7, FUN = mean, fill = NA, align = "center"),
      upper = zoo::rollapply(upper_raw, width = 7, FUN = mean, fill = NA, align = "center"),
      Dataset = dataset_name
    )
}

# Load empirical data
raw <- fread("data/fig4.csv")
groups <- list(
  Mammals = raw %>% filter(Taxa == "Large mammal herbivores"),
  Seagrass = raw %>% filter(Taxa == "Seagrass"),
  Fish = raw %>% filter(Taxa == "Fish"),
  Broadleaf = raw %>% filter(Taxa == "Broadleaf"),
  Coniferous = raw %>% filter(Taxa == "Needleleaf"),
  Grassland = raw %>% filter(Taxa %in% c("Pinales: Cupressaceae", "Pinales: Pinaceae")),
  MarineBacteria = raw %>% filter(Taxa == "Marine sediment bacteria")
)

# Build curve data with ribbon from all datasets
curve_with_ribbon <- map_dfr(names(results), function(name) {
  model <- results[[name]]
  if (is.null(model)) return(NULL)
  
  x_vals <- groups[[name]]$x_gm %>%
    as.numeric() %>%
    .[is.finite(.) & . > 0]
  
  if (length(x_vals) == 0) return(NULL)
  x_range <- range(x_vals)
  
  get_curve_draws(model, name, x_range)
})

# Extract raw points
point_data <- map_dfr(names(groups), function(name) {
  df <- groups[[name]]
  df %>%
    mutate(
      log10_x = log10(as.numeric(x_gm)),
      log10_y = log10(as.numeric(y_gm)),
      Dataset = name
    ) %>%
    filter(is.finite(log10_x), is.finite(log10_y)) %>%
    select(Dataset, log10_x, log10_y)
})

# Prepare transformed versions of curve data
curve_plot_data <- curve_with_ribbon %>%
  mutate(
    N = 10^log10_N,
    dNdt_mean = 10^mean,
    dNdt_lower = 10^lower,
    dNdt_upper = 10^upper
  )

point_plot_data <- point_data %>%
  mutate(
    N = 10^log10_x,
    dNdt = 10^log10_y
  )

# Log-log plot
ggplot() +
  geom_ribbon(
    data = curve_plot_data,
    aes(x = N, ymin = dNdt_lower, ymax = dNdt_upper, fill = Dataset),
    alpha = 0.2
  ) +
  geom_line(
    data = curve_plot_data,
    aes(x = N, y = dNdt_mean, color = Dataset),
    linewidth = 1
  ) +
  geom_point(
    data = point_plot_data,
    aes(x = N, y = dNdt, color = Dataset),
    size = 1.5, stroke = 0, alpha = 0.25
  ) +
  labs(
    x = "Biomass",
    y = "Production ",
    title = "Refitting Hatton et al. data"
  ) +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    minor_breaks = log10(as.vector(outer(1:9, 10^(-2:4), `*`)))
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    minor_breaks = log10(as.vector(outer(1:9, 10^(-2:4), `*`)))
  ) +
  annotation_logticks(sides = "bl") +  # adds log tick marks on bottom and left
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  )

#
version_tag <- tools::file_path_sans_ext(basename(rds_path))
ggsave(paste0("plots/fit_", version_tag, ".jpg"), width = 6, height = 5)
ggsave(paste0("plots/fit_", version_tag, ".pdf"), width = 6, height = 5)
