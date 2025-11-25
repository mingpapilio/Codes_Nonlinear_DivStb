# Load required libraries
library(ggplot2)
library(readr)
library(dplyr)

# Read the dataset
df <- read_csv("data/mean_microbial_v1.2.csv")

# Add a 'spacer' row for Oscar header
gap_row <- tibble(
  Dataset = "Oscar_header",
  mean_a = NA_real_,
  lower_a = NA_real_,
  upper_a = NA_real_,
  mean_r = NA_real_,
  lower_r = NA_real_,
  upper_r = NA_real_,
  label = "Oscar et al. 1999"
)

# Clean labels for display
df <- df %>%
  mutate(label = case_when(
    Dataset == "FE32" ~ "Our experiment",
    Dataset == "Monod_1941" ~ "Monod 1941",
    Dataset == "Varma_1994" ~ "Varma et al. 1994",
    Dataset == "Deng_2021" ~ "Deng et al. 2021",
    grepl("Oscar", Dataset) ~ gsub("Oscar_1999_batchBo", "Batch ", Dataset),
    TRUE ~ Dataset
  ))

# Combine with gap row
df_plot <- bind_rows(
  df %>% filter(Dataset %in% c("FE32", "Monod_1941", "Varma_1994", "Deng_2021")),
  gap_row,
  df %>% filter(grepl("Oscar_1999_batchBo", Dataset)) %>%
    arrange(match(Dataset, c(
      "Oscar_1999_batchBo54", "Oscar_1999_batchBo62", "Oscar_1999_batchBo63",
      "Oscar_1999_batchBo92", "Oscar_1999_batchBo93", "Oscar_1999_batchBo99",
      "Oscar_1999_batchBo105", "Oscar_1999_batchBo119"
    )))
)

# Explicit order: top to bottom
label_order <- c(
  "Our experiment",
  "Monod 1941",
  "Varma et al. 1994",
  "Deng et al. 2021",
  "Oscar et al. 1999",
  "Batch 54",
  "Batch 62",
  "Batch 63",
  "Batch 92",
  "Batch 93",
  "Batch 99",
  "Batch 105",
  "Batch 119"
)

df_plot$label <- factor(df_plot$label, levels = rev(label_order))
# Calculate global mean of alpha
global_mean <- mean(df$mean_a, na.rm = TRUE)
arrow_y <- length(levels(df$label)) + 1.0

# Plot
ggplot(df_plot, aes(x = mean_a, y = label)) +
  geom_errorbarh(aes(xmin = lower_a, xmax = upper_a), height = 0.3, na.rm = TRUE) +
  geom_point(shape = 15, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 0) +
  annotate("segment",
           x = global_mean, xend = global_mean,
           y = arrow_y, yend = arrow_y - 0.5,
           colour = "red", arrow = arrow(length = unit(0.2, "cm")),
           linewidth = 0.8) +
  labs(
    x = expression(alpha ~ "estimates (95% CI)"),
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold")
  )
#
ggsave("plots/errbar_microbial_v2.pdf", width= 7, height= 7)
