# Load libraries
library(ggplot2)
library(readr)
library(dplyr)

# Read data
df <- read_csv("data/mean_Hatton_v1.2.csv")

# Add sample sizes
df$n <- c(118, 152, 76, 469, 330, 204, 207)

# Create label with n
df <- df %>%
  mutate(label = paste0(Dataset, " (n=", n, ")"))

# Preserve original order (top to bottom)
df$label <- factor(df$label, levels = rev(df$label))

# Calculate global mean of alpha
global_mean <- mean(df$mean_a, na.rm = TRUE)
arrow_y <- length(levels(df$label)) + 0.5

# Plot
ggplot(df, aes(x = mean_a, y = label)) +
  geom_errorbarh(aes(xmin = lower_a, xmax = upper_a), height = 0.2) +
  geom_point(shape = 15, size = 2) +
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

# Save
ggsave("plots/errbar_Hatton_v2.pdf", width=5, height=5)
