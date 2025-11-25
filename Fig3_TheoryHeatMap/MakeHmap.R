library(data.table)
library(ggplot2)
library(grid)
library(scales)  
library(ggpubr)
theme_set(theme_minimal())

# Read in the file
results<- fread("result_array.txt")

# Subset the data for N_species = 25 and N_species = 100
div_1 <- 10
div_2 <- 100
data_1  <- subset(results, N_species == div_1)
data_2 <- subset(results, N_species == div_2)

# Merge the two subsets by the common parameters (alpha and w)
merged <- merge(data_2, data_1, by = c("alpha", "w"), 
                suffixes = c("_2", "_1"))

# Compute the difference: (eigenvalue for N_species = 100) minus (for N_species = 25)
merged$diff <- merged$leading_eigenvalue_2 - merged$leading_eigenvalue_1

# Create a new data frame (hmap) with the renamed columns:
#    b = alpha, c = w, h = diff
hmap <- merged[, c("alpha", "w", "diff")]
names(hmap) <- c("x", "y", "z")

# 1. Define the overall range and the region you want as pure white
hmin <- -0.5
hmax <-  0.5
white_min <- -0.01
white_max <-  0.01

# 2. Create vectors of 'breakpoints' (xvals) and colors (colorvals).
#    We'll map the segment [-0.01, 0.01] entirely to white, 
#    while distributing the rest of the colors around it.

xvals <- c(hmin,   -0.05,  white_min, white_min, white_max, white_max,  0.05,   hmax)
#         ^        ^       ^          ^           ^          ^           ^       ^
#         |        |       identical  identical   identical  identical   |       |
#         lower    small   left edge  left edge   right edge right edge  small   upper
#         limit    negative of white  of white    of white   of white    positive limit

colorvals <- c('#2166ac',    # deep blue
               '#67a9cf',    # lighter blue
               '#d1e5f0',    # very light blue
               'white',      # from -0.01 ...
               'white',      # ... to +0.01 all white
               '#fddbc7',    # light salmon
               '#ef8a62',    # salmon
               '#b2182b')    # deep red

# Note: length(xvals) == length(colorvals).
#       The pairs of identical xvals around -0.01 and 0.01
#       ensure that we get a "plateau" of white in that region.

# 3. Convert xvals into [0..1] scale for scale_fill_gradientn via rescale()
rescaled_vals <- rescale(xvals, from = c(hmin, hmax))
# This ensures that the left-most breakpoint is mapped to 0, and the right-most to 1.

# 4. Plot with geom_tile() and the custom gradient
ggplot(hmap, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = colorvals,
    values = rescaled_vals,
    limits = c(hmin, hmax),
    oob = squish,         # squish values outside the limits
    na.value = "grey",    # color cells grey if h is NA
    name = ""
  ) +
  ggtitle(paste0("Leading EV differences (", div_2, "–", div_1, " spcs)"))+
  xlab("Nonlinear exponent, alpha") +
  ylab("Interaction similarity, w") +
  theme(
    plot.title = element_text(size = 10, face = "italic"),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.key.height = unit(2, "cm"),
    legend.key.width  = unit(0.5, "cm")
  )

ggsave("wa_draft.pdf", width = 5.5, height = 5)
