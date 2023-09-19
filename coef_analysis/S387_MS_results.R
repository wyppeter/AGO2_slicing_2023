library(tidyverse)
library(ggforce)
library(ggnewscale)

# Set ggplot2 theme ----
theme0 = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  axis.text = element_text(color = "black", size = 12),
  axis.title = element_text(color = "black", size = 12),
  legend.background = element_blank(),
  legend.key = element_blank(),
  strip.background = element_blank(),
  strip.text.x = element_text(size = 12, colour = "black"),
  strip.text.y = element_text(size = 12, colour = "black", angle = 0)
)

# Initial data wrangling ----
df.k = read.csv(paste0(
  "./",
  "slicing_ktable.csv")) %>%
  left_join(namedict, by = c("miR", "rxnID"))

df.TIC = read.csv(paste0(
  "./",
  "peptide_TIC_phosphorylation.csv"))

df = df.TIC %>%
  left_join(df.k %>%
              select(rxnID,
                     kcat, kcat.lo, kcat.hi,
                     Fa, Fa.lo, Fa.hi), by = "rxnID")

# Get regression stats ----
lm(
  formula = log10(kcat) ~ fracphos,
  data = df
) %>%
  summary()
lm(
  formula = Fa ~ fracphos,
  data = df
) %>%
  summary()

# Plot ----
df %>%
  ggplot(aes(
    x = fracphos,
    
    # Toggle:
    y = kcat
    # y = Fa
    
  )) +
  geom_errorbar(
    
    # Toggle:
    aes(ymin = kcat.lo, ymax = kcat.hi),
    # aes(ymin = Fa.lo, ymax = Fa.hi),
    
    width = 0,
    color = "springgreen4", alpha = 0.5) +
  geom_point(shape = 16,
             color = "springgreen4") +
  scale_y_continuous(trans = "log2") +
  scale_x_continuous(breaks = c(0,0.02,0.04)) +
  theme0
