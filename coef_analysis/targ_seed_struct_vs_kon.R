library(tidyverse)

## Setup and data input ----

# Set ggplot2 theme
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

difflim = 60

# Initial data wrangling
df.k = read.csv(paste0(
  "./",
  "slicing_ktable",
  ".csv")) %>%
  left_join(namedict, by = c("miR", "rxnID"))

df.dG = read.csv(paste(
  "./",
  "targ-seqs+seedStruct.csv",
  sep = ""))

df = df.k %>%
  select(rxnID, kon, kon.lo, kon.hi, kon.pP) %>%
  left_join(df.dG, by = "rxnID") %>%
  filter(!is.na(kon))

lmOut = df %>%
  lm(data = .,
     formula = log10(kon) ~ seed.dG)
coefs = coef(lmOut)
lm.plot = data.frame(
  seed.dG = seq(-7,0,0.01)
) %>%
  mutate(kon = 10^(seed.dG * coefs[2] + coefs[1]))

taub = cor.test(
  x = df$seed.dG,
  y = log10(df$kon),
  method = "kendall")
taub.pval = taub$p.value

df %>%
  ggplot(aes(x = seed.dG, y = kon)) +
  geom_hline(yintercept = difflim, linetype = "dashed",
             color = "black", linewidth = 0.5) +
  geom_line(data = lm.plot, color = "gray20", linewidth = 0.5) +
  geom_errorbar(aes(ymin = kon.lo, ymax = kon.hi),
                linewidth = 0.5, width = 0.1,
                color = "firebrick", alpha = 0.25) +
  geom_point(shape = 16, color = "firebrick", alpha = 0.5, size = 2) +
  annotate("text",
    label = paste0(
      "p_tau-b = ",
      signif(taub.pval, 2), "\n",
      "N = ",
      nrow(df)
    ),
    size = 3,
    x = -7, y = 30, hjust = 0,
    color = "gray20") +
  scale_y_continuous(trans = "log10") +
  labs(y = "kon (nM-1 min-1)",
       x = "Predicted dG of target RNA structure at the seed site (kcal mol-1)") +
  theme0
