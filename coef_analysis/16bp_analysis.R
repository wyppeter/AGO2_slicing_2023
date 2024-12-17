library(tidyverse)
library(ggrepel)
library(minpack.lm)

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

# Read list
df.kcat = read.csv(paste0(
  "./",
  "slicing_ktable.csv")) %>%
  select(miR, rxnID, kcat, kcat.lo, kcat.hi)

# List of miRs
IDs_measured = c(
  "lsy-6",
  "let-7a",
  "miR-122",
  "miR-124.M1",
  "miR-124.M2",
  "miR-124",
  "miR-196a.M1",
  "miR-196a",
  "miR-451a",
  "miR-7.24mer",
  "miR-7.M2",
  "miR-7"
)

df.kcat.full = df.kcat %>%
  filter(rxnID %in% IDs_measured)
df.kcat.16bp = df.kcat %>%
  filter(rxnID %in% paste0(IDs_measured, "_16bp")) %>%
  transmute(rxnID = sub("_16bp", "", rxnID),
            kcat.16bp = kcat,
            kcat.lo.16bp = kcat.lo,
            kcat.hi.16bp = kcat.hi) %>%
  full_join(df.kcat.full, by = "rxnID") %>%
  mutate(rel.kcat = kcat.16bp/kcat,

         log.CI = log(kcat.hi/kcat.lo),
         log.CI.16bp = log(kcat.hi.16bp/kcat.lo.16bp),
         log.CI.rel = sqrt(log.CI^2 + log.CI.16bp^2),

         rel.kcat.lo = rel.kcat / exp(log.CI.rel),
         rel.kcat.hi = rel.kcat * exp(log.CI.rel)
  )

df.dG = read.delim(paste0(
  "./",
  "16bp_analysis-dGvals.tsv"),
  stringsAsFactors = T
) %>%
  arrange(miR) %>%
  mutate(rxnID = gsub("#.+|_.+", "", miR)) %>%
  filter(rxnID %in% IDs_measured)

df = df.kcat.16bp %>%
  left_join(df.dG %>% select(-miR), by = c("rxnID"))

# Thermodynamic fitting
RT = 1.9872042586 * 310.15 / 1000

nlsFit = df %>%
  ungroup() %>%
  nlsLM(
    formula = log(rel.kcat) ~ -log((exp((cent_dG+dGthres)/RT) + 1)),
    data = .,
    start = c(dGthres = 10),
    algorithm = "LM"
  )
summary(nlsFit)
dGthres.fit = unname(coef(nlsFit)["dGthres"])
CI   = summary(nlsFit)$coefficients[[2]] * 1.96
pval = summary(nlsFit)$coefficients[[4]]

plot.vals = data.frame(
  dG = seq(-100,20,0.1)
) %>%
  mutate(rel.kcat = 1/(1/exp(-(dG+dGthres.fit)/RT) + 1))

df %>%
  ggplot(aes(
    x = cent_dG,
    y = rel.kcat)) +

  # 1:1 ref line
  geom_hline(yintercept = 1, color = "black", linetype = "dashed") +

  # Fit curve
  geom_line(data = plot.vals,
            aes(x = dG),
            color = "black", alpha = 1.0, linewidth = 0.5) +
  annotate(geom = "text", x = -Inf, y = 0.001,
           hjust = 0, vjust = 0, color = "black",
           size = 4.25, label = paste0(
             "   dGthres = \n   ",
             signif(dGthres.fit, 2),
             " (", signif(dGthres.fit-CI, 2), " - ", signif(dGthres.fit+CI, 2), ")",
             " kcal/mol",
             "\n   p = ",
             signif(pval, 2)
           )) +

  # Data points
  geom_errorbar(aes(ymin = rel.kcat.lo, ymax = rel.kcat.hi),
                linewidth = 0.5, width = 0.15,
                color = "blue4", alpha = 0.5) +
  geom_point(shape = 18,
             size = 3,
             alpha = 1.0, color = "blue4") +

  # Labels
  geom_text_repel(aes(label = rxnID),
                  size = 4.25, alpha = 0.5,
                  min.segment.length = 0,
                  box.padding = 0.5,
                  show.legend = F) +

  scale_y_continuous(trans = "log10", breaks = 10^seq(-10,10,1)) +
  scale_x_continuous(breaks = seq(-20,20,2)) +
  coord_cartesian(
    ylim = c(0.001,3),
    xlim = c(-10.5,-1)
  ) +
  labs(y = "kcat_16bp/kcat_perfect",
       x = "Predicted dG of base-pairing at positions 9-12 (kcal mol-1)") +
  theme0

#############

# SAVE
write.csv(df %>%
            transmute(rxnID = rxnID,
                      Rkcat.16bp.perf = rel.kcat,
                      Rkcat.16bp.perf.lo = rel.kcat.lo,
                      Rkcat.16bp.perf.hi = rel.kcat.hi
                      ),
          paste0(
            "./",
            "16bp_Rkcat.csv"
          ),
          row.names = F, quote = F)
