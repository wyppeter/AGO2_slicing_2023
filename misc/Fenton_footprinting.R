# Processing and plotting data from Fenton probing
library(tidyverse)

# Defining lists and mappings ----
miRlvls = c(
  "miR-7",
  "let-7a",
  "miR-196a",
  "miR-7.M4",
  "miR-430a",
  "miR-124",
  "miR-124.M2"
)
spacebuffer = "___________________________"
sampleTimeList = c(
  "notarget_15",
  "notarget_30",
  "seed_30",
  "seedsupp_30",
  "tdmd_30",
  "partial_4",
  "partial_10",
  "partial_30",
  "slicing_1",
  "slicing_4",
  "slicing_10",
  "slicing_15",
  "slicing_30",
  spacebuffer
)
sampleTimeList.names = c(
  "No target (15 min)",
  "No target (30 min)",
  "Seed (30 min)",
  "Seed + suppl. (30 min)",
  "Seed + ext. 3' (TDMD) (30 min)",
  "16bp (4 min)",
  "16bp (10 min)",
  "16bp (30 min)",
  "Perfect (1 min)",
  "Perfect (4 min)",
  "Perfect (10 min)",
  "Perfect (15 min)",
  "Perfect (30 min)",
  spacebuffer
)
sampleTimeList.colors = c(
  "black",
  "black",
  "dodgerblue1",
  "slateblue",
  "springgreen3",
  "goldenrod1",
  "darkgoldenrod3",
  "darkgoldenrod4",
  "lightsalmon1",
  "indianred2",
  "red",
  "red3",
  "red4",
  "white"
)
considered.range = c(4, 20)

# Set ggplot2 theme ----
theme0 = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black", linewidth = 0.75),
  axis.ticks = element_line(colour = "black", linewidth = 0.75),
  axis.text = element_text(color = "black", size = 12),
  axis.title = element_text(color = "black", size = 12),
  legend.background = element_blank(),
  legend.key = element_blank(),
  strip.background = element_blank(),
  strip.text.x = element_text(size = 12, colour = "black"),
  strip.text.y = element_text(size = 12, colour = "black", angle = 0)
)

# Read and process data ----
## Raw input and initial processing
df.raw = read.csv("fenton_values_upto21.csv") %>%
  mutate(rxnID = miR) %>%
  filter(!(miR %in% c("miR-124", "miR-124.M2") & sample == "tdmd")) %>%
  
  ### Data QC filters
  filter(!(miR == "miR-430a" & expe == "2")) %>%  # Outlier set; uniform & elevated signal, confirmed by code below
  filter(!(miR == "miR-430a"   & pos >= 15)) %>%  # Interference from trimmed isoforms
  filter(!(miR == "miR-124"    & pos >= 15)) %>%  # Interference from trimmed isoforms
  filter(!(miR == "miR-124.M2" & pos >= 15))      # Interference from trimmed isoforms
  

## Linear-rescaling normalization to the range between quenched and naked guide ----
df = left_join(
    df.raw %>%
      filter(!sample %in% c("naked","quenched")),
    df.raw %>%
      filter(sample %in% c("naked","quenched")) %>%
      select(-time) %>%
      pivot_wider(names_from = "sample", values_from = "fracfenton"),
    by = c("expe","pos","miR","rxnID")  # by each miR, each experiment, and each nt size band
  ) %>%
  rename(signal.raw = fracfenton) %>%
  filter(!is.na(signal.raw)) %>%
  mutate(
    signal.max = if_else(naked > quenched, naked - quenched, NA_real_),
    signal.subtracted = signal.raw - quenched,
    signal = signal.subtracted / signal.max) %>%  # value is now in the relative range
  mutate(sample_time = paste(sample, time, sep = "_")) %>%
  arrange(rxnID, time, sample, -pos) %>%
  filter(pos >= considered.range[1], pos <= considered.range[2]) %>%
  group_by(expe, miR, sample, time, sample_time, pos) %>%
  mutate(sample.enum = row_number())

## Calculate mean and SEM ----
sem.safe = function(x) {
  if_else(length(x) > 1,
          sqrt(var(x)/length(x)),
          0)
}
# Get mean and spread of values over replicates, calculated at a per-sample level
df.summ = df %>%
  group_by(miR, rxnID, time, sample, sample_time, pos) %>%  # by miR, time point, treatment, nt size band
  filter(!is.na(signal)) %>%
  summarize(data.frame(
    signal.mean = mean(signal),
    signal.sem  = sem.safe(signal),
    signal.sd   = sd(signal) # !! for outlier detection
    ),
  .groups = "drop_last") %>%
  mutate(
    signal.max = signal.mean + signal.sem*1.96,
    signal.min = signal.mean - signal.sem*1.96  # raw value distribution
  )

#@#@#@#@#@#@#@#@
# # Special section: checking that miR-430a, seed+supp, exp 2 is an outlier: ----
# # Calculate how many SDs are the values away from the mean, per pos
# avg.sd = sqrt(mean((ungroup(df.summ)$signal.sd)**2))
# 
# df.dev = df %>%
#   left_join(df.summ %>%
#               ungroup() %>%
#               select(miR, time, pos, sample,
#                      signal.mean),
#             by = c("miR", "time", "pos", "sample")) %>%
#   ungroup() %>%
#   mutate(devia = (signal - signal.mean)/avg.sd)
# df.dev.summ = df.dev %>%
#   group_by(expe, sample.enum, time, rxnID, sample) %>%
#   mutate(Q.lo = quantile(devia, probs = c(0.25), na.rm = T),
#          Q.me = quantile(devia, probs = c(0.50), na.rm = T),
#          Q.hi = quantile(devia, probs = c(0.75), na.rm = T)
#          ) %>%
#   ungroup() %>%
#   select(expe, sample.enum, rxnID, sample, time,
#          Q.lo, Q.me, Q.hi) %>%
#   distinct() %>%
#   arrange(-abs(Q.me))
# df.dev %>%
#   ggplot(aes(y = factor(paste(expe, sample, sample.enum)),
#              x = devia,
#              group = factor(paste(expe, sample, sample.enum)),
#              color = sample)) +
#   facet_wrap("rxnID", scales = "free_y") +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_vline(xintercept = c(-2, 2), color = "salmon1", linetype = "dashed") +
#   geom_vline(xintercept = c(-2.5, 2.5), color = "salmon3", linetype = "dashed") +
#   geom_boxplot() +
#   geom_point(shape = 1, size = 1) +
#   labs(y = "Experiment", x = "SEMs away from mean") +
#   theme0
# #### >>>> dre-miR-430a, B, expe 2 is an outlier set. Median > 2 SDs out!
#@#@#@#@#@#@#@#@

# Plots ----
# Reactivity for each target type overlaid

df.ND = data.frame(rxnID = c(
  "miR-430a",
  "miR-124",
  "miR-124.M2"
  )) %>%
  mutate(rxnID = factor(rxnID, levels = miRlvls))

df.summ %>%
  
  ### ***Filter which results to show
  # filter(time == 30) %>%
  # filter(time != 30) %>%
  filter(!(time == 15 & sample == "notarget" & miR == "miR-196a")) %>%
  filter(rxnID %in% c(
    "miR-7",
    # "let-7a",
    "miR-196a",
    "miR-7.M4",
    "miR-430a",
    "miR-124",
    "miR-124.M2"
  )) %>%
  filter(sample %in% c(
    "notarget",
    # "seed","seedsupp",
    # "partial",
    "slicing",
    # "tdmd",
    spacebuffer
  )) %>%
  
  bind_rows(data.frame(
    sample_time = spacebuffer,
    sample = spacebuffer,
    # miR = "miR-124", rxnID = "miR-124"
    miR = "miR-7-L22_v2", rxnID = "miR-7"
    )) %>%  # makes display more consistent between plots
  mutate(rxnID = factor(rxnID, levels = miRlvls)) %>%
  
  ### Plot
  ggplot(aes(x = pos, y = signal.mean, color = sample_time)) +
  facet_wrap("rxnID", ncol = 6, scales = "free_x") +
  
  # ***Highlight pos 9 to 11 for slicing data plots
  annotate("rect", fill = "lightcoral", alpha = 0.15,
           xmin = 8.5, xmax = 11.5, ymin = -Inf, ymax = Inf) +
  
  geom_hline(yintercept = c(0, 1), color = "black",
             linewidth = 0.75, linetype = "dashed") +
  
  geom_errorbar(aes(ymax = signal.max, ymin = signal.min),
                width = 0.4, linewidth = 0.5, alpha = 0.3) +
  geom_line(linewidth = 0.75, alpha = 0.75) +
  
  # ***Indicate indeterminable data range
  geom_text(data = df.ND, inherit.aes = F,
            aes(x = 18, y = 0.5, label = "N.D."),
            size = 4, color = "gray50",
            hjust = 0.5, vjust = 0.5) +
  
  scale_color_manual(
    breaks = sampleTimeList,
    labels = sampleTimeList.names,
    values = sampleTimeList.colors) +
  
  coord_cartesian(ylim = c(-0.05,1.05), xlim = considered.range) +
  scale_x_continuous(breaks = seq(0,25,2)) +
  labs(x = "Guide nucleotide position", y = "*OH reactivity",
       color = "Target bound") +
  theme0

#===================
df.countN = df %>% 
  filter(pos == 10) %>% 
  group_by(rxnID, sample, time) %>%
  count()
