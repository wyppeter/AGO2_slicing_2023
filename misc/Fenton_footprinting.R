# Processing and plotting data from Fenton probing

library(tidyverse)

# Defining lists and mappings ----
miRlvls = c(
  "miR-7",
  "let-7a",
  "miR-196a",
  "miR-430a"
)
miRcols = c(
  "deeppink3",
  "royalblue3",
  "green4",
  "purple4"
)
spacebuffer = "_______________"
samplelanesList = c(
  "N","A","B",
  "E","C","D","C_k",
  spacebuffer
)
samplelanesList.names = c(
  "No target", "Seed", "Seed+Suppl.",
  "16-bp", "Perfect", "TDMD", "Perfect (1 min)",
  spacebuffer
)
samplelanesList.cols = c(
  "gray50",
  "dodgerblue2",
  "royalblue4",
  "goldenrod1",
  "red4",
  "springgreen4",
  "tomato1",
  "white"
)

kineticexpes = c(8, 9)
considered.range = c(4,20)

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

# Read and process data ----
## Raw input and initial processing
df.raw = read.csv(paste0(
    "./",
    "fenton-data.csv"
  )) %>%
  mutate(
    miR = if_else(miR == "miR-7-L22", "miR-7-L22_v2", miR),
    rxnID = idLookUp(miR),
    expetype = if_else(expe %in% kineticexpes, "ki", "eq")) %>%
  filter(rxnID %in% miRlvls) %>%
  select(-signalzeroed) %>%

  ### Data QC filters
  # filter(!(miR == "miR-430a" & band >= 15)) %>%  # Interference from trimmed isoforms
  filter(!(miR == "miR-430a" & expe == 2 & lane == "B")) %>%  # Outlier set; B sample had uniform & elevated signal, confirmed by code below. Toggle when checking outlier code.
  
  mutate()  # dummy

## Normalization to the range between Q (quenched) and G (naked guide) ----
df = left_join(
    df.raw %>%
      filter(!lane %in% c("G","Q")),
    df.raw %>%
      filter(lane %in% c("G","Q")) %>%
      select(-lane.expe) %>%
      pivot_wider(names_from = "lane", values_from = "signalnorm"),
    by = c("expe","expetype","band","miR","rxnID")  # by each miR, each experiment, and each nt size band
  ) %>%
  rename(sample = lane, signal.raw = signalnorm) %>%
  filter(!is.na(signal.raw)) %>%
  mutate(
    max.inter = if_else(G > Q, G - Q, NA_real_),
    signal.Qsub = signal.raw - Q,
    signal = signal.Qsub / max.inter) %>%  # value is now in the relative range starting from Q to G
  rename(pos = band) %>%
  arrange(miR, expetype, sample, lane.expe, pos) %>%
  filter(pos >= considered.range[1], pos <= considered.range[2])

## Calculate mean and SEM ----
sem.safe = function(x) {
  if_else(length(x) > 1,
          sqrt(var(x)/length(x)),
          0)
}
df.summ = df %>%
  group_by(miR, rxnID, expetype, sample, pos) %>%  # by miR, time point, treatment, nt size band
  filter(!is.na(signal)) %>%
  summarize(data.frame(
    signal.mean = mean(signal),
    signal.sem = sem.safe(signal),
    signal.sd = sd(signal) # !! for outlier detection
    ),
  .groups = "drop_last") %>%
  mutate(
    signal.max = signal.mean + signal.sem*1.96,
    signal.min = signal.mean - signal.sem*1.96  # raw value distribution
  )

## Calculate delta reactivity from no-target samples ----
# i.e. increased reactivity upon binding to substrate
df.delta = df %>%
  left_join(df.summ %>%  # using mean values!
              ungroup() %>%
              filter(sample == "N", expetype == "eq") %>%
              transmute(miR = miR, pos = pos, signal.N = signal.mean) %>%
              distinct(),
            by = c("miR","pos") # by miR, time point, nt size band; all treatments relative to the N treatment
  ) %>%
  mutate(
    signal.d = signal - signal.N
  )
df.delta.summ = df.delta %>%
  group_by(miR, rxnID, expetype, sample, pos) %>%  # by miR, time point, treatment, nt size band
  filter(!is.na(signal.d)) %>%
  summarize(data.frame(
    signal.d.mean = mean(signal.d),
    signal.d.sem  = sem.safe(signal.d)
  ),
  .groups = "drop_last") %>%
  mutate(
    signal.d.max = signal.d.mean + signal.d.sem*1.96,
    signal.d.min = signal.d.mean - signal.d.sem*1.96  # delta value distribution
  )

## Incorporate kslice values for relationship, calc fraction occupancy of duplex (slicing) state ----
df.k = read.csv(
  "./slicing_ktable.csv"
)

# Quantify the bridge reactivity of the C state as a proxy of proportion in duplex
duplex.peak.range = c(9, 11)  # region used as proxy for occupancy
df.fracduplex.raw = df.delta %>%
  # Bind in kslices
  left_join(df.k %>%
              select(rxnID, kcat, kcat.lo, kcat.hi),
            by = c("rxnID")) %>%
  filter(pos %in% seq(
    duplex.peak.range[1],
    duplex.peak.range[2]
    ),  # focus on the bridge region that gets exposed upon duplexing
         sample == "C") %>% 
  ungroup() %>%
  arrange(
    miR, rxnID,  # these separate different guides
    pos,  # nt position
    expetype,  # this separates 1 min and 30 min time points
    expe, lane.expe  # these separate by lane (samples)
    )
df.fracduplex = df.fracduplex.raw %>%
  left_join(
    # Get the ceiling of steady state, per guide, per pos
    df.fracduplex.raw %>%
      filter(expetype == "eq") %>%
      group_by(miR, pos) %>%
      transmute(
        miR = miR,
        bridge.ss = mean(signal.d)/(1 - exp(- kcat * 30))  # math approx of where the max really is based on 30 min time point
      ) %>%
      distinct(),
    by = c("miR", "pos")
  ) %>%
  # Calculate fraction occupancy
  mutate(time = if_else(expetype == "ki", 1, 30),
         bridge.norm = signal.d / bridge.ss) %>%  # relative proportion of duplex state
  select(miR, rxnID, expetype, time, pos,
         expe, lane.expe,
         kcat, kcat.lo, kcat.hi,
         bridge.ss, bridge.norm) %>%
  # Now we average out the positions per rep. Diff positions are not reps!
  group_by(miR, rxnID, expetype, time,
           expe, lane.expe) %>%
  select(-pos) %>%
  summarize_all(mean) %>%
  group_by(miR, rxnID, expetype)

df.fracduplex.summ = df.fracduplex %>%
  # Make summarizing stats averaging across reps
  group_by(miR, rxnID, kcat, kcat.lo, kcat.hi, expetype, time) %>%
  summarize(data.frame(
    bridge.norm.mean = mean(bridge.norm),
    bridge.norm.sem  = sem.safe(bridge.norm)
  ),
  .groups = "drop_last") %>%
  mutate(
    bridge.norm.max = bridge.norm.mean + bridge.norm.sem*1.96,
    bridge.norm.min = bridge.norm.mean - bridge.norm.sem*1.96,  # duplex state distribution
    fraccleaved     = 1 - exp(- kcat * 1),
    fraccleaved.max = 1 - exp(- kcat.hi * 1),
    fraccleaved.min = 1 - exp(- kcat.lo * 1)  # fraction sliced calculation
  )

# #@#@#@#@#@#@#@#@
# ## Special section: checking that miR-430a, seed+supp, exp 2 is an outlier: ----
# # Calculate how many avg SDs are the values away from the mean, per pos
# avg.sd = sqrt(mean((ungroup(df.summ)$signal.sd)**2))
# 
# df.dev = df %>%
#   left_join(df.summ %>%
#               ungroup() %>%
#               select(miR, expetype, pos, sample,
#                      signal.mean),
#             by = c("miR", "expetype", "pos", "sample")) %>%
#   ungroup() %>%
#   filter(sample != "Co") %>%
#   mutate(devia = (signal - signal.mean)/avg.sd) %>%
#   mutate(sample.wtype = if_else(expetype == "ki", paste(sample, "k", sep = "_"), sample))
# df.dev.summ = df.dev %>%
#   group_by(expe, lane.expe, expetype, rxnID, sample) %>%
#   mutate(Q.lo = quantile(devia, probs = c(0.25), na.rm = T),
#          Q.me = quantile(devia, probs = c(0.50), na.rm = T),
#          Q.hi = quantile(devia, probs = c(0.75), na.rm = T)
#          ) %>%
#   ungroup() %>%
#   select(expe, lane.expe, rxnID, sample, expetype,
#          Q.lo, Q.me, Q.hi) %>%
#   distinct() %>%
#   arrange(-abs(Q.me)) %>%
#   mutate(sample.wtype = if_else(expetype == "ki", paste(sample, "k", sep = "_"), sample))
# df.dev %>%
#   ggplot(aes(y = factor(interaction(expe, lane.expe)),
#              x = devia,
#              group = factor(interaction(expe, lane.expe)),
#              color = sample.wtype)) +
#   facet_wrap("rxnID", scales = "free_y") +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_vline(xintercept = c(-2, 2), color = "salmon1", linetype = "dashed") +
#   geom_vline(xintercept = c(-2.5, 2.5), color = "salmon3", linetype = "dashed") +
#   geom_boxplot() +
#   geom_point(shape = 1, size = 1) +
#   labs(y = "Experiment", x = "Avg SDs away from mean") +
#   theme0
# #### >>>> dre-miR-430a, B, expe 2 is an outlier set. Median > 2 SDs out.
# #@#@#@#@#@#@#@#@

# All plots ----
# Reactivity for each target type overlaid
df.reactivity = df.summ %>%
  mutate(sample.wtype = if_else(expetype == "ki", paste(sample, "k", sep = "_"), sample)) %>%
  bind_rows(data.frame(
    sample.wtype = spacebuffer,
    rxnID = "miR-7")) %>%
  mutate(rxnID = factor(rxnID, levels = miRlvls))
df.430atrim = data.frame(
  rxnID = "miR-430a"
) %>%
  mutate(rxnID = factor(rxnID, levels = miRlvls))
f1 = df.reactivity %>%
  # Filter which target results to show
  filter(sample.wtype %in% c(
    "N",
    # "A","B",
    # "E",
    "C",
    # "D",
    "C_k",
    spacebuffer
  )) %>%
  # Plot
  ggplot(aes(x = pos, y = signal.mean, color = sample.wtype)) +
  facet_wrap("rxnID", ncol = 4, scales = "free_x") +

  # Highlight pos 9 to 11 for slicing ones
  annotate("rect", fill = "lightcoral", alpha = 0.15,
           xmin = 8.5, xmax = 11.5, ymin = -Inf, ymax = Inf) +

  geom_errorbar(aes(ymax = signal.max, ymin = signal.min),
                width = 0.25, linewidth = 0.5, alpha = 0.33) +
  geom_line(linewidth = 0.75, alpha = 0.75) +

  # Block out interference from trimmed isoforms of miR-430a
  geom_rect(data = df.430atrim,
            inherit.aes = FALSE,
            fill = "white", alpha = 1,
            xmin = 14.5, xmax = 22.5,
            ymin = -Inf, ymax = Inf) +

  geom_hline(yintercept = c(0,1), color = "black",
             linewidth = 0.5, linetype = "dashed") +
  
  scale_color_manual(
    breaks = samplelanesList,
    labels = samplelanesList.names,
    values = samplelanesList.cols) +
  coord_cartesian(ylim = c(-0.05,1.05), xlim = considered.range) +
  scale_x_continuous(breaks = seq(0,25,2)) +
  labs(x = "Guide nucleotide position", y = "OH reactivity") +
  theme0
f1

# Mech graph ----
f2 = df.fracduplex.summ %>%
  filter(expetype == "ki") %>%
  ggplot(aes(y = bridge.norm.mean,
             x = fraccleaved,
             color = rxnID)) +

  geom_hline(yintercept = c(0,1),
             color = "gray70", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = c(0,1),
             color = "gray70", linewidth = 0.5, linetype = "dashed") +
  geom_abline(
    slope = 1,
    intercept = 0,
    linewidth = 0.5,
    color = "gray70"
  ) +

  geom_point(data = df.fracduplex %>%
               filter(expetype == "ki"),
             aes(y = bridge.norm,
                 x = 1 - exp(- kcat * 1)),
             shape = 16, alpha = 0.6, size = 3) +

  geom_point(shape = 16, alpha = 1, size = 1.5) +

  geom_errorbar(aes(ymin = bridge.norm.min, ymax = bridge.norm.max),
                width = 0.05,
                linewidth = 0.5,
                alpha = 0.8,
                # color = "black"
                ) +
  geom_errorbarh(aes(xmin = fraccleaved.min, xmax = fraccleaved.max),
                 height = 0.05,
                 linewidth = 0.5,
                 alpha = 0.8,
                 # color = "black"
                 ) +

  scale_color_manual(
    breaks = miRlvls,
    values = miRcols) +
  scale_y_continuous(breaks = seq(-10,10,0.5)) +

  coord_fixed(ylim = c(-0.25, 1.25), xlim = c(-0.25, 1.25)) +
  labs(x = "Expected fraction sliced based on kslice, at 1 min",
       y = "Fraction Fenton reactivity at g9-11, at 1 min") +
  theme0
f2

