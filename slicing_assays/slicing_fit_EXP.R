# Slicing assay fitting and graphing
# Fitting using exponential approximation for very slow reaction trajectories:
# - log-transformed fit and CI calculation
# - no koff and koffP considerations
# - no kon considerations; assumes kslice << kon
# - no second phase kinetics; assumes known Fa
# Requires helper file lookup-helpers.R
# Bartel Lab
# Peter Y. Wang 2023

# LIBRARIES ----
# tidyverse
library(tidyverse)
# graphing
library(ggnewscale)
library(viridis)
library(ggpubr)
# NLS
library(nlstools)
library(minpack.lm)

# SAVE OR NOT? ----
SAVE_PLOT = T
SAVE_TABLE = T

# DISPLAY SETUP
options(scipen = 5)
formatCoefs = function(x, dp = 3){ as.character(format(round(as.numeric(x), dp), nsmall = dp)) }

# CONSTANTS ----
difflim = log(60/60)  # nM-1 s-1 (10-9 M-1 s-1); diffusion limit; 60 in nM-1 min-1, 1 in nM-1 s-1
inf.conc = 1e6  # nM; a large number to graph limit behavior (no kon effects)
Rconc.given = 0.02  # nM; used for making graphs

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

# READ DATA ----
# [A] >= [R]*4
df.raw = read.csv("./slicing-assay-data_collated-STO.csv", stringsAsFactors = T)

# PROCESS AND COMBINE DATA
# Note suffixes: miR-000a#XXX_XXX --> #XXX = RISC notes, _XXX = targ notes
df = df.raw %>%
  # Map back to legacy headings
  transmute(miR = reaction.symbol,
            assay = assayID,
            conc = RISC.conc,
            Rconc = RNA.conc,
            time = time..s.,
            fraccleaved = fracsliced
  ) %>%
  
  ## Filter to specific data points
  filter(miR %in% c(
    # Slow 16bp rxns
    "miR-124.M2_16bp",
    "miR-196a_16bp",
    "miR-196a.M1_16bp",
    # Slow mm17 rxns
    "miR-196a_mm17GA",
    "miR-196a.M1_mm17AA"
  )) %>%
  mutate(miR = factor(miR)) %>%

  ## MISC CLEANUP
  mutate(assay = factor(as.character(assay)), conc_d = factor(conc)) %>%  # factor-ize some variables
  arrange(miR, conc, assay)

# Note down the guide names used here
guides = as.character(unname(levels(df$miR)))

##############################
# Get ODE data for Fa ----

suff = "CTS"
fit.coefs = read.csv(paste0("./",
                            "slicing_ktable_",
                            suff,
                            ".csv"
                            ), stringsAsFactors = T) %>%
  select(!ends_with(c("\\.lo","\\.hi","\\.p")), -sampleN)

fit.coefs.Fa = fit.coefs %>%
  transmute(miR.given = as.character(miR), Fa.given = Fa)

df.dt = df %>%
  mutate(miR_parent = sub("_16bp|_mm17[AUGC]{2}", "", as.character(miR))) %>%
  left_join(., fit.coefs.Fa, by = c("miR_parent" = "miR.given")) %>%
  select(miR, Fa.given, conc, Rconc, time, fraccleaved)  # Simplify df for fitting

Fa.get = function(y) {
  y$Fa.given %>% unique() %>% unname()
}

# Get number of data points
df.dt.N = df.dt %>%
  group_by(miR) %>%
  mutate(sampleN = n()) %>%
  select(miR, sampleN) %>%
  ungroup() %>%
  distinct()

# NLS fit ----
nlsout = df.dt %>%
  group_by(miR) %>%
  do(
    data.frame(model = broom::tidy(
      nlsLM(  # Levenberg-Marquardt
        data = .,
        fraccleaved ~ Fa.get(.) * (1-exp(-time * exp(kcat))),
        algorithm = "LM",
        start = c(kcat = log(0.01/60)),
        control = nls.lm.control(maxiter = 100)
      )
    ))
  ) %>%
  transmute(
    miR = miR,
    miR_parent = sub("_16bp|_mm17[AUGC]{2}", "", as.character(miR)),
    kcat = exp(model.estimate),
    kcat.lo = exp(model.estimate-model.std.error*1.96),
    kcat.hi = exp(model.estimate+model.std.error*1.96),
    kcat.pP = model.p.value
  ) %>%
  left_join(., fit.coefs.Fa, by = c("miR_parent" = "miR.given")) %>%
  rename(Fa = Fa.given) %>%
  left_join(df.dt.N, by = "miR")

###############################################################

# Get fit graphs ----
expandTime = 1/3
timeInt = 5  # seconds
timeEnd = 4 * 60 * 60  # hrs
t = seq(0, timeEnd, timeInt)
calcFC = function(kcat, Fa, time){
  Fa * (1-exp(-time * kcat))
}

n.guides = length(guides)

nlsFit.graph = data.frame(
  time =        rep(t, times = n.guides),
  miR  = factor(rep(guides, each = timeEnd/timeInt+1), levels = guides)
) %>%
  left_join(nlsout, by = "miR") %>%
  mutate(fraccleaved = calcFC(kcat, Fa, time)) %>%
  filter(time > 0)

##############################

# PLOT ----
FACETS = as.integer(round(sqrt(n.guides*3/2)))
plt = df.dt %>%
  ggplot(aes(x = time/60, y = fraccleaved,
             color = conc)) +

  geom_line(data = nlsFit.graph,
            linewidth = 0.5, alpha = 0.75, color = "black") +
  geom_point(size = 0.75, stroke = 1, shape = 1, alpha = 0.75) +

  # Fitted parameters print
  # print in min-1
  geom_text(data = nlsout,
            inherit.aes = FALSE,
            size = 3,
            y = Inf,
            x = -Inf,
            vjust = 1,
            hjust = -0.1,
            color = "gray25",
            aes(
              label = paste(miR, ":\n",
                            "k_slice (min-1) = \n",
                            formatCoefs(kcat*60*1000),
                            " (",
                            formatCoefs(kcat.lo*60*1000),
                            "-",
                            formatCoefs(kcat.hi*60*1000),
                            ") * 10^-3",
                            ",\nF_a = ",
                            formatCoefs(Fa),
                            ",\nN = ",
                            sampleN,
                            sep = ""))) +

  facet_wrap("miR", ncol = FACETS, scales = "free") +
  scale_color_viridis(direction = -1, option= "D", discrete = F, na.value = "gray80",
                      end = 0.97, limits = c(0.05, 50),
                      trans = "log10",
                      name = "[RISC] (nM)") +

  # scale_x_continuous(expand = expansion(mult = c(0, expandTime), add = 0)) +
  scale_x_continuous(expand = expansion(mult = c(0, expandTime), add = 0), trans = "log10") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) +

  xlab("Time (min)") +
  ylab("Fraction sliced") +
  theme0
plt

##############

## Save plot ----
if(SAVE_PLOT){
  ggsave(
    filename = paste0(
      "./",
      "slicingplot_EXP",
      ".pdf"),
    plot = plt,

    width = 8, height = 5,

    units = "in"
  )
}

# Save output values ----
if(SAVE_TABLE){
  nlsout.write = nlsout %>%
    ungroup() %>%
    transmute(
      rxnID     = miR,
      miR       = miR,

      kon       = exp(difflim)*60,
      kon.lo    = NA_real_,
      kon.hi    = NA_real_,
      kon.pP    = NA_real_,

      kcat      = kcat*60,
      kcat.lo   = kcat.lo*60,
      kcat.hi   = kcat.hi*60,
      kcat.pP   = -log10(kcat.pP),

      kph2      = 0,
      kph2.lo   = NA_real_,
      kph2.hi   = NA_real_,
      kph2.pP   = NA_real_,

      Fa        = Fa,
      Fa.lo     = NA_real_,
      Fa.hi     = NA_real_,
      Fa.pP     = NA_real_,

      sampleN = sampleN
    )
  write.csv(nlsout.write, paste0(
    "./",
    "slicing_ktable_EXP.csv"),
    quote = F, row.names = F
  )
}

#@#@#
print("Done!")
