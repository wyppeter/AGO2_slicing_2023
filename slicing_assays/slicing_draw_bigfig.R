# Redraw slicing kinetics plots from saved fit results
# Requires helper file ODElib.R
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
library(lemon)
# ODEs
library(deSolve)

# SELECT MODEL MODE ----
modelmode = "CTS"

# SAVE OR NOT? ----
SAVE_PLOT = T

# DISPLAY SETUP ----
options(scipen = 5)
formatCoefs = function(x, dp = 3){
  as.character(format(round(as.numeric(x), dp), nsmall = dp))
}
formatCoefs.print = function(x){
  paste0("", sub("\\.$", "", sprintf("%#.2g", signif(x, 2))), "")
}
sci10 = function(x){
  parse(text = gsub("\\+", "",
                    gsub("1e", "10^",
                         scales::scientific_format()(x))))
}

# CONSTANTS ----
difflim = log(60/60)  # nM-1 s-1 (10-9 M-1 s-1); diffusion limit; 60 in nM-1 min-1, 1 in nM-1 s-1
inf.conc = 1e6  # nM; a large number to graph limit behavior (no kon effects)
Rconc.given = 0.02  # nM; used for making graphs

# Set ggplot2 theme
theme0.megaplot = theme(
  text = element_text(family = "sans"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  axis.text = element_text(color = "black", size = 12, family = "sans"),
  axis.title = element_text(color = "black", size = 12, family = "sans"),
  legend.background = element_blank(),
  legend.key = element_blank(),
  strip.background = element_blank(),
  strip.text.x = element_text(size = 12, colour = "black", 
                              family = "sans", hjust = 0.5),
  strip.text.y = element_text(size = 12, colour = "black", 
                              family = "sans", hjust = 0.5, angle = 0)
)

# READ DATA ----
# [A] >= [R]*4
df.raw = read.csv(paste0(
  "./",
  "slicing_data_STO.csv"),
  stringsAsFactors = T
)

# READ COEFS ----
ODEfit = read.csv(paste(
    "./",
    "slicing_ktable.csv",
    sep = "")) %>%
  select(miR, 
         kon, kon.lo, kon.hi,
         kcat, kcat.lo, kcat.hi,
         kph2, kph2.lo, kph2.hi,
         Fa, Fa.lo, Fa.hi,
         sampleN) %>%
  mutate(
    kon = kon/60,
    kcat = kcat/60,
    kph2 = kph2/60
  )

# PROCESS AND COMBINE DATA ----
df = df.raw %>%
  
  ### Note suffixes:
  ###
  ### miR-000a#XXX_XXX --> #XXX = RISC notes, _XXX = targ notes
  ###
  
  # Add note to non-PW prep AGO rxns
  mutate(
    miR = as.character(miR) %>%
      if_else(prep == "PW", ., paste0(., "#", prep)) %>%
      factor()
  ) %>%

  ## Remove special data points
  filter(!miR %in% c(
    # Single-conc rxns that cannot be used
    "miR-122#TP",
    "miR-124#TP",
    "miR-7#TP"
  )) %>%

  ## Concentration correction (guide* quant)
  mutate(miR_c = sub("_.+", "", as.character(miR)),  # get AGO prep name for quant lookup
         conc_raw = conc,
         conc = if_else(
           qtype == "guide",
           conc * qLookUp(miR_c),
           conc),  # correct for quant ratio
         approx.quant = (qtype == "guide" & qLookUp.bool(miR_c)), 
         
         ## MISC CLEANUP
         time = as.integer(round(time*60))) %>%  # rounding to the time resolution of ODE, and convert to seconds
  filter(conc != 0, time != 0) %>%  # discard control points
  mutate(assay = factor(as.character(assay)), conc_d = factor(conc)) %>%  # factor-ize some variables
  arrange(miR, conc, assay)

# Note down the guide names used here
guides = as.character(unname(levels(df$miR)))

# Note down quant approximations from guide*
df.q.note = df %>%
  select(miR, approx.quant) %>%
  distinct()

# Simplify df
df.dt = df %>%
  select(miR, conc, Rconc, time, fraccleaved)

#############################################################
# ODE SETUP ----

# Asymptotic floating point protection
eps = 1e-25 # 1 molecule in 10 uL is 1/(Avo*10^-5) M ~ 10^-19 M = 10^-10 nM.
# Event triggered if state variable <= eps
rootfun = function (t, y, pars) y - eps
# Sets state variable to zero
eventfun = function (t, y, pars) {
  y[y <= eps] = 0
  return(y)
  }

# Generate time points for ODE
# First non-zero time point is one second
Tend = 10 * 60 * 60  # 10 hours
Tthres = c(10, 30, 60) * 60  # 10 min, 0.5 hr, then 1 hr
allTime = c(
  seq(0,         Tthres[1]-1, 1),
  seq(Tthres[1], Tthres[2]-1, 5),
  seq(Tthres[2], Tthres[3]-1, 60),
  seq(Tthres[3], Tend,        5*60)
)
timeGen = function(Tmax) {
  if(Tmax>Tend){stop("time exceeds expected upper bound")}
  return(allTime[allTime <= Tmax])
}

############################################################
# ODE FUNCTIONS ----

initConc = function(d){
  # Generate conc setup as a named vector
  c(conc = d$conc[1], Rconc = d$Rconc[1])
}

ODEcalc = function(state.setup, times, parameters, params.static){
  # Given initial concs, time points, and current parameters,
  # provide the entire simulation state trajectory over time points

  # Define A_T and R_T
  A_T = state.setup["conc"]  %>% unname()
  R_T = state.setup["Rconc"] %>% unname()

  # Include static params
  paramsplus = c(parameters, params.static)

  # Generate the initial conditions
  state = ODEstateinit(modelmode, A_T, R_T, paramsplus)

  # Apply ODE
  out.model = ode(
    y = state,
    times = times,
    func = ODElib(modelmode),
    parms = paramsplus,
    method = "lsode",
    mf = 22,
    rtol = 1e-6,
    atol = 1e-10,
    # these two lines provide a floating point protection at asymptotes
    rootfun = rootfun,
    events = list(func = eventfun, root = T)
  ) %>%
    as.data.frame() %>%
    mutate(conc = A_T, Rconc = R_T,
           fraccleaved.ode = P/R_T)
  out.model # return value
}

####################################################################

# SETUP FOR PLOTS ----

## Get fit graphs ----
expandTime = 1/3

ODEcalc.graph = function(d){
  # Given a one-row df of setup parameters, generate an ODE trajectory

  t.plot = timeGen(d$Tmax*(1+expandTime))

  params.static = c(
    kon = log(d$kon),
    kcat = log(d$kcat),
    kph2 = log(d$kph2),
    Fa = d$Fa
    )

  return(
    ODEcalc(
      state.setup = c(conc = d$conc, Rconc = d$Rconc),
      times = t.plot,
      parameters = c(),
      params.static = params.static
    ) %>%
      filter(time > 0)  # remove 0 time point from sim; we are plotting in log
  )
}

# Get Tmax for plot
df.times = df.dt %>%
  ungroup() %>%
  select(miR, time) %>%
  group_by(miR) %>%
  summarize_all(max) %>%
  rename(Tmax = time)

ODEfit.plot = ODEfit %>%
  inner_join(., df.times, by = c("miR"))

# Expe concs
df.concs = df.dt %>%
  ungroup() %>%
  select(miR, conc) %>%
  distinct() %>%
  mutate(Rconc = Rconc.given)
df.graph.setup = ODEfit.plot %>%
  inner_join(df.concs, ., by = c("miR"))
df.graph = df.graph.setup %>%
  group_by(miR, conc, Rconc) %>%
  do(ODEcalc.graph(.))

# Limit line (no kon delay, ~ Inf conc)
n.guides    = length(guides)
df.concs.inf = data.frame(
  miR   = guides,
  conc  = rep(inf.conc,   times = n.guides),
  Rconc = rep(Rconc.given, times = n.guides)
)
df.graph.setup.inf = ODEfit.plot %>%
  inner_join(df.concs.inf, ., by = c("miR"))
df.graph.inf = df.graph.setup.inf %>%
  group_by(miR, conc, Rconc) %>%
  do(ODEcalc.graph(.)) %>%
  mutate(conc = Inf)

###########################################################################

# PLOT ----
FACETS = 7
plt = df.dt %>%
  left_join(namedict.print.pub, by = "miR") %>%
  ggplot(aes(x = time/60, y = fraccleaved,
             color = conc, group = conc)) +

  geom_line(data = df.graph %>%
              left_join(namedict.print.pub, by = "miR"),
            aes(y = fraccleaved.ode), linewidth = 0.5, alpha = 0.75) +
  geom_line(data = df.graph.inf %>%
              left_join(namedict.print.pub, by = "miR"),
            aes(y = fraccleaved.ode), linewidth = 0.5, alpha = 0.75, color = "black") +
  geom_point(size = 0.85, stroke = 1, shape = 1, alpha = 0.75) +

  geom_text(data = ODEfit.plot %>%
              left_join(namedict.print.pub, by = "miR") %>%
              left_join(df.q.note, by = "miR"),
            inherit.aes = FALSE, size = 4.25, color = "black",
            y = 0.95, vjust = 1,
            x = -Inf, hjust = -0.15,

            parse = TRUE,
            aes(
              label = paste0("N == ",
                             sampleN
              ))) +

  facet_wrap("rxnID.print", ncol = FACETS, scales = "free") +
  scale_color_viridis(direction = -1, option= "D", discrete = F, na.value = "gray80",
                      end = 0.97, limits = c(0.05, 50),
                      trans = "log10",
                      name = "[RISC] (nM)") +

  coord_capped_cart(
    ylim = c(0, 1.2),
    left = "top", bottom = "none") +

  scale_x_continuous(expand = expansion(mult = c(0, expandTime), add = 0),
                     trans = "log10", labels = sci10) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = seq(0, 1, 0.2)) +

  xlab("Time (min)") +
  ylab("Fraction sliced") +
  theme0.megaplot
# plt

##############################################################
# SAVING ----

suff = "COMBINED"

## Save plot ----
if(SAVE_PLOT){
  ggsave(
    filename = paste0(
      "./",
      "slicingplot_",
      suff,
      "_DRAW",

      ".pdf"),
    plot = plt,

    width = 14, height = 19,

    units = "in"
  )
}

#@#@#
print("Done!")
