# Plot residuals of fits to data
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
# ODEs
library(deSolve)
# Levene's test
library(car)

# DISPLAY SETUP ----
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
df.raw = read.csv(paste0(
  "./",
  "slicing_data_STO.csv"),
  stringsAsFactors = T
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
    # Slow 16bp rxns
    "miR-124-X2_16bp",
    "miR-196a_16bp",
    "miR-196a-X1_16bp",
    
    # Slow mm17 rxns
    "miR-196a_mm17GA",
    "miR-196a-X1_mm17AA",
    
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

ODEcalc = function(state.setup, times, parameters, params.static, modelmode){
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

ODEcompare.mod = function(df.this, modelmode){
  # Given the data of a miR at a particular conc set, calculate the ODE outcome,
  # and give a vector of residuals

  # Generate time sequence
  Tmax = max(df.this$time)
  t = timeGen(Tmax)

  # Get params
  kph2.check = log(unname(unique(df.this$kph2)))
  params = c(
    kon = log(unname(unique(df.this$kon))),
    kcat = log(unname(unique(df.this$kcat))),
    kph2 = if_else(is.na(kph2.check), -Inf, kph2.check),
    Fa = unname(unique(df.this$Fa))
  )

  # Calculate ODE simulation and compare
  overlaid = ODEcalc(
    state.setup = initConc(df.this),
    times = t,
    parameters = params,
    params.static = c(),
    modelmode = modelmode
  ) %>%
    left_join(
      df.this, .,
      by = c("conc","Rconc","time"),
      keep = F
    ) %>%
    ungroup() %>%
    mutate(residuals = fraccleaved - fraccleaved.ode)

  if(any(is.na(overlaid$fraccleaved.ode))){
    print(overlaid)
    stop("Crit. error at ODEcompare!")
  }

  overlaid
}

####################################################################

expandTime = 1/3

getResi = function(modelmode){
  # Get residuals for fit of given modelmode

  # READ COEFS ----
  ODEfit = read.csv(paste(
    "./",
    "slicing_ktable_", modelmode, ".csv",
    sep = "")) %>%
    select(miR, kon, kcat, kph2, Fa, sampleN) %>%
    mutate(
      kon = kon/60,
      kcat = kcat/60,
      kph2 = kph2/60
    )

  # Get Tmax for plot
  df.times = df.dt %>%
    ungroup() %>%
    select(miR, time) %>%
    group_by(miR) %>%
    summarize_all(max) %>%
    rename(Tmax = time)

  df.params = ODEfit %>%
    inner_join(., df.times, by = c("miR"))

  # GET RESIDUALS ----
  df.res = df.params %>%
    left_join(df.dt %>%
                ungroup(),
              .,
              by = "miR") %>%
    group_by(miR, conc, Rconc) %>%
    do(ODEcompare.mod(., modelmode)) %>%
    filter(time / (log(2)/kcat) > 3) %>%  # take late time points only
    ungroup() %>%
    select(residuals) %>%
    mutate(model = modelmode)

  df.res
}

# Run here ----
res.CTS = getResi("CTS")
res.MAX = getResi("MAX")

df.res = bind_rows(
  res.CTS,
  res.MAX
  ) %>%
  mutate(model = factor(model))

# Levene's test ----
leveneOut = leveneTest(
  residuals~model,
  center = function(x){0},
  data = df.res
  )

# Plot ----
df.res %>%
  ggplot(aes(x = residuals, fill = model)) +
  geom_histogram(stat = "bin", position = "identity",
                 bins = 50,
                 alpha = 0.5) +
  geom_vline(xintercept = 0, color = "gray20", linetype = "dashed") +
  scale_fill_manual(values = c(
    "turquoise3",
    "goldenrod3")
    ) +
  coord_cartesian(xlim = c(-0.2, 0.2)) +
  labs(y = "Data points count", x = "Fit residual") +
  theme0

print(paste0("p = ", leveneOut[1,3] %>% signif(4)))
print("Sample sizes (N) = ")
nrow(res.CTS)
nrow(res.MAX)
