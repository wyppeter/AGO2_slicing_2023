# Slicing assay fitting and graphing; multiple-turnover to fit koffP
# Fitting using ODEs:
# - log-transformed fit and CI calculation
# - no koff considerations
# - uses known constant values from STO fitting
# - universal adaptive architecture to consider different ODE models
# - includes a separate Fmax for MTO, as inactive RNA
# Requires helper file ODElib-MTO.R
# Requires helper file lookUp-helpers.R
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
library(FME)

# SELECT MODEL MODE ----
modelmode = "CTS"

# modelmode = "MAX"
# modelmode = "CTR"
# modelmode = "UNL"
# modelmode = "CFI"
# modelmode = "CFR"
# modelmode = "CFS"
# modelmode = "BNS"
# modelmode = "BNR"

# SAVE OR NOT? ----
SAVE_TABLE = T

# DISPLAY SETUP ----
options(scipen = 5)
formatCoefs = function(x, dp = 3){ as.character(format(round(as.numeric(x), dp), nsmall = dp)) }

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
# Who needs fitting for Fmax?
fit.Fmax.and.QQ = c(
  "lsy-6",
  "lsy-6#S387A",
  "lsy-6#S387D",

  # "miR-7-X3",
  # "miR-7-X3#S387A",
  # "miR-7-X3#S387D"

  "miR-196a-X3"
  # "miR-7-X5",
  # "miR-M7"
)
# All fit for QQ

df.raw = read.csv(paste0(
  "./",
  "slicing_data_MTO.csv"
  ),
  stringsAsFactors = T
)

# PROCESS AND COMBINE DATA ----
df.mto = df.raw %>%
  mutate(expe_set = paste0("expe",str_extract(as.character(expe), "^[0-9]+"))) %>%  # experiment set

  ## Concentration correction (guide* quant)
  mutate(miR_c = sub("_.+", "", as.character(miR)),  # get AGO prep name for quant lookup
         conc_raw = conc,
         conc = if_else(
           qtype == "guide",
           conc * qLookUp(miR_c),
           conc),  # correct for quant ratio
         approx.quant = (qtype == "guide" & qLookUp.bool(miR_c))
  ) %>%

  ## MISC CLEANUP
  mutate(conc.c = formatCoefs(conc))

# Note down the guide names used here
mto.guides = as.character(unname(levels(df.mto$miR)))

## Get k values from STO data ----
df.k = read.csv(paste0(
  "./",
  "slicing_ktable_", modelmode, ".csv"),
  stringsAsFactors = T
)

## Add k values ----
df = df.mto %>%
  left_join(df.k %>% select(miR, kon, kcat, kph2, Fa), by = "miR")

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

# Generate time points for ODE (in minutes)
Tend = 10 * 60  # 10 hours
Tthres = c(20, 100)
allTime = c(
  seq(0,         Tthres[1]-1, 0.05),
  seq(Tthres[1], Tthres[2]-1, 0.25),
  seq(Tthres[2], Tend,        1)
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

  # Include static params
  paramsplus = c(parameters, params.static)

  ### Special params ----
  # QQ is quantification deviation of RISC (basically fitting size of burst)
  QQ = unname(paramsplus["QQ"])

  # Fmax is max fraction of RNA that can be sliced (basically fitting plateau specific to MTO)
  Fmax = unname(paramsplus["Fmax"])

  # Define A_T and R_T
  A_T = (state.setup["conc"]  %>% unname()) * QQ
  R_T = (state.setup["Rconc"] %>% unname()) * Fmax

  # Generate the initial conditions
  state = ODEstateinit(modelmode, A_T, R_T, params.static)  # params.static contains the Fa value we need

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
    mutate(conc = A_T/QQ, Rconc = R_T/Fmax,
           fraccleaved.ode = (P+AP+IP)/Rconc,
           conc.c = formatCoefs(conc))
  out.model # return value
}

ODEcompare = function(df.this, params, params.static){
  # Given the data of a miR at a particular conc set, calculate the ODE outcome,
  # and give a vector of residuals

  # Generate time sequence
  Tmax = max(df.this$time)
  t = timeGen(Tmax)

  # Calculate ODE simulation and compare
  overlaid = ODEcalc(
    state.setup = initConc(df.this),
    times = t,
    parameters = params,
    params.static = params.static
  ) %>%
    left_join(
      df.this, .,
      by = c("conc.c","Rconc","time"),
      keep = F
    ) %>%
    ungroup() %>%
    mutate(residuals = fraccleaved.ode - fraccleaved)

  if(any(is.na(overlaid$fraccleaved.ode))){
    View(overlaid)
    stop("Crit. error at ODEcompare!")
  }

  overlaid$residuals # return value
}

ODEcost = function(params, df.data, params.static) {
  # Cost function wrapper
  # For every conc-Rconc pair,
  # given parameters, generate simulations, and calculate how much it deviates from the data
  # Returns a compiled vector of residuals (per modFit requirements)

  ### Extract given details ----
  # log transform for ODE
  kon   = log(unique(df.data$kon))
  kcat  = log(unique(df.data$kcat))
  kph2  = log(unique(df.data$kph2))
  Fa    =     unique(df.data$Fa)
  params.static = c(c(kon = kon, kcat = kcat, kph2 = kph2, Fa = Fa), params.static)

  ### Calculate cost ----
  df.cost = df.data %>%
    group_by(conc, Rconc) %>%
    do(res = ODEcompare(
      df.this = .,
      params = params,
      params.static = params.static
    )
    )

  unlist(df.cost$res) # return value
}

###############################

# Coefs and goodness-of-fit ----
# CI is 95%, assuming normal distribution around model
CIcalc = function(param, fitpars, val){
  if( !is.na(val) & any(!is.na(fitpars)) ) {
    tryCatch(unname(unlist(fitpars[param,2])*1.96), error = function(e){NA_real_})
  } else {NA_real_}
}
Pcalc  = function(param, fitpars, val){
  if( !is.na(val) & any(!is.na(fitpars)) ) {
    tryCatch(unname(fitpars[param,4]), error = function(e){NA_real_})
  } else {NA_real_}
}
runaway = 1e3 # run away error margin (when no data to fit)
ci_protect = function(ci, runaway.thres) {
  if_else(is.na(ci), NA_real_, if_else(ci > runaway.thres, Inf, ci))
}
extractCoef = function(df.here, fitO){
  # Extract coefs and their CI and p values

  # Check which miR this is first
  miR.this = unique(as.character(df.here$miR))

  # Get number of data points for statistics (sampleN)
  sampleN = nrow(df.here)

  print(
    paste0(
      "Now fitting: ",
      miR.this,
      " (N = ",
      sampleN,
      ")"
    )
  )

  # print(fitO)

  coefs   = coef(fitO)
  fitpars = tryCatch(
    summary(fitO)$par,
    warning = function(w){
      warnMsg = conditionMessage(w)
      message(warnMsg)
      return(NA)
    }
  )
  # print(fitpars)

  koffP = coefs["koffP"]
  QQ    = coefs["QQ"]
  QQ    = if_else(is.na(QQ), 1, QQ)
  Fmax  = coefs["Fmax"]
  Fmax  = if_else(is.na(Fmax), 1, Fmax)

  koffP.CI = ci_protect(CIcalc("koffP", fitpars, koffP), log(runaway))
  QQ.CI    = ci_protect(CIcalc("QQ",    fitpars, QQ),        runaway)
  Fmax.CI  = ci_protect(CIcalc("Fmax",  fitpars, Fmax),      runaway)

  koffP.p = Pcalc("koffP", fitpars, koffP)
  QQ.p    = Pcalc("QQ",    fitpars, QQ)
  Fmax.p  = Pcalc("Fmax",  fitpars, Fmax)

  print("Fit output:")
  print(suppressWarnings(as.numeric(formatCoefs(c(
    koffP = exp(koffP),
    QQ    = QQ,
    Fmax  = Fmax
  )))))

  # return in unit of min-1, linear
  return(data.frame(
    koffP    = exp(koffP),
    koffP.lo = exp(koffP  - koffP.CI),
    koffP.hi = exp(koffP  + koffP.CI),
    koffP.p  = koffP.p,
    QQ       = QQ,
    QQ.lo    = QQ    - QQ.CI,
    QQ.hi    = QQ    + QQ.CI,
    QQ.p     = QQ.p,
    Fmax     = Fmax,
    Fmax.lo  = Fmax  - Fmax.CI,
    Fmax.hi  = Fmax  + Fmax.CI,
    Fmax.p   = Fmax.p,

    sampleN = sampleN,

    row.names = miR.this
  ))

}

###############################################################

# FITTING TIME ----

MAXITER = 200

fitmax.any = any( mto.guides %in% fit.Fmax.and.QQ) > 0
fitqq.any  = any(!mto.guides %in% fit.Fmax.and.QQ) > 0

# Make guesses
koffP.guess.fast = log(10)
koffP.guess.slow = log(1)
QQ.guess   = 1.0
QQ.lim   = c(1/2, 2)
Fmax.guess = 0.95
Fmax.lim = c(0.0, 1.0)
# First, fit for ones that do not need to fit Fmax (Fmax = 1) ----
MTOfit.default = data.frame()
if(fitqq.any){
  MTOfit.default = df %>%
    filter(!miR %in% fit.Fmax.and.QQ) %>%
    group_by(expe_set, miR) %>%  # separate experiments
    do(extractCoef(.,
                   modFit(
                     f = ODEcost,
                     p = c(koffP = koffP.guess.slow, QQ = QQ.guess),
                     df.data = .,
                     params.static = c(Fmax = 1.0),
                     lower = c(-Inf, QQ.lim[1]),
                     upper = c( Inf, QQ.lim[2]),
                     method = "Marq",
                     jac = NULL,
                     control = list(
                       nprint = 0,
                       maxiter = MAXITER
                     ),
                     hessian = T
                   )
    ))
}

# Then, fit for ones that do need to fit Fmax ----
MTOfit.Fmax = data.frame()
if(fitmax.any){
  MTOfit.Fmax = df %>%
    filter(miR %in% fit.Fmax.and.QQ) %>%
    group_by(expe_set, miR) %>%  # separate experiments
    do(extractCoef(.,
                   modFit(
                     f = ODEcost,
                     p = c(koffP = koffP.guess.fast, QQ = QQ.guess, Fmax = Fmax.guess),
                     df.data = .,
                     params.static = c(QQ = 1.0),
                     lower = c(-Inf, QQ.lim[1], Fmax.lim[1]),
                     upper = c( Inf, QQ.lim[2], Fmax.lim[2]),
                     method = "Marq",
                     jac = NULL,
                     control = list(
                       nprint = 0,
                       maxiter = MAXITER
                     ),
                     hessian = T
                   )
    ))
}

# Combine results ----
MTOfit = rbind(
  MTOfit.default,
  MTOfit.Fmax
)

# Add in the fit outputs
df.result = df %>%
  full_join(MTOfit, by = c("expe_set", "miR")) %>%
  mutate(burst.lvl = conc*QQ/Rconc)

####################################################################

# SETUP FOR PLOTS ----

## Get fit graphs ----
expandTime = 1/3

ODEcalc.graph = function(d){
  # Given a one-row df of setup parameters, generate an ODE trajectory

  # Generate time sequence
  Tmax = max(d$time)
  t.plot = timeGen(Tmax*(1+expandTime))

  params.static = c(kon   = log(unique(d$kon)),
                    kcat  = log(unique(d$kcat)),
                    kph2  = log(unique(d$kph2)),
                    Fa    =     unique(d$Fa),
                    koffP = log(unique(d$koffP)),
                    Fmax  =     unique(d$Fmax))

  return(
    ODEcalc(
      state.setup = initConc(d),
      times = t.plot,
      parameters = c(QQ = unique(d$QQ)),
      params.static = params.static
    )
  )
}


# Now we create the graphs
df.graph = df.result %>%
  group_by(expe, miR) %>%
  do(ODEcalc.graph(.)) %>%
  mutate(expe_set = paste0("expe",str_extract(as.character(expe), "^[0-9]+")))  # experiment

df.coefs = df.result %>% select(-time, -fraccleaved) %>% distinct()

###########################################################################

# PLOT ----
MTOgraph = df.result %>%
  left_join(namedict.print, by = "miR") %>%
  mutate(plot.group = paste(expe_set, rxnID.print, sep = ": ")) %>%
  ggplot(aes(x = time,
             group = expe)) +

  geom_hline(yintercept = 0, linewidth = 0.5, color = "gray85") +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "gray85") +

  geom_hline(
    data = df.coefs %>%
      left_join(namedict.print, by = "miR") %>%
      mutate(plot.group = paste(expe_set, rxnID.print, sep = ": ")),
    aes(
      yintercept = burst.lvl * Fa
    ),
    linewidth = 0.5, linetype = "dashed",
    alpha = 0.3, color = "dodgerblue") +

  geom_line(data = df.graph %>%
              left_join(namedict.print, by = "miR") %>%
              mutate(plot.group = paste(expe_set, rxnID.print, sep = ": ")),
            linewidth = 0.5, alpha = 0.75, color = "gray25",
            aes(y = fraccleaved.ode)) +

  geom_point(size = 0.75, stroke = 1, shape = 1, alpha = 0.75,
             color = "black",
             aes(y = fraccleaved)) +

  # Fitted parameters print
  # print in min-1
  geom_text(data = df.coefs %>%
              left_join(namedict.print, by = "miR") %>%
              mutate(plot.group = paste(expe_set, rxnID.print, sep = ": ")),
            x = Inf, y = 0,
            vjust = 0, hjust = 1,
            size = 3,
            aes(
              label = paste0(
                "koffP (min-1) = \n",
                formatCoefs(koffP),
                " (",
                formatCoefs(koffP.lo),
                "-",
                formatCoefs(koffP.hi),
                "\nQQ = \n",
                formatCoefs(QQ),
                " (",
                formatCoefs(QQ.lo),
                "-",
                formatCoefs(QQ.hi),
                ")",
                "\nFmax = \n",
                formatCoefs(Fmax),
                " (",
                formatCoefs(Fmax.lo),
                "-",
                formatCoefs(Fmax.hi),
                ")",
                "\nN = ",
                sampleN
              )
            )
            ) +

  facet_wrap("plot.group",
             ncol = 3, scales = "free") +

  scale_y_continuous(limits = c(0, NA)) +

  xlab("Time (min)") +
  ylab("Fraction sliced") +
  theme0
MTOgraph

###############
# Save output values

if(SAVE_TABLE){
  MTOfit.write = MTOfit %>%
    ungroup() %>%
    left_join(namedict.print, by = "miR") %>%
    transmute(
      expe_set  = expe_set,
      miR       = miR,
      rxnID     = rxnID,

      koffP       = koffP,
      koffP.lo    = koffP.lo,
      koffP.hi    = koffP.hi,
      koffP.pP    = -log10(koffP.p),

      QQ        = QQ,
      QQ.lo     = QQ.lo,
      QQ.hi     = QQ.hi,
      QQ.pP     = -log10(QQ.p),

      Fmax        = Fmax,
      Fmax.lo     = Fmax.lo,
      Fmax.hi     = Fmax.hi,
      Fmax.pP     = -log10(Fmax.p),

      sampleN.MTO = sampleN
    )
  write.csv(MTOfit.write, paste0(
    "./",
    "slicing_ktable_", modelmode, "-MTO",
    ".csv"),
    quote = F, row.names = F
  )
}

#@#@#
print("Done!")
