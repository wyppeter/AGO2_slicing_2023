library(tidyverse)
library(minpack.lm)

quantoligoconc = 1000  # pM

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
df.stock.approx = read.csv(paste0(
  "./",
  "AGOconc-approx.csv"
), stringsAsFactors = T) %>%
  mutate(conc.stock.approx = conc.stock.approx * 1000) # convert nM to pM

df = read.csv(paste0(
  "./",
  "FB_quant-data.csv"
  ), stringsAsFactors = T) %>%
  left_join(df.stock.approx, by = "miR") %>%
  group_by(miR) %>%
  mutate(dil = conc/conc.stock.approx)

# Threshold for saturation fit ----
LINthres = 0.6
df.thres = df %>%
  filter(fracbound <= LINthres)
miR.sat = as.character(unname(filter(df, fracbound > LINthres)$miR))
miR.F = df$miR %>% as.character() %>% unique()

# Equation for fraction bound ----
calcFB = function(At, dil, Rt, Kd, yint){
  ((At*dil+Rt+10^Kd) - sqrt((At*dil+Rt+10^Kd)^2 - 4*At*dil*Rt))/(2*Rt) * (1-yint) + yint
}

# Linear modeling ----
LMtable = df.thres %>%
  summarise(
    beta =    coef(lm(fracbound ~ dil))[[2]],
    yint =    coef(lm(fracbound ~ dil))[[1]],
    r.sq = summary(lm(fracbound ~ dil))$r.squared,
    pval = summary(lm(fracbound ~ dil))$coefficients[2,4],
    lwCI = confint(lm(fracbound ~ dil), 'dil', level=0.95)[[1]],
    hiCI = confint(lm(fracbound ~ dil), 'dil', level=0.95)[[2]]
  ) %>%
  mutate(At    = beta * quantoligoconc,
         CI.lo = lwCI * quantoligoconc,
         CI.hi = hiCI * quantoligoconc
  )

# NLS modeling ----
if(length(miR.sat) > 0){

  ## NLS fitting ----
  nlsFit = df %>%
    mutate(Rt = quantoligoconc) %>%
    filter(miR %in% miR.sat) %>%
    mutate(miR = factor(miR, levels = unique(miR.sat))) %>%
    group_by(miR) %>%
    do(data.frame(model = broom::tidy(
      nlsLM(
        formula = fracbound ~ ((At*dil+Rt+10^Kd) - sqrt((At*dil+Rt+10^Kd)^2 - 4*At*dil*Rt))/(2*Rt) * (1-yint) + yint,
        start   = c(At = 30 * 1000,   Kd = 2,  yint = 0.1),
        upper   = c(At = 1000 * 1000, Kd = 4,  yint = 1),
        lower   = c(At = 0,           Kd = -2, yint = 0),
        control = c(maxiter = 1000),
        data = .
      )
    ))) %>%
    ungroup() %>%
    transmute(miR = miR,
              coeff = model.term,
              v  = model.estimate,
              hi = model.estimate+model.std.error*1.96,
              lo = model.estimate-model.std.error*1.96,
              p = model.p.value) %>%
    pivot_wider(id_cols = miR, names_from = coeff, names_sep = "_", values_from = c(v, hi, lo, p))

  ## Create NLS fit model graph ----
  df.maxconc = df %>%
    group_by(miR) %>%
    summarize_at("dil", max) %>%
    rename(max.dil = dil)
  xmax = 0.15
  xstep = 0.001
  nlsFit.graph = data.frame(
    dil = rep(seq(0, xmax, xstep), times = length(miR.sat)),
    miR = factor(rep(miR.sat, each = xmax/xstep+1), levels = unique(miR.sat))
  ) %>%
    left_join(df.maxconc, by = "miR") %>%
    filter(dil < max.dil * 1.5)
  At.V   = as.numeric(nlsFit$v_At)
  Kd.V   = as.numeric(nlsFit$v_Kd)
  yint.V = as.numeric(nlsFit$v_yint)
  nlsFit.graph$fracbound = with(nlsFit.graph,
                                  calcFB(
                                    At.V[as.numeric(miR)],
                                    dil,
                                    quantoligoconc,
                                    Kd.V[as.numeric(miR)],
                                    yint.V[as.numeric(miR)]
                                  )
  )
}

nameconvert = function(miR){
  sub("\\.2", ", 2",
      sub("\\(", " (",
          sub("#", ", ",
              sub("_", ", ",
                  gsub("\\|", ", ",
                       idLookUp(miR, T, T))))))
}

# Plotting time ----
df %>%
  mutate(rxnID.cleaned = nameconvert(miR)) %>%
  ggplot(aes(x = dil, y = fracbound, color = miR)) +
  geom_abline(data = LMtable %>%
                mutate(rxnID.cleaned = nameconvert(miR)),
              aes(slope = beta,
              intercept = yint,
              color = miR
              ),
              linetype = "solid",
              linewidth = 0.5, alpha = 0.5,
              show.legend = F
              ) +

  geom_line(data = nlsFit.graph %>%
              mutate(rxnID.cleaned = nameconvert(miR)),
            linetype = "dashed",
            linewidth = 0.5, alpha = 0.5,
            show.legend = F) +

  geom_point(shape = 1, size = 1.2, stroke = 0.75, alpha = 0.8,
             show.legend = F) +

  facet_wrap("rxnID.cleaned", scales = "free", ncol = 7) +

  scale_x_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0.05, 0.2), add = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2), add = 0)) +

  labs(x = "Dilution from stock", y = "Fraction bound of target RNA with a seed site") +

  theme0


# Save output ----
out.table = LMtable %>%
  full_join(nlsFit, by = "miR") %>%
  mutate(lm.conc.nM = At/1000,
         nls.conc.nM = v_At/1000,
         adopted.conc.nM = if_else(!is.na(p_At) & p_At < 0.05, nls.conc.nM, lm.conc.nM)
  ) %>%
  mutate(rxnID.cleaned = nameconvert(miR) %>%
           gsub("_.+$", "", .)) %>%
  select(miR, rxnID.cleaned, adopted.conc.nM, lm.conc.nM, nls.conc.nM, everything())

out.table %>%
  write.csv(paste0(
    "./",
    "FB_quant_results.csv"
  ),
  row.names = F,
  quote = T)
