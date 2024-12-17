library(tidyverse)
library(viridis)
library(ggpubr)
library(nlstools)
library(minpack.lm)

options(scipen = 5)

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

difflim = 60 # nM-1 min-1 (10-9 M-1 s-1)

# READ DATA
parentDir = "./"
df.raw = read.csv(paste0(
  parentDir,
  "kon-FB-data.csv"), stringsAsFactors = T
) %>%
  rename(miR = miR.targ)

guides = unique(df.raw$miR)
df.dt = df.raw %>%
  mutate(time = time/60) %>%
  arrange(miR, time) %>%
  group_by(miR)

# NLS fitting (Levenberg-Marquardt)
nlsFit = df.dt %>%
  group_by(miR) %>%
  do(data.frame(model = broom::tidy(nlsLM(
        fracbound ~ A * (1-exp(-time * exp(kon) * conc)),
        algorithm = "LM",
        start = c(A = 0.8, kon = log(2)),
        upper = c(A = 1.0, kon = log(difflim)),
        lower = c(A = 0.0, kon = -Inf),
        control = nls.lm.control(maxiter = 100),
        data = .)
    ))) %>%
  ungroup() %>%
  transmute(miR   = miR,
            coef = model.term,
            val   = model.estimate,
            CI.hi = model.estimate+model.std.error*1.96,
            CI.lo = model.estimate-model.std.error*1.96) %>%
  pivot_wider(id_cols = miR,
              names_from = coef,
              names_sep = "_",
              values_from = c(val, CI.hi, CI.lo)) %>%
  transmute(miR = miR,
            A    = val_A,
            A.hi = CI.hi_A,
            A.lo = CI.lo_A,
            kon    = val_kon,
            kon.hi = CI.hi_kon,
            kon.lo = CI.lo_kon) %>%
  left_join(df.dt %>% filter(!is.na(conc)) %>%
              select(miR, conc) %>% distinct(),
            by = "miR")

# Make graph
t.plot = seq(0,4,0.01)
fitGraph = data.frame(
  time = rep(t.plot, times = length(guides)),
  miR  = rep(guides, each =  length(t.plot))
  ) %>%
  left_join(nlsFit, by = "miR") %>%
  mutate(fracbound = A * (1-exp(-time * exp(kon) * conc)))

formatCoefs = function(x){
  as.character(format(round(as.numeric(x), 3), nsmall = 3))
}

df.dt %>%
  filter(!is.na(conc)) %>%
  ggplot(aes(x = time, y = fracbound)) +
  geom_line(data = fitGraph, linewidth = 0.5, color = "gray60") +
  geom_point(shape = 1, size = 1, stroke = 1, color = "gray20") +
  geom_text(data = nlsFit, inherit.aes = FALSE,
            size = 3, x = 2.5, y = 0.4, color = "gray25",
            aes(label = paste(miR, " (", conc, " nM) :\n",
                              "A = ",
                              formatCoefs(A),
                              " (",
                              formatCoefs(A.lo),
                              "-",
                              formatCoefs(A.hi),
                              ")",
                              ",\nk_on (nM-1 min-1) = \n",
                              formatCoefs(exp(kon)),
                              " (",
                              formatCoefs(exp(kon.lo)),
                              "-",
                              formatCoefs(exp(kon.hi)),
                              ")",
                              sep = "")
            )
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  facet_wrap("miR", scales = "free") +
  labs(x = "Time (min)", y = "Fraction bound") +
  theme0
