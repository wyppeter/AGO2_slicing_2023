# Plot reporter results

library(tidyverse)
library(ggnewscale)

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

# Get data
df.hela.1 = read.csv(paste0(
  "./",
  "rep_expe_I_export.csv"
  )
) %>%
  mutate(expe = "1",
         cells = "HeLa")

df.hela.2 = read.csv(paste0(
  "./",
  "rep_expe_II_export.csv"
  )
) %>%
  mutate(expe = "2",
         cells = "HeLa")

df.hek.1 = read.csv(paste0(
  "./",
  "rep_expe_III_export.csv"
)
) %>%
  mutate(expe = "3",
         cells = "HEK293T")

df.hek.2 = read.csv(paste0(
  "./",
  "rep_expe_IV_export.csv"
)
) %>%
  mutate(expe = "4",
         cells = "HEK293T")

df = bind_rows(df.hela.1, df.hela.2,
               df.hek.1,  df.hek.2)

# Stat functions
geomean = function(x){
  exp(mean(log(x)))
}
SEM = function(x) {
  sd(x)/sqrt(length(x))
}

# Processing
df.l = df %>%
  pivot_longer(
    cols = -c(cells, expe, rxnID, rxnID.explicit),
    names_sep = "\\.",
    names_to = c("context", "iron", "rep"),
    values_to = "KD"
  ) %>%
  mutate(miR.par = gsub("\\.X[0-9]$", "", rxnID),
         miR.mut = case_when(
           rxnID %in% c("miR-7","miR-196a.M3") ~ "Fast",
           rxnID %in% c("miR-7.M4","miR-196a") ~ "Slow",
           TRUE ~ "NA"
         )
         ) %>%

  # Convert to fold repression values
  filter(!is.na(KD)) %>%
  mutate(FR = 1/KD,
         logFR = log10(FR)) %>%

  mutate(cells = factor(cells, levels = c("HeLa", "HEK293T")),
         miR.par = factor(miR.par, levels = c("miR-7", "miR-196a")),
         iron = factor(iron, levels = c("DFO", "FAC")),
         miR.mut = factor(miR.mut, levels = c("Slow","Fast"))) %>%
  arrange(cells, miR.par, iron, miR.mut, expe) %>%
  mutate(plot.group = paste(cells, miR.par, iron),
         plot.group = fct_inorder(plot.group))

df.l.mean = df.l %>%
  group_by(cells,
           rxnID, rxnID.explicit, miR.par, miR.mut,
           context, iron,
           plot.group) %>%
  summarize_at("logFR", c(mean, SEM)) %>%
  rename(mean.FR = fn1, sem.FR = fn2) %>%

  mutate(cells = factor(cells, levels = c("HeLa", "HEK293T")),
         miR.par = factor(miR.par, levels = c("miR-7", "miR-196a")),
         iron = factor(iron, levels = c("DFO", "FAC")),
         miR.mut = factor(miR.mut, levels = c("Slow","Fast")))

# Stats
df.t = df.l %>%
  group_by(plot.group,
           cells,
           miR.par,
           context, iron) %>%
  do(., {
    df.this = .
    t.out = t.test(logFR ~ miR.mut, data = df.this)
    delt = t.out$estimate[2] - t.out$estimate[1]
    pval = unname(t.out$p.value)
    df.this %>%
      select(miR.par, context, iron) %>%
      distinct() %>%
      mutate(pval = pval,
             delt = delt,
             delt.fold = 10^(delt))
  })

###############################
# Plot!
psignif = function(p){
  if_else(p < 0.05,
          if_else(p < 0.01,
                  if_else(p < 0.001,
                          if_else(p < 0.0001,
                                  "****",
                                  "***"),
                          "**"),
                  "*"),
          "n.s.")
}
df.l.mean %>%
  ggplot(aes(
    x = miR.mut, y = 10^mean.FR,
    color = miR.mut, group = interaction(context, iron)
  )) +

  facet_grid(context~plot.group) +

  geom_hline(yintercept = 1, linewidth = 0.5, color = "gray10",
             linetype = "dashed") +


  geom_line(alpha = 0.75, linewidth = 0.5, color = "black") +
  geom_errorbar(aes(ymin = 10^(mean.FR-1.96*sem.FR),
                    ymax = 10^(mean.FR+1.96*sem.FR)),
                alpha = 0.75, color = "black",
                linewidth = 0.5, width = 0.3) +
  geom_point(data = df.l,
             aes(y = FR),
             shape = 16, alpha = 0.5, size = 1.5, show.legend = F) +

  scale_color_manual(breaks = c("Slow","Fast"),
                     values = c("tomato", "forestgreen")) +


  geom_text(data = df.t,
            aes(x = 1.5, y = Inf,
                label = format(round(delt.fold, 1), nsmall = 1),
                color = pval < 0.05),
            color = "black",
            vjust = 1, hjust = 0.5, size = 3.2,
            show.legend = FALSE) +
  geom_text(data = df.t,
            aes(x = 1.5, y = 1.1,
                label = psignif(pval),
                color = pval < 0.05),
            color = "black",
            vjust = 0, hjust = 0.5, size = 3.2,
            show.legend = FALSE) +

  scale_y_continuous(trans = "log10") +
  coord_cartesian(ylim = c(1, 45)) +
  labs(x = "Guide sequence variant",
       y = "Fold knockdown in NanoLuc, normalized by firefly luciferase") +
  theme0

