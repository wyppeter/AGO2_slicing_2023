# Analyses and plotting of k coefficients
# Bartel Lab
# Peter Y Wang 2023

# Table of contents ----
#  k1: kslice values with seq dets indicated
#  k2: Seq det fold changes
# k2b: 16-bp fold changes upon miR-124 9~12 subst
#  k3: Fold changes of mm vs seq dets
# k3b: Fold change magnitudes of k4 in barplot form
#  k4: kon for miR-7 v1 vs v2 targets
# k4b: kon for lsy-6/miR-124 chimera
# k5a: Correlation of kslice and kon
# k5b: Correlation of kslice and kphase2
#  k6: Head-to-head of k coefs for each miR
#  k7: S387mut coefs

## Libraries ----
library(tidyverse)
library(ggforce)
library(ggnewscale)

## Setup and data input ----
workDir = "./"

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
formatCoefs.print = function(x){
  sub("\\.$", "", sprintf("%#.2g", signif(x, 2)))
}

# Initial data wrangling
df.offP = read.csv(paste0(
  workDir,
  "slicing_ktable_CTS-MTO.csv")) %>%
  select(miR, starts_with(c("koffP")))

df.k = read.csv(paste0(
  workDir,
  "slicing_ktable.csv")) %>%
  select(-rxnID) %>%
  left_join(namedict, by = "miR") %>%
  mutate(mainrxn = rxnID == guidename)

k_TO_PLOT = "kcat"
df.k.l = df.k %>%
  pivot_longer(cols = starts_with(c("kon","kcat","kph2","koffP","Fa")),
               names_sep = "\\.",
               names_to = c("coeff", "piece")
               ) %>%
  mutate(piece = case_when(
    is.na(piece) ~ "val",
    piece %in% c("lo", "hi") ~ paste0("CI.", piece),
    TRUE ~ piece
  )
  ) %>%
  pivot_wider(names_from = piece, values_from = value)
df.k.this = df.k.l %>%
  filter(coeff == k_TO_PLOT) %>%
  arrange(val) %>%
  mutate(rxnID = fct_inorder(rxnID),
         rxnID.explicit = fct_inorder(rxnID.explicit))

############### K PLOTS ########################################################

# k + seq dets plot ----
otherkcat = data.frame(
  rxnID = c("hsa-let-7a.21mer", "hsa-miR-21.21mer"),
  rxnID.explicit = c("hsa-let-7a.21mer", "hsa-miR-21.21mer"),
  val   = c(3.30, 5.22),
  mainrxn = c(NA, NA),
  W7 = c(T,T),
  R10 = c(T,T),
  W17 = c(T,T),
  short = c(T,T)
)
df.k.this.dets = df.k.l %>%
  filter(coeff == "kcat", mainrxn) %>%
  mutate(seq = seqLookUp(miR.base)) %>%
  mutate(
    W7 = detLookUp(seq, "W7") == 1,
    R10 = detLookUp(seq, "R10") == 1,
    W17 = detLookUp(seq, "W17") == 1,
    short = detLookUp(seq, "len") <= 22
  ) %>%
  arrange(val) %>%
  bind_rows(otherkcat, .) %>%  # add Becker data
  mutate(rxnID = fct_inorder(rxnID),
         rxnID.explicit = fct_inorder(rxnID.explicit))
spacer = 1.25
sqs.start = 30
k1 = df.k.this.dets %>%
  ggplot(aes(x = val, y = rxnID)) +
  geom_errorbarh(aes(xmin = CI.lo, xmax = CI.hi),
                 color = "gray50", linewidth = 0.5, height = 0.5) +
  geom_point(shape = 18, color = "orangered3", size = 2.25) +
  geom_text(aes(label = formatexplicitID(rxnID.explicit),
                x = if_else(val > 0.85,
                            0.015,
                            15),
                hjust = if_else(val > 0.85, 0, 1)),
            vjust = 0.5, size = 2.5, color = "black") +

  geom_point(
    aes(x = sqs.start*spacer^0, color = R10),
    shape = 15, size = 2.5, show.legend = F) +
  geom_point(
    aes(x = sqs.start*spacer^1, color = W17),
    shape = 15, size = 2.5, show.legend = F) +
  geom_point(
    aes(x = sqs.start*spacer^2, color = W7),
    shape = 15, size = 2.5, show.legend = F) +
  geom_point(
    aes(x = sqs.start*spacer^3, color = short),
    shape = 15, size = 2.5, show.legend = F) +

  scale_color_manual(values = c("salmon1","forestgreen")) +
  scale_x_continuous(trans = "log10",
                     breaks = 10^seq(-10,10,1)) +
  coord_cartesian(xlim = c(0.02, 80)) +
  xlab(k_TO_PLOT) +
  theme0 +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(size = 12, colour = "black", angle = 90),
        panel.grid.major.x = element_line(colour = "gray90"))
k1

############### INITIAL LIST ####################################

originalSet = c(
  "let-7a",
  "miR-451a",
  "lsy-6",
  "miR-1",
  "miR-7",
  "miR-122",
  "miR-124",
  "miR-20a.23mer",
  "miR-430a",
  "miR-671.23mer",
  "miR-196a",
  "miR-302a.23mer",
  "miR-155.23mer"
)
df.seqVSkcat = df.k.this.dets %>%
  filter(rxnID %in% originalSet) %>%
  select(rxnID, seq, val) %>%
  mutate(seq.s = seq) %>%
  rename(kcat = val) %>%
  arrange(-kcat)

# Write list
write.csv(df.seqVSkcat,
          file = paste0(
            workDir,
            "WT_miRs_kslice_rank",
            ".csv"
          ), row.names = F)

############### COMPARISONS ####################################################

# Perturbation wrangling ----
suffixregex = "[\\._#]([a-zA-Z0-9\\|\\(\\)~]+)$"
df.kcomp.perturb = df.k.this %>%
  mutate(
    ID.base    = sub(suffixregex, "", rxnID),
    ID.perturb = str_match(rxnID.explicit, suffixregex)[,2]
  ) %>%
  mutate(perturb.group = case_when(

    grepl("^alt", ID.perturb) ~ "alt.targ",
    ID.perturb == "UUCflank" ~ "alt.targ.flank",
    grepl("[A-Z]{2}prep$", ID.perturb) ~ "alt.prot",

    grepl("^S387[A-Z]$", ID.perturb) ~ "S387mut",

    grepl("^mm6[AUGC]{2}$", ID.perturb) ~ "mm6",
    grepl("^mm7[AUGC]{2}$", ID.perturb) ~ "mm7",
    grepl("^mm17[AUGC]{2}$", ID.perturb) ~ "mm17",

    ID.perturb == "16bp" ~ "16bp",

    grepl("\\([AUGC]7[AUGC]\\)", ID.perturb) ~ "W7",
    grepl("\\([AUGC]10[AUGC]\\)", ID.perturb) ~ "R10",
    grepl("\\([AUGC]17[AUGC]\\)", ID.perturb) ~ "W17",
    grepl("\\([AUGC]7[AUGC]\\|[AUGC]10[AUGC]\\|[AUGC]17[AUGC]\\)", ID.perturb) ~ "W7|R10|W17",
    grepl("^[0-9]+mer$", ID.perturb) ~ "len",
    grepl("\\([AUGC]10[AUGC]\\|[AUGC]11[AUGC]\\)", ID.perturb) ~ "R10|Y11",

    grepl("\\([SW][0-9~]+[SW]\\)", ID.perturb) ~ "dom.mut",

    T ~ ID.perturb
  )) %>%
  mutate(
    ID.perturb.print = if_else(substr(ID.perturb,1,1) == "M",
                               gsub("\\(", " (",
                                    gsub("\\|", ", ",
                                         ID.perturb)),
                               ID.perturb)

  )

perturb.groups = c(
  "R10", "W17", "W7",
  "W7|R10|W17",
  "len",
  "R10|Y11", "dom.mut",
  "mm6", "mm7", "mm17",
  "16bp",
  "S387mut",
  "alt.prot", "alt.targ", "alt.targ.flank"
)
comp.base.dict = data.frame(
  perturb.group = perturb.groups,
  ID.perturb = c(
    ".WT", ".WT", ".WT",
    ".WT",
    ".22mer",
    ".WT", ".WT",
    ".perfect", ".perfect", ".perfect",
    ".perfect",
    ".S387S",
    "...", "...", "..."
  )
)

df.kcomp.base = df.kcomp.perturb %>%
  select(-ID.base, -ID.perturb, -ID.perturb.print, -perturb.group) %>%
  crossing(comp.base.dict) %>%
  mutate(ID.base = rxnID,
         ID.perturb.print = ID.perturb)

df.kcomp.raw = bind_rows(
  df.kcomp.base,
  df.kcomp.perturb
) %>%
  group_by(perturb.group, ID.base) %>%
  filter(n() > 1 |
           (ID.base == "miR-196a.M2" & perturb.group == "mm7")) %>%

  left_join(
    df.kcomp.base %>%
      transmute(
        ID.base = ID.base,
        kcat.base = val
      ) %>%
      distinct(),
    by = "ID.base") %>%
  mutate(k.ratio = val/kcat.base)

df.kcomp.control = df.kcomp.raw %>%
  filter(perturb.group %in% c(
    "alt.prot",
    "alt.targ",
    "alt.targ.flank"
  ),
  !mainrxn) %>%
  ungroup() %>%
  select(rxnID, k.ratio) %>%
  mutate(lr = log10(k.ratio))
control.CI = 10^(mad(df.kcomp.control$lr, center = 0) * 1.96)

df.kcomp = df.kcomp.raw %>%
  mutate(FC.dir = factor(
           if_else(
             abs(log(k.ratio)) > abs(log(control.CI)),
               if_else(k.ratio > 1, +1, -1),
             0
             )
           ),
         FC = if_else(k.ratio > 1, k.ratio, 1/k.ratio),
         FC.print =
           if_else(k.ratio == 1, "",
             paste(
             if_else(k.ratio > 1, "^", "v"),
             format(round(FC, 1), nsmall = 1)
             )
           )
         ) %>%
  mutate(perturb.group = factor(perturb.group, levels = perturb.groups)) %>%
  mutate(plot.group = paste(perturb.group, ID.base,
                                  sep = "   ")) %>%
  arrange(perturb.group, ID.base) %>%
  mutate(plot.group = fct_inorder(plot.group))

df.kcomp.subset = df.kcomp %>%
  filter(!perturb.group %in% c(
    "S387mut",
    "16bp",
    "dom.mut"
  ),
  !grepl("^mm[0-9]+$", as.character(perturb.group))
  )

# Plot here!
k2 = df.kcomp.subset %>%
  ggplot(aes(x = val, y = ID.perturb.print)) +
  facet_col("plot.group",
            space = "free",
            scales = "free_y",
            strip.position = "left") +

  geom_hline(yintercept = -Inf, linewidth = 0.5, color = "black") +

  geom_errorbarh(aes(xmin = CI.lo, xmax = CI.hi),
                 color = "gray50", linewidth = 0.5, height = 0.5) +
  geom_point(aes(color = k.ratio == 1), shape = 18, size = 2.25, show.legend = F) +
  scale_color_manual(breaks = c(T,F),
                     values = c("gray60","darkslategray")) +

  new_scale_color() +

  geom_text(aes(label = FC.print, color = FC.dir, x = Inf), hjust = 1, vjust = 0.5, size = 3, show.legend = F) +
  scale_color_manual(breaks = c(+1, 0, -1),
                     values = c("forestgreen", "gray70", "tomato")) +

  scale_x_continuous(trans = "log10", expand = c(0.1, 0), breaks = 10^seq(-10,10,1)) +
  scale_y_discrete(limits = rev) +
  theme0 +
  theme(
    strip.placement = "outside",
    axis.title.y = element_blank(),
    strip.text.y.left = element_text(size = 9, colour = "black", angle = 0),
    panel.grid.major.x = element_line(colour = "gray90"))
k2

# 9-12 mut of miR-124
df.kcomp.dom.mut.16bp = df.kcomp %>%
  filter(perturb.group == "16bp",
         rxnID %in% c(
           "miR-124",
           "miR-124_16bp",
           "miR-124.M2",
           "miR-124.M2_16bp"
         ))
k2b = df.kcomp.dom.mut.16bp %>%
  ggplot(aes(x = val, y = ID.perturb.print)) +
  facet_col("plot.group",
            space = "free",
            scales = "free_y",
            strip.position = "left") +

  geom_hline(yintercept = -Inf, linewidth = 0.5, color = "black") +

  geom_errorbarh(aes(xmin = CI.lo, xmax = CI.hi),
                 color = "gray50", linewidth = 0.5, height = 0.5) +
  geom_point(aes(color = k.ratio == 1), shape = 18, size = 3.5, show.legend = F) +
  scale_color_manual(breaks = c(T,F),
                     values = c("gray60","darkslategray")) +

  new_scale_color() +

  geom_text(aes(label = FC.print, color = FC.dir, x = Inf), hjust = 1, vjust = 0.5, size = 3, show.legend = F) +
  scale_color_manual(breaks = c(+1, 0, -1),
                     values = c("forestgreen", "gray70", "tomato")) +

  scale_x_continuous(trans = "log10", expand = c(0.1, 0), breaks = 10^seq(-10,10,1)) +
  scale_y_discrete(limits = rev) +
  theme0 +
  theme(
    strip.placement = "outside",
    axis.title.y = element_blank(),
    strip.text.y.left = element_text(size = 9, colour = "black", angle = 0),
    panel.grid.major.x = element_line(colour = "gray90"))
k2b


############### SPECIFIC COMPARISONS: mmVdets ###########################################

# mm vs dets ----
df.mmVdets.perturb = df.kcomp %>%
  ungroup() %>%
  filter(
    grepl("^mm[0-9]+$", perturb.group)
    ) %>%
  left_join(
    df.k.this.dets %>%
      select(rxnID, seq, W7, R10, W17) %>%
      rename(ID.base = rxnID),
    by = "ID.base"
  ) %>%
  mutate(perturb.supergroup =
           if_else(perturb.group == "mm17",
                   "mm17",
                   "mm6/7"
                   )) %>%
  select(-ID.perturb.print, -plot.group, -mainrxn,
         -kcat.base, -k.ratio) %>%
  rename(bp.pattern = ID.perturb,
         FC_vsPerf.dir = FC.dir,
         FC_vsPerf = FC,
         FC_vsPerf.print = FC.print
         ) %>%
  mutate(ID.basefam =
           gsub("\\.[A-Za-z0-9]+$", "", ID.base)) %>%
  distinct() %>%
  mutate(det.here =
           if_else(
             perturb.supergroup == "mm6/7",
             W7,
             W17
           )) %>%
  mutate(ID.wt = if_else(substr(bp.pattern, 1, 1) == ".",
                         ID.basefam,
                         paste(ID.basefam, perturb.group, sep = "_"))
    ) %>%
  arrange(perturb.supergroup, ID.basefam, ID.base, bp.pattern)
df.mmVdets.base = df.mmVdets.perturb  %>%
  filter(ID.wt == gsub("[AUGC]{2}$", "", rxnID)) %>%
  select(perturb.supergroup, val, ID.wt) %>%
  distinct() %>%
  rename(kcat.wt = val)
df.mmVdets = df.mmVdets.perturb %>%
  left_join(df.mmVdets.base,
            by = c("perturb.supergroup", "ID.wt")) %>%
  mutate(k.vsWT.ratio = val/kcat.wt,
         FC_vsWT.dir = factor(
           if_else(
             abs(log(k.vsWT.ratio)) > abs(log(1.5)),
             if_else(k.vsWT.ratio > 1, +1, -1),
             0
           )
         ),
         FC_vsWT = if_else(k.vsWT.ratio > 1, k.vsWT.ratio, 1/k.vsWT.ratio),
         FC_vsWT.print =
           if_else(k.vsWT.ratio == 1, "",
                   paste(
                     if_else(k.vsWT.ratio > 1, "^", "v"),
                     format(round(FC_vsWT, 1), nsmall = 1)
                   )
           )
  ) %>%
  mutate(targ.group = str_match(rxnID, "_(mm[0-9]{1,2})")[,2],
         targ.group = if_else(is.na(targ.group), ".perfect", targ.group)) %>%
  select(-perturb.group) %>%
  select(perturb.supergroup, targ.group,
         ID.basefam, ID.base, bp.pattern, det.here,
         rxnID, ID.wt,
         FC_vsPerf.dir, FC_vsPerf, FC_vsPerf.print,
         FC_vsWT.dir, FC_vsWT, FC_vsWT.print,
         everything()
         ) %>%
  distinct() %>%
  mutate(plot.group = paste(
    perturb.supergroup, ID.basefam, sep = " "
  )) %>%
  mutate(ID.base = if_else(ID.base == "miR-7.24mer", "miR-7.X24mer", ID.base)) %>%
  arrange(perturb.supergroup, ID.basefam, ID.base, bp.pattern) %>%
  mutate(rxnID.explicit = gsub("\\.|_", "    ", as.character(rxnID.explicit))) %>%
  mutate(rxnID.explicit = fct_inorder(as.character(rxnID.explicit)))


k3 = df.mmVdets %>%
  ggplot(aes(x = val, y = rxnID.explicit)) +
  facet_col("plot.group",
            space = "free",
            scales = "free_y",
            strip.position = "left") +

  geom_hline(yintercept = -Inf, linewidth = 0.5, color = "black") +

  geom_errorbarh(aes(xmin = CI.lo, xmax = CI.hi),
                 color = "gray50", linewidth = 0.5, height = 0.5) +
  geom_point(aes(color = targ.group), shape = 18, size = 2.25) +

  scale_color_manual(
    values = c("gray60", "magenta4", "violetred1", "cyan3"),
    breaks = c(".perfect", "mm6", "mm7", "mm17")) +

  new_scale_color() +

  geom_text(aes(label = FC_vsPerf.print, color = FC_vsPerf.dir, x = 50), hjust = 0, vjust = 0.5, size = 3, show.legend = F) +
  
  scale_color_manual(breaks = c(+1, 0, -1),
                     values = c("forestgreen", "gray70", "tomato")) +


  scale_x_continuous(trans = "log10", expand = c(0.1, 0),
                     breaks = 10^seq(-10,10,1)) +
  scale_y_discrete(limits = rev) +
  theme0 +
  theme(
    strip.placement = "outside",
    axis.title.y = element_blank(),
    axis.text.y = element_text(hjust = 0),
    strip.text.y.left = element_text(size = 9, colour = "black", angle = 0),
    panel.grid.major.x = element_line(colour = "gray90"))
k3

# Bar plots for fold changes in the mm vs dets analysis for 6/7
df.mmVdets.bars.base = df.mmVdets %>%
  filter(perturb.supergroup == "mm6/7") %>%
  select(rxnID, rxnID.explicit,
         ID.basefam, ID.base, ID.wt, bp.pattern,
         det.here,
         FC_vsPerf.dir, FC_vsPerf,
         FC_vsWT.dir, FC_vsWT,
         val, CI.lo, CI.hi
         ) %>%
  pivot_longer(cols = ends_with(".dir"),
               names_prefix = "FC_vs",
               names_to = "compref.dir",
               values_to = "dirFC") %>%
  mutate(compref.dir = sub("\\.dir$", "", compref.dir)) %>%
  pivot_longer(cols = starts_with("FC"),
               names_prefix = "FC_vs",
               names_to = "compref",
               values_to = "valFC") %>%
  filter(compref == compref.dir) %>%
  mutate(CI.rad = (log2(CI.hi) - log2(CI.lo)) / 2)

df.mmVdets.bars = df.mmVdets.bars.base %>%
  filter(valFC != 1) %>%
  select(-compref.dir) %>%
  mutate(merger = if_else(compref == "Perf", ID.base, ID.wt)) %>%
  left_join(df.mmVdets.bars.base %>%
              mutate(
                merger = if_else(compref == "Perf", rxnID, sub("GG$", "", rxnID))
                ) %>%
              transmute(
                merger = merger,
                compref = compref,
                CI.rad.ref = CI.rad),
            by = c("compref", "merger")) %>%
  mutate(logFC = log2(valFC ^ if_else(dirFC == -1, -1, 1)),
         CI.rad.total = sqrt(CI.rad^2 + CI.rad.ref^2),
         CI.FC.lo = logFC - CI.rad.total,
         CI.FC.hi = logFC + CI.rad.total) %>%
  mutate(comparison.here = if_else(
    compref == "Perf",
    bp.pattern,
    str_extract(rxnID.explicit, "M[0-9][0-9AUGC\\(\\|]+\\)")
  )) %>%
  mutate(plotgroup = interaction(compref, ID.basefam),
         mmVdet = paste(
           gsub("\\(", " (",
                gsub("\\|", ", ",
                str_extract(rxnID.explicit, "M[0-9][0-9AUGC\\(\\|]+\\)")
                )),
           sub("^\\.p", "P", bp.pattern),
           det.here,
           sep = "\n") %>%
           sub("^NA", "WT", .)
         ) %>%
  group_by(plotgroup) %>%
  mutate(order = 1:n()) %>%
  complete(order = 1:3, fill = list(Value = 0, mmVdet = "")) %>%
  ungroup() %>%
  mutate(order = 1:n())

k3b = df.mmVdets.bars %>%
  ggplot(aes(y = 2^logFC,
             x = order,
             fill = dirFC)) +
  facet_wrap("plotgroup", ncol = 2, scales = "free") +
  geom_rect(ymin = log2(control.CI^(-1)), ymax = log2(control.CI),
            xmin = -Inf, xmax = Inf,
            fill = "gray80") +
  geom_col(position = "dodge", width = 0.75,
           show.legend = F) +
  geom_errorbar(aes(ymin = 2^CI.FC.lo, ymax = 2^CI.FC.hi),
                color = "black", width = 0.2) +
  geom_text(aes(label = formatCoefs.print(valFC), color = dirFC,
                vjust = if_else(dirFC == +1, -0.5, +1.5)),
            size = 4.25, show.legend = F
            ) +
  geom_hline(yintercept = 2^0, color = "black", linetype = "dashed") +
  coord_cartesian(ylim = 2^c(-3.2, 3.2)) +
  scale_y_continuous(
    trans = "log2",
    breaks = scales::trans_breaks("log2", function(x) 2^x),
    labels = scales::trans_format("log2", scales::math_format(2^.x))
  ) +
  scale_x_continuous(
    breaks = df.mmVdets.bars$order,
    labels = df.mmVdets.bars$mmVdet) +
  scale_fill_manual(breaks = c(+1, 0, -1),
                     values = c("forestgreen", "gray20", "tomato")) +
  scale_color_manual(breaks = c(+1, 0, -1),
                    values = c("forestgreen", "gray20", "tomato")) +
  labs(y = "Fold change (log2)", x = "Perturbation") +
  theme0
k3b


# kon for miR-7 v1 vs v2 target with diff structure
df.miR7.kon = df.k.l %>%
  filter(coeff == "kon", rxnID %in% c(
    "miR-7", "miR-7_altTarg",
    "miR-7.24mer", "miR-7.24mer_altTarg"
  )) %>%
  mutate(targ.vers = if_else(
    grepl("_altTarg$", rxnID),
    "v1",
    "v2"
  ))
k4 = df.miR7.kon %>%
  ggplot(aes(x = val, y = targ.vers)) +
  facet_col("guidename",
            space = "free",
            scales = "free_y",
            strip.position = "left") +

  geom_hline(yintercept = -Inf, linewidth = 0.5, color = "black") +

  geom_errorbarh(aes(xmin = CI.lo, xmax = CI.hi),
                 color = "gray50", linewidth = 0.5, height = 0.5) +
  geom_point(aes(color = targ.vers),
             shape = 18, size = 2.25,
             show.legend = F) +

  scale_color_manual(
    values = c("gray60", "olivedrab4"),
    breaks = c("v2", "v1")) +
  scale_x_continuous(trans = "log10", expand = c(0.1, 0),
                     breaks = c(1,3,10)) +
  scale_y_discrete(limits = rev) +
  theme0 +
  theme(
    strip.placement = "outside",
    axis.title.y = element_blank(),
    axis.text.y = element_text(hjust = 0),
    strip.text.y.left = element_text(size = 9, colour = "black", angle = 0),
    panel.grid.major.x = element_line(colour = "gray90"))
k4

# Looking at lsy-6/miR-124 chimera
df.chim.kon = df.k.l %>%
  filter(coeff == "kon", rxnID %in% c(
    "lsy-6", "lsy-6/miR-124", "miR-124"
  ))
k4b = df.chim.kon %>%
  ggplot(aes(x = val, y = rxnID.explicit)) +

  annotate("rect",
           xmin = difflim, xmax = Inf,
           ymin = -Inf, ymax = Inf,
           fill = "black", color = NA, alpha = 0.1) +

  geom_errorbarh(aes(xmin = CI.lo, xmax = CI.hi),
                 color = "gray50", linewidth = 0.5, height = 0.5) +
  geom_point(shape = 18, size = 3,
             show.legend = F) +

  scale_x_continuous(trans = "log10", expand = c(0.1, 0),
                     breaks = c(1,3,10,30)) +
  scale_y_discrete(limits = rev) +
  theme0 +
  theme(
    strip.placement = "outside",
    axis.title.y = element_blank(),
    axis.text.y = element_text(hjust = 0),
    strip.text.y.left = element_text(size = 9, colour = "black", angle = 0),
    panel.grid.major.x = element_line(colour = "gray90"))
k4b

################### k-k CORRELATIONS #######################

# Do k values correlate? ----

### kcat <-> kon PREP ----
kcatVkon = df.k %>%
  filter(mainrxn) %>%
  filter(kcat.pP > -log10(0.05),
         kon.pP > -log10(0.05) | kon == difflim)
kcatVkon.lm = lm(
  formula = log10(kon) ~ log10(kcat),
  data = kcatVkon
)

### kcat <-> kph2 PREP ----
kcatVkph2 = df.k %>%
  filter(mainrxn) %>%
  filter(kcat.pP > -log10(0.05),
         kph2.pP > -log10(0.05))
kcatVkph2.lm = lm(
  formula = log10(kph2) ~ log10(kcat),
  data = kcatVkph2
)

kcatVk.line = data.frame(
  kcat = 10^seq(-2,1.5,0.01)) %>%
  mutate(kon = 10^(log10(kcat) * coef(kcatVkon.lm)[2] +
                     coef(kcatVkon.lm)[1])) %>%
  mutate(kph2 = 10^(log10(kcat) * coef(kcatVkph2.lm)[2] +
                      coef(kcatVkph2.lm)[1]))

### Plot corr ----
k5a = kcatVkon %>%
  ggplot(aes(
    x = kcat,
    y = kon
    )) +
  geom_hline(yintercept = difflim, linetype = "dashed") +
  geom_line(
    data = kcatVk.line,

    color = "black", alpha = 0.6
    ) +
  geom_errorbar(aes(
    ymin = kon.lo, ymax = kon.hi
  ),
  color = "midnightblue", alpha = 0.4, width = 0
  ) +
  geom_errorbarh(aes(
    xmin = kcat.lo, xmax = kcat.hi
  ),
  color = "midnightblue", alpha = 0.4, height = 0
  ) +
  geom_point(color = "midnightblue", alpha = 0.8,
             shape = 16, size = 1) +
  annotate("text", size = 2, x = 0, y = 0, hjust = 0, vjust = 0,
           label = paste0("R2 = ",
                          summary(kcatVkon.lm)$r.squared %>% signif(2),
                          "\np = ",
                          summary(kcatVkon.lm)$coefficients[2,4] %>% signif(2),
                          "\nN = ",
                          nrow(kcatVkon)
                          )) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  theme0
k5a


k5b = kcatVkph2 %>%
  ggplot(aes(
    x = kcat,
    y = kph2
  )) +
  geom_line(
    data = kcatVk.line,
    
    color = "black", alpha = 0.6
  ) +
  geom_errorbar(aes(
    ymin = kph2.lo, ymax = kph2.hi
  ),
  color = "hotpink4", alpha = 0.4, width = 0
  ) +
  geom_errorbarh(aes(
    xmin = kcat.lo, xmax = kcat.hi
  ),
  color = "hotpink4", alpha = 0.4, height = 0
  ) +
  geom_point(color = "hotpink4", alpha = 0.8,
             shape = 16, size = 1) +
  annotate("text", size = 2, x = 0, y = 0, hjust = 0, vjust = 0,
           label = paste0("R2 = ",
                          summary(kcatVkph2.lm)$r.squared %>% signif(2),
                          "\np = ",
                          summary(kcatVkph2.lm)$coefficients[2,4] %>% signif(2),
                          "\nN = ",
                          nrow(kcatVkph2)
           )) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  theme0
k5b


#################### COMPARING kslice and koffP ####################

df.k.battle = df.k.l %>%
  filter(mainrxn,
         !is.na(val),
         coeff %in% c("kcat", "koffP")) %>%
  ungroup() %>%
  mutate(coeff = factor(coeff, levels = c("kcat", "koffP"))) %>%
  mutate(withoffP = miR %in% df.offP$miR) %>%
  arrange(withoffP, coeff, rxnID)

k6 = df.k.battle %>%
  ggplot(aes(y = val, x = coeff,
             alpha = withoffP,
             group = rxnID)) +
  geom_errorbar(aes(ymin = CI.lo, ymax = CI.hi),
                 linewidth = 0.5, width = 0,
                 color = "black",
                 show.legend = F) +
  geom_point(shape = 5, size = 1, stroke = 1,
             color = "black",
             show.legend = F) +
  geom_path(linewidth = 0.5,
            color = "black",
            show.legend = F) +
  scale_alpha_manual(values = c(0.3, 1.0)) +
  scale_y_continuous(trans = "log10", expand = c(0.1, 0)) +
  labs(x = "Kinetic parameter",
       y = "Value (min-1)") +
  theme0
k6

#################### S387X analyses #######################

df.S387 = df.k.l %>%
  filter(rxnID %in% c(
    "miR-7.M1",
    "miR-7.M1#S387A",
    "miR-7.M1#S387D",
    "lsy-6",
    "lsy-6#S387A",
    "lsy-6#S387D"
  ),
  coeff %in% c(
    "kon",
    "kcat",
    "koffP"
  )) %>%
  mutate(S387 = factor(if_else(mainrxn,
                        "S",
                        str_sub(rxnID, -1, -1)),
                       levels = c("D", "A", "S")),
         coeff = factor(coeff, levels = c(
           "kon",
           "kcat",
           "koffP"
         ))) %>%
  arrange(coeff, guidename, S387)

k7 = df.S387 %>%
  ggplot(aes(x = val, y = S387)) +
  geom_errorbarh(aes(xmin = CI.lo, xmax = CI.hi),
                 color = "gray50", linewidth = 0.5, height = 0.5) +
  geom_point(aes(color = S387), shape = 18, size = 2.25,
             show.legend = F) +
  scale_x_continuous(trans = "log10", expand = c(0.2, 0)) +
  scale_color_manual(
    values = c("gray60", "royalblue", "darkgoldenrod"),
    breaks = c("S", "A", "D")) +
  facet_grid(guidename~coeff, scales = "free_x") +
  labs(y = "Residue 387",
       x = "Value (nM-1 min-1 or min-1)") +
  theme0
k7
