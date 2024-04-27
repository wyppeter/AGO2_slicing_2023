library(tidyverse)

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

# Read data ----
PDBcodes = c("7SWF", "4NCB")
PDBdict = data.frame(
  PDB = PDBcodes,
  ago.state = c("athAGO10", "TtAgo")
)
blank = PDBdict %>%
  mutate(filler = 1) %>%
  crossing(data.frame(
    nt = seq(1,22,1)
  ))

df = read.csv(
  "pAgo_torsion_values+7swf.csv",
  stringsAsFactors = T
) %>%
  filter(PDB %in% c("4NCB", "7SWF")) %>%
  left_join(PDBdict, by = "PDB") %>%
  filter(bp.state) %>%
  # rename(nt = pos) %>%
  mutate(PDB = factor(PDB, levels = PDBcodes)) %>%
  full_join(blank, by = c("PDB", "ago.state", "nt")) %>%
  arrange(PDB, nt)

# All-angle analysis ----
anglesOI = c(
  "delta", "epsilon", "v2", "zeta", "alpha", "beta", "gamma"
)

df.angles = df %>%
  select(PDB, ago.state, nt, all_of(anglesOI))

df.angles.l = df.angles %>%
  pivot_longer(cols = -c("PDB", "ago.state", "nt"),
               names_to = "angle") %>%
  mutate(PDB = factor(PDB, levels = PDBdict$PDB),
         ago.state = factor(ago.state, levels = PDBdict$ago.state)) %>%
  mutate(angle = factor(angle, levels = anglesOI)) %>%
  arrange(PDB, angle)

df.angles.l %>%
  ggplot(aes(x = nt, y = ago.state,
             fill = value)) +
  geom_tile() +
  facet_wrap("angle",
             ncol = 1,
             strip.position = "left") +
  scale_fill_gradientn(
    colors = c("gray10","dodgerblue3","gray80","tomato","gray10"),
    breaks = seq(-180,180,90),
    limits = c(-180,180),
    na.value = "white"
  ) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(1,21),
                     breaks = c(2,6,7,10,13,16,21)) +
  theme0 +
  theme(strip.placement = "outside")

# Conformer analysis ----
df.conf = df %>%
  select(PDB, ago.state, nt, gamma, delta) %>%
  pivot_longer(cols = -c("PDB", "ago.state", "nt"),
               names_to = "angle") %>%
  mutate(PDB = factor(PDB, levels = PDBdict$PDB),
         ago.state = factor(ago.state, levels = PDBdict$ago.state)) %>%

  mutate(nt = if_else(angle == "delta", nt, nt - 0.5)) %>%  # gamma is a suite angle
  mutate(value = if_else(nt < 2, NA_real_, value)) %>%

  mutate(value = if_else(value < 0, value + 360, value)) %>%  # formatting
  arrange(angle, rev(PDB), nt)

df.refs = read.csv(
  "ref-angles.csv",
  stringsAsFactors = T
) %>%
  mutate(value = if_else(value < 0, value + 360, value),
         lower = if_else(lower < 0, lower + 360, lower),
         upper = if_else(upper < 0, upper + 360, upper))

df.conf %>%
  arrange(PDB) %>%
  ggplot(aes(x = nt, y = value,
             color = ago.state)) +
  facet_wrap("angle",
             scales = "free",
             ncol = 1,
             strip.position = "left") +

  # Reference angles
  geom_rect(data = df.refs %>% filter(ref.type == "conformer"),
            inherit.aes = F,
            aes(xmin = -Inf, xmax = Inf,
                ymin = lower, ymax = upper),
            fill = "gray95"
            ) +
  geom_hline(data = df.refs %>% filter(ref.type == "conformer"),
             aes(yintercept = value),
             color = "gray85", linewidth = 0.5
  ) +

  # Normal A-form angles
  geom_errorbar(data = df.refs %>% filter(ref.type == "Aform"),
                aes(x = 22.5, ymin = lower, ymax = upper),
                color = "turquoise", width = 0, linewidth = 2
  ) +
  geom_point(data = df.refs %>% filter(ref.type == "Aform"),
             aes(x = 22.5),
             color = "turquoise4", size = 1, shape = 3
  ) +

  # Torsion angles plots
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_point(aes(shape = ago.state), stroke = 1.2, size = 1.2, alpha = 0.8) +

  scale_color_manual(breaks = PDBdict$ago.state,
                     values = c("#BB0022", "#4411FF")) +
  scale_shape_manual(breaks = PDBdict$ago.state,
                     values = c(5,4)) +

  scale_x_continuous(limits = c(1,23),
                     breaks = c(2,6,7,10,13,16,21)) +
  scale_y_continuous(breaks = seq(-360, 360, 60),
                     expand = c(0,30)
                     # limits = c(0, 360)
                     ) +
  labs(x = "Guide nucleotide position", y = "Torsion angle") +
  theme0 +
  theme(strip.placement = "outside")

