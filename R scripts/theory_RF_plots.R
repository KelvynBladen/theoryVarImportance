library(tidyr)
library(tidyverse)
library(GGally)
library(patchwork)
library(latex2exp)

mc <- 100
sd <- 0.01

Delta <- seq(0, 0.99, by = .11)
n_set <- round(2 * 10^(seq(from = 1, to = 3, length = 5)))
p_set <- 1:4 * 3

############################ Load Data ############################
####### These files can be obtained from running sum_RF.R #########

loco <- read.csv("../data/loco.csv")
pap <- read.csv("../data/pap.csv")
loco_rf <- read.csv("../data/loco_rf.csv")
pap_rf <- read.csv("../data/pap_rf.csv")

d_vec <- rep(Delta, times = length(p_set) * length(n_set))
p_vec <- rep(rep(p_set, each = length(Delta)), times = length(n_set))
n_vec <- rep(n_set, each = length(Delta) * length(p_set))

pap_theory <- 1 * sqrt(2 * (1 + (p_vec - 1) * d_vec^2))
loco_theory <- 1 * (1 - d_vec) *
  sqrt((1 + (p_vec - 1) * d_vec)^2 /
         ((1 + (p_vec - 2) * d_vec)^2 + (p_vec - 1) * d_vec^2))

t_theory <- (1 / sd) *
  sqrt(((n_vec - 1) * ((1 - d_vec)^2) * (1 + (p_vec - 1) * d_vec)^2) /
         ((1 + (p_vec - 1) * d_vec) *
            (1 + (p_vec - 3) * d_vec) + p_vec * d_vec^2))

df <- data.frame(
  Delta = rep(d_vec, mc),
  p = rep(p_vec, mc),
  n = rep(n_vec, mc),
  PaP = as.vector(t(pap)),
  PaP_rf = as.vector(t(pap_rf)),
  pap_theory = rep(pap_theory, mc),
  pap_theory_fix = rep(pap_theory, mc) *
    (rep(n_vec, mc) - 1) / rep(n_vec, mc),
  LOCO = as.vector(t(loco)),
  LOCO_rf = as.vector(t(loco_rf)),
  loco_theory = rep(loco_theory, mc),
  loco_theory_fix = rep(loco_theory, mc) *
    sqrt(rep(n_vec, mc) / (rep(n_vec, mc) - rep(p_vec, mc))),
  t_theory = rep(t_theory, mc) * sd / sqrt(rep(n_vec, mc) - 1)
)

piv_df <- df |>
  tidyr::pivot_longer(
    cols = c(PaP, LOCO), names_to = "Metric",
    values_to = "linear"
  ) |>
  dplyr::select(Delta, p, n, t_theory, Metric, linear)

piv_df$rf <- df |>
  tidyr::pivot_longer(
    cols = c(PaP_rf, LOCO_rf), names_to = "Metric",
    values_to = "rf"
  ) |>
  pull(rf)

piv_df$theory <- df |>
  tidyr::pivot_longer(
    cols = c(pap_theory, loco_theory), names_to = "metric",
    values_to = "theory"
  ) |>
  pull(theory)

piv_df$theory_fixed <- df |>
  tidyr::pivot_longer(
    cols = c(pap_theory_fix, loco_theory_fix),
    names_to = "metric", values_to = "theory_fixed"
  ) |>
  pull(theory_fixed)

piv_df$Metric <- factor(piv_df$Metric, levels = c("PaP", "LOCO"))

##############################################################
############# PaP and LOCO dependence on Delta ###############

summary_df <- piv_df |>
  dplyr::group_by(Delta, p, n, Metric) |> # maybe Delta
  dplyr::summarise(
    y = mean(rf, na.rm = TRUE),
    sd = sd(rf, na.rm = TRUE),
    ymin = y - 2 * sd,
    ymax = y + 2 * sd,
    theory = unique(theory),
    yc = mean(rf / sqrt(1 + (p - 1) * Delta^2), na.rm = TRUE),
    sdc = sd(rf / sqrt(1 + (p - 1) * Delta^2), na.rm = TRUE),
    yminc = yc - 2 * sdc,
    ymaxc = yc + 2 * sdc,
    theoryc = unique(theory / sqrt(1 + (p - 1) * Delta^2)),
    .groups = "drop"
  )

############### Wide #########################

df_long <- summary_df |>
  dplyr::filter(p %in% c(6, 12), n %in% c(200, 2000)) |>
  tidyr::pivot_longer(
    cols = c(y, yc, sd, sdc, ymin, yminc, ymax, ymaxc, theory, theoryc),
    names_to = c(".value", "type"),
    names_pattern = "^(y|ymin|ymax|sd|theory)(c?)$"
  ) %>%
  dplyr::mutate(type = ifelse(type == "", "Covariance", "Correlation"))

df_long$type <- factor(df_long$type, levels = c("Covariance", "Correlation"))

ggplot(
  df_long |> dplyr::filter(type == "Covariance"),
  aes(x = Delta, y = y, group = Metric)
) +
  ylab("Importance") +
  geom_line(aes(x = Delta, y = theory), linewidth = 0.7, color = "red") +
  geom_point(size = 2, aes(color = Metric)) +
  geom_errorbar(aes(
    ymin = ymin,
    ymax = ymax,
    color = Metric
  ), width = 0.07, linewidth = 0.9) +
  ggeasy::easy_remove_legend_title() +
  ggh4x::facet_grid2(p ~ n,
    labeller = "label_both",
    scales = "free_y"
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      p == 6 ~ scale_y_continuous(limits = c(-0.25, 4.12), breaks = 0:2 * 2),
      p == 12 ~ scale_y_continuous(limits = c(-0.25, 6.82), breaks = 0:3 * 2)
    )
  ) +
  scale_x_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = c("blue", "black")) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

ggsave("../plots/rf_cov_vs_delta.jpg", dpi = 1200, width = 6, height = 4)

############### Tall #########################

df_long <- summary_df |>
  dplyr::filter(p %in% c(6, 12), n %in% c(200, 2000)) |>
  tidyr::pivot_longer(
    cols = c(y, yc, sd, sdc, ymin, yminc, ymax, ymaxc, theory, theoryc),
    names_to = c(".value", "type"),
    names_pattern = "^(y|ymin|ymax|sd|theory)(c?)$"
  ) %>%
  dplyr::mutate(type = ifelse(type == "", "Covariance", "Correlation"))

df_long$type <- factor(df_long$type, levels = c("Covariance", "Correlation"))

ggplot(
  df_long |> dplyr::filter(type == "Covariance"),
  aes(x = Delta, y = y, group = Metric)
) +
  ylab("Importance") +
  geom_line(aes(x = Delta, y = theory), linewidth = 0.7, color = "red") +
  geom_point(size = 2, aes(color = Metric)) +
  geom_errorbar(aes(
    ymin = ymin,
    ymax = ymax,
    color = Metric
  ), width = 0.07, linewidth = 0.9) +
  ggeasy::easy_remove_legend_title() +
  ggh4x::facet_grid2(p ~ n,
    labeller = "label_both",
    scales = "free_y"
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      p == 6 ~ scale_y_continuous(limits = c(-0.25, 4.12), breaks = 0:2 * 2),
      p == 12 ~ scale_y_continuous(limits = c(-0.25, 6.82), breaks = 0:3 * 2)
    )
  ) +
  scale_x_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = c("blue", "black")) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# ggsave("../plots/rf_cov_vs_delta_tall.jpg",
#        dpi = 1200, width = 4.6, height = 4.6)
