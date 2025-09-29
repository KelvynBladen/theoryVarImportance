library(tidyr)
library(tidyverse)
library(GGally)
library(patchwork)
library(latex2exp)

# Takes about 2.5 minutes to run 100 mc simulations.
# Adjusting 'mc' should scale time linearly, while memory is more than linear.
mc <- 100
sd <- 0.01
Delta <- seq(0, 0.99, by = .11)
n_set <- round(2 * 10^(seq(from = 1, to = 3, length = 5)))
p_set <- 1:4 * 3

beta <- mat.or.vec(mc, length(Delta) * length(p_set) * length(n_set))
loco <- mat.or.vec(mc, length(Delta) * length(p_set) * length(n_set))
pap <- mat.or.vec(mc, length(Delta) * length(p_set) * length(n_set))
t_stat <- mat.or.vec(mc, length(Delta) * length(p_set) * length(n_set))

set.seed(123)

s <- Sys.time()
for (m in 1:mc) {
  k <- 1
  for (n in n_set) {
    for (p in p_set) {
      for (i in seq_len(length(Delta))) {
        A <- matrix(Delta[i], p, p)
        diag(A) <- 1
        # A = A / sqrt(1 + (p-1)*Delta[i]^2) # Cor instead of Cov
        sig <- diag(1, p, p)
        Z <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = sig)
        X <- Z %*% A
        # X = scale(X1, center = T, scale = T)
        # scale = T yields t_stat = t_theory, correlation based
        # scale = F yields t_stat = loco_theory, covariance based

        Zv <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = sig)
        Xv <- Zv %*% A
        y <- X %*% rep(1, p) + rnorm(n, sd = sd)
        df <- data.frame(X, y)
        yv <- Xv %*% rep(1, p) + rnorm(n, sd = sd)
        dfv <- data.frame(Xv, yv)

        mod <- lm(y ~ ., data = df)
        t_stat[m, k] <- summary(mod)$coefficients[2, 3]
        pr <- predict(mod, dfv)
        mv <- mean((pr - dfv$yv)^2)

        df_new <- df

        df_new[, 1] <- df_new[sample(1:n), 1]
        reg <- lm(y ~ ., data = df_new)
        beta[m, k] <- mean(reg$coefficients[-c(1:2)])

        prl <- predict(reg, dfv)
        new_mvl <- mean((prl - dfv$yv)^2)
        loco[m, k] <- sqrt(max((new_mvl - mv), 0))

        dfv_new <- dfv

        dfv_new[, 1] <- dfv_new[sample(1:n), 1]
        prp <- predict(mod, dfv_new)
        new_mvp <- mean((prp - dfv$yv)^2)
        pap[m, k] <- sqrt(max((new_mvp - mv), 0))

        k <- k + 1
      }
    }
  }
  print(m)
}
ss <- Sys.time() - s
ss

############## Data ############################################################

d_vec <- rep(Delta, times = length(p_set) * length(n_set))
p_vec <- rep(rep(p_set, each = length(Delta)), times = length(n_set))
n_vec <- rep(n_set, each = length(Delta) * length(p_set))
c_mat <- beta - 1
c_vec <- as.vector(t(c_mat))

c_theory <- (2 * d_vec + (p_vec - 2) * d_vec^2) /
  ((1 - d_vec)^2 + (2 * d_vec + (p_vec - 2) * d_vec^2) * (p_vec - 1))
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
  c_theory = rep(c_theory, mc),
  c = c_vec,
  PaP = as.vector(t(pap)),
  pap_theory = rep(pap_theory, mc),
  pap_theory_fix = rep(pap_theory, mc) *
    (rep(n_vec, mc) - 1) / rep(n_vec, mc),
  LOCO = as.vector(t(loco)),
  loco_theory = rep(loco_theory, mc),
  loco_theory_fix = rep(loco_theory, mc) *
    sqrt(rep(n_vec, mc) / (rep(n_vec, mc) - rep(p_vec, mc))),
  t_stat = as.vector(t(t_stat)) * sd / sqrt(rep(n_vec, mc) - 1),
  t_theory = rep(t_theory, mc) * sd / sqrt(rep(n_vec, mc) - 1),
  t_theory_fixed = (rep(t_theory, mc) * sd / sqrt(rep(n_vec, mc) - 1)) /
    (sqrt((rep(n_vec, mc) - 1) / (rep(n_vec, mc) - rep(p_vec, mc) + 1)))
)

piv_df <- df |>
  tidyr::pivot_longer(
    cols = c(PaP, LOCO),
    names_to = "Metric", values_to = "empirical"
  ) |>
  dplyr::select(Delta, p, n, t_stat, t_theory, Metric, empirical)

piv_df$theory <- df |>
  tidyr::pivot_longer(
    cols = c(pap_theory, loco_theory),
    names_to = "metric", values_to = "theory"
  ) |>
  pull(theory)

piv_df$theory_fixed <- df |>
  tidyr::pivot_longer(
    cols = c(pap_theory_fix, loco_theory_fix),
    names_to = "metric", values_to = "theory_fixed"
  ) |>
  pull(theory_fixed)

piv_df$Metric <- factor(piv_df$Metric, levels = c("PaP", "LOCO"))

######### C ##################################################################

summary_cdf <- df |>
  #filter(n > 20) |>
  dplyr::group_by(Delta, p, n) |>
  dplyr::summarise(
    y = mean((c - c_theory) / c_theory, na.rm = TRUE) * 100,
    sd = sd((c - c_theory) / c_theory, na.rm = TRUE) * 100,
    ymin = y - 2 * sd,
    ymax = y + 2 * sd,
    yc = mean(c, na.rm = TRUE),
    sdc = sd(c, na.rm = TRUE),
    yminc = yc - 2 * sdc,
    ymaxc = yc + 2 * sdc,
    theory = unique(c_theory),
    .groups = "drop"
  )

ggplot(
  summary_cdf |> filter(p %in% c(6, 12), n %in% c(200, 2000)),
  aes(x = Delta, y = yc)
) +
  geom_line(aes(x = Delta, y = theory, group = p),
    color = "red", linewidth = 0.7
  ) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = yminc, ymax = ymaxc),
                width = 0.07, linewidth = 0.9) +
  ylab(latex2exp::TeX("$\\textit{c}$: Empirical Distributions")) +
  scale_x_continuous(limits = c(-0.04, 1.03), breaks = 0:5 / 5) +
  ggh4x::facet_grid2(p ~ n, scales = "free_y", labeller = "label_both") +
  # maybe change to fixed
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(breaks = c(0, 0.1, 0.2)),
      scale_y_continuous(breaks = c(-0.03, 0, 0.03, 0.06, 0.09))
    )
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

ggsave("../plots/c_vs_delta.jpg", dpi = 800, width = 6, height = 4)


ggplot(
  summary_cdf,
  aes(x = Delta, y = yc)
) +
  geom_line(aes(x = Delta, y = theory, group = p),
            color = "red", linewidth = 0.7
  ) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = yminc, ymax = ymaxc),
                width = 0.07, linewidth = 0.9) +
  ylab(latex2exp::TeX("$\\textit{c}$: Empirical Distributions")) +
  scale_x_continuous(limits = c(-0.04, 1.03), breaks = 0:5 / 5) +
  ggh4x::facet_grid2(p ~ n, scales = "free_y", labeller = "label_both") +
  # maybe change to fixed
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6)),
      scale_y_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2)),
      scale_y_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2)),
      scale_y_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2))
    )
  ) +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

ggsave("../plots/c_vs_delta_full1.jpg", dpi = 800, width = 10, height = 9)


######## Combo Plots ###########################################################

summary_df <- piv_df |>
  dplyr::group_by(p, n, Metric) |> # maybe Delta
  dplyr::summarise(
    y = mean((empirical - theory) / theory, na.rm = TRUE) * 100,
    sd = sd((empirical - theory) / theory, na.rm = TRUE) * 100,
    ymin = y - 2 * sd,
    ymax = y + 2 * sd,
    yf = mean((empirical - theory_fixed) / theory_fixed, na.rm = TRUE) * 100,
    sdf = sd((empirical - theory_fixed) / theory_fixed, na.rm = TRUE) * 100,
    yminf = yf - 2 * sdf,
    ymaxf = yf + 2 * sdf,
    .groups = "drop"
  )

ss <- summary_df |> tidyr::pivot_longer(
  cols = c(y, yf), names_to = "Type",
  values_to = "Mean"
)
ss$Type <- rep(c("Raw", "Adjusted"), times = 40)

ss$min <- summary_df |>
  tidyr::pivot_longer(
    cols = c(ymin, yminf),
    names_to = "Type", values_to = "min"
  ) |>
  pull(min)
ss$max <- summary_df |>
  tidyr::pivot_longer(
    cols = c(ymax, ymaxf),
    names_to = "Type", values_to = "max"
  ) |>
  pull(max)

ss <- ss |> dplyr::select(p, n, Metric, Type, Mean, min, max)
ss$Type <- factor(ss$Type, levels = c("Raw", "Adjusted"))


ggplot(ss |> filter(p %in% c(6, 12)), aes(x = n, y = Mean)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = min, ymax = max), width = 0.3, linewidth = 0.9) +
  xlab(latex2exp::TeX("Sample Size ($\\textit{n}$)")) +
  ylab("Importance: % Change Empirical vs. Theory") +
  scale_x_log10(breaks = c(20, 200, 2000)) +
  ggh4x::facet_nested(
    rows = vars(p), cols = vars(Metric, Type),
    labeller = labeller(
      p = label_both,
      Type = label_value,
      Metric = label_value
    ),
    scales = "free_y"
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(-100, 165), breaks = -2:3 * 50),
      scale_y_continuous(limits = c(-100, 165), breaks = -2:3 * 50)
    )
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

ggsave("../plots/combo_parity_w.jpg", dpi = 800, width = 7.2, height = 4.5)

ggplot(ss |> filter(p %in% c(6, 12)), aes(x = n, y = Mean)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = min, ymax = max), width = 0.3, linewidth = 0.9) +
  xlab(latex2exp::TeX("Sample Size ($\\textit{n}$)")) +
  ylab("Importance: % Change Empirical vs. Theory") +
  scale_x_log10(breaks = c(20, 200, 2000)) +
  ggh4x::facet_nested(
    rows = vars(Metric, p), cols = vars(Type),
    labeller = labeller(
      p = label_both,
      Type = label_value,
      Metric = label_value
    ),
    scales = "free_y"
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(-103, 103)),
      scale_y_continuous(limits = c(-103, 103)),
      scale_y_continuous(limits = c(-103, 103)),
      scale_y_continuous(limits = c(-100, 168), breaks = -2:3 * 50)
    )
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

# ggsave("../plots/combo_parity_t.jpg", dpi = 1200, height = 5.8, width = 4.4)


ggplot(ss, aes(x = n, y = Mean)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = min, ymax = max), width = 0.3, linewidth = 0.9) +
  xlab(latex2exp::TeX("Sample Size ($\\textit{n}$)")) +
  ylab("Importance: % Change Empirical vs. Theory") +
  scale_x_log10(breaks = c(20, 200, 2000)) +
  ggh4x::facet_nested(
    rows = vars(p), cols = vars(Metric, Type),
    labeller = labeller(
      p = label_both,
      Type = label_value,
      Metric = label_value
    ),
    scales = "free_y"
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(-100, 165), breaks = -2:3 * 50),
      scale_y_continuous(limits = c(-100, 165), breaks = -2:3 * 50),
      scale_y_continuous(limits = c(-100, 165), breaks = -2:3 * 50),
      scale_y_continuous(limits = c(-100, 165), breaks = -2:3 * 50)
    )
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

ggsave("../plots/combo_parity_w_full.jpg", dpi = 800, width = 7.2,
       height = 8)

######################### LOCO vs T ##########################

ggplot(
  df |> filter(n %in% c(200, 2000), p %in% c(6, 12)),
  aes(x = t_stat, y = LOCO)
) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
  xlab("t-stat") +
  geom_abline(slope = 1, col = "red", linewidth = 2) +
  geom_point(size = 2, alpha = 0.3) +
  ggh4x::facet_grid2(p ~ n, scales = "fixed", labeller = "label_both") +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12)
  )

ggsave("../plots/loco_vs_t.jpg", dpi = 800, width = 5, height = 4)


ggplot(
  df,
  aes(x = t_stat, y = LOCO)
) +
  xlab("t-stat") +
  geom_abline(slope = 1, col = "red", linewidth = 2) +
  geom_point(size = 2, alpha = 0.3) +
  ggh4x::facet_grid2(p ~ n, scales = "fixed", independent = F, labeller = "label_both") +
  scale_x_continuous(limits = c(-0.05, max(df$LOCO)), breaks = 0:3) +
  # ggh4x::facetted_pos_scales(
  #   x = list(
  #     n == 63 ~ scale_x_continuous(limits = c(-0.01, 1.5), breaks = 0:3 / 2),
  #     n == 200 ~ scale_x_continuous(limits = c(-0.01, 1.2), breaks = 0:2 / 2),
  #     n == 632 ~ scale_x_continuous(limits = c(-0.01, 1.1), breaks = 0:2 / 2),
  #     n == 2000 ~ scale_x_continuous(limits = c(-0.01, 1.05), breaks = 0:2 / 2)),
  #   y = list(
  #     p == 3 ~ scale_y_continuous(breaks = 0:2),
  #     p == 6 ~ scale_y_continuous(breaks = 0:2))
  #   ) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 16)
  )

ggsave("../plots/loco_vs_t_full1.jpg", dpi = 800, width = 12, height = 8)

##############################################################
############# PaP and LOCO dependence on Delta ###############

summaryd_df <- piv_df |>
  group_by(Delta, p, n, Metric) |> # maybe Delta
  summarise(
    y = mean(empirical, na.rm = TRUE),
    sd = sd(empirical, na.rm = TRUE),
    ymin = y - 2 * sd,
    ymax = y + 2 * sd,
    theory = unique(theory),
    yc = mean(empirical / sqrt(1 + (p - 1) * Delta^2), na.rm = TRUE),
    sdc = sd(empirical / sqrt(1 + (p - 1) * Delta^2), na.rm = TRUE),
    yminc = yc - 2 * sdc,
    ymaxc = yc + 2 * sdc,
    theoryc = unique(theory / sqrt(1 + (p - 1) * Delta^2)),
    .groups = "drop"
  )

df_long <- summaryd_df |>
  dplyr::filter(p %in% c(6, 12), n %in% c(200, 2000)) |>
  tidyr::pivot_longer(
    cols = c(y, yc, sd, sdc, ymin, yminc, ymax, ymaxc, theory, theoryc),
    names_to = c(".value", "type"),
    names_pattern = "^(y|ymin|ymax|sd|theory)(c?)$"
  ) %>%
  dplyr::mutate(type = ifelse(type == "", "Covariance", "Correlation"))

df_long$type <- factor(df_long$type, levels = c("Covariance", "Correlation"))

df_longf <- summaryd_df |>
  tidyr::pivot_longer(
    cols = c(y, yc, sd, sdc, ymin, yminc, ymax, ymaxc, theory, theoryc),
    names_to = c(".value", "type"),
    names_pattern = "^(y|ymin|ymax|sd|theory)(c?)$"
  ) %>%
  dplyr::mutate(type = ifelse(type == "", "Covariance", "Correlation"))

df_longf$type <- factor(df_longf$type, levels = c("Covariance", "Correlation"))

######## Wide ##########
ggplot(
  df_long |> filter(type == "Covariance"),
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
      p == 6 ~ scale_y_continuous(limits = c(-0.1, 4.1), breaks = 0:2 * 2),
      p == 12 ~ scale_y_continuous(limits = c(-0.1, 6), breaks = 0:3 * 2)
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

ggsave("../plots/cov_vs_delta.jpg", dpi = 800, width = 6, height = 4)


ggplot(
  df_longf |> filter(type == "Covariance"),
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
      p == 3 ~ scale_y_continuous(limits = c(-0.2, 3.7)),
      p == 6 ~ scale_y_continuous(limits = c(-0.2, 6), breaks = 0:3 * 2),
      p == 9 ~ scale_y_continuous(limits = c(-0.2, 6.8), breaks = 0:3 * 2),
      p == 12 ~ scale_y_continuous(limits = c(-0.2, 8), breaks = 0:4 * 2)
    )
  ) +
  scale_x_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = c("blue", "black")) +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 16)
  )

ggsave("../plots/cov_vs_delta_full1.jpg", dpi = 800, width = 10, height = 8)

###### Tall #######

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
      p == 6 ~ scale_y_continuous(limits = c(-0.1, 4.1), breaks = 0:2 * 2),
      p == 12 ~ scale_y_continuous(limits = c(-0.1, 6), breaks = 0:3 * 2)
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

# ggsave("../plots/cov_vs_delta_tall.jpg",
#        dpi = 1200, width = 4.6, height = 4.6)

######## correlation ##########

ggplot(
  df_long |> filter(type == "Correlation"),
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
  ggh4x::facet_grid2(n ~ p,
    labeller = "label_both",
    scales = "free_y"
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

# ggsave("../plots/cor_vs_delta.jpg", dpi = 800, width = 6, height = 4.5)
