library(tidyr)
library(tidyverse)
library(randomForest)

mc <- 100
sd <- 0.01
Delta <- seq(0, 0.99, by = .11)
n_set <- round(2 * 10^(seq(from = 1, to = 3, length = 5)))
p_set <- 1:4 * 3

beta <- mat.or.vec(length(Delta) * length(p_set) * length(n_set), mc)
loco <- mat.or.vec(length(Delta) * length(p_set) * length(n_set), mc)
pap <- mat.or.vec(length(Delta) * length(p_set) * length(n_set), mc)
loco_rf <- loco
pap_rf <- pap

set.seed(123)

s <- Sys.time()
for (m in 1:mc) {
  k <- 1
  ts <- Sys.time()

  for (n in n_set) {
    for (p in p_set) {
      for (i in seq_len(length(Delta))) {
        A <- matrix(Delta[i], p, p)
        diag(A) <- 1
        sig <- diag(1, p, p)
        Z <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = sig)
        X <- Z %*% A

        Zv <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = sig)
        Xv <- Zv %*% A

        y <- X %*% rep(1, p) + rnorm(n, sd = sd)
        df <- data.frame(X, y)
        yv <- Xv %*% rep(1, p) + rnorm(n, sd = sd)
        dfv <- data.frame(Xv, yv)

        mod <- lm(y ~ ., data = df)
        pr <- predict(mod, dfv)
        mv <- mean((pr - dfv$yv)^2)

        rf <- randomForest(y ~ ., data = df, mtry = max(2, floor(p / 3)))
        prr <- predict(rf, dfv)
        mvr <- mean((prr - dfv$yv)^2)

        df_new <- df

        df_new[, 1] <- df_new[sample(1:n), 1]
        reg <- lm(y ~ ., data = df_new)
        reg$coefficients
        beta[k, m] <- mean(reg$coefficients[-c(1:2)])

        prl <- predict(reg, dfv)
        new_mvl <- mean((prl - dfv$yv)^2)
        loco[k, m] <- sqrt(max((new_mvl - mv), 0))

        rf_new <- randomForest(y ~ ., data = df_new,
                               mtry = max(2, floor(p / 3)))
        prrl <- predict(rf_new, dfv)
        new_mvr <- mean((prrl - dfv$yv)^2)
        loco_rf[k, m] <- sqrt(max((new_mvr - mvr), 0))

        dfv_new <- dfv

        dfv_new[, 1] <- dfv_new[sample(1:n), 1]
        prp <- predict(mod, dfv_new)
        new_mvp <- mean((prp - dfv$yv)^2)
        pap[k, m] <- sqrt(max((new_mvp - mv), 0))

        prpr <- predict(rf, dfv_new)
        new_mvpr <- mean((prpr - dfv$yv)^2)
        pap_rf[k, m] <- sqrt(max((new_mvpr - mvr), 0))

        k <- k + 1
      }
    }
  }
  print(m)
  print(Sys.time() - ts)
}

ss <- Sys.time() - s
ss

beta <- t(beta)
loco <- t(loco)
pap <- t(pap)
loco_rf <- t(loco_rf)
pap_rf <- t(pap_rf)

write.csv(beta, "../data/beta.csv", row.names = FALSE)
write.csv(loco, "../data/loco.csv", row.names = FALSE)
write.csv(pap, "../data/pap.csv", row.names = FALSE)
write.csv(loco_rf, "../data/loco_rf.csv", row.names = FALSE)
write.csv(pap_rf, "../data/pap_rf.csv", row.names = FALSE)
