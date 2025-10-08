#DL論文
library(dplyr)
library(minpack.lm)

#15種

data15 <- read.csv("data_15sp.csv")
data15

plot(F1 ~ traval, data15)



# ---- ② 非線形モデル（飽和型）----
# 2-1) ミカエリス–メンテン: F1 = A * N / (N + N50)
m_mm <- nlsLM(
  F1 ~ A * traval / (traval + N50),
  data = data15,
  start = list(A = min(0.99, max(data15$F1) * 1.05), N50 = median(data15$traval, na.rm = TRUE)),
  lower = c(0.5, 1), upper = c(1.0, Inf),
  control = nls.lm.control(maxiter = 500)
)
summary(m_mm)

AIC(m_mm)


# ---- ②-2) 指数飽和: F1 = A * (1 - exp(-k * N)) ----

m_exp <- nlsLM(
  F1 ~ A * (1 - exp(-k * traval)),
  data = data15,
  start = list(A = min(0.99, max(data15$F1) * 1.05), k = 0.01),
  lower = c(0.5, 0),      # A ∈ [0.5, 1.0], k >= 0
  upper = c(1.0, Inf),
  control = nls.lm.control(maxiter = 500)
)
summary(m_exp)
AIC(m_exp)

# ---- ②-3) ロジスティック（x=log10N）: F1 = A / (1 + exp(-k * (log10N - c))) ----
data15$log10N <- log10(pmax(data15$traval, 1))  # 念のための安全策

m_logit <- nlsLM(
  F1 ~ A / (1 + exp(-k * (log10N - c))),
  data = data15,
  start = list(A = min(0.99, max(data15$F1) * 1.05), k = 2, c = mean(data15$log10N)),
  lower = c(0.5, 0, -Inf),
  upper = c(1.0, Inf, Inf),
  control = nls.lm.control(maxiter = 500)
)
summary(m_logit)
AIC(m_logit)

#AICが低いミカエリスメンテンでプロット

par(mar=c(4.5,4.5,2.5,1))
plot(
  log10(data15$traval), data15$F1,
  pch = 19, col = "gray30",
  xlab = expression(log[10]*"(Training sample size, N)"),
  ylab = expression("F"[1]*"-score"),
  main = "Michaelis–Menten fit (log-scaled)",
  cex.lab = 1.2, cex.main = 1.2
)
lines(log10(xseq), yhat, col = "blue", lwd = 2)
#abline(v = log10(coef(m_mm)["N50"]), col = "red", lty = 2)


