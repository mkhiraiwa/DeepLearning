#DL論文
library(dplyr)
library(tidyverse)
library(minpack.lm)

library(ggplot2)
library(patchwork) 
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




# ---- モデルに基づく予測値生成 ----
xseq <- seq(min(data15$traval), max(data15$traval), length.out = 200)

# Michaelis–Mentenモデルの推定係数を使って予測
A_mm <- coef(m_mm)["A"]
N50_mm <- coef(m_mm)["N50"]
yhat <- A_mm * xseq / (xseq + N50_mm)

# ---- プロット ----
par(mar=c(4.5,4.5,2.5,1))
plot(
  log10(data15$traval), data15$F1,
  pch = 19, col = "gray30",
  xlab = expression(log[10]*"(Training sample size, N)"),
  ylab = expression("F"[1]*"-score"),
  main = "Michaelis–Menten fit (log-scaled)",
  cex.lab = 1.2, cex.main = 1.2
)

# フィット曲線（logスケールに変換して描く）
lines(log10(xseq), yhat, col = "blue", lwd = 2)

# N50 の位置に縦線を追加（任意）
abline(v = log10(N50_mm), col = "red", lty = 2)



#recall

# 2-1) ミカエリス–メンテン: F1 = A * N / (N + N50)
m_mm <- nlsLM(
  Recall ~ A * traval / (traval + N50),
  data = data15,
  start = list(A = min(0.99, max(data15$F1) * 1.05), N50 = median(data15$traval, na.rm = TRUE)),
  lower = c(0.5, 1), upper = c(1.0, Inf),
  control = nls.lm.control(maxiter = 500)
)
summary(m_mm)

library(ggplot2)

# ---- 予測データフレーム ----
pred_df <- data.frame(
  traval = xseq,
  Recall = yhat
)

# ---- 作図 ----
p <- ggplot(data15, aes(x = log10(traval), y = Recall)) +
  geom_point(
    size = 2.8,
    alpha = 0.75,
    color = "#D55E00"
  ) +
  geom_line(
    data = pred_df,
    aes(x = log10(traval), y = Recall),
    linewidth = 1.4,
    color = "#D55E00"
  ) +
  labs(
    x = expression(log[10] * "(Training sample size, N)"),
    y = "Recall",
    title = "Michaelis–Menten fit"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 17, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

p
ggsave("recall_vs_total.pdf", width = 4, height = 4)



#Precision

# 2-1) ミカエリス–メンテン: F1 = A * N / (N + N50)
m_mm <- nlsLM(
  Precision ~ A * traval / (traval + N50),
  data = data15,
  start = list(A = min(0.99, max(data15$F1) * 1.05), N50 = median(data15$traval, na.rm = TRUE)),
  lower = c(0.5, 1), upper = c(1.0, Inf),
  control = nls.lm.control(maxiter = 500)
)
summary(m_mm)

library(ggplot2)

# ---- 予測データフレーム ----
pred_df <- data.frame(
  traval = xseq,
  Precision = yhat
)

# ---- 作図 ----
p <- ggplot(data15, aes(x = log10(traval), y = Precision)) +
  geom_point(
    size = 2.8,
    alpha = 0.75,
    color = "#D55E00"
  ) +
  geom_line(
    data = pred_df,
    aes(x = log10(traval), y = Precision),
    linewidth = 1.4,
    color = "#D55E00"
  ) +
  coord_cartesian(ylim = c(0, 1))+
  labs(
    x = expression(log[10] * "(Training sample size, N)"),
    y = "Precision",
    title = "Michaelis–Menten fit"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 17, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

p
ggsave("Precision_vs_total.pdf", width = 4, height = 4)




#############################
#　精度比較
###############################
#写真idとタンクの対応
mesopic <- read_csv("meso_pic_20260317.csv", locale = readr::locale(encoding = "CP932"))%>%
  mutate(photo_id = unlist(strsplit(photo_name, ".JPG")), id = paste(block, tank, net, week, sep = "_"), use = 揃ってるデータ, time = アノテーション時間, anno = 自力アノテーション)%>%
  dplyr::filter(use == 1)%>%
  dplyr::select(id, photo_id, week, block, tank, net, use, time, anno)
mesopic

#Ground Truth
GT <- read_csv("annotation_summary.csv", locale = readr::locale(encoding = "CP932"))%>%
  mutate(photo_id = unlist(strsplit(xml_name, ".xml")))%>%
  right_join(mesopic, by = "photo_id")%>%
  mutate(method = "GT", total9 = boufura + hiru + itotonbo + makigai + meiga + nimaigai + tonbo + yanma + yusurika)%>%
  dplyr::select(id, photo_id, total9, boufura, hiru, itotonbo, makigai, meiga, nimaigai, tonbo, yanma, yusurika)
  #dplyr::select(id, photo_id, week, block, tank, net, method, total9, boufura, hiru, itotonbo, makigai, meiga, nimaigai, tonbo, yanma, yusurika)
GT

#目視
direct <- read_csv("manual.csv", locale = readr::locale(encoding = "CP932"))%>% 
  mutate(boufura = ハエ目+ボウフラ+ボウフラさなぎ+ボウフラ緑さなぎ+オニボウフラ, hiru = ヒル+ウズムシ, itotonbo = イトトンボ亜目+イトトンボ科+アオイトトンボ科+キイトトンボ+アオモンイトトンボ+クロイトトンボ, makigai = サカマキガイ+モノアラ+稚貝, meiga = ミズメイガ+メイガ蛹+メイガ蛹葉っぱなし, nimaigai = ドブシジミ, tonbo = トンボ科+アカネ属+シオカラ属+オオシオカラ+ショウジョウ+ハラビロトンボ, yanma = ギンヤンマ属+クロスジギンヤンマ+ギンヤンマ, yusurika = ユスリカ赤+ユスリカ+ユスリカ赤さなぎ+ユスリカ褐色+ユスリカ褐色さなぎ+ユスリカ緑+ユスリカ緑さなぎ)%>%
  mutate(id = paste(block, tank, net, week, sep = "_"), total9 = boufura + hiru + itotonbo + makigai + meiga + nimaigai + tonbo + yanma + yusurika, method = "direct")%>%
  dplyr::select(id, week, block, tank, net, method, time, counter, total, total9, boufura, hiru, itotonbo, makigai, meiga, nimaigai, tonbo, yanma, yusurika)%>%
  merge(GT, by = "id", suffixes = c("", "_GT"))
direct

#write.csv(direct,"direct_GT.csv")


#plot(hiru ~ hiru_GT, direct)


#photo
photo <- read_csv("annotation_summary_photo.csv", locale = readr::locale(encoding = "CP932"))%>%
  mutate(photo_id = unlist(strsplit(xml_name, ".xml")))%>%
  right_join(mesopic, by = "photo_id")%>%
  mutate(method = "photo", total9 = boufura + hiru + itotonbo + makigai + meiga + nimaigai + tonbo + yanma + yusurika)%>%
  mutate(total = total9 + kaiebi + kagerou + matsumomushi + gamushi_youchu + amagaeru + gamushi_seichu + maruhananomi + hanaabu)%>%
  dplyr::select(id, photo_id, week, block, tank, net, method, time, total, total9, boufura, hiru, itotonbo, makigai, meiga, nimaigai, tonbo, yanma, yusurika)%>%
  #merge(GT, by = "photo_id", suffixes = c("", "_GT"))
  left_join(GT, by = "photo_id", suffix = c("", "_GT"))
photo


DLtime <- read_csv("split_time_summary.csv", locale = readr::locale(encoding = "CP932"))%>%
  mutate(photo_id = image_name, time = split_time_sec + predict_time_sec + merge_count_time_sec)%>%
  dplyr::select(photo_id, time)
DLtime

#DL
DL <- read_csv("predict_9sp.csv", locale = readr::locale(encoding = "CP932"))%>%
  mutate(photo_id = unlist(strsplit(file_name, ".jpg")))%>%
  right_join(mesopic, by = "photo_id")%>%
  mutate(method = "DL", total9 = boufura + hiru + itotonbo + makigai + meiga + nimaigai + tonbo + yanma + yusurika)%>%
  mutate(total = total9)%>%
  dplyr::select(id, photo_id, week, block, tank, net, method, total9, boufura, hiru, itotonbo, makigai, meiga, nimaigai, tonbo, yanma, yusurika, total)%>%
  merge(GT, by = "photo_id", suffixes = c("", "_GT"))%>%
  merge(DLtime, by = "photo_id", suffixes = c("", "_GT"))
DL


df <- direct%>%
  dplyr::bind_rows(photo)%>%
  dplyr::bind_rows(DL)%>%
  mutate(method = factor(method, levels = c("direct", "photo","DL")))
df

plot(total9 ~ total9_GT, df, col=c(1:3)[as.factor(df$method)])





species_list <- c("total9", "boufura", "hiru", "itotonbo",
                  "makigai", "meiga", "nimaigai", "tonbo", "yanma", "yusurika")

plots <- list()  # 空リストを作る

for (sp in species_list) {
  axis_max <- max(df[[paste0(sp, "_GT")]], df[[sp]], na.rm = TRUE)
  p <- ggplot(df, aes_string(x = paste0(sp, "_GT"), y = sp)) +
    facet_wrap(~method, nrow = 1) +
    geom_point(alpha = 0.7, color = "black") +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 0.8) +
    coord_fixed(ratio = 1, xlim = c(0, axis_max), ylim = c(0, axis_max)) +
    theme_bw(base_size = 10) +
    labs(title = sp, x = NULL, y = NULL) +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  plots[[sp]] <- p
}

# patchworkで 10行 × 3列に配置
final_plot <- wrap_plots(plots, ncol = 3)
final_plot



#
#
#
library(ggplot2)
library(dplyr)

species_list <- c("total9", "boufura", "hiru", "itotonbo",
                  "makigai", "meiga", "nimaigai", "tonbo", "yanma", "yusurika")

method_cols <- c(
  "DL"     = "#D55E00",
  "photo"  = "#0072B2",
  "direct" = "#009E73"
)

# 作業ディレクトリに保存
outdir <- getwd()

# R2を保存するための入れ物
r2_all <- list()

for (sp in species_list) {
  
  xvar <- paste0(sp, "_GT")
  yvar <- sp
  
  df_sp <- df %>%
    select(method, all_of(xvar), all_of(yvar)) %>%
    filter(!is.na(.data[[xvar]]), !is.na(.data[[yvar]]))
  
  df_sp$method <- factor(df_sp$method, levels = c("direct", "photo", "DL"))
  
  axis_max <- max(df_sp[[xvar]], df_sp[[yvar]], na.rm = TRUE)
  
  # methodごとのR2
  r2_df <- df_sp %>%
    group_by(method) %>%
    summarise(
      r2 = summary(lm(as.formula(paste(yvar, "~", xvar)), data = cur_data()))$r.squared,
      .groups = "drop"
    ) %>%
    mutate(
      species = sp,
      x = axis_max * 0.80,
      y = axis_max * 0.05,
      label = paste0("R^2 == ", sprintf("%.2f", r2))
    )
  
  # 平均R2計算用に保存
  r2_all[[sp]] <- r2_df %>% select(species, method, r2)
  
  p <- ggplot(df_sp, aes(x = .data[[xvar]], y = .data[[yvar]], color = method)) +
    facet_wrap(~method, nrow = 1) +
    geom_point(alpha = 0.8, size = 3) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.5) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", color = "gray30", linewidth = 0.8
    ) +
    geom_text(
      data = r2_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      parse = TRUE,
      size = 5,
      color = "black"
    ) +
    scale_color_manual(values = method_cols, drop = FALSE) +
    coord_fixed(
      ratio = 1,
      xlim = c(0, axis_max),
      ylim = c(0, axis_max),
      expand = FALSE
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  outfile <- file.path(outdir, paste0(sp, ".pdf"))
  
  pdf(outfile, width = 9, height = 3.5)
  print(p)
  dev.off()
}

# ---------------------------
# methodごとの平均R2を計算
# total9は除外
# ---------------------------
r2_summary <- bind_rows(r2_all) %>%
  filter(species != "total9") %>%
  group_by(method) %>%
  summarise(
    mean_r2 = mean(r2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean_r2 = round(mean_r2, 3))

print(r2_summary)

# 必要ならcsv出力
#write.csv(r2_summary, file.path(outdir, "mean_r2_by_method_excluding_total9.csv"), row.names = FALSE)


#R2
library(dplyr)
library(purrr)
library(ggplot2)

species_list <- c("total9", "boufura", "hiru", "itotonbo",
                  "makigai", "meiga", "nimaigai", "tonbo", "yanma", "yusurika")

# 種×methodごとに R2 を計算
r2_df <- map_dfr(species_list, function(sp) {
  
  gt_col  <- paste0(sp, "_GT")
  obs_col <- sp
  
  df %>%
    select(method, all_of(gt_col), all_of(obs_col)) %>%
    rename(GT = all_of(gt_col),
           OBS = all_of(obs_col)) %>%
    filter(!is.na(GT), !is.na(OBS)) %>%
    group_by(method) %>%
    summarise(
      species = sp,
      n = n(),
      R2 = if (n() >= 2 && var(GT) > 0 && var(OBS) > 0) {
        summary(lm(OBS ~ GT))$r.squared
      } else {
        NA_real_
      },
      .groups = "drop"
    )
})

r2_df


ggplot(r2_df, aes(x = method, y = R2, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  theme_bw(base_size = 14) +
  xlab("Method") +
  ylab(expression(R^2)) +
  coord_cartesian(ylim = c(0, 1))


glm_r2 <- glmmTMB(R2 ~ method+(1|species), r2_df, family = "gaussian")
summary(glm_r2)


#種数 NA処理はまだできていない
df$method

df$spno <- apply(ifelse(df[,c(species_list[-1])] > 0, 1, 0), 1, sum) 

paste0(species_list[-1], "_GT")
df$spno_GT <- apply(ifelse(df[,c(paste0(species_list[-1], "_GT"))] > 0, 1, 0), 1, sum) 
p <- ggplot(df, aes_string(x = "spno_GT", y = "spno")) +
  facet_wrap(~method, nrow = 1) +
  geom_point(alpha = 0.7, color = "black") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 0.8) +
  coord_fixed(ratio = 1, xlim = c(0, 9), ylim = c(0, 9)) +
  theme_bw(base_size = 10) +
  labs(title = sp, x = NULL, y = NULL) +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.minor = element_blank()
  )
p



p <- ggplot(df, aes(x = spno_GT, y = spno)) +
  facet_wrap(~method, nrow = 1) +
  geom_count(alpha = 0.7, color = "black") +  # ←ここが変更点
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 0.8) +
  coord_fixed(ratio = 1, xlim = c(0, 7), ylim = c(0, 7)) +
  theme_bw(base_size = 10) +
  labs(title = "sp_no", x = NULL, y = NULL) +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.minor = element_blank()
  )
p


#誤差～個体数

df$dif_total9 <- (df$total9-df$total9_GT)#/df$total9_GT

p <- ggplot(df, aes(x = log10(total9_GT), y = dif_total9)) +
  facet_wrap(~method, nrow = 1) +
  geom_count(alpha = 0.7, color = "black") +  # ←ここが変更点
  geom_smooth(method = NULL, color = "blue", se = FALSE) #+
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 0.8) #+
  #coord_fixed(ratio = 1, xlim = c(0, max(df$total9_GT)), ylim = c(min(df$dif_total9), max(df$dif_total9))) +
  #theme_bw(base_size = 10) +
  #labs(title = "sp_no", x = NULL, y = NULL) +
  #theme(
  #  plot.title = element_text(size = 10, hjust = 0.5),
  #  panel.grid.minor = element_blank()
  #)
p



# 種リスト
species_list <- c("total9", "boufura", "hiru", "itotonbo",
                  "makigai", "meiga", "nimaigai", "tonbo", "yanma", "yusurika")

# 散布図
axis_max <- max(df$total9_GT, df$total9, na.rm = TRUE)
p <- ggplot(df, aes(x = total9_GT, y = total9)) +
  facet_wrap(~method)+
  geom_point(alpha = 0.7, color = "black") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 0.8) +  # 1:1線
  coord_fixed(ratio = 1, xlim = c(0, axis_max), ylim = c(0, axis_max)) +  # 軸の比率を1:1に固定
  theme_bw(base_size = 12) +
  #labs(
  #  x = NULL, y = NULL,
  #  title = paste0(sp, "\n(r = ", sprintf("%.2f", r), ")")
  #) +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    panel.grid.minor = element_blank()
  )
p

p <- ggplot(df, aes(x = log10(total9_GT), y = log10(total9), color = method)) +
  #facet_wrap(~method)+
  geom_point(alpha = 0.7)+#, color = "black") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 0.8) +  # 1:1線
  #coord_fixed(ratio = 1, xlim = c(0, axis_max), ylim = c(0, axis_max)) +  # 軸の比率を1:1に固定
  theme_bw(base_size = 12) +
  #labs(
  #  x = NULL, y = NULL,
  #  title = paste0(sp, "\n(r = ", sprintf("%.2f", r), ")")
  #) +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    panel.grid.minor = element_blank()
  )
p


#多様度
# install.packages("vegan")  # 未インストールなら
library(vegan)
library(dplyr)

# 9種の列名
sp_cols <- c("boufura", "hiru", "itotonbo", "makigai", "meiga",
             "nimaigai", "tonbo", "yanma", "yusurika")

# GT 側の9種列名
gt_cols <- paste0(sp_cols, "_GT")

# Shannon と Simpson を計算
library(vegan)
library(dplyr)

df2 <- df %>%
  mutate(
    # 多様度指数
    shannon_obs = diversity(select(., all_of(sp_cols)), index = "shannon"),
    simpson_obs = diversity(select(., all_of(sp_cols)), index = "simpson"),
    shannon_GT  = diversity(select(., all_of(gt_cols)), index = "shannon"),
    simpson_GT  = diversity(select(., all_of(gt_cols)), index = "simpson"),
    
    # 種数（richness）
    richness_obs = specnumber(select(., all_of(sp_cols))),
    richness_GT  = specnumber(select(., all_of(gt_cols)))
  )

head(df2[, c("id", "method",
             "richness_obs", "richness_GT",
             "shannon_obs", "shannon_GT",
             "simpson_obs", "simpson_GT")])



library(vegan)
library(dplyr)

df2 <- df %>%
  rowwise() %>%
  mutate(
    # NAチェック
    has_na_obs = any(is.na(c_across(all_of(sp_cols)))),
    has_na_GT  = any(is.na(c_across(all_of(gt_cols)))),
    
    # 多様度指数
    shannon_obs = if (has_na_obs) NA_real_ else diversity(c_across(all_of(sp_cols)), index = "shannon"),
    simpson_obs = if (has_na_obs) NA_real_ else diversity(c_across(all_of(sp_cols)), index = "simpson"),
    shannon_GT  = if (has_na_GT)  NA_real_ else diversity(c_across(all_of(gt_cols)), index = "shannon"),
    simpson_GT  = if (has_na_GT)  NA_real_ else diversity(c_across(all_of(gt_cols)), index = "simpson"),
    
    # 種数
    richness_obs = if (has_na_obs) NA_real_ else specnumber(c_across(all_of(sp_cols))),
    richness_GT  = if (has_na_GT)  NA_real_ else specnumber(c_across(all_of(gt_cols)))
  ) %>%
  ungroup()




library(ggplot2)

ggplot(df2, aes(x = shannon_GT, y = shannon_obs)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() +
  xlab("GT Shannon") +
  ylab("Observed Shannon")

ggplot(df2, aes(x = simpson_GT, y = simpson_obs)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() +
  xlab("GT Simpson") +
  ylab("Observed Simpson")


ggplot(df2, aes(x = shannon_GT, y = shannon_obs, color = method)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() +
  xlab("GT Shannon") +
  ylab("Observed Shannon")

ggplot(df2, aes(x = simpson_GT, y = simpson_obs, color = method)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() +
  xlab("GT Simpson") +
  ylab("Observed Simpson")


ggplot(df2, aes(x = shannon_GT, y = shannon_obs)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~ method) +
  theme_bw() +
  xlab("GT Shannon") +
  ylab("Observed Shannon")

ggplot(df2, aes(x = simpson_GT, y = simpson_obs)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~ method) +
  theme_bw() +
  xlab("GT Simpson") +
  ylab("Observed Simpson")

ggplot(df2, aes(x = richness_GT, y = richness_obs)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~ method) +
  theme_bw() +
  xlab("GT richness") +
  ylab("Observed richness")



ggplot(df2, aes(x = richness_obs, y = simpson_obs)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~ method) +
  theme_bw() +
  xlab("GT Simpson") +
  ylab("Observed Simpson")


#時間
plot(time ~ log10(total), photo)
plot(time ~ total, direct, xlim = c(0,1200))
max(photo$time, na.rm = T)
max(photo$total, na.rm = T)
hist(photo$time, na.rm = T)


direct$log10_total <- log10(direct$total)
photo$log10_total  <- log10(photo$total)
DL$log10_total     <- log10(DL$total)
df$log10_total     <- log10(df$total)

library(glmmTMB)
glm_time_direct <- glmmTMB(log(time) ~ log10_total + (1|counter) + (1|week), direct, family = gaussian)
summary(glm_time_direct)

glm_time_photo <- glmmTMB(log(time) ~ log10_total + (1|week), photo, family = gaussian)
summary(glm_time_photo)

glm_time_DL <- glmmTMB(log(time) ~ log10_total + (1|week), DL, family = gaussian)
summary(glm_time_DL)


library(ggeffects)
library(ggplot2)
library(dplyr)

pred_direct <- ggpredict(glm_time_direct, terms = "log10_total")
pred_direct$method <- "direct"

pred_photo <- ggpredict(glm_time_photo, terms = "log10_total")
pred_photo$method <- "photo"

pred_DL <- ggpredict(glm_time_DL, terms = "log10_total")
pred_DL$method <- "DL"

pred_all <- bind_rows(pred_direct, pred_photo, pred_DL)
pred_all$method <- factor(pred_all$method, levels = c("direct", "photo","DL"))

method_cols <- c(
  "DL"     = "#D55E00",
  "photo"  = "#0072B2",
  "direct" = "#009E73"
)

ggplot() +
  geom_point(
    data = df,
    aes(x = log10_total, y = time, color = method),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = pred_all,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = method),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(
    data = pred_all,
    aes(x = x, y = predicted, color = method),
    linewidth = 1.2
  ) +
  
  # ← ここが重要
  scale_color_manual(values = method_cols) +
  scale_fill_manual(values = method_cols) +
  
  theme_bw() +
  xlab("log10 total individuals") +
  ylab("Time (s)") +
  labs(color = "Method", fill = "Method") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )
ggsave("figure_time_vs_total.pdf", width = 6, height = 4)

tapply(df$time, df$method, mean, na.rm = T)


library(ggplot2)

ggplot(direct, aes(x = log10(total), y = time, color = counter)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  xlab("log10 total individuals") +
  ylab("Time") +
  labs(color = "Counter")


library(ggplot2)

ggplot(df, aes(x = log10(total), y = time, color = method)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "glm",
              method.args = list(family = Gamma(link = "log")),
              se = F) +
  theme_bw() +
  labs(x = "Total", y = "Time", color = "Method")


ggplot(df, aes(x = total, y = time, color = method)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm",se = F) +
  theme_bw() +
  labs(x = "Total", y = "Time", color = "Method")


ggplot(df, aes(x = log10(total), y = log10(time), color = method)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm",
              se = FALSE) +
  theme_bw() +
  labs(x = "Total", y = "Time", color = "Method")



ggplot(df, aes(x = total, y = time, color = method)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_bw() +
  labs(x = "Total", y = "Time", color = "Method")



#検出率　サイズ

library(dplyr)
library(tidyr)
library(stringr)

# GT列の名前から種名を取得
gt_species <- c("boufura", "hiru", "itotonbo",
                "makigai", "meiga", "nimaigai", "tonbo", "yanma", "yusurika")

# 必要列だけ取り出してlong形式へ
det_rate_df <- df %>%
  select(method, photo_id, all_of(gt_species), all_of(paste0(gt_species, "_GT"))) %>%
  pivot_longer(
    cols = -c(method, photo_id),
    names_to = "name",
    values_to = "count"
  ) %>%
  mutate(
    species = str_remove(name, "_GT$"),
    type = if_else(str_detect(name, "_GT$"), "GT", "det")
  ) %>%
  select(-name) %>%
  pivot_wider(
    names_from = type,
    values_from = count
  ) %>%
  mutate(
    detection_rate = ifelse(GT == 0, NA, det / GT)
  ) %>%
  select(method, photo_id, species, detection_rate)%>%
  mutate(image_name = paste0(photo_id, ".jpg"))

det_rate_df


df_size <- read_csv("bbox_size_summary_by_image_species.csv", locale = readr::locale(encoding = "CP932"))


merged_df <- det_rate_df %>%
  left_join(df_size, by = c("image_name", "species"))


library(ggplot2)
library(dplyr)

plot_df <- merged_df %>%
  filter(!is.na(detection_rate), !is.na(mean_area_px))

ggplot(plot_df, aes(x = log10(mean_area_px), y = detection_rate)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ method) +
  theme_bw() +
  xlab("Mean bbox area (px)") +
  ylab("Detection rate")






#photoとDLだけ比較


df_size <- read_csv("bbox_level_detection_01_common_only.csv", locale = readr::locale(encoding = "CP932"))
df_size
df_yosoku <- read_csv("yosoku15.csv", locale = readr::locale(encoding = "CP932"))%>%
  mutate(image_id = photo_id)%>%
  dplyr::select(image_id,noise)
df_yosoku

df_size <- df_size%>%
  left_join(df_yosoku, by = "image_id")
df_size


library(ggplot2)
library(readr)
library(dplyr)

#df_size <- read_csv("D:/deeplearning_meso/bbox_level_detection_01_common_only.csv")

# NA除去（念のため）
plot_df <- df_size %>%
  filter(!is.na(log10_gt_area_px), !is.na(detected))


method_cols <- c(
  "DL"     = "#D55E00",
  "photo"  = "#0072B2",
  "direct" = "#009E73",
  "DLphoto" = "black"
)

ggplot(plot_df, aes(x = log10_gt_area_px, y = detected,
                    color = method, fill = method)) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = TRUE,
              linewidth = 1.2) +
  scale_color_manual(values = method_cols) +
  scale_fill_manual(values = method_cols) +
  theme_bw() +
  xlab("log10 bbox area (px)") +
  ylab("Detection probability") +
  labs(color = "Method", fill = "Method")+
  theme(
    axis.text = element_text(size = 14),   # ← 目盛り
  )
ggsave("figure_detect_vs_size.pdf", width = 6, height = 4)

plot_df$method <- factor(plot_df$method, levels = c("photo", "DL", "DLphoto"))

glm_detect <- glm(detected ~ log10_gt_area_px * method + noise * method, family = binomial, plot_df)
summary(glm_detect)


glm_detect <- glmmTMB(detected ~ log10_gt_area_px * method + noise * method + (1|image_id) + (1|species), family = binomial, plot_df)
summary(glm_detect)


#######kakunin
library(readr)
library(dplyr)
library(glmmTMB)
library(ggplot2)
library(ggeffects)

df_size <- read_csv("bbox_level_detection_01_common_only.csv",
                    locale = readr::locale(encoding = "CP932"))

df_yosoku <- read_csv("yosoku15.csv",
                      locale = readr::locale(encoding = "CP932")) %>%
  mutate(image_id = photo_id) %>%
  dplyr::select(image_id, noise)

df_size <- df_size %>%
  left_join(df_yosoku, by = "image_id")

df_size$method <- factor(df_size$method, levels = c("photo", "DL", "DLphoto"))

method_cols <- c(
  "DL"      = "#D55E00",
  "photo"   = "#0072B2",
  "direct"  = "#009E73",
  "DLphoto" = "black"
)

glm_detect <- glmmTMB(
  detected ~ log10_gt_area_px * method + noise * method + (1|image_id) + (1|species),
  family = binomial,
  data = df_size
)

summary(glm_detect)

pred_size <- ggpredict(
  glm_detect,
  terms = c("log10_gt_area_px", "method"),
  condition = c(noise = mean(df_size$noise, na.rm = TRUE))
)

ggplot(pred_size, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = method_cols) +
  scale_fill_manual(values = method_cols) +
  theme_bw() +
  xlab("log10 bbox area (px)") +
  ylab("Detection probability") +
  labs(color = "Method", fill = "Method") +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

ggsave("figure_detect_vs_size_ggpredict.pdf", width = 6, height = 4)




#######Fig 5 作図!!!!!
library(tidyverse)
library(glmmTMB)
library(ggplot2)
library(dplyr)

target_species <- c(
  "tonbo", "makigai", "yanma", "yusurika",
  "boufura", "itotonbo", "hiru", "meiga", "nimaigai"
)

method_cols <- c(
  "direct"  = "#009E73",  # 緑（人間・現場）
  "photo"   = "#0072B2",  # 青（人間・画像）
  "DL"      = "#D55E00",  # オレンジ（機械）
  "DLphoto" = "#E69F00"   # 明るいオレンジ（ハイブリッド）
)

df_yosoku <- read_csv("yosoku15.csv", locale = readr::locale(encoding = "CP932"))%>%
  mutate(image_id = photo_id)%>%
  dplyr::select(image_id,noise)
df_yosoku


df_recall <- read_csv("bbox_level_detection_01_common_only.csv", locale = readr::locale(encoding = "CP932")) %>%
  filter(species %in% target_species)%>%
  left_join(df_yosoku, by = "image_id")%>%
  filter(!is.na(log10_gt_area_px), !is.na(detected))%>%
  mutate(method = factor(plot_df$method, levels = c("photo", "DL", "DLphoto")))
df_recall


tapply(df_size$detected, df_size$method, sum) / tapply(df_size$detected, df_size$method, length)

glm_recall_species <- glmmTMB(
  detected ~ log10_gt_area_px * method + noise * method + (1|image_id) + (1|species),
  #detected ~ log10_gt_area_px * method + noise * method + (1|image_id),
  family = binomial,
  data = df_recall
)
summary(glm_recall_species)


glm_recall <- glmmTMB(
  #detected ~ log10_gt_area_px * method + noise * method + (1|image_id) + (1|species),
  detected ~ log10_gt_area_px * method + noise * method + (1|image_id),
  family = binomial,
  data = df_recall
)

summary(glm_recall)

# -------------------------------
# 予測用データ
# -------------------------------
newdat <- expand.grid(
  log10_gt_area_px = seq(min(df_recall$log10_gt_area_px, na.rm = TRUE),
                         max(df_recall$log10_gt_area_px, na.rm = TRUE),
                         length.out = 200),
  method = levels(df_recall$method),
  noise = mean(df_yosoku$noise, na.rm = TRUE)
)
newdat
# glmmTMBから予測値（母平均）
pred <- predict(glm_recall,
                newdata = newdat,
                type = "link",
                se.fit = TRUE,
                re.form = NA)

newdat$fit_link <- pred$fit
newdat$se_link  <- pred$se.fit

# ロジット→確率に変換
newdat$fit <- plogis(newdat$fit_link)
newdat$lwr <- plogis(newdat$fit_link - 1.96 * newdat$se_link)
newdat$upr <- plogis(newdat$fit_link + 1.96 * newdat$se_link)

# -------------------------------
# 作図
# -------------------------------
f5_1 <- ggplot(newdat, aes(x = log10_gt_area_px, y = fit,
                   color = method, fill = method)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = method_cols) +
  scale_fill_manual(values = method_cols) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  theme_bw() +
  xlab("body size (log10 bbox area (px))") +
  ylab("Recall") +
  labs(color = "Method", fill = "Method") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
f5_1

#ggsave("figure_detect_vs_size_glmmTMB.pdf", width = 6, height = 4)


#noise
# -------------------------------
# 予測用データ（noiseを変化させる）
# -------------------------------
newdat_noise <- expand.grid(
  noise = seq(min(plot_df$noise, na.rm = TRUE),
              max(plot_df$noise, na.rm = TRUE),
              length.out = 200),
  method = levels(df_recall$method),
  log10_gt_area_px = mean(plot_df$log10_gt_area_px, na.rm = TRUE)
)

# 予測
pred_noise <- predict(glm_recall,
                      newdata = newdat_noise,
                      type = "link",
                      se.fit = TRUE,
                      re.form = NA)

newdat_noise$fit_link <- pred_noise$fit
newdat_noise$se_link  <- pred_noise$se.fit

# ロジット→確率
newdat_noise$fit <- plogis(newdat_noise$fit_link)
newdat_noise$lwr <- plogis(newdat_noise$fit_link - 1.96 * newdat_noise$se_link)
newdat_noise$upr <- plogis(newdat_noise$fit_link + 1.96 * newdat_noise$se_link)

# -------------------------------
# 作図
# -------------------------------
f5_2 <- ggplot(newdat_noise, aes(x = noise, y = fit,
                         color = method, fill = method)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = method_cols) +
  scale_fill_manual(values = method_cols) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  theme_bw() +
  xlab("Noise") +
  ylab("Recall") +
  labs(color = "Method", fill = "Method") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
f5_2
#ggsave("figure_detect_vs_noise_glmmTMB.pdf", width = 6, height = 4)



#Presicion

# -------------------------------
# データ読み込み
# -------------------------------
df_precision <- read_csv("bbox_level_precision_01_common_only.csv", locale = readr::locale(encoding = "CP932"))%>%
  filter(species %in% target_species)%>%
  left_join(df_yosoku, by = "image_id")%>%
  filter(!is.na(log10_pred_area_px), !is.na(correct))%>%
  mutate(method = factor(df_precision$method, levels = c("photo", "DL", "DLphoto")))
df_precision

# -------------------------------
# モデル
# -------------------------------
glm_precision_species <- glmmTMB(
  correct ~ log10_pred_area_px * method + noise * method + (1|image_id) + (1|species),
  family = binomial,
  data = df_precision
)

summary(glm_precision)


glm_precision <- glmmTMB(
  correct ~ log10_pred_area_px * method + noise * method + (1|image_id),
  family = binomial,
  data = df_precision
)

summary(glm_precision)

# -------------------------------
# サイズを横軸にした予測図
# noiseは平均値に固定
# -------------------------------
newdat_precision_size <- expand.grid(
  log10_pred_area_px = seq(
    min(df_precision$log10_pred_area_px, na.rm = TRUE),
    max(df_precision$log10_pred_area_px, na.rm = TRUE),
    length.out = 200
  ),
  method = levels(df_precision$method),
  noise = mean(df_precision$noise, na.rm = TRUE)
)

pred_precision_size <- predict(
  glm_precision,
  newdata = newdat_precision_size,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

newdat_precision_size$fit_link <- pred_precision_size$fit
newdat_precision_size$se_link  <- pred_precision_size$se.fit

newdat_precision_size$fit <- plogis(newdat_precision_size$fit_link)
newdat_precision_size$lwr <- plogis(newdat_precision_size$fit_link - 1.96 * newdat_precision_size$se_link)
newdat_precision_size$upr <- plogis(newdat_precision_size$fit_link + 1.96 * newdat_precision_size$se_link)

f5_3 <- ggplot(newdat_precision_size, aes(x = log10_pred_area_px, y = fit,
                                  color = method, fill = method)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = method_cols) +
  scale_fill_manual(values = method_cols) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  theme_bw() +
  xlab("log10 bbox area (px)") +
  ylab("Precision") +
  labs(color = "Method", fill = "Method") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
f5_3
#ggsave("figure_precision_vs_size_glmmTMB.pdf", width = 6, height = 4)

# -------------------------------
# noiseを横軸にした予測図
# bbox sizeは平均値に固定
# -------------------------------
newdat_precision_noise <- expand.grid(
  noise = seq(
    min(df_precision$noise, na.rm = TRUE),
    max(df_precision$noise, na.rm = TRUE),
    length.out = 200
  ),
  method = levels(df_precision$method),
  log10_pred_area_px = mean(df_precision$log10_pred_area_px, na.rm = TRUE)
)

pred_precision_noise <- predict(
  glm_precision,
  newdata = newdat_precision_noise,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

newdat_precision_noise$fit_link <- pred_precision_noise$fit
newdat_precision_noise$se_link  <- pred_precision_noise$se.fit

newdat_precision_noise$fit <- plogis(newdat_precision_noise$fit_link)
newdat_precision_noise$lwr <- plogis(newdat_precision_noise$fit_link - 1.96 * newdat_precision_noise$se_link)
newdat_precision_noise$upr <- plogis(newdat_precision_noise$fit_link + 1.96 * newdat_precision_noise$se_link)

f5_4 <- ggplot(newdat_precision_noise, aes(x = noise, y = fit,
                                   color = method, fill = method)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = method_cols) +
  scale_fill_manual(values = method_cols) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  theme_bw() +
  xlab("Noise") +
  ylab("Precision") +
  labs(color = "Method", fill = "Method") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
f5_4


#ggsave("figure_precision_vs_noise_glmmTMB.pdf", width = 6, height = 4)
library(patchwork)
(f5_1 | f5_2) / (f5_3 | f5_4) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.justification = "top"
  )
ggsave("figure_5.pdf", width=8, height=6)



###############
#############
#############




plot(detected ~ log10(gt_area_px), df_size)


plot(detected ~ log10(gt_area_px), df_size[df_size$method == "photo",])
plot(detected ~ log10(gt_area_px), df_size[df_size$method == "DL",])


plots <- list()
cor_results <- data.frame(species = character(), cor_value = numeric(), stringsAsFactors = FALSE)

# 各種ループ
for (sp in species_list) {
  df_wide <- df %>%
    group_by(id, method) %>%
    summarise(value = mean(.data[[sp]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = method, values_from = value)
  
  # 相関係数
  r <- cor(df_wide$GT, df_wide$direct, use = "pairwise.complete.obs")
  cor_results <- rbind(cor_results, data.frame(species = sp, cor_value = r))
  
  # 軸の最大値を共通に（GTとdirectの最大値のうち大きい方）
  axis_max <- max(df_wide$GT, df_wide$direct, na.rm = TRUE)
  
  # 散布図
  p <- ggplot(df_wide, aes(x = GT, y = direct)) +
    geom_point(alpha = 0.7, color = "black") +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 0.8) +  # 1:1線
    coord_fixed(ratio = 1, xlim = c(0, axis_max), ylim = c(0, axis_max)) +  # 軸の比率を1:1に固定
    theme_bw(base_size = 12) +
    labs(
      x = NULL, y = NULL,
      title = paste0(sp, "\n(r = ", sprintf("%.2f", r), ")")
    ) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  plots[[sp]] <- p
}






# df はすでに読み込まれていると仮定
df_wide <- df %>%
  group_by(id, method) %>%
  summarise(total9 = mean(total9, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = total9)
df_wide
cor(df_wide$GT, df_wide$direct, use = "pairwise.complete.obs")

# df はすでに読み込まれていると仮定
df_wide <- df %>%
  group_by(id, method) %>%
  summarise(yusurika = mean(yusurika, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = yusurika)
df_wide
cor(df_wide$GT, df_wide$direct, use = "pairwise.complete.obs")

ggplot(df_wide, aes(x = GT, y = direct)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", color = "blue") +
  theme_bw() +
  labs(x = "Ground truth (total9)", 
       y = "Direct count (total9)",
       title = "Correlation between direct and ground truth")



library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(patchwork)

# 種リスト
species_list <- c("total9", "boufura", "hiru", "itotonbo",
                  "makigai", "meiga", "nimaigai", "tonbo", "yanma", "yusurika")

plots <- list()
cor_results <- data.frame(species = character(), cor_value = numeric(), stringsAsFactors = FALSE)

# 各種ループ
for (sp in species_list) {
  df_wide <- df %>%
    group_by(id, method) %>%
    summarise(value = mean(.data[[sp]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = method, values_from = value)
  
  # 相関係数
  r <- cor(df_wide$GT, df_wide$direct, use = "pairwise.complete.obs")
  cor_results <- rbind(cor_results, data.frame(species = sp, cor_value = r))
  
  # 軸の最大値を共通に（GTとdirectの最大値のうち大きい方）
  axis_max <- max(df_wide$GT, df_wide$direct, na.rm = TRUE)
  
  # 散布図
  p <- ggplot(df_wide, aes(x = GT, y = direct)) +
    geom_point(alpha = 0.7, color = "black") +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 0.8) +  # 1:1線
    coord_fixed(ratio = 1, xlim = c(0, axis_max), ylim = c(0, axis_max)) +  # 軸の比率を1:1に固定
    theme_bw(base_size = 12) +
    labs(
      x = NULL, y = NULL,
      title = paste0(sp, "\n(r = ", sprintf("%.2f", r), ")")
    ) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  plots[[sp]] <- p
}

# ---- 2×5配置 ----
combined_plot <- wrap_plots(plots, ncol = 5, nrow = 2) +
  plot_annotation(
    title = "Direct Count vs Ground Truth (1:1 reference line in red)",
    subtitle = "Dashed red line = 1:1 relationship, Blue = linear regression fit",
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5))
  )

# 表示
combined_plot

# 保存（任意）
ggsave("D:/deeplearning_meso/GT_direct_correlation_2x5_1to1.png",
       combined_plot, width = 12, height = 5, dpi = 300)

# 相関結果
cor_results







photoGT <- photo %>%
  left_join(GT, by = "photo_id")

photoGT 


ggplot(photoGT, aes(x = total9.y, y = total9.x)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", color = "blue") +
  theme_bw() +
  labs(x = "Ground truth (total9)", 
       y = "Photo count (total9)",
       title = "Correlation between direct and ground truth")




DLGT <- DL %>%
  left_join(GT, by = "photo_id")

DLGT 


ggplot(DLGT, aes(x = total9.y, y = total9.x)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", color = "blue") +
  theme_bw() +
  labs(x = "Ground truth (total9)", 
       y = "Photo count (total9)",
       title = "Correlation between direct and ground truth")





photodirect <- photo %>%
  left_join(GT, by = "photo_id")

photoGT 


ggplot(photodirect, aes(x = total9.y, y = total9.x)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", color = "blue") +
  theme_bw() +
  labs(x = "Ground truth (total9)", 
       y = "Photo count (total9)",
       title = "Correlation between direct and ground truth")









