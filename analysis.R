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
  dplyr::select(id, week, block, tank, net, method, time, total, total9, boufura, hiru, itotonbo, makigai, meiga, nimaigai, tonbo, yanma, yusurika)%>%
  merge(GT, by = "id", suffixes = c("", "_GT"))
direct

write.csv(direct,"direct_GT.csv")


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
  dplyr::bind_rows(DL)
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



#時間
plot(time ~ log10(total), photo)
plot(time ~ total, direct, xlim = c(0,1200))
max(photo$time, na.rm = T)
max(photo$total, na.rm = T)
hist(photo$time, na.rm = T)

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









