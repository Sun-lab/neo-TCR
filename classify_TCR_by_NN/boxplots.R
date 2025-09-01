library(ggplot2)
library(stringr)


df = read.csv("results/ensemble_auc.csv", header=TRUE)
dim(df)

df$n_CNN = as.character(df$n_CNN)


df$V = rep("with V", nrow(df))
df$V[which(df$feature%in%c("alpha_cdr3", "beta_cdr3", "both_cdr3"))] = "without V"
df$V = factor(df$V, levels=c("without V", "with V"))

df$chain = str_extract(df$feature, "[^_]+")
df$layer_sizes = factor(df$n_units, levels=c("32_16", "64_32", "128_64", "256_128"))

df

pdf(file = "./figures/auc_chain_v_boxplot.pdf", width = 4, height = 2.5)
ggplot(df, aes(x=chain, y=auc, fill=V)) + 
  geom_boxplot() +
  geom_point(position=position_jitterdodge(), size=0.1) + 
  theme_classic()
dev.off()

pdf(file = "./figures/auc_chain_v_stripchart.pdf", width = 4, height = 2.5)
ggplot(df, aes(x=chain, y=auc, color=V)) +
  geom_point(position = position_jitterdodge(), size=1) + theme_classic()
dev.off()

pdf(file = "./figures/auc_chain_v_encoding_stripchart.pdf", width = 5, height = 2.6)
ggplot(df, aes(x=chain, y=auc, color=V, shape=encoding)) +
  geom_jitter(position=position_dodge(0.6)) + theme_classic()
dev.off()

pdf(file = "./figures/auc_chain_v_layer_sizes_stripchart.pdf", width = 5, height = 2.6)
ggplot(df, aes(x=chain, y=auc, color=V, shape=layer_sizes)) +
  geom_jitter(position=position_dodge(0.8)) + theme_classic()
dev.off()

pdf(file = "./figures/auc_chain_v_n_cnn_stripchart.pdf", width = 4, height = 2.6)
ggplot(df, aes(x=chain, y=auc, color=V, shape=n_CNN)) +
  geom_jitter(position=position_dodge(0.8)) + theme_classic()
dev.off()



sessionInfo()

q(save="no")
