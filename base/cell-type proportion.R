setwd("Renv")
library(EpiDISH)
library(tidyr)
library(ggplot2)
library(reshape2)
library(dplyr)



load("/mnt/local-disk/data/guoxiaolong/data/BBC/dataBBC.Rd")
load("/mnt/local-disk/data/guoxiaolong/data/BBC/PhenoTypesBBC.Rd")
load("/mnt/local-disk/data/guoxiaolong/data/BBC/annoEPICv1B2.Rd")

# predict fractions of each cell-type using EpiDISH
BloodFrac.m <- epidish(beta.m = bmiqBBC.m, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
save(BloodFrac.m, file = "BFrac.rda")
boxplot(BloodFrac.m)
pheatmap::pheatmap(BloodFrac.m)

### arrange format
# unify col names
facsBBC.m <- data.frame(facsBBC.m)
BloodFrac.m <- data.frame(BloodFrac.m)
colnames(facsBBC.m)[colnames(facsBBC.m) %in% c('Neu', 'Mon', 'Eos/Bas')] <- c('Neutro', 'Mono', 'Eosino')
# unify row names
namelist <- data.frame(PhenoTypesBBC.lv$ChIP, PhenoTypesBBC.lv$Pos, PhenoTypesBBC.lv$PID, PhenoTypesBBC.lv$Visit)
colnames(namelist) <- c('ChIP', 'Pos', 'PID', 'Visit')
namelist <- unite(namelist, "ChIP_Pos", ChIP, Pos) %>% unite("PID_Visit", PID, Visit, sep = "-V")

for (i in 1:nrow(namelist)) {
  rownames(BloodFrac.m)[i] <- namelist$PID_Visit[namelist$ChIP_Pos %in% rownames(BloodFrac.m)[i]]
} # replace
#save(BloodFrac.m, file = "BFrac.rda")
### bar plot
# average
dish_aver <- data.frame(apply(BloodFrac.m, 2, mean))
colnames(dish_aver)[1] <- "EpiDISH"

facs_aver <- data.frame(apply(facsBBC.m, 2, mean, na.rm = TRUE))
colnames(facs_aver)[1] <- "FACS"
# merge
dish_vs_facs <- merge(dish_aver, facs_aver, by = "row.names", all = FALSE)
rownames(dish_vs_facs) <- dish_vs_facs$Row.names
dish_vs_facs <- subset(dish_vs_facs, select = -c(Row.names))
dish_vs_facs$cell_type <- rownames(dish_vs_facs)
d.melt <- melt(dish_vs_facs, id.vars = "cell_type")

# mean RMSE
options(digits = 2)
RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}
mean_rmse <- RMSE(dish_vs_facs$EpiDISH, dish_vs_facs$FACS)

# draw bar plot

ggplot(d.melt, aes(x = cell_type, y = value, fill = variable)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("EpiDISH" = "blue", "FACS" = "#56A0D3")) +
  labs(y = "AvFrac") +
  theme_classic()+
  theme(axis.title.x = element_text(size = 0), legend.title = element_blank())+
  scale_y_continuous(limits = c(0, 0.8), expand = c(0,0), breaks = seq(0, 0.8, by = 0.2))+
  annotate("text", x = max(d.melt$cell_type), y = max(d.melt$value), 
           label = paste0("RMSE = ", format(mean_rmse, nsmall = 3)), 
           hjust = 2, vjust = 0, size = 5, color = "black")

### scatter plot

# melt two sets

facsBBC.m$Sample <- rownames(facsBBC.m)
BloodFrac.m$Sample <- rownames(BloodFrac.m)
## melt the data frames
facsBBC_melt <- melt(facsBBC.m, id.vars = c("Sample"))
BloodFrac_melt <- melt(BloodFrac.m, id.vars = c("Sample"))
## merge the melted dataframes
merged_data <- merge(facsBBC_melt, BloodFrac_melt, by = c("Sample", "variable"))
colnames(merged_data) <- c("Sample", "Cell_type", "facsBBC.m", "BloodFrac.m")
facsBBC.m <- subset(facsBBC.m, select = -c(Sample))
BloodFrac.m <- subset(BloodFrac.m, select = -c(Sample))
# Select a random sample n points for each cell type
merged_data_sample <- merged_data %>% 
  group_by(Cell_type) %>% 
  sample_frac(size = 0.06)

# delete NA.
merged_data_sample_clean <- merged_data_sample %>% filter(!is.na(BloodFrac.m), !is.na(facsBBC.m))

# group the data by cell type
grouped_data <- group_by(merged_data, Cell_type) %>% filter(!is.na(BloodFrac.m), !is.na(facsBBC.m))

# calculate R-squared for each cell type
r_squared_values <- grouped_data %>% 
  do(mod = lm(BloodFrac.m ~ facsBBC.m, data = .)) %>% 
  mutate(R_squared = summary(mod)$r.squared)

AvR2 <- round(mean(r_squared_values$R_squared), digits = 2) # average R2
# draw plot

ggplot(merged_data_sample_clean, aes(x=BloodFrac.m, y=facsBBC.m, color=Cell_type)) +
  geom_point(shape = 2) +
  geom_abline(intercept = 0, slope = 1, color = "black", size = 0.3, linetype = 2) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange", "black", "pink")) + 
  theme_classic() +
  theme(legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        legend.title = element_blank(), )+
  labs(x = "Frac(EpiDISH)", y = "Frac(FACS)")+
  annotate("text",
           label = paste0("R2(", r_squared_values$Cell_type, ") = ", round(r_squared_values$R_squared, digits = 2))
           , y = c(0.4,0.43,0.46,0.49,0.52,0.55,0.58), x = 0.05, color = "black", hjust = 0, vjust = 0)+ 
  annotate("text", label = paste0("AvR2=", AvR2), x= 0.4, y = 0.2, size = 6, color = "black", hjust = 0, vjust = 0)







