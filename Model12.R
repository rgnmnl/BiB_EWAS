#########################################################################
## Title: Model 12 EWAS
## Version: 2 
## Note: Rerun to exclude duplicate samples identified in Feb 2020.
##		See EWAS_analysis_planv5annotated
## Author: Regina Manansala
## Date: 21-February-2020
#########################################################################

library(car)
library(data.table)
library(dplyr)
library(tidyverse)
library(qvalue)

load("M_val_IDsub.RData")
load("CellCountsBloodgse35069Complete.Rdata")
ph <- read.csv("bibgwas_16JAN2020.csv") %>% subset(., !(sentrix %in% c("201057150183_R01C01", "201096090153_R05C01", "200992330056_R05C01", "201046290095_R08C01")))
#pcs <- fread("PC1-10.txt")

#### Make sure the methylation and phenotype files are sorted exactly the same way.
counts <- data.frame(counts)
counts$sentrix <- rownames(counts)
ph2 <- merge(ph, counts) %>% 
  mutate(., R = as.factor(R), alpreg = as.factor(alpreg), bn_mtben = as.factor(bn_mtben), edufa = as.factor(edufa), 
  edufa2 = as.factor(edufa2), edufa4c = as.factor(edufa4c), eduma = as.factor(eduma), eduma2 = as.factor(eduma2), 
  eduma_di = as.factor(eduma_di), eduma4c = as.factor(eduma4c), empstma = as.factor(empstma), empma = as.factor(empma), 
  Eth9 = as.factor(Eth9), gender = as.factor(gender), MatSD = as.factor(MatSD), Marstat = as.factor(Marstat), 
  quin07 = as.factor(quin07), quin07di = as.factor(quin07di), quin073c = as.factor(quin073c), smkpreg = as.factor(smkpreg), 
  hsten = as.factor(hsten), tenure = as.factor(tenure), hsten2 = as.factor(hsten2), tenown2c = as.factor(tenown2c), 
  tenown3c = as.factor(tenown3c), finsec = as.factor(finsec))
M_val <- M_val_IDsub[, match(ph2$sentrix, colnames(M_val_IDsub))]

#########################################################################
#########################################################################
#########################################################################

## Model 12 - With DNA Methylation data

# Initialize results list
results <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){ #
	x <- M_val[i,]
	mod <- lm(x ~ relevel(ph2$empma, ref = "1") + ph2$AgeMom + relevel(ph2$Eth9, ref = "1") + relevel(ph2$eduma, ref = "4") + ph2$Bcell + ph2$CD4T + ph2$CD8T + ph2$Eos + ph2$Mono + ph2$Neu + ph2$NK)
	results[[i]] <- summary(mod)$coef[2:3,]
	if(i %in% c(1, nrow(M_val))) {
		print(Sys.time())
		print("Model 12")
	}
}

# Compile results list into data frame
names(results) <- rownames(M_val)
for(i in 1:length(results)){
	results[[i]] <- results[[i]] %>% as.data.frame() %>% rownames_to_column()
}
big_data <- do.call(rbind, results) %>% as.data.frame()
names(big_data)[1:5] <- c("coefficients", "beta", "se", "t", "pvalue")

# Calculate qvalue and filter by pvalue
qval <- big_data %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq <- big_data %>% rownames_to_column('dnam_id') %>% left_join(., qval, by = "dnam_id") %>% subset(., qval <= 0.05)
# write.table(resq, "Model12_New_Out/V5/model12_qfilt_v5.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

top50 <- big_data %>% rownames_to_column('dnam_id') %>% left_join(., qval, by = "dnam_id") %>% arrange(pvalue)
# write.table(top50[1:50,], "Model12_New_Out/V5/model12_top50_v5.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Export model results to RData file
save(list = ls(pattern = "results|big_data|qval|resq"), file = "Model12_New_Out/V5/Model12_unstrat_v5.RData")

#########################################################################
#########################################################################
#########################################################################

# Model 12a - Stratified by Offspring Sex, With DNA Methylation Data

# Gender 1 - Initialize results matrices, Run Models
results_gen1 <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){
  x <- M_val[i, colnames(M_val) %in% ph2[ph2$gender == 1, "sentrix"]]
  df <- ph2 %>% subset(., gender == 1)
  mod <- lm(x ~ relevel(df$empma, ref = "1") + df$AgeMom + relevel(df$Eth9, ref = "1") + relevel(df$eduma, ref = "4") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_gen1[[i]] <- summary(mod)$coef[2:3,]
  if(i %in% c(1, nrow(M_val))) {
    print(Sys.time())
    print("Model 12a - Gender == 1")
  }
}

# Compile results list into data frame
names(results_gen1) <- rownames(M_val)
for(i in 1:length(results_gen1)){
	results_gen1[[i]] <- results_gen1[[i]] %>% as.data.frame() %>% rownames_to_column()
}
big_data_gen1 <- do.call(rbind, results_gen1) %>% as.data.frame()
names(big_data_gen1)[1:5] <- c("coefficients", "beta", "se", "t", "pvalue")

# Calculate qvalue and filter by pvalue
qval_gen1 <- big_data_gen1 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_gen1 <- big_data_gen1 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen1, by = "dnam_id") %>% subset(., qval <= 0.05)
# write.table(resq_gen1, "Model12_New_Out/V5/model12a_gen1_qfilt_v5.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

top50_gen1 <- big_data_gen1 %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen1, by = "dnam_id") %>% arrange(pvalue)
# write.table(top50_gen1[1:50,], "Model12_New_Out/V5/model12_gen1_top50_v5.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Export model results to RData file
save(list = ls(pattern = "gen1"), file = "Model12_New_Out/V5/Model12_gen1_v5.RData")

#########################################################################

# Gender 2 - Initialize results matrices, Run Models
results_gen2 <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){
  x <- M_val[i, colnames(M_val) %in% ph2[ph2$gender == 2, "sentrix"]]
  df <- ph2 %>% subset(., gender == 2)
  mod <- lm(x ~ relevel(df$empma, ref = "1") + df$AgeMom + relevel(df$Eth9, ref = "1") + relevel(df$eduma, ref = "4") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_gen2[[i]] <- summary(mod)$coef[2:3,]
  if(i %in% c(1, nrow(M_val))) {
    print(Sys.time())
    print("Model 12a - Gender == 2")
  }
}

# Compile results list into data frame
names(results_gen2) <- rownames(M_val)
for(i in 1:length(results_gen2)){
	results_gen2[[i]] <- results_gen2[[i]] %>% as.data.frame() %>% rownames_to_column()
}
big_data_gen2 <- do.call(rbind, results_gen2) %>% as.data.frame()
names(big_data_gen2)[1:5] <- c("coefficients", "beta", "se", "t", "pvalue")

# Calculate qvalue and filter by pvalue
qval_gen2 <- big_data_gen2 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_gen2 <- big_data_gen2 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen2, by = "dnam_id") %>% subset(., qval <= 0.05)
# write.table(resq_gen2, "Model12_New_Out/V5/model12a_gen2_qfilt_v5.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

top50_gen2 <- big_data_gen2 %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen2, by = "dnam_id") %>% arrange(pvalue)
# write.table(top50_gen2[1:50,], "Model12_New_Out/V5/model12_gen2_top50_v5.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Export model results to RData file
save(list = ls(pattern = "gen2"), file = "Model12_New_Out/V5/Model12_gen2_v5.RData")

#########################################################################
#########################################################################
#########################################################################

## Model 12b - Stratified by Mother's Ethnicity, With DNA Methylation Data

# Ethnicity 1 - Initialize results matrices, Run Models
results_eth1 <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){
  x <- M_val[i, colnames(M_val) %in% ph2[ph2$Eth9 == 1, "sentrix"]]
  df <- ph2 %>% subset(., Eth9 == 1)
  mod <- lm(x ~ relevel(df$empma, ref = "1") + df$AgeMom + relevel(df$eduma, ref = "4") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_eth1[[i]] <- summary(mod)$coef[2:3,]
  if(i %in% c(1, nrow(M_val))) {
    print(Sys.time())
    print("Model 12b - Ethnicity == 1")
  }
}
#results_eth1 <- results_gen1 %>% as.data.frame() %>% subset(beta != 0)

# Compile results list into data frame
names(results_eth1) <- rownames(M_val)
for(i in 1:length(results_eth1)){
	results_eth1[[i]] <- results_eth1[[i]] %>% as.data.frame() %>% rownames_to_column()
}
big_data_eth1 <- do.call(rbind, results_eth1) %>% as.data.frame()
names(big_data_eth1)[1:5] <- c("coefficients", "beta", "se", "t", "pvalue")

# Calculate qvalue and filter by pvalue
qval_eth1 <- big_data_eth1 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_eth1 <- big_data_eth1 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth1, by = "dnam_id") %>% subset(., qval <= 0.05)
# write.table(resq_eth1, "Model12_New_Out/V5/model12b_eth1_qfilt_v5.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

top50_eth1 <- big_data_eth1 %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth1, by = "dnam_id") %>% arrange(pvalue)
# write.table(top50_eth1[1:50,], "Model12_New_Out/V5/model12_eth1_top50_v5.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Export model results to RData file
save(list = ls(pattern = "eth1"), file = "Model12_New_Out/V5/Model12_eth1_v5.RData")

#########################################################################

# Ethnicity 7 - Initialize results matrices, Run Models
results_eth7 <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){
  x <- M_val[i, colnames(M_val) %in% ph2[ph2$Eth9 == 7, "sentrix"]]
  df <- ph2 %>% subset(., Eth9 == 7)
  mod <- lm(x ~ relevel(df$empma, ref = "1") + df$AgeMom + relevel(df$eduma, ref = "4") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_eth7[[i]] <- summary(mod)$coef[2:3,]
  if(i %in% c(1, nrow(M_val))) {
    print(Sys.time())
    print("Model 12b - Ethnicity == 7")
  }
}

# Compile results list into data frame
names(results_eth7) <- rownames(M_val)
for(i in 1:length(results_eth7)){
	results_eth7[[i]] <- results_eth7[[i]] %>% as.data.frame() %>% rownames_to_column()
}
big_data_eth7 <- do.call(rbind, results_eth7) %>% as.data.frame()
names(big_data_eth7)[1:5] <- c("coefficients", "beta", "se", "t", "pvalue")

# Calculate qvalue and filter by pvalue
qval_eth7 <- big_data_eth7 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_eth7 <- big_data_eth7 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth7, by = "dnam_id") %>% subset(., qval <= 0.05)
# write.table(resq_eth7, "Model12_New_Out/V5/model12b_eth7_qfilt_v5.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

top50_eth7 <- big_data_eth7 %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth7, by = "dnam_id") %>% arrange(pvalue)
# write.table(top50_eth7[1:50,], "Model12_New_Out/V5/model12_eth7_top50_v5.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Export model results to RData file
save(list = ls(pattern = "eth7"), file = "Model12_New_Out/V5/Model12_eth7_v5.RData")





