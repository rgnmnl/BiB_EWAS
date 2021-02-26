#########################################################################
## Title: Model 3 EWAS
## Version: 3 - Binomial Maternal Education
## Author: Regina Manansala
## Date: 21-November-2019
#########################################################################

library(car)
library(data.table)
library(dplyr)
library(tidyverse)
library(qvalue)

load("normalized_betas_kids_for_analysis.RData")
load("CellCountsBloodgse35069Complete.Rdata")
ph <- read.csv("bibgwas_07NOV2019.csv")
#pcs <- fread("PC1-10.txt")

#### Make sure the methylation and phenotype files are sorted exactly the same way.
counts <- data.frame(counts)
counts$sentrix <- rownames(counts)
ph2 <- merge(ph, counts) %>% 
  mutate(., Eth9 = as.factor(Eth9), edufa = as.factor(edufa), edufa2 = as.factor(edufa2), MatSD = as.factor(MatSD), 
         eduma = as.factor(eduma), eduma2 = as.factor(eduma2), bn_mtben = as.factor(bn_mtben), quin07 = as.factor(quin07), 
         empstma = as.factor(empstma), smkpreg = as.factor(smkpreg), alpreg = as.factor(alpreg), gender = as.factor(gender)) 
         #%>% left_join(., pcs, by = c("sentrix" = "V1"))
betas <- norm.beta3 + 0.0001 
M_val <- log2(betas/(1-betas))
M_val <- M_val[, match(ph2$sentrix, colnames(M_val))]

#########################################################################
#########################################################################
#########################################################################

# Model 3 - With DNA Methylation data

#Initialize results list
results <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){
	x <- M_val[i,]
	mod <- lm(x ~ relevel(ph2$eduma2, ref = "2") + ph2$AgeMom + relevel(ph2$Eth9, ref = "1") + ph2$Bcell + ph2$CD4T + ph2$CD8T + ph2$Eos + ph2$Mono + ph2$Neu + ph2$NK)
	results[[i]] <- summary(mod)$coef[2,]
	if(i %in% c(1, nrow(M_val))) {
		print(Sys.time())
		print("Model 3")
	}
}

# Compile results list into data frame
names(results) <- rownames(M_val)
big_data <- do.call(rbind, results) %>% as.data.frame()
names(big_data)[1:4] <- c("beta", "se", "t", "pvalue")

# Calculate qvalue and filter by pvalue
qval <- big_data %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq <- big_data %>% rownames_to_column('dnam_id') %>% left_join(., qval, by = "dnam_id") %>% subset(., pvalue <= 5e-06)
#write.table(resq, "Model3_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Export model results to RData file
save(list = ls(pattern = "big_data|qval|resq"), file = "Model3_unstrat.RData")

#########################################################################
#########################################################################
#########################################################################

# Model 3a - Stratified by Offspring Sex, With DNA Methylation Data

# Gender 1 - Initialize results matrices, Run Models
results_gen1 <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){
  x <- M_val[i, colnames(M_val) %in% ph2[ph2$gender == 1, "sentrix"]]
  df <- ph2 %>% subset(., gender == 1)
  mod <- lm(x ~ relevel(df$eduma2, ref = "2") + df$AgeMom + relevel(df$Eth9, ref = "1") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_gen1[[i]] <- summary(mod)$coef[2,]
  if(i %in% c(1, nrow(M_val))) {
    print(Sys.time())
    print("Model 3a - Gender == 1")
  }
}

# Compile results list into data frame
names(results_gen1) <- rownames(M_val)
big_data_gen1 <- do.call(rbind, results_gen1) %>% as.data.frame()
names(big_data_gen1)[1:4] <- c("beta", "se", "t", "pvalue")

# Calculate qvalue and filter by pvalue
qval_gen1 <- big_data_gen1 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_gen1 <- big_data_gen1 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen1, by = "dnam_id") %>% subset(., pvalue <= 5e-06)
#write.table(resq_gen1, "Model3a_gen1_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Export model results to RData file
save(list = ls(pattern = "gen1"), file = "Model3_ed_gen1.RData")

#########################################################################

# Gender 2 - Initialize results matrices, Run Models
results_gen2 <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){
  x <- M_val[i, colnames(M_val) %in% ph2[ph2$gender == 2, "sentrix"]]
  df <- ph2 %>% subset(., gender == 2)
  mod <- lm(x ~ relevel(df$eduma2, ref = "2") + df$AgeMom + relevel(df$Eth9, ref = "1") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_gen2[[i]] <- summary(mod)$coef[2,]
  if(i %in% c(1, nrow(M_val))) {
    print(Sys.time())
    print("Model 3a - Gender == 2")
  }
}

# Compile results list into data frame
names(results_gen2) <- rownames(M_val)
big_data_gen2 <- do.call(rbind, results_gen2) %>% as.data.frame()
names(big_data_gen2)[1:4] <- c("beta", "se", "t", "pvalue")

# Calculate qvalue and filter by pvalue
qval_gen2 <- big_data_gen2 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_gen2 <- big_data_gen2 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen2, by = "dnam_id") %>% subset(., pvalue <= 5e-06)
#write.table(resq_gen2, "Model3a_gen2_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Export model results to RData file
save(list = ls(pattern = "gen2"), file = "Model3_ed_gen2.RData")

#########################################################################
#########################################################################
#########################################################################

## Model 3b - Stratified by Mother's Ethnicity, With DNA Methylation Data

# Ethnicity 1 - Initialize results matrices, Run Models
results_eth1 <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){
  x <- M_val[i, colnames(M_val) %in% ph2[ph2$Eth9 == 1, "sentrix"]]
  df <- ph2 %>% subset(., Eth9 == 1)
  mod <- lm(x ~ relevel(df$eduma2, ref = "2") + df$AgeMom + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_eth1[[i]] <- summary(mod)$coef[2,]
  if(i %in% c(1, nrow(M_val))) {
    print(Sys.time())
    print("Model 3b - Ethnicity == 1")
  }
}
#results_eth1 <- results_gen1 %>% as.data.frame() %>% subset(beta != 0)

# Compile results list into data frame
names(results_eth1) <- rownames(M_val)
big_data_eth1 <- do.call(rbind, results_eth1) %>% as.data.frame()
names(big_data_eth1)[1:4] <- c("beta", "se", "t", "pvalue")

# Calculate qvalue and filter by pvalue
qval_eth1 <- big_data_eth1 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_eth1 <- big_data_eth1 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth1, by = "dnam_id") %>% subset(., pvalue <= 5e-06)
#write.table(resq_eth1, "Model3b_eth1_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Export model results to RData file
save(list = ls(pattern = "eth1"), file = "Model3_ed_eth1.RData")

#########################################################################

# Ethnicity 7 - Initialize results matrices, Run Models
results_eth7 <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){
  x <- M_val[i, colnames(M_val) %in% ph2[ph2$Eth9 == 7, "sentrix"]]
  df <- ph2 %>% subset(., Eth9 == 7)
  mod <- lm(x ~ relevel(df$eduma2, ref = "2") + df$AgeMom + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_eth7[[i]] <- summary(mod)$coef[2,]
  if(i %in% c(1, nrow(M_val))) {
    print(Sys.time())
    print("Model 3b - Ethnicity == 7")
  }
}

# Compile results list into data frame
names(results_eth7) <- rownames(M_val)
big_data_eth7 <- do.call(rbind, results_eth7) %>% as.data.frame()
names(big_data_eth7)[1:4] <- c("beta", "se", "t", "pvalue")

# Calculate qvalue and filter by pvalue
qval_eth7 <- big_data_eth7 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_eth7 <- big_data_eth7 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth7, by = "dnam_id") %>% subset(., pvalue <= 5e-06)
#write.table(resq_eth7, "Model3b_eth7_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Export model results to RData file
save(list = ls(pattern = "eth7"), file = "Model3_ed_eth7.RData")

