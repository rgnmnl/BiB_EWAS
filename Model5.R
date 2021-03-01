#########################################################################
## Title: Model 5 EWAS
## Version: 2 
## Note: Reran Model 5 on Feb 2020 after finding duplicate IDs in 
##		original dataset. Model covariates did not change from previous
##		version.
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
  mutate(., Eth9 = as.factor(Eth9), edufa = as.factor(edufa), MatSD = as.factor(MatSD), 
         eduma = as.factor(eduma), bn_mtben = as.factor(bn_mtben), quin07 = as.factor(quin07), 
         empstma = as.factor(empstma), smkpreg = as.factor(smkpreg), alpreg = as.factor(alpreg), 
         gender = as.factor(gender)) #%>% left_join(., pcs, by = c("sentrix" = "V1"))
M_val <- M_val_IDsub[, match(ph2$sentrix, colnames(M_val_IDsub))]

#########################################################################
#########################################################################
#########################################################################

# Model 4: Index of multiple deprivation 2007 -> DNAm (also adjusted for mother’s age and ethnicity as well as cell type composition)

# Initialize results list
results <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){ 
	x <- M_val[i,]
	mod <- lm(x ~ relevel(ph2$quin07, ref = "1") + ph2$AgeMom + relevel(ph2$Eth9, ref = "1") + relevel(ph2$eduma, ref = "4") + ph2$Bcell + ph2$CD4T + ph2$CD8T + ph2$Eos + ph2$Mono + ph2$Neu + ph2$NK)
	results[[i]] <- summary(mod)$coef[2:5,]
	if(i %in% c(1, nrow(M_val))) {
		print(Sys.time())
		print("Model 4")
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
resq <- big_data %>% rownames_to_column('dnam_id') %>% left_join(., qval, by = "dnam_id") %>% subset(., pvalue <= 5e-06)
top50 <- big_data %>% rownames_to_column('dnam_id') %>% left_join(., qval, by = "dnam_id") %>% arrange(pvalue)

# Export model results to RData file
save(list = ls(pattern = "results|big_data|qval|resq"), file = "Model5_New_Out/Model5_Feb_rerun.RData")

#########################################################################
#########################################################################
#########################################################################

## Model 5a - Stratified by Offspring Sex, With DNA Methylation Data

# Gender 1 - Initialize results matrices, Run Models
results_gen1 <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){
  x <- M_val[i, colnames(M_val) %in% ph2[ph2$gender == 1, "sentrix"]]
  df <- ph2 %>% subset(., gender == 1)
  mod <- lm(x ~ relevel(df$quin07, ref = "1") + df$AgeMom + relevel(df$Eth9, ref = "1") + relevel(df$eduma, ref = "4") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_gen1[[i]] <- summary(mod)$coef[2:5,]
  if(i %in% c(1, nrow(M_val))) {
    print(Sys.time())
    print("Model 5a - Gender == 1")
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
resq_gen1 <- big_data_gen1 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen1, by = "dnam_id") %>% subset(., pvalue <= 5e-06)
top50_gen1 <- big_data_gen1 %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen1, by = "dnam_id") %>% arrange(pvalue)

# Export model results to RData file
save(list = ls(pattern = "gen1"), file = "Model5_New_Out/Model5_gen1.RData")

#########################################################################

# Gender 2 - Initialize results matrices, Run Models
results_gen2 <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){ #
  x <- M_val[i, colnames(M_val) %in% ph2[ph2$gender == 2, "sentrix"]]
  df <- ph2 %>% subset(., gender == 2)
  mod <- lm(x ~ relevel(df$quin07, ref = "1") + df$AgeMom + relevel(df$Eth9, ref = "1") + relevel(df$eduma, ref = "4") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_gen2[[i]] <- summary(mod)$coef[2:5,]
  if(i %in% c(1, nrow(M_val))) {
    print(Sys.time())
    print("Model 5a - Gender == 2")
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
resq_gen2 <- big_data_gen2 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen2, by = "dnam_id") %>% subset(., pvalue <= 5e-06)
top50_gen2 <- big_data_gen2 %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen2, by = "dnam_id") %>% arrange(pvalue)

# Export model results to RData file
save(list = ls(pattern = "gen2"), file = "Model5_New_Out/Model5_gen2.RData")

#########################################################################
#########################################################################
#########################################################################

## Model 4b - Stratified by Mother's Ethnicity, With DNA Methylation Data

# Ethnicity 1 - Initialize results matrices, Run Models
results_eth1 <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){
  x <- M_val[i, colnames(M_val) %in% ph2[ph2$Eth9 == 1, "sentrix"]]
  df <- ph2 %>% subset(., Eth9 == 1)
  mod <- lm(x ~ relevel(df$quin07, ref = "1") + df$AgeMom + relevel(df$eduma, ref = "4") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_eth1[[i]] <- summary(mod)$coef[2:5,]
  if(i %in% c(1, nrow(M_val))) {
    print(Sys.time())
    print("Model 5b - Ethnicity == 1")
  }
}

# Compile results list into data frame
names(results_eth1) <- rownames(M_val)
for(i in 1:length(results_eth1)){
	results_eth1[[i]] <- results_eth1[[i]] %>% as.data.frame() %>% rownames_to_column()
}
big_data_eth1 <- do.call(rbind, results_eth1) %>% as.data.frame()
names(big_data_eth1)[1:5] <- c("coefficients", "beta", "se", "t", "pvalue")

# Calculate qvalue and filter by pvalue
qval_eth1 <- big_data_eth1 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_eth1 <- big_data_eth1 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth1, by = "dnam_id") %>% subset(., pvalue <= 5e-06)
top50_eth1 <- big_data_eth1 %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth1, by = "dnam_id") %>% arrange(pvalue)

# Export model results to RData file
save(list = ls(pattern = "eth1"), file = "Model5_eth1.RData")

#########################################################################

# Ethnicity 7 - Initialize results matrices, Run Models
results_eth7 <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){
  x <- M_val[i, colnames(M_val) %in% ph2[ph2$Eth9 == 7, "sentrix"]]
  df <- ph2 %>% subset(., Eth9 == 7)
  mod <- lm(x ~ relevel(df$quin07, ref = "1") + df$AgeMom + relevel(df$eduma, ref = "4") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_eth7[[i]] <- summary(mod)$coef[2:5,]
  if(i %in% c(1, nrow(M_val))) {
    print(Sys.time())
    print("Model 5b - Ethnicity == 7")
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
resq_eth7 <- big_data_eth7 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth7, by = "dnam_id") %>% subset(., pvalue <= 5e-06)
top50_eth7 <- big_data_eth7 %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth7, by = "dnam_id") %>% arrange(pvalue)

# Export model results to RData file
save(list = ls(pattern = "eth7"), file = "Model5_New_Out/Model5_eth7.RData")

