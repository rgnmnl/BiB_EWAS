#########################################################################
## Title: Model 1 EWAS
## Version: 1
## Author: Regina Manansala
## Date: 29-October-2019
#########################################################################

library(car)
library(data.table)
library(dplyr)
library(tidyverse)
library(qvalue)

load("normalized_betas_kids_for_analysis.RData")
load("CellCountsBloodgse35069Complete.Rdata")
ph <- read.csv("bibgwas.csv")
#pcs <- fread("PC1-10.txt")

#### Make sure the methylation and phenotype files are sorted exactly the same way.
counts <- data.frame(counts)
counts$sentrix <- rownames(counts)
ph2 <- merge(ph, counts) %>% 
  mutate(., Eth9 = as.factor(Eth9), edufa = as.factor(edufa), MatSD = as.factor(MatSD), 
         eduma = as.factor(eduma), bn_mtben = as.factor(bn_mtben), quin07 = as.factor(quin07), 
         empstma = as.factor(empstma), smkpreg = as.factor(smkpreg), alpreg = as.factor(alpreg), 
         gender = as.factor(gender)) #%>% left_join(., pcs, by = c("sentrix" = "V1"))
betas <- norm.beta3 + 0.0001 
M_val <- log2(betas/(1-betas))
M_val <- M_val[, match(ph2$sentrix, colnames(M_val))]

#########################################################################


### INCORRECT MODELS ###
# Model 1 - With PCs, Without DNA Methylation DataLossWarning
# 
# mod_pcs <- lm(MatSD ~ AgeMom + Eth9 + eduma + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ph2)
# vif(mod_pcs)
# 
# Model 1 - Without DNA Methylation data
# base_mod1 <- lm(MatSD ~ AgeMom + relevel(Eth9, ref = "1") + relevel(eduma, ref = "4") + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK, data = ph2)
# 
# results_base_mod1 <- summary(base_mod1)$coef[2:14,c(1,2,4)] %>% as.data.frame()
# write.table(results_base_mod1, "model1_results_base.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
###########################

# Model 1 - With DNA Methylation data

### THIS IS NO LONGER NEEDED ###
## Initialize results matrix
# results <- matrix(0, nrow=dim(betas)[1], ncol=4)
# rownames(results) <- rownames(betas)
# colnames(results) <- c("beta", "se", "t", "pvalue")
#################################

results <- list()

# results <- data.frame(matrix(ncol = 5, nrow = 0))
# colnames(results) <- c("beta", "se", "t", "pvalue", "i")

for(i in 1:10){
	x <- M_val[i,]
	#mod <- lm(ph2$MatSD ~ x + ph2$AgeMom + ph2$Eth9 + ph2$eduma + ph2$Bcell + ph2$CD4T + ph2$CD8T + ph2$Eos + ph2$Mono + ph2$Neu + ph2$NK)
	mod <- lm(x ~ relevel(ph2$MatSD, ref = "1") + ph2$AgeMom + relevel(ph2$Eth9, ref = "1") + relevel(ph2$eduma, ref = "4") + ph2$Bcell + ph2$CD4T + ph2$CD8T + ph2$Eos + ph2$Mono + ph2$Neu + ph2$NK)
	results[[i]] <- summary(mod)$coef[2:5,]
	#names(datalist)[[i]] <- rownames(M_val[i,])
	#results <- summary(mod)$coef[2:5,] %>% as.data.frame() %>% mutate(i = rownames(M_val)[i]) %>% rbind(results)
	#print(rownames(M_val)[i])
# 	results[i,] <- summary(mod)$coef[2,]
	if(i %in% c(1, nrow(M_val))) {
		print(Sys.time())
		print("Model 1 - With DNA Methylation")
	}
}


# which(rownames(M_val) == "cg26709118", arr.ind = TRUE)
# # 1051, 1992, 5753, 8715
# foo <- M_val[27587,] %>% as.data.frame() %>% rownames_to_column('sentrix')
# names(foo)[2] <- "cpg"
# head(foo[order(foo$cpg),])
# head(foo[foo$cpg == "-Inf",])
# 
# M_val[8715, "200992320006_R05C01"]
# norm.beta3[8715, "200992320006_R05C01"]
# 
# do.call(rbind, lapply(foo, function(x) summary(x)))

# bh <- results %>% as.data.frame() %>% select(., pvalue) %>% p.adjust(., method = "BH")
qval <- results %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq <- results %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval, by = "dnam_id") %>% subset(., qval <= .05)

write.table(resq, "model1_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#########################################################################
#########################################################################
#########################################################################

## Model 1a - Stratified by Offspring Sex, Without DNA Methylation Data

base_mod1a_gen1 <- lm(MatSD ~ AgeMom + relevel(Eth9, ref = "1") + relevel(eduma, ref = "4") + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK, data = ph2[ph2$gender == 1,])
base_mod1a_gen2 <- lm(MatSD ~ AgeMom + relevel(Eth9, ref = "1") + relevel(eduma, ref = "4") + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK, data = ph2[ph2$gender == 2,])

results_base_mod1a_gen1 <- summary(base_mod1a_gen1)$coef[2:14,c(1,2,4)] %>% as.data.frame() %>% rownames_to_column('covariate') %>% mutate(gender = 1)
results_base_mod1a_gen2 <- summary(base_mod1a_gen2)$coef[2:14,c(1,2,4)] %>% as.data.frame() %>% rownames_to_column('covariate') %>% mutate(gender = 2)

results_base_mod1a <- rbind(results_base_mod1a_gen1, results_base_mod1a_gen2)
write.table(results_base_mod1a, "model1a_results_base.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

## Model 1a - Stratified by Offspring Sex, With DNA Methylation Data

# Gender 1 - Initialize results matrices, Run Models
results_gen1 <- matrix(0, nrow=dim(betas)[1], ncol=4)
rownames(results_gen1) <- rownames(betas)
colnames(results_gen1) <- c("beta", "se", "t", "pvalue")

for(i in 1:nrow(betas)){
  x <- betas[i, colnames(betas) %in% ph2[ph2$gender == 1, "sentrix"]]
  df <- ph2 %>% subset(., gender == 1)
  mod <- lm(df$MatSD ~ x + df$AgeMom + relevel(df$Eth9, ref = "1") + relevel(df$eduma, ref = "4") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_gen1[i,] <- summary(mod)$coef[2,]
  if(i %in% c(1, nrow(betas))) {
    print(Sys.time())
    print("Model 1a - Gender == 1")
  }
}
#results_gen1 <- results_gen1 %>% as.data.frame() %>% subset(beta != 0)

qval_gen1 <- results_gen1 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_gen1 <- results_gen1 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen1, by = "dnam_id") %>% subset(., qval_gen1 <= .05)
write.table(resq_gen1, "model1a_gen1_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


# Gender 2 - Initialize results matrices, Run Models
results_gen2 <- matrix(0, nrow=dim(betas)[1], ncol=4)
rownames(results_gen2) <- rownames(betas)
colnames(results_gen2) <- c("beta", "se", "t", "pvalue")

for(i in 1:nrow(betas)){
  x <- betas[i, colnames(betas) %in% ph2[ph2$gender == 2, "sentrix"]]
  df <- ph2 %>% subset(., gender == 2)
  mod <- lm(df$MatSD ~ x + df$AgeMom + relevel(df$Eth9, ref = "1") + relevel(df$eduma, ref = "4") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_gen2[i,] <- summary(mod)$coef[2,]
  if(i %in% c(1, nrow(betas))) {
    print(Sys.time())
    print("Model 1a - Gender == 2")
  }
}

qval_gen2 <- results_gen2 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_gen2 <- results_gen2 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen2, by = "dnam_id") %>% subset(., qval_gen2 <= .05)
write.table(resq_gen2, "model1a_gen2_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#########################################################################
#########################################################################
#########################################################################

## Model 1b - Stratified by Mother's Ethnicity, Without DNA Methylation Data

base_mod1b_eth1 <- lm(MatSD ~ AgeMom + relevel(eduma, ref = "4") + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK, data = ph2[ph2$Eth9 == 1,])
base_mod1b_eth7 <- lm(MatSD ~ AgeMom + relevel(eduma, ref = "4") + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK, data = ph2[ph2$Eth9 == 7,])

results_base_mod1b_eth1 <- summary(base_mod1b_eth1)$coef[2:13,c(1,2,4)] %>% as.data.frame() %>% rownames_to_column('covariate') %>% mutate(Eth9 = 1)
results_base_mod1b_eth7 <- summary(base_mod1b_eth7)$coef[2:13,c(1,2,4)] %>% as.data.frame() %>% rownames_to_column('covariate') %>% mutate(Eth9 = 7)

results_base_mod1b <- rbind(results_base_mod1b_eth1, results_base_mod1b_eth7)
write.table(results_base_mod1b, "model1b_results_base.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

## Model 1b - Stratified by Mother's Ethnicity, With DNA Methylation Data

# Ethnicity 1 - Initialize results matrices, Run Models

results_eth1 <- matrix(0, nrow=dim(betas)[1], ncol=4)
rownames(results_eth1) <- rownames(betas)
colnames(results_eth1) <- c("beta", "se", "t", "pvalue")

for(i in 1:nrow(betas)){
  x <- betas[i, colnames(betas) %in% ph2[ph2$Eth9 == 1, "sentrix"]]
  df <- ph2 %>% subset(., Eth9 == 1)
  mod <- lm(df$MatSD ~ x + df$AgeMom + relevel(df$eduma, ref = "4") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_eth1[i,] <- summary(mod)$coef[2,]
  if(i %in% c(1, nrow(betas))) {
    print(Sys.time())
    print("Model 1b - Ethnicity == 1")
  }
}
#results_eth1 <- results_gen1 %>% as.data.frame() %>% subset(beta != 0)

qval_eth1 <- results_eth1 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_eth1 <- results_eth1 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth1, by = "dnam_id") %>% subset(., qval_eth1 <= .05)
write.table(resq_eth1, "model1b_eth1_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Ethnicity 7 - Initialize results matrices, Run Models

results_eth7 <- matrix(0, nrow=dim(betas)[1], ncol=4)
rownames(results_eth7) <- rownames(betas)
colnames(results_eth7) <- c("beta", "se", "t", "pvalue")

for(i in 1:nrow(betas)){
  x <- betas[i, colnames(betas) %in% ph2[ph2$Eth9 == 7, "sentrix"]]
  df <- ph2 %>% subset(., Eth9 == 7)
  mod <- lm(df$MatSD ~ x + df$AgeMom + relevel(df$eduma, ref = "4") + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_eth7[i,] <- summary(mod)$coef[2,]
  if(i %in% c(1, nrow(betas))) {
    print(Sys.time())
    print("Model 1b - Ethnicity == 7")
  }
}

qval_eth7 <- results_eth7 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_eth7 <- results_eth7 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth7, by = "dnam_id") %>% subset(., qval_eth7 <= .05)
write.table(resq_eth7, "model1b_eth7_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# results_gen1 <- qval()
#write.table(results_gen1, "~/Data/Simanek_BiB/DNAm/m1a_gen1.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

# results_gen2 <- qval()
#write.table(results_gen2, "~/Data/Simanek_BiB/DNAm/m1a_gen2.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

save.image("Model1.RData")


# results2 <- results %>% as.data.frame() %>% subset(pvalue < 5e-06) %>% rownames_to_column('dnam_id') %>% left_join(., qval, by = "dnam_id")
# results2_gen1 <- results_gen1 %>% as.data.frame() %>% subset(pvalue < 5e-06) %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen1, by = "dnam_id")
# results2_gen2 <- results_gen2 %>% as.data.frame() %>% subset(pvalue < 5e-06) %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen2, by = "dnam_id")
# results2_eth1 <- results_eth1 %>% as.data.frame() %>% subset(pvalue < 5e-06) %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth1, by = "dnam_id")
# results2_eth7 <- results_eth7 %>% as.data.frame() %>% subset(pvalue < 5e-06) %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth7, by = "dnam_id")
# 
# results2$model  <- "1"
# results2_gen1$model <- "1a; Gender = 1"
# results2_eth1$model <- "1b; Ethnicity = 1"
# results2_eth7$model <- "1b; Ethnicity = 7"
# 
# model1_out <- rbind(results2, results2_gen1, results2_eth1, results2_eth7)
# write.table(model1_out, "model1_pfilt.txt", sep = "\t", row.names = F, col.names = T, quote = F)

