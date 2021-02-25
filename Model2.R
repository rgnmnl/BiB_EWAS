#########################################################################
## Title: Model 2 EWAS
## Version: 1
## Author: Regina Manansala
## Date: 31-October-2019
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
ph2 <- merge(ph, counts) #%>% left_join(., pcs, by = c("sentrix" = "V1"))
betas <- norm.beta3[, match(ph2$sentrix, colnames(norm.beta3))]

#Model 2: Mother’s educational attainment -> DNAm (also adjusted for mother’s age and ethnicity as well as cell type composition)

## Model 2 - Without DNA Methylation data
base_mod2 <- lm(eduma ~ AgeMom + Eth9 + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK, data = ph2)

results_base_mod2 <- summary(base_mod2)$coef[2:10,c(1,2,4)] %>% as.data.frame()
write.table(results_base_mod2, "Model2_Out/model2_results_base.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

## Model 2 - With DNA Methylation data

results <- matrix(0, nrow=dim(betas)[1], ncol=4)
rownames(results) <- rownames(betas)
colnames(results) <- c("beta", "se", "t", "pvalue")

for(i in 1:nrow(betas)){
	x <- betas[i,]
	mod <- lm(ph2$eduma ~ x + ph2$AgeMom + ph2$Eth9 + ph2$Bcell + ph2$CD4T + ph2$CD8T + ph2$Eos + ph2$Mono + ph2$Neu + ph2$NK)
	results[i,] <- summary(mod)$coef[2,]
	if(i %in% c(1, dim(results)[1])) {
		print(Sys.time())
		print("Model 2 - With DNA Methylation")
	}
}

qval <- results %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq <- results %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval, by = "dnam_id") %>% subset(., qval <= .05)
write.table(resq, "Model2_Out/model2_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


#Models 2a: Also stratify by offspring sex

## Model 2a - Stratified by Offspring Sex, Without DNA Methylation Data

base_mod2a_gen1 <- lm(eduma ~ AgeMom + Eth9 + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK, data = ph2[ph2$gender == 1,])
base_mod2a_gen2 <- lm(eduma ~ AgeMom + Eth9 + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK, data = ph2[ph2$gender == 2,])

results_base_mod2a_gen1 <- summary(base_mod2a_gen1)$coef[2:10,c(1,2,4)] %>% as.data.frame() %>% rownames_to_column('covariate') %>% mutate(gender = 1)
results_base_mod2a_gen2 <- summary(base_mod2a_gen2)$coef[2:10,c(1,2,4)] %>% as.data.frame() %>% rownames_to_column('covariate') %>% mutate(gender = 2)

results_base_mod2a <- rbind(results_base_mod2a_gen1, results_base_mod2a_gen2)
write.table(results_base_mod2a, "Model2_Out/model2a_results_base.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

## Model 2a - Stratified by Offspring Sex, With DNA Methylation Data

# Gender 1 - Initialize results matrices, Run Models
results_gen1 <- matrix(0, nrow=dim(betas)[1], ncol=4)
rownames(results_gen1) <- rownames(betas)
colnames(results_gen1) <- c("beta", "se", "t", "pvalue")

for(i in 1:nrow(betas)){
  x <- betas[i, colnames(betas) %in% ph2[ph2$gender == 1, "sentrix"]]
  df <- ph2 %>% subset(., gender == 1)
  mod <- lm(df$eduma ~ x + df$AgeMom + df$Eth9 + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_gen1[i,] <- summary(mod)$coef[2,]
  if(i %in% c(1, nrow(betas))) {
    print(Sys.time())
    print("Model 2a - Gender == 1")
  }
}
#results_gen1 <- results_gen1 %>% as.data.frame() %>% subset(beta != 0)

qval_gen1 <- results_gen1 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_gen1 <- results_gen1 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen1, by = "dnam_id") %>% subset(., qval_gen1 <= .05)
write.table(resq_gen1, "Model2_Out/model2a_gen1_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


# Gender 2 - Initialize results matrices, Run Models
results_gen2 <- matrix(0, nrow=dim(betas)[1], ncol=4)
rownames(results_gen2) <- rownames(betas)
colnames(results_gen2) <- c("beta", "se", "t", "pvalue")

for(i in 1:nrow(betas)){
  x <- betas[i, colnames(betas) %in% ph2[ph2$gender == 2, "sentrix"]]
  df <- ph2 %>% subset(., gender == 2)
  mod <- lm(df$eduma ~ x + df$AgeMom + df$Eth9 + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_gen2[i,] <- summary(mod)$coef[2,]
  if(i %in% c(1, nrow(betas))) {
    print(Sys.time())
    print("Model 2a - Gender == 2")
  }
}

qval_gen2 <- results_gen2 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_gen2 <- results_gen2 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_gen2, by = "dnam_id") %>% subset(., qval_gen2 <= .05)
write.table(resq_gen2, "Model2_Out/model2a_gen2_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


#Models 2b: Also stratify by mother’s ethnicity (instead of adjusted for)

base_mod2b_eth1 <- lm(eduma ~ AgeMom + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK, data = ph2[ph2$Eth9 == 1,])
base_mod2b_eth7 <- lm(eduma ~ AgeMom + Bcell + CD4T + CD8T + Eos + Mono + Neu + NK, data = ph2[ph2$Eth9 == 7,])

results_base_mod2b_eth1 <- summary(base_mod2b_eth1)$coef[2:9,c(1,2,4)] %>% as.data.frame() %>% rownames_to_column('covariate') %>% mutate(Eth9 = 1)
results_base_mod2b_eth7 <- summary(base_mod2b_eth7)$coef[2:9,c(1,2,4)] %>% as.data.frame() %>% rownames_to_column('covariate') %>% mutate(Eth9 = 7)

results_base_mod2b <- rbind(results_base_mod2b_eth1, results_base_mod2b_eth7)
write.table(results_base_mod2b, "Model2_Out/model2b_results_base.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

## Model 2b - Stratified by Mother's Ethnicity, With DNA Methylation Data

# Ethnicity 1 - Initialize results matrices, Run Models

results_eth1 <- matrix(0, nrow=dim(betas)[1], ncol=4)
rownames(results_eth1) <- rownames(betas)
colnames(results_eth1) <- c("beta", "se", "t", "pvalue")

for(i in 1:nrow(betas)){
  x <- betas[i, colnames(betas) %in% ph2[ph2$Eth9 == 1, "sentrix"]]
  df <- ph2 %>% subset(., Eth9 == 1)
  mod <- lm(df$eduma ~ x + df$AgeMom + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_eth1[i,] <- summary(mod)$coef[2,]
  if(i %in% c(1, nrow(betas))) {
    print(Sys.time())
    print("Model 2b - Ethnicity == 1")
  }
}
#results_eth1 <- results_gen1 %>% as.data.frame() %>% subset(beta != 0)

qval_eth1 <- results_eth1 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_eth1 <- results_eth1 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth1, by = "dnam_id") %>% subset(., qval_eth1 <= .05)
write.table(resq_eth1, "Model2_Out/model2b_eth1_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Ethnicity 7 - Initialize results matrices, Run Models

results_eth7 <- matrix(0, nrow=dim(betas)[1], ncol=4)
rownames(results_eth7) <- rownames(betas)
colnames(results_eth7) <- c("beta", "se", "t", "pvalue")

for(i in 1:nrow(betas)){
  x <- betas[i, colnames(betas) %in% ph2[ph2$Eth9 == 7, "sentrix"]]
  df <- ph2 %>% subset(., Eth9 == 7)
  mod <- lm(df$eduma ~ x + df$AgeMom + df$Bcell + df$CD4T + df$CD8T + df$Eos + df$Mono + df$Neu + df$NK)
  results_eth7[i,] <- summary(mod)$coef[2,]
  if(i %in% c(1, nrow(betas))) {
    print(Sys.time())
    print("Model 2b - Ethnicity == 7")
  }
}

qval_eth7 <- results_eth7 %>% as.data.frame() %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq_eth7 <- results_eth7 %>% as.data.frame() %>% rownames_to_column('dnam_id') %>% left_join(., qval_eth7, by = "dnam_id") %>% subset(., qval_eth7 <= .05)
write.table(resq_eth7, "Model2_Out/model2b_eth7_results_qfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# results_gen1 <- qval()
#write.table(results_gen1, "~/Data/Simanek_BiB/DNAm/m1a_gen1.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

# results_gen2 <- qval()
#write.table(results_gen2, "~/Data/Simanek_BiB/DNAm/m1a_gen2.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

save.image("Model2_Out/Model2.RData")




























