#########################################################################
## Title: Model 4 EWAS
## Version: 1.1 
## Note: Originally Model 3 - Version 1, but shifted to Model 4
##		see EWAS_analysis_plan_v3
## Author: Regina Manansala
## Date: 07-November-2019
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

# Model 4c 

#Initialize results list
results <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){ #
	x <- M_val[i,]
	mod <- lm(x ~ relevel(ph2$bn_mtben, ref = "1") + ph2$AgeMom + relevel(ph2$Eth9, ref = "1") + relevel(ph2$eduma, ref = "4") + relevel(ph2$smkpreg, ref = "2") + relevel(ph2$alpreg, ref = "0") + ph2$BMImom + ph2$Bcell + ph2$CD4T + ph2$CD8T + ph2$Eos + ph2$Mono + ph2$Neu + ph2$NK)
	results[[i]] <- summary(mod)$coef[2,]
	if(i %in% c(1, nrow(M_val))) {
		print(Sys.time())
		print("Model 4c")
	}
}

# Compile results list into data frame
names(results) <- rownames(M_val)
big_data <- do.call(rbind, results) %>% as.data.frame()
names(big_data)[1:4] <- c("beta", "se", "t", "pvalue")

# Calculate qvalue and filter by pvalue
qval <- big_data %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq <- big_data %>% rownames_to_column('dnam_id') %>% left_join(., qval, by = "dnam_id") %>% subset(., pvalue <= 5e-06)
top50 <- big_data %>% rownames_to_column('dnam_id') %>% left_join(., qval, by = "dnam_id") %>% arrange(pvalue)
# write.table(resq, "Model4c_results_pfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
# write.table(top50, "Model4c_top50.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Export model results to RData file
save(list = ls(pattern = "results|big_data|qval|resq"), file = "Model4c_unstrat.RData")

#########################################################################
#########################################################################
#########################################################################

# Model 4d

#Initialize results list
results <- list()

# Run model for each cpg site
for(i in 1:nrow(M_val)){ #nrow(M_val)
	x <- M_val[i,]
	mod <- lm(x ~ relevel(ph2$bn_mtben, ref = "1") + ph2$AgeMom + relevel(ph2$Eth9, ref = "1") + relevel(ph2$eduma, ref = "4") + relevel(ph2$smkpreg, ref = "2") + relevel(ph2$alpreg, ref = "0") + ph2$BMImom + ph2$pregwks + ph2$Bcell + ph2$CD4T + ph2$CD8T + ph2$Eos + ph2$Mono + ph2$Neu + ph2$NK)
	results[[i]] <- summary(mod)$coef[2,]
	if(i %in% c(1, nrow(M_val))) {
		print(Sys.time())
		print("Model 4d")
	}
}

# Compile results list into data frame
names(results) <- rownames(M_val)
big_data <- do.call(rbind, results) %>% as.data.frame()
names(big_data)[1:4] <- c("beta", "se", "t", "pvalue")

# Calculate qvalue and filter by pvalue
qval <- big_data %>% select(., pvalue) %>% qvalue() %>% `[[`("qvalues") %>% rename(qval = pvalue) %>% rownames_to_column('dnam_id')
resq <- big_data %>% rownames_to_column('dnam_id') %>% left_join(., qval, by = "dnam_id") %>% subset(., pvalue <= 5e-06)
top50 <- big_data %>% rownames_to_column('dnam_id') %>% left_join(., qval, by = "dnam_id") %>% arrange(pvalue)
# write.table(resq, "Model4d_results_pfilt.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
# write.table(top50, "Model4d_top50.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Export model results to RData file
save(list = ls(pattern = "results|big_data|qval|resq"), file = "Model4d_unstrat.RData")
