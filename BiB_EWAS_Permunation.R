args <- commandArgs(trailingOnly = TRUE)
model.cl <-  as.numeric(args[1])
us.cl <- as.logical(args[2])
gender.cl <- as.numeric(args[3])
ethnicity.cl <-  as.numeric(args[4])
filename.cl <- as.character(args[5])

### BiB EWAS Permutation Model 3

library(car)
library(data.table)
library(dplyr)
library(tidyverse)
library(qvalue)

dyn.load("VT_perms.so")
source("Perm_Function.R")

load("betas_cpg_sub.RData")
load("CellCountsBloodgse35069Complete.Rdata")
ph <- read.csv("bibgwas_16JAN2020.csv")  %>% subset(., !(sentrix %in% c("201057150183_R01C01", "201096090153_R05C01", "200992330056_R05C01", "201046290095_R08C01")))

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
M_val_sub <- M_val_ID_cpg[, match(ph2$sentrix, colnames(M_val_ID_cpg))]

perm_model3 <- function(model, us.strata = "FALSE", gender.strata = NA, ethnicity.strata = NA){
	
	filepath <- grep(paste0("/Model", model, "/CpGs$"), list.dirs(), value = TRUE, ignore.case = TRUE)
  
	if(us.strata == "TRUE"){
  		file <- grep("Total", list.files(path = filepath), value = TRUE)
  	}
  
  	if(gender.strata %in% c(1,2)){
  		file <- grep(paste0("Gen", gender.strata), list.files(path = filepath), value = TRUE)
  	}
  
  	if(ethnicity.strata %in% c(1,7)){
  		file <- grep(paste0("Race", ethnicity.strata), list.files(path = filepath), value = TRUE)
  	}
  
  	for(i in file){
  		cpg_list <- read.table(paste(filepath, i, sep = "/"), header = TRUE) #%>% mutate(., DNAm_ID = as.character(DNAm_ID))
  		cpg_list$level <- gsub("^.*\\.", "", cpg_list$DNAm_ID) %>% as.numeric()
  		cpg_list$DNAm_ID <- gsub("(.*)[.](.*)", "\\1", cpg_list$DNAm_ID)
  	}

  #########
  
  	if(model == 5){
    	main <- "quin07"
    	ref.level <- 1
  	}
  	if(model == 10){
    	main <- "tenown3c"
    	ref.level <- 3
    	ph2 <- ph2 %>% subset(., !is.na(tenown3c)) %>% mutate(., tenown3c = relevel(tenown3c, ref = "3"))
  	}
  
  	p.ref <- matrix(nrow = length(cpg_list$DNAm_ID), ncol = 6)
  	index.temp <- which(!is.na(ph2$eduma) & !is.na(ph2$BMImom))
  	
  	if(us.strata == "TRUE"){
  		ph3 <-  ph2[index.temp,]
  	}
  	if(gender.strata %in% c(1,2)){
  		ph3 <-  ph2[index.temp,] %>% subset(., gender == gender.strata)
  	}
  	if(ethnicity.strata %in% c(1,7)){
  		ph3 <-  ph2[index.temp,] %>% subset(., Eth9 == ethnicity.strata)
  	}
  	
  	
	start <- proc.time()
  	for(b in 1:length(cpg_list$DNAm_ID)){
  
		if(us.strata == "TRUE" | is.na(gender.strata) & is.na(ethnicity.strata)){
			y.perm <- M_val_sub[cpg_list$DNAm_ID[b],]
			y.perm <- y.perm[names(y.perm) %in% ph3$sentrix]
			y.perm <- y.perm[match(ph3$sentrix, names(y.perm))]
  		} else {
  			y.perm <- M_val_sub[cpg_list$DNAm_ID[b],]
			y.perm <- y.perm[names(y.perm) %in% ph3$sentrix]
			y.perm <- y.perm[match(ph3$sentrix, names(y.perm))]
			
  		}
		
		if(ethnicity.strata %in% c(1,7)){
			mod_ref <- lm(y.perm ~ ph3[[main]] + ph3$AgeMom + relevel(ph3$eduma, ref = "4") + relevel(ph3$smkpreg, ref = "2") + relevel(ph3$alpreg, ref = "0") + ph3$BMImom + ph3$pregwks + ph3$Bwt + ph3$Bcell + ph3$CD4T + ph3$CD8T + ph3$Eos + ph3$Mono + ph3$Neu + ph3$NK)
			p.ref[b, 1:4] <- summary(mod_ref)$coef[cpg_list$level[b] + 1,]
		} else {
		
			mod_ref <- lm(y.perm ~ ph3[[main]] + ph3$AgeMom + relevel(ph3$Eth9, ref = "1") + relevel(ph3$eduma, ref = "4") + relevel(ph3$smkpreg, ref = "2") + relevel(ph3$alpreg, ref = "0") + ph3$BMImom + ph3$pregwks + ph3$Bwt + ph3$Bcell + ph3$CD4T + ph3$CD8T + ph3$Eos + ph3$Mono + ph3$Neu + ph3$NK)
			p.ref[b, 1:4] <- summary(mod_ref)$coef[cpg_list$level[b] + 1,]
		}
		
		temp <- model.matrix(mod_ref)
		x.new <- temp[,cpg_list$level[b] + 1]
		z.new <- temp[, -c(1, cpg_list$level[b]+1)]
	
		mod_perm <- perm.func(xvar = x.new, y = y.perm, z = z.new, max.perms=10000000, adaptive.perms.cutoff = 10000000)
	
		p.ref[b, 5] <- mod_perm$p.value
		
		p.ref[b, 6] <- cpg_list$level[b]
		print(b)
	}
	
	end <- proc.time()
	dimnames(p.ref) <- list(cpg_list$DNAm_ID, c("Beta", "Std.Error", "T.value", "P.ref", "P.perm", "Contrast.Level"))
  	return(p.ref)
}


temp <- perm_model3(model.cl, us.cl, gender.cl, ethnicity.cl)
write.table(temp, paste0("Permutations/Model", model.cl, "/Output/", filename.cl, ".txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

dyn.unload("VT_perms.so")















