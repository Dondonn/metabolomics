library(magrittr)
library(ggplot2)
library(stringr)
library(dplyr)
library(ggbeeswarm)


setwd("D:\\Study thingie\\FACH\\S5 WS2021\\GOBI\\metabolomics")
meta_data <- read.table("MetabolomicsAndGenetics_4Praktikum.data.tsv",h=T)

meta_data$GENDER <- ifelse(meta_data$GENDER == "male",1,2)
meta_concentration <- log2(meta_data[4:154])
meta_concentration[,-1] <- lapply(meta_concentration[,-1],
                             function(x) replace(x,abs(scale(x))>3,NA))

meta_data_v2 <- data.frame(meta_data)
meta_data_v2[4:154] <- meta_concentration
ncol(meta_data_v2)
SNPs <- meta_data_v2[(155:164)]

head(meta_data_v2[, which(sum(is.na(SNPs)) > 500)])
SNPs <- SNPs[, which(sum(is.na(SNPs)) > 500)]

meta_names <- colnames(meta_data_v2[4:154])
snp_names <- colnames(meta_data_v2[155:164])

number_pairs = length(meta_names)*length(snp_names)
number_rows = nrow(meta_data_v2)

mat = matrix(ncol = number_pairs, nrow = number_rows)
df=data.frame(mat)


#create colnames with pairs and meta

list_names <- matrix(ncol = number_pairs, nrow = number_rows)
meta_snp_data=data.frame(list_names)
meta_snp_summary =  data.frame(list_names)

big_name <- c()
for(meta in meta_names){
  for (snp in snp_names){
  
    big_name <- c(big_name,paste(meta,snp,sep = "-"))
  }
}

#create table to store est sd and p-value for SNP-met effect
snp_met_effect <- matrix(ncol = 5, nrow = 0)
meta_snp_effects=data.frame(snp_met_effect)


#create colnames
colnames(meta_snp_data) <- big_name


meta_data_v2$AGE <- as.numeric(meta_data_v2$AGE)


  for(meta in 1:length(meta_names)){
    for (snp in 1:length(snp_names)){
      
      col_name <- paste(meta_names[meta],snp_names[snp],sep = "-")
      
      meta_snp_summary <- summary(lm(meta_data_v2[,meta_names[meta]]~meta_data_v2[,"GENDER"]+
                                        meta_data_v2[,"AGE"]+
                                        meta_data_v2[,snp_names[snp]]))
      non_missing <- sum(!is.na(meta_data_v2[,meta_names[meta]]))
      result_snp <- meta_snp_summary$coefficients[4,] #1 for estimate, 2 for std error, 4 for p-value  
      meta_snp_effects <- rbind(meta_snp_effects, c("pairs" = col_name,"estimate"= result_snp[1], "std_error"= result_snp[2],"p_value" =result_snp[4], "N" = non_missing))

    }
  }


colnames(meta_snp_effects) <- c("pairs","estimate","std_error","p_value", "N")



meta_snp_effects$p_value <- as.numeric(meta_snp_effects$p_value)
meta_snp_effects$N <- as.numeric(meta_snp_effects$N)


sum((meta_snp_effects$p_value < (0.05/meta_snp_effects$N)) == TRUE )

cleaned_meta_snp_effects <- meta_snp_effects[meta_snp_effects$p_value < (0.05/meta_snp_effects$N),]



cleaned_meta_snp_effects$snp_id <-as.numeric(gsub("^[^r]*rs", "\\1", cleaned_meta_snp_effects$pairs))

ordered_meta_snp_effects <- cleaned_meta_snp_effects[order(cleaned_meta_snp_effects$snp_id, cleaned_meta_snp_effects$p_value),]

nrow(ordered_meta_snp_effects)
#cleaned_meta_snp_effects %>% group_by(snp_id) %>% arrange(snp_id,p_value) 

#3. Plotting if results
gsub("\\-.*", "", cleaned_meta_snp_effects$pairs)
ordered_meta_snp_effects$meta <- gsub("\\-.*", "", ordered_meta_snp_effects$pairs)
View(ordered_meta_snp_effects)

best_related_pairs <- ordered_meta_snp_effects %>% group_by(snp_id) %>% top_n(-1,p_value)
the_best <- best_related_pairs[which.min(best_related_pairs$p_value),]


gsub(".*-", "", ordered_meta_snp_effects$pairs)
ordered_meta_snp_effects$snp_id <- gsub(".*-", "", ordered_meta_snp_effects$pairs)
best_pair <- meta_data_v2 %>% select(GENDER,SNP_rs2066938,C4)
best_pair$SNP_rs2066938 <- as.character(best_pair$SNP_rs2066938)
best_pair$GENDER <- as.character(best_pair$GENDER)

ggplot(best_pair %>% filter(!is.na(SNP_rs2066938)), aes(SNP_rs2066938, C4))+
  geom_boxplot(na.rm = T)+
  facet_wrap(~GENDER)

ggplot(best_pair %>% filter(!is.na(SNP_rs2066938)), aes(SNP_rs2066938, C4))+
  geom_quasirandom(na.rm = T, aes(color = GENDER))

#4. Biological context of results/ resporing of results
#C4 (Biocrates, P = 2.9×10???44) = butyrylcarnitine (Metabolon, P = 1.75×10???114)

