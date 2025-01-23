library(tidyverse)

#Test using random data of right type, but no basis in biology
pilot_vcf<-tibble(SNP_ID = c(paste0("SNP",1:50)),
                  Mex_AF = runif(n=50,min=0,max=1),
                  IND_1 = round(runif(n=50, min=0,max=2)),
                  IND_2 = round(runif(n=50, min=0,max=2)),
                  IND_3 = round(runif(n=50, min=0,max=2)),
                  IND_4 = round(runif(n=50, min=0,max=2)),
                  IND_5 = round(runif(n=50, min=0,max=2)),
                  IND_6 = round(runif(n=50, min=0,max=2)),
                  IND_7 = round(runif(n=50, min=0,max=2)),
                  IND_8 = round(runif(n=50, min=0,max=2)),
                  IND_9 = round(runif(n=50, min=0,max=2)),
                  IND_10 = round(runif(n=50, min=0,max=2)))
pilot_admix<-tibble(Individual = c("IND_1","IND_2","IND_3","IND_4","IND_5","IND_6","IND_7","IND_8","IND_9","IND_10"),
                    alpha = runif(n=10,min=0,max=1))

(pilot_vcf[,3] - 2*pilot_vcf[,2]*pull(pilot_admix[which(pilot_admix[,1] == "IND_1"),2]))/(1-pull(pilot_admix[which(pilot_admix[,1] == "IND_1"),2]))
#test first entry:
(1 - 2*0.172*0.00514)/(1-0.00514)

pilot_vcf_transformed<-select(pilot_vcf, SNP_ID)
for(i in pilot_admix$Individual){
  alpha<-pilot_admix %>% filter(Individual == i) %>% select(alpha) %>% pull()
  pilot_vcf_transformed<-add_column(pilot_vcf_transformed, !!i := pull((pilot_vcf[,i] - 2*pilot_vcf[,2]*alpha)/(1-alpha)))
}
remove(ls(pattern = "pilot"))

# Test on a subset of the real data

#Final draft
