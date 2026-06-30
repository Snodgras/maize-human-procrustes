library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)

cores = detectCores()
print(cores)
registerDoParallel(cores - 1)

args<-commandArgs(T)
language = args[1]
dataset = args[2]
# language = 'otomanguean'
# dataset = 'Haynie'

language_dataset = paste0(language, dataset)
print(language_dataset)

if (grepl('Haynie', language_dataset)) {
  elevation_mat = read.table(sprintf('../results/Haynie/%s_elevationSigned_distances.csv', language_dataset), sep = ',', row.names = 1, header = T)
  admix_mat = read.table(sprintf('../results/Haynie/%s_admixProp_distances.csv', language_dataset), row.names = 1, sep = ',', header = T)
  geographic_mat = read.table(sprintf('../results/Haynie/%s_haversine_distances.csv', language_dataset), row.names = 1, sep = ',', header = T)
  language_mat = read.table(sprintf('../results/Haynie/%s_levenshtein_language_distances.csv', language_dataset), row.names = 1, sep = ',', header = T)
  # branch_mat = read.table(sprintf('../results/Haynie/%s_branch_language_distances.csv', language), row.names = 1, sep = ',', header = T)
  # climate_mat = read.table(sprintf('../results/Haynie/%s_climate_distances.csv', language), row.names = 1, sep = ',', header = T)
} else {
  elevation_mat = fread(sprintf('../results/NL/%s_elevationSigned_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
  admix_mat = fread(sprintf('../results/NL/%s_admixProp_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
  geographic_mat = fread(sprintf('../results/NL/%s_haversine_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
  language_mat = fread(sprintf('../results/NL/%s_levenshtein_language_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
  # branch_mat = fread(sprintf('../results/NL/%s_branch_language_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
  # climate_mat = fread(sprintf('../results/NL/%s_climate_distances.csv', language)) %>% remove_rownames %>% column_to_rownames(var="V1")
}

language_mat = as.matrix(language_mat)
elevation_mat = abs(as.matrix(elevation_mat))
geographic_mat = as.matrix(geographic_mat)
admix_mat = abs(as.matrix(admix_mat))

print('imported matrices')

language_pca = prcomp(language_mat, rank = 5)
elevation_pca = prcomp(elevation_mat, rank = 5)
admix_pca = prcomp(admix_mat, rank = 5)
geographic_pca = prcomp(geographic_mat, rank = 5)
pca_predictors = cbind(geographic_pca$rotation, elevation_pca$rotation, admix_pca$rotation, language_pca$rotation)
colnames(pca_predictors) = outer(c('PC1', 'PC2', 'PC3', 'PC4', 'PC5'), c('geographic', 'elevation', 'admix', 'language'), paste, sep = '_')


print("Loaded PCs")

genotypes = fread(sprintf('../data/%s/%s.%s.LD-pruned.dosage.regression.vcf', language, language, dataset), sep = '\t', header = T)
genotypes_mat = as.matrix(genotypes)

print("Read in raw genotypes")

# ## for test
# test_genotypes = genotypes_mat[1:10, ]

snp_results = foreach(snp = c(1:nrow(genotypes_mat)), .combine = bind_rows, .errorhandling = 'remove') %dopar% {
		result = tryCatch({

			if(snp %% 1000 == 0) {
				print(snp)
			}

			genotype = genotypes_mat[snp, ]

			lm_results = lm(genotype ~ pca_predictors)
			lm.summary <-summary(lm_results)

			
			},
			warning = function(war) {
				result = c(snp = snp, r.squared = NA, adj.r.squared = NA)
				},
			error = function(err) {
				results = c(snp = snp, r.squared = NA, adj.r.squared = NA)
				},
			finally = {
				results = c(snp = snp, r.squared = lm.summary$r.squared, adj.r.squared = lm.summary$adj.r.squared)
				})
		return(results)
}

print("Finished regression")


write.csv(snp_results, sprintf('../results/regression_%s.csv', language_dataset), quote = F, row.names = F)

print("Wrote results")

# maize_genotypes<-read_tsv(args[2])
# maize_genotypes<-select(maize_genotypes, ID, contains("SEEDGWAS"))
# #transpose the genotypes so the SNPs are columns and individuals are rows and conserve names
# print(paste("Transposing genotypes...", Sys.time()))
# t_maize_genotypes<-t(select(maize_genotypes, contains("GWAS")))  %>% as_tibble(rownames = NA) %>% rownames_to_column()
# colnames(t_maize_genotypes)<-c("sample_id",maize_genotypes$ID)
# print(paste("Joining transposed genotypes to PCs...", Sys.time()))
# t_maize_genotypes<-mutate(t_maize_genotypes, sample_id = str_split(sample_id, "\\.", simplify = T)[,1]) %>% left_join(., maize_wHumanPCcentroids, by="sample_id")


