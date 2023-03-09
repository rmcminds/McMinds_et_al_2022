
nthreads <- 7
raw_data_prefix <- path.expand('raw_data/20221215_primate_allometry/')
output_prefix <- path.expand('outputs/primates_20230304/')

data_species <- c('Callithrix_jacchus', 'Daubentonia_madagascariensis', 'Homo_sapiens', 'Lemur_catta', 'Macaca_mulatta', 'Microcebus_murinus', 'Papio_anubis', 'Pongo_abelii', 'Sapajus_appella')

species_strings <- c(Callithrix_jacchus='ENSCJA', Homo_sapiens='ENS', Microcebus_murinus='ENSMIC', Macaca_mulatta='ENSMMU', Papio_anubis='ENSPAN', Pongo_abelii='ENSPPY')

cat('Importing cladogram\n')
species_tree <- ape::read.tree(file.path(raw_data_prefix, 'primates.newick'))
species_tree_norm <- species_tree
species_maxtime <- max(phytools::nodeHeights(species_tree_norm))
species_tree_norm$edge.length <- species_tree_norm$edge.length / species_maxtime
species_tree_norm <- reorder(species_tree_norm, order='postorder')
##

## import and wrangle sample metadata
sample_data <- read.table(file.path(raw_data_prefix, 'sample_data.txt'), sep='\t', header=T)
# summarize pooled Microcebus murinus samples into one pseudo-sample using geometric mean for numeric values (all should be positive and this makes sense for size and time)
sample_data <- rbind(sample_data, apply(sample_data[sample_data$Genus == 'Microcebus',], 2, function(x) {temp <- unique(x); if(length(temp)==1) return(temp) else return(exp(mean(log(as.numeric(temp)))))}))
sample_data$Animal.ID[nrow(sample_data)] <- 'mmurPool'
#
sample_data$genus_species <- paste(sample_data$Genus, sample_data$Species, sep='_')
##

## import counts
countfiles <- Sys.glob(file.path(output_prefix, '/03_generate_counts/*/quant.sf'))
samples <- sapply(basename(dirname(countfiles)), \(x) paste(strsplit(x,'_')[[1]][3:4], collapse='_'))
names(countfiles) <- names(samples)
sample_data_filt <- sample_data[sample_data$Animal.ID %in% sub('_.*', '', samples),]

cat('Retrieving ensembl genes\n')
human_globin_gene_ids <- c('ENSG00000206172', 'ENSG00000188536', 'ENSG00000244734', 'ENSG00000229988', 'ENSG00000223609', 'ENSG00000213931', 'ENSG00000213934', 'ENSG00000196565', 'ENSG00000206177', 'ENSG00000086506', 'ENSG00000130656', 'ENSG00000206178') # https://doi.org/10.1038/s41598-020-62801-6 Table S1
numtries <- 5
for(i in 1:numtries) {
  try({
    human_globin_peptide_ids <- biomaRt::getBM(attributes = 'ensembl_peptide_id', 
                                               filters = 'ensembl_gene_id', 
                                               values = human_globin_gene_ids, 
                                               mart = biomaRt::useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', version='Ensembl Genes 109'), 
                                               uniqueRows = TRUE)[,1]
    break
  }, silent = FALSE)
}

orthologs <- read.table(file.path(output_prefix,'02_find_orthologs','of_out','Results_Mar08','Orthogroups','Orthogroups.tsv'), header=TRUE, row.names = 1, sep='\t')
orthologs_1_1 <- orthologs[!apply(orthologs, 1, \(x) any(grepl(',',x) | (nchar(x) == 0))),]
orthologs_1_1[,colnames(orthologs_1_1) != 'Homo_ensembl'] <- apply(orthologs_1_1[,colnames(orthologs_1_1) != 'Homo_ensembl'], 2, \(x) sapply(x, \(y) paste(strsplit(y,'.',fixed=TRUE)[[1]][1:2],collapse='.')))
orthologs_1_1[,'Homo_ensembl'] <- sapply(orthologs_1_1[,'Homo_ensembl'], \(x) strsplit(x,'.',fixed=TRUE)[[1]][1])
orthologs_1_1 <- orthologs_1_1[!orthologs_1_1[,'Homo_ensembl'] %in% human_globin_peptide_ids,]

for(i in 1:numtries) {
  try({
  ensembl2ext <- biomaRt::getBM(attributes = c('ensembl_peptide_id', 'ensembl_gene_id', 'external_gene_name', 'go_id'), 
                                filters = 'ensembl_peptide_id', 
                                values = orthologs_1_1[,'Homo_ensembl'], 
                                mart = biomaRt::useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', version='Ensembl Genes 109'), 
                                uniqueRows = TRUE)
  break
  }, silent = FALSE)
}

## feels weird to round decimals for poisson error, but since data were counts at one point, zeros are possible, so doesn't make sense to just log-transform and use gaussian error; and error probably still scales like poisson such that low counts are less meaningful (and the decimal values rounded off would contribute negligible information) (cite Z1000 tximport paper)
raw_txi <- lapply(data_species, \(x) {
  cat(paste0('Importing ', x, ' counts\n'))
  cfiles <- countfiles[grep(x, names(countfiles), ignore.case = TRUE)]
  tnames <- unique(do.call(rbind, lapply(cfiles, read.table, header=TRUE, sep='\t'))$Name)
  tx2gene <- data.frame(transcript=tnames, gene=sapply(tnames, \(y) paste(strsplit(y,'.',fixed=TRUE)[[1]][1:2],collapse='.')))
  txi <- tximport::tximport(cfiles, type = 'salmon', tx2gene = tx2gene)
  return(txi)
})
names(raw_txi) <- data_species
##

## import reference data for species average body sizes
body_size_ref <- read.csv(file.path(raw_data_prefix, 'gyz043_suppl_Supplement_Data.csv'))
body_size_ref$genus_species <- paste(body_size_ref$genus, body_size_ref$species, sep='_')
body_size_ref$genus_species[body_size_ref$genus_species == 'Cebus_apella'] <- 'Sapajus_apella' # misspelled in some versions of previous analysis; careful 
##

## calculate differences of individual sizes from species means and log-transform
sample_data_filt$body_mass_log <- log(as.numeric(sample_data_filt$Body.Mass..g.))

body_mass_log_sp_means <- sapply(unique(sample_data_filt$genus_species), \(x) log(as.numeric(body_size_ref$Mean_body_mass_g[body_size_ref$genus_species == x])))

sample_data_filt$body_mass_log_sp <- body_mass_log_sp_means[sample_data_filt$genus_species]
sample_data_filt$body_mass_log_diff <- sample_data_filt$body_mass_log - sample_data_filt$body_mass_log_sp
##

## calculate a theoretical optimum body size on which to center the model (somewhat arbitrary but could help interpret 'main effects')
penalized_likelihood_ou <- function(x) {
  fit <- ape::compar.ou(body_mass_log_sp_means, species_tree_norm, alpha = exp(x))
  return(fit$deviance - 2 * dlnorm(fit$para[2,1], mean(body_mass_log_sp_means), sd(body_mass_log_sp_means), log = TRUE))
} ## optimizing based on the compar.ou deviance alone produced crazy estimates for the optimum (>15000), so this essentially contains it with an 'empirical bayes'-like prior
opt_ou_alpha <- exp(optim(0, penalized_likelihood_ou, method = 'BFGS')$par)
body_mass_opt <- ape::compar.ou(body_mass_log_sp_means, species_tree_norm, alpha = opt_ou_alpha) ## by using an estimate of the evolutionary 'optimum', all main effects (such as LPS) can be interpreted as effects at this optimum, rather than at the arbitrary mean of all samples
body_mass_log_center <- body_mass_opt$para[2,1]
##

## standardize masses for model input
body_mass_log_sd <- sd(sample_data_filt$body_mass_log_sp)
body_mass_log_diff_sd <- sd(sample_data_filt$body_mass_log_diff)

sample_data_filt$body_mass_log_sp_std <- (sample_data_filt$body_mass_log_sp - body_mass_log_center) / body_mass_log_sd
sample_data_filt$body_mass_log_diff_std <- sample_data_filt$body_mass_log_diff / body_mass_log_diff_sd
##

ortholog_genes <- unlist(parallel::mclapply(orthologs_1_1[,'Homo_ensembl'], \(x) ensembl2ext$ensembl_gene_id[match(x, ensembl2ext$ensembl_peptide_id)], mc.cores = nthreads))

ortholog_mat <- orthologs_1_1[,names(raw_txi)]
txi_ortho <- list(abundance = raw_txi[[1]]$abundance[ortholog_mat[,1],], 
                  counts    = raw_txi[[1]]$counts[ortholog_mat[,1],],
                  length    = raw_txi[[1]]$length[ortholog_mat[,1],],
                  countsFromAbundance = raw_txi[[1]]$countsFromAbundance)
for(i in 2:length(data_species)) {
  txi_ortho$abundance <- cbind(txi_ortho$abundance, raw_txi[[i]]$abundance[ortholog_mat[,i],])
  txi_ortho$counts <- cbind(txi_ortho$counts, raw_txi[[i]]$counts[ortholog_mat[,i],])
  txi_ortho$length <- cbind(txi_ortho$length, raw_txi[[i]]$length[ortholog_mat[,i],])
} 
rownames(txi_ortho$abundance) <- rownames(txi_ortho$counts) <- rownames(txi_ortho$length) <- ensembl2ext$ensembl_gene_id[match(orthologs_1_1[match(rownames(txi_ortho$length), orthologs_1_1[,1]), 'Homo_ensembl'], ensembl2ext$ensembl_peptide_id)]

## import immune annotations

## all ensembl identifiers corresponding to immune annotations
modules <- read.table(file.path(raw_data_prefix, 'immune_modules.txt'), sep='\t', quote='', header=T)
ensembl2ext$module <- modules$IIG_class[match(ensembl2ext$external_gene_name, modules$HGNC_symbol)]
ensembl2ext$module[is.na(ensembl2ext$module)] <- 'no_module_annotation'

immuneGenes <- unique(ensembl2ext$ensembl_gene_id[ensembl2ext$module != 'no_module_annotation'])

## pool immune-annotated 1:1 orthologous genes
geneGroups <- rownames(txi_ortho$abundance)
geneGroups[geneGroups %in% immuneGenes] <- 'immune'
txi_ortho_immune_pool <- txi_ortho
txi_ortho_immune_pool$abundance <- rowsum(txi_ortho$abundance, geneGroups)
txi_ortho_immune_pool$counts <- rowsum(txi_ortho$counts, geneGroups)
txi_ortho_immune_pool$length <- rowsum(txi_ortho$abundance * txi_ortho$length, geneGroups) / txi_ortho_immune_pool$abundance ## weighted arithmetic mean length
txi_ortho_immune_pool$length[rownames(txi_ortho_immune_pool$length) != 'immune',] <- txi_ortho$length[rownames(txi_ortho_immune_pool$length)[rownames(txi_ortho_immune_pool$length) != 'immune'],]

## get normalization factors
des_ortho_pool <- DESeq2:::DESeqDataSetFromTximport(txi_ortho_immune_pool, data.frame(Int=rep(1,ncol(txi_ortho_immune_pool$counts))), ~1)
des_ortho_pool <- DESeq2::estimateSizeFactors(des_ortho_pool)

## model total differential expression of immune genes
dat_ortho <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des_ortho_pool)), '_'), \(x) x[[3]]), 
                        species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des_ortho_pool)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                        treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des_ortho_pool)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                        norm_factor = log(DESeq2::normalizationFactors(des_ortho_pool)['immune',]),
                        count       = DESeq2::counts(des_ortho_pool)['immune',])
levels(dat_ortho$species)[levels(dat_ortho$species) == 'Sapajus_appella'] <- 'Sapajus_apella'
dat_ortho <- dat_ortho[dat_ortho$individual %in% sample_data_filt$Animal.ID,]
dat_ortho$body_mass_log_sp_std <- sapply(dat_ortho$individual, function(z) sample_data_filt$body_mass_log_sp_std[sample_data_filt$Animal.ID==z])
dat_ortho$body_mass_log_diff_std <- sapply(dat_ortho$individual, function(z) sample_data_filt$body_mass_log_diff_std[sample_data_filt$Animal.ID==z])

## define phyr formula
phy_formula <- count ~ offset(norm_factor) + body_mass_log_sp_std * treatment + body_mass_log_diff_std + body_mass_log_diff_std:treatment + (1 | individual) + (1 | species__) + (1 | species@treatment) + (1 | species__@treatment)

## define linear combinations of effects of interest
lc1 <- INLA::inla.make.lincomb('body_mass_log_sp_std'              = body_mass_log_sd, 
                               'body_mass_log_sp_std:treatmentLPS' = body_mass_log_sd) 
names(lc1) = "induced_evo_allometry"

lc2 <- INLA::inla.make.lincomb('body_mass_log_diff_std'              = body_mass_log_diff_sd, 
                               'treatmentLPS:body_mass_log_diff_std' = body_mass_log_diff_sd)
names(lc2) = "induced_intra_allometry"

## scaling by first sd to get both on same scale, then weighting the combination by the variance of each (so if intra-specific size variation is larger than inter-specific size variation, then the intra-specific coefficient dominates)
lc3 <- INLA::inla.make.lincomb('body_mass_log_diff_std' = body_mass_log_diff_sd^3 / (body_mass_log_diff_sd^2 + body_mass_log_sd^2), 
                               'body_mass_log_sp_std'   = body_mass_log_sd^3 / (body_mass_log_diff_sd^2 + body_mass_log_sd^2))
names(lc3) = "overall_allometry"

lc4 <- INLA::inla.make.lincomb('treatmentLPS:body_mass_log_diff_std' = body_mass_log_diff_sd^3 / (body_mass_log_diff_sd^2 + body_mass_log_sd^2), 
                               'body_mass_log_sp_std:treatmentLPS'   = body_mass_log_sd^3 / (body_mass_log_diff_sd^2 + body_mass_log_sd^2))
names(lc4) = "overall_allometric_response"

lc5 <- INLA::inla.make.lincomb('body_mass_log_diff_std'              = body_mass_log_diff_sd^3 / (body_mass_log_diff_sd^2 + body_mass_log_sd^2), 
                               'treatmentLPS:body_mass_log_diff_std' = body_mass_log_diff_sd^3 / (body_mass_log_diff_sd^2 + body_mass_log_sd^2), 
                               'body_mass_log_sp_std'                = body_mass_log_sd^3 / (body_mass_log_diff_sd^2 + body_mass_log_sd^2), 
                               'body_mass_log_sp_std:treatmentLPS'   = body_mass_log_sd^3 / (body_mass_log_diff_sd^2 + body_mass_log_sd^2))
names(lc5) = "induced_overall_allometry"

cat('Modeling immune orthologs\n')
ortho_immune_res <- phyr::pglmm(formula       = phy_formula, 
                                family        = "poisson", 
                                cov_ranef     = list(species = species_tree_norm), 
                                data          = dat_ortho, 
                                bayes         = TRUE, 
                                bayes_options = list(lincomb = c(lc1,lc2,lc3,lc4,lc5)))
ortho_immune_restab <- ortho_immune_res$inla.model$summary.lincomb.derived
##

## add method to pool ALL immune genes - not just 1:1 orthologs. will need to use full orthologs table to find any rows that contain an immune annotated human gene

save.image(file.path(output_prefix,'02_phyr_int_results.RData'))

## fit per-gene models

## get per-observation normalization factors
des <- DESeq2:::DESeqDataSetFromTximport(txi_ortho, data.frame(Int=rep(1,ncol(txi_ortho$counts))), ~1)
des <- DESeq2::estimateSizeFactors(des)

cat('Fitting models\n')
fits <- parallel::mclapply(rownames(des), function(x) {
  cat(round(which(x == rownames(des)) / nrow(des),2) * 100, '%: ', x,'\n')
  dat <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des)), '_'), \(x) x[[3]]), 
                    species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                    treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                    norm_factor = log(DESeq2::normalizationFactors(des)[x,]),
                    count       = DESeq2::counts(des)[x,])
  levels(dat$species)[levels(dat$species) == 'Sapajus_appella'] <- 'Sapajus_apella'
  dat <- dat[dat$individual %in% sample_data_filt$Animal.ID,]
  dat$body_mass_log_sp_std <- sapply(dat$individual, function(z) sample_data_filt$body_mass_log_sp_std[sample_data_filt$Animal.ID==z])
  dat$body_mass_log_diff_std <- sapply(dat$individual, function(z) sample_data_filt$body_mass_log_diff_std[sample_data_filt$Animal.ID==z])
  
  if(sd(DESeq2::counts(des)[x,]) > 0) {

    fit <- tryCatch(phyr::pglmm(formula       = phy_formula, 
                                family        = "poisson", 
                                cov_ranef     = list(species = species_tree_norm), 
                                data          = dat, 
                                bayes         = TRUE, 
                                bayes_options = list(lincomb = c(lc1,lc2,lc3,lc4,lc5))),
                    error = function(e) NA)
    if(!all(is.na(fit))) {
      coefs <- c(phyr::fixef(fit)$Value, fit$inla.model$summary.lincomb.derived[1:5,'mean'])
      coef_cis <- rbind(fit$B.ci, LPS_beta = fit$inla.model$summary.lincomb.derived[1:5, c('0.025quant','0.975quant')])
      fit$inla.model <- NULL
    } else {
      coefs <- NA
      coef_cis <- NA
    }
  } else {
    fit <- NA
    coefs <- NA
    coef_cis <- NA
  }
  return(list(data=dat, fit=fit, coefs=coefs, coef_cis=coef_cis))
}, mc.cores = nthreads)
names(fits) <- rownames(des)

save.image(file.path(output_prefix,'02_phyr_results.RData'))

