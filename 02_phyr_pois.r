## modify using logic tested in 05b_run_batches_sparseOU_gp.r to ensure species means are distinct from individual sizes, and other logic therein

nthreads <- 6
raw_data_prefix <- path.expand('raw_data/20221215_primate_allometry/')
output_prefix <- path.expand('outputs/primates_20230224/')

species_strings <- c(Callithrix_jacchus='ENSCJA', Homo_sapiens='ENS', Microcebus_murinus='ENSMIC', Macaca_mulatta='ENSMMU', Papio_anubis='ENSPAN', Pongo_abelii='ENSPPY')

## download and normalize species chronogram
cat('Retrieving cladogram\n')
species_tree <- datelife::summarize_datelife_result(datelife::get_datelife_result(c('Callithrix jacchus', 'Homo sapiens', 'Macaca mulatta', 'Microcebus murinus', 'Papio anubis', 'Pongo abelii')), summary_format='phylo_biggest')

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
countfiles <- Sys.glob(file.path(output_prefix, '/01_generate_counts/*/quant.sf'))
samples <- sapply(basename(dirname(countfiles)), \(x) paste(strsplit(x,'_')[[1]][3:4], collapse='_'))
names(countfiles) <- names(samples)
sample_data_filt <- sample_data[sample_data$Animal.ID %in% sub('_.*', '', samples),]

martlist <- list(ENSCJA='cjacchus_gene_ensembl', ENS='hsapiens_gene_ensembl', ENSMIC='mmurinus_gene_ensembl', ENSMMU='mmulatta_gene_ensembl', ENSPAN='panubis_gene_ensembl', ENSPPY='pabelii_gene_ensembl')

genetrees <- ape::read.tree(file.path(output_prefix, '00_references/Compara.109.protein_default.newick'))
names(genetrees) <- paste('family', 1:length(genetrees), sep='_') ## arbitrary gene family names
genetrees <- lapply(genetrees, \(x) ape::keep.tip(x, grep(paste(paste0(species_strings,'P[[:digit:]]'),collapse='|'), x$tip.label)))
genetrees <- genetrees[!sapply(genetrees,is.null)]

cat('Retrieving ensembl genes\n')
human_globin_gene_ids <- c('ENSG00000206172', 'ENSG00000188536', 'ENSG00000244734', 'ENSG00000229988', 'ENSG00000223609', 'ENSG00000213931', 'ENSG00000213934', 'ENSG00000196565', 'ENSG00000206177', 'ENSG00000086506', 'ENSG00000130656', 'ENSG00000206178') # https://doi.org/10.1038/s41598-020-62801-6 Table S1
human_globin_peptide_ids <- biomaRt::getBM(attributes = 'ensembl_peptide_id', 
                                           filters = 'ensembl_gene_id', 
                                           values = human_globin_gene_ids, 
                                           mart = biomaRt::useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', version='Ensembl Genes 109'), 
                                           uniqueRows = TRUE)[,1]

globin_trees <- names(genetrees)[sapply(genetrees, \(x) any(human_globin_peptide_ids %in% x$tip.label))]
tax_globin_peptide_ids <- lapply(names(species_strings), \(x) {
  grep(paste0(species_strings[[x]],'P[[:digit:]]'), unlist(sapply(genetrees[globin_trees], \(y) y$tip.label)), value=TRUE)
})
names(tax_globin_peptide_ids) <- names(species_strings)

## feels weird to round decimals for poisson error, but since data were counts at one point, zeros are possible, so doesn't make sense to just log-transform and use gaussian error; and error probably still scales like poisson such that low counts are less meaningful (and the decimal values rounded off would contribute negligible information) (cite Z1000 tximport paper)
raw_txi <- lapply(names(species_strings), \(x) {
  cat(paste0('Retrieving ensembl genes for ', x, ' \n'))
  cfiles <- countfiles[grep(x, names(countfiles), ignore.case = TRUE)]
  tnames <- unique(do.call(rbind, lapply(cfiles, read.table, header=TRUE, sep='\t'))$Name)
  ensembl2ext <- biomaRt::getBM(attributes = c('ensembl_transcript_id_version','ensembl_gene_id','external_gene_name','ensembl_peptide_id','go_id'), 
                                filters = 'ensembl_transcript_id_version', 
                                values = tnames, 
                                mart = biomaRt::useMart('ENSEMBL_MART_ENSEMBL', martlist[[species_strings[[x]]]], version='Ensembl Genes 109'), 
                                uniqueRows = TRUE)
  tx2gene <- ensembl2ext[, c('ensembl_transcript_id_version', 'ensembl_gene_id')]
  globins <- ensembl2ext$ensembl_gene_id[ensembl2ext$ensembl_peptide_id %in% tax_globin_peptide_ids[[x]]]
  txi <- tximport::tximport(cfiles, type = 'salmon', tx2gene = tx2gene)
  txi$abundance <- txi$abundance[!rownames(txi$abundance) %in% globins,]
  txi$counts <- txi$counts[!rownames(txi$counts) %in% globins,]
  txi$length <- txi$length[!rownames(txi$length) %in% globins,]
  txi$ensembl2ext <- ensembl2ext
  return(txi)
})
names(raw_txi) <- names(species_strings)
##

## import reference data for species average body sizes
body_size_ref <- read.csv(file.path(raw_data_prefix, 'gyz043_suppl_Supplement_Data.csv'))
body_size_ref$genus_species <- paste(body_size_ref$genus, body_size_ref$species, sep='_')
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

## find 1:1 orthologs for size factor calcs and for later per-gene analyses 
cat('Finding orthologs\n')
ortholog_peptides <- parallel::mclapply(genetrees, function(tree) {
  nodes <- unique(tree$edge[,1])
  alldescs <- lapply(phangorn::Descendants(tree, nodes), function(y) tree$tip.label[y])
  names(alldescs) <- as.character(nodes)
  is.ortho <- sapply(alldescs, function(y) {
    all(sapply(species_strings, function(z) any(grepl(paste0(z,'P[[:digit:]]'), y)))) & length(y) == length(species_strings)
  })
  return(alldescs[is.ortho])
}, mc.cores = nthreads)
names(ortholog_peptides) <- names(genetrees)
ortholog_peptides <- lapply(unlist(ortholog_peptides,recursive=F), sort)

ensembl2ext_full <- do.call(rbind, lapply(raw_txi, \(x) x$ensembl2ext))
ortholog_genes <- parallel::mclapply(ortholog_peptides, \(x) ensembl2ext_full$ensembl_gene_id[match(x, ensembl2ext_full$ensembl_peptide_id)], mc.cores = nthreads)

ortholog_mat <- na.omit(do.call(rbind, ortholog_genes)[,c(1,4,2,3,5,6)]) ## order is hard-coded to match the order of the 'species_strings' variable for easier indexing
colnames(ortholog_mat) <- names(species_strings)

txi_ortho <- list(abundance = raw_txi[[1]]$abundance[ortholog_mat[,1],], 
                  counts    = raw_txi[[1]]$counts[ortholog_mat[,1],],
                  length    = raw_txi[[1]]$length[ortholog_mat[,1],],
                  countsFromAbundance = raw_txi[[1]]$countsFromAbundance)
for(i in 2:length(species_strings)) {
  txi_ortho$abundance <- cbind(txi_ortho$abundance, raw_txi[[i]]$abundance[ortholog_mat[,i],])
  txi_ortho$counts <- cbind(txi_ortho$counts, raw_txi[[i]]$counts[ortholog_mat[,i],])
  txi_ortho$length <- cbind(txi_ortho$length, raw_txi[[i]]$length[ortholog_mat[,i],])
} 
rownames(txi_ortho$abundance) <- rownames(txi_ortho$counts) <- rownames(txi_ortho$length) <- ortholog_mat[match(rownames(txi_ortho$length), ortholog_mat[,1]), 'Homo_sapiens']

## import immune annotations

## all ensembl identifiers corresponding to immune annotations
modules <- read.table(file.path(raw_data_prefix, 'immune_modules.txt'), sep='\t', quote='', header=T)
ensembl2ext_full$module <- modules$IIG_class[match(ensembl2ext_full$external_gene_name, modules$HGNC_symbol)]
ensembl2ext_full$module[is.na(ensembl2ext_full$module)] <- 'no_module_annotation'

immuneGenes <- unique(ensembl2ext_full$ensembl_gene_id[ensembl2ext_full$module != 'no_module_annotation'])

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
dat_ortho <- dat_ortho[dat_ortho$individual %in% sample_data_filt$Animal.ID,]
dat_ortho$body_mass_log_sp_std <- sapply(dat_ortho$individual, function(z) sample_data_filt$body_mass_log_sp_std[sample_data_filt$Animal.ID==z])
dat_ortho$body_mass_log_diff_std <- sapply(dat_ortho$individual, function(z) sample_data_filt$body_mass_log_diff_std[sample_data_filt$Animal.ID==z])

## define linear combinations of effects of interest
lc1 <- INLA::inla.make.lincomb(body_mass_log_sp_std = 1, 'body_mass_log_sp_std:treatmentLPS' = 1)
names(lc1) = "induced_allometry"
lc2 <- INLA::inla.make.lincomb(body_mass_log_diff_std = 1, 'treatmentLPS:body_mass_log_diff_std' = 1)
names(lc2) = "induced_ontology"

cat('Modeling immune orthologs\n')
ortho_immune_res <- phyr::pglmm(count ~ offset(norm_factor) + body_mass_log_sp_std * treatment + body_mass_log_diff_std + body_mass_log_diff_std:treatment + (1 | individual) + (1 | species__) + (1 | species@treatment) + (1 | species__@treatment), family = "poisson", cov_ranef = list(species = species_tree_norm), data=dat_ortho, bayes=TRUE, bayes_options=list(lincomb=c(lc1,lc2)))
ortho_immune_restab <- ortho_immune_res$inla.model$summary.lincomb.derived
##

## pool ALL immune-annotated genes
txi_all_immune_pool <- txi_ortho_immune_pool
txi_all_immune_pool$abundance['immune',] <- unlist(lapply(raw_txi, \(x) colSums(x$abundance[rownames(x$abundance) %in% immuneGenes, , drop=FALSE])))
txi_all_immune_pool$counts['immune',] <- unlist(lapply(raw_txi, \(x) colSums(x$counts[rownames(x$counts) %in% immuneGenes, , drop=FALSE])))
txi_all_immune_pool$length['immune',] <- unlist(lapply(raw_txi, \(x) colSums(x$length[rownames(x$length) %in% immuneGenes, , drop=FALSE] * x$abundance[rownames(x$abundance) %in% immuneGenes, , drop=FALSE]))) / txi_all_immune_pool$abundance['immune',]  ## weighted arithmetic mean length

## get normalization factors
des_all_pool <- DESeq2:::DESeqDataSetFromTximport(txi_all_immune_pool, data.frame(Int=rep(1,ncol(txi_all_immune_pool$counts))), ~1)
des_all_pool <- DESeq2::estimateSizeFactors(des_all_pool)

## model total differential expression of immune genes
dat_all <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des_all_pool)), '_'), \(x) x[[3]]), 
                      species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des_all_pool)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                      treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des_all_pool)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                      norm_factor = log(DESeq2::normalizationFactors(des_all_pool)['immune',]),
                      count       = DESeq2::counts(des_all_pool)['immune',])
dat_all <- dat_all[dat_all$individual %in% sample_data_filt$Animal.ID,]
dat_all$body_mass_log_sp_std <- sapply(dat_all$individual, function(z) sample_data_filt$body_mass_log_sp_std[sample_data_filt$Animal.ID==z])
dat_all$body_mass_log_diff_std <- sapply(dat_all$individual, function(z) sample_data_filt$body_mass_log_diff_std[sample_data_filt$Animal.ID==z])

cat('Modeling all immune genes\n')
all_immune_res <- phyr::pglmm(count ~ offset(norm_factor) + body_mass_log_sp_std * treatment + body_mass_log_diff_std + body_mass_log_diff_std:treatment + (1 | individual) + (1 | species__) + (1 | species@treatment) + (1 | species__@treatment), family = "poisson", cov_ranef = list(species = species_tree_norm), data=dat_all, bayes=TRUE, bayes_options=list(lincomb=c(lc1,lc2)))
all_immune_restab <- all_immune_res$inla.model$summary.lincomb.derived
##


## fit per-gene models

## get per-observation normalization factors
des <- DESeq2:::DESeqDataSetFromTximport(txi_ortho, data.frame(Int=rep(1,ncol(txi_ortho$counts))), ~1)
des <- DESeq2::estimateSizeFactors(des)

cat('Fitting models\n')
fits <- parallel::mclapply(rownames(des), function(x) {
  cat(round(which(x==rownames(des)) / nrow(des),2) * 100, '%: ', x,'\n')
  dat <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des)), '_'), \(x) x[[3]]), 
                    species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                    treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                    norm_factor = log(DESeq2::normalizationFactors(des)[x,]),
                    count       = DESeq2::counts(des)[x,])
  dat <- dat[dat$individual %in% sample_data_filt$Animal.ID,]
  dat$body_mass_log_sp_std <- sapply(dat$individual, function(z) sample_data_filt$body_mass_log_sp_std[sample_data_filt$Animal.ID==z])
  dat$body_mass_log_diff_std <- sapply(dat$individual, function(z) sample_data_filt$body_mass_log_diff_std[sample_data_filt$Animal.ID==z])
  
  if(sd(DESeq2::counts(des)[x,])>0) {

    fit <- tryCatch(phyr::pglmm(count ~ offset(norm_factor) + body_mass_log_sp_std * treatment + body_mass_log_diff_std + body_mass_log_diff_std:treatment + (1 | individual) + (1 | species__) + (1 | species@treatment) + (1 | species__@treatment), family = "poisson", cov_ranef = list(species = species_tree_norm), data=dat, bayes=TRUE, bayes_options=list(lincomb=c(lc1,lc2))),
                    error = function(e) NA)
    if(!all(is.na(fit))) {
      coefs <- c(phyr::fixef(fit)$Value, fit$inla.model$summary.lincomb.derived[1:2,'mean'])
      coef_cis <- rbind(fit$B.ci, LPS_beta=fit$inla.model$summary.lincomb.derived[1:2,c('0.025quant','0.975quant')])
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

