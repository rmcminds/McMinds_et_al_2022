## incorporate 'download_GO' here directly
## modify using logic tested in 05b_run_batches_sparseOU_gp.r to ensure species means are distinct from individual sizes, and other logic therein

nthreads <- 6
raw_data_prefix <- path.expand('raw_data/20221215_primate_allometry/')
output_prefix <- path.expand('outputs/primates_20230224/')

species_strings <- c(Callithrix_jacchus='ENSCJAG', Homo_sapiens='ENSG', Microcebus_murinus='ENSMICG', Macaca_mulatta='ENSMMUG', Papio_anubis='ENSPANG', Pongo_abelii='ENSPPYG')

## download and normalize species chronogram
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
countfiles <- Sys.glob(file.path(output_prefix, '/01_generate_counts/*/t_data.ctab'))
samples <- basename(dirname(countfiles))
names(countfiles) <- samples
sample_data_filt <- sample_data[sample_data$Animal.ID %in% sub('_.*', '', samples),]

raw_txi <- lapply(names(species_strings), \(x) {
  cfiles <- countfiles[sub('_.*', '', samples) %in% sample_data_filt$Animal.ID[sample_data_filt$genus_species == x]]
  tmp <- read.table(cfiles[[1]], header = TRUE, sep = '\t')
  tmp$gene_id <- sub('.*:', '', tmp$gene_id)
  ensembl2ext <- unique(tmp[, c('gene_id', 'gene_name')])
  tx2gene <- tmp[, c('t_name', 'gene_id')]
  return(c(tximport::tximport(cfiles, type = 'stringtie', tx2gene = tx2gene, readLength = 100), ensembl2ext = ensembl2ext))
})
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
body_mass_log_center <- body_mass_opt$para[3,1]
##

## standardize masses for model input
body_mass_log_sd <- sd(sample_data_filt$body_mass_log_sp)
body_mass_log_diff_sd <- sd(sample_data_filt$body_mass_log_diff)

sample_data_filt$body_mass_log_sp_std <- (sample_data_filt$body_mass_log_sp - body_mass_log_center) / body_mass_log_sd
sample_data_filt$body_mass_log_diff_std <- sample_data_filt$body_mass_log_diff / body_mass_log_diff_sd
##

## pool all immune gene counts and analyze univariate binomial model (percent of expressed transcripts annotated with immune functions)
## what is the most interesting contrast here? immune annotated vs literally everything else including ribosomes and noncoding? filter out ribosomes? just look at protein coding? ribosomes and globin were targeted for removal during library prep so could have odd patterns and still dominate the contrast. probably should just look at protein coding genes, and then also filter out globin reads by finding anything in the same gene family as the harrington gene list
## feels weird to round decimals for poisson error, but since data were counts at one point, zeros are possible, so doesn't make sense to just log-transform and use gaussian error; and error probably still scales like poisson such that low counts are less meaningful (and the decimal values rounded off would contribute negligible information)
## how to pool normalization factors? arithmetic, geometric, harmonic mean, or something else? in terms of energy use, maybe shouldn't normalize by length at all...? three levels: total rna produced / total amino acids, total transcripts by number pooling gene duplications, and total transcripts per gene. each tells different
## refs (maybe not for pooled calc but just for later?):
## DESeq2:::DESeqDataSetFromTximport
## getMethod(DESeq2::estimateSizeFactors, 'DESeqDataSet')
## DESeq2:::estimateNormFactors
percent_immune_res <- phyr::pglmm(count ~ offset(size_factor) + body_mass_log_sp_std * treatment + body_mass_log_diff_std * treatment + (1 | individual) + (1 | species__) + (1 | species@treatment) + (1 | species__@treatment), family = "poisson", cov_ranef = list(species = species_tree_norm), data=dat, bayes=TRUE, bayes_options=list(lincomb=lc))
##

##
genetrees <- read.tree(file.path(output_prefix, '00_references/Compara.109.protein_default.newick'))
names(genetrees) <- paste('family', 1:length(genetrees), sep='_') ## arbitrary gene family names

orthologs <- lapply(genetrees, function(tree) {
  nodes <- unique(tree$edge[,1])
  alldescs <- lapply(phangorn::Descendants(tree, nodes), function(y) tree$tip.label[y])
  names(alldescs) <- as.character(nodes)
  is.ortho <- sapply(alldescs, function(y) {
    all(sapply(species_strings, function(z) any(grepl(z, y)))) & length(y) == length(species_strings)
  })
  return(alldescs[is.ortho])
})
names(orthologs) <- names(genetrees)
orthologs <- unlist(orthologs,recursive=F)

cat('Filtering\n')
filtered <- lapply(names(species_strings), function(x) raw[[x]][grep(species_strings[[x]],unlist(orthologs),value=T),] )
allcounts <- do.call(cbind,filtered)
rownames(allcounts) <- names(orthologs)

individuals <- unique(gsub('_.*','',colnames(allcounts)))
N_individuals <- length(individuals)


## fit per-gene models
cat('Fitting models\n')
fits <- parallel::mclapply(rownames(allcounts), function(x) {
  cat(round(which(x==rownames(allcounts)) / nrow(allcounts),2) * 100, '%: ', x,'\n')
  dat <- data.frame(individual=sub('_.*','',colnames(allcounts)), treatment=factor(sub('.*_', '', colnames(allcounts)), levels=c('Null','LPS')), count = unname(t(allcounts[x,])))
  
  ## change this to be more deseq-like and incorporate gene lengths!!
  dat$size_factor <- apply(allcounts,2,function(x) sum(log(x[x>0])) / length(x))
  ##
  
  dat$species <- as.factor(sapply(dat$individual, function(z) sample_data_filt$genus_species[sample_data_filt$Animal.ID==z]))
  dat$body_mass_log_sp_std <- sapply(dat$individual, function(z) sample_data_filt$body_mass_log_sp_std[sample_data_filt$Animal.ID==z])
  dat$body_mass_log_diff_std <- sapply(dat$individual, function(z) sample_data_filt$body_mass_log_diff_std[sample_data_filt$Animal.ID==z])

  if(sd(allcounts[x,])>0) {
    lc <- INLA::inla.make.lincomb(body_mass_log_sp_std = 1, 'body_mass_log_sp_std:treatmentLPS' = 1)
    ##add lcomb for diffs
    fit <- tryCatch(phyr::pglmm(count ~ offset(size_factor) + body_mass_log_sp_std * treatment + body_mass_log_diff_std * treatment + (1 | individual) + (1 | species__) + (1 | species@treatment) + (1 | species__@treatment), family = "poisson", cov_ranef = list(species = species_tree_norm), data=dat, bayes=TRUE, bayes_options=list(lincomb=lc)),
                    error = function(e) NA)
    if(!all(is.na(fit))) {
      coefs <- c(phyr::fixef(fit)$Value, fit$inla.model$summary.lincomb.derived[1,'mean'])
      coef_cis <- rbind(fit$B.ci, LPS_beta=as.numeric(fit$inla.model$summary.lincomb.derived[1,c('0.025quant','0.975quant')]))
      body_mass_constitutive_significant <- (prod(sign(fit$B.ci[2,])) == 1)
      conserved_LPS_significant <- (prod(sign(fit$B.ci[3,])) == 1)
      body_mass_LPS_significant <- (prod(sign(fit$B.ci[4,])) == 1)
      fit$inla.model <- NULL
    } else {
      coefs <- NA
      coef_cis <- NA
      body_mass_constitutive_significant <- NA
      conserved_LPS_significant <- NA
      body_mass_LPS_significant <- NA
    }
  } else {
    fit <- NA
    coefs <- NA
    coef_cis <- NA
    body_mass_constitutive_significant <- NA
    conserved_LPS_significant <- NA
    body_mass_LPS_significant <- NA
  }
  return(list(data=dat, fit=fit, coefs=coefs, coef_cis=coef_cis, body_mass_constitutive_significant=body_mass_constitutive_significant, conserved_LPS_significant=conserved_LPS_significant, body_mass_LPS_significant=body_mass_LPS_significant))
}, mc.cores = nthreads)
names(fits) <- names(orthologs)

sigs_body_mass_constitutive <- which(as.logical(sapply(fits, function(x) x$body_mass_constitutive_significant)))
sigs_conserved_LPS_significant <- which(as.logical(sapply(fits, function(x) x$conserved_LPS_significant)))
sigs_body_mass_LPS <- which(as.logical(sapply(fits, function(x) x$body_mass_LPS_significant)))

load('outputs/primates/ensembl2GO.Rdata')
gene2GO <- do.call(rbind, gene_data)
gene_names <- sapply(orthologs,function(x) {
  y <- unique(gene2GO[gene2GO[,1] %in% x,2])
  y <- y[y != '']
  y <- y[!is.na(y)]
  if(length(y) == 0) y <- 'Unknown gene'
  return(paste(y,collapse='/'))
})

orthoinv <- setNames(rep(names(orthologs), each=length(orthologs[[1]])), unlist(orthologs))

save.image(file.path(output_prefix,'03a_orthologs_phyr_results_noINLA.RData'))

