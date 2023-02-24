## incorporate 'filter_genetrees' logic here directly (for orthologs, might be redundant, but for summed annotations with paralogs, want to think)
## incorporate 'download_GO' here directly
## modify using logic tested in 05b_run_batches_sparseOU_gp.r to ensure species means are distinct from individual sizes, and other logic therein


load('outputs/primates/filter_genetrees.RData')

nthreads <- 6
output_prefix <- path.expand('outputs/primates/')

genetree_paths <- list.files(path='outputs/primates/filtered_ensembl_trees', pattern='*.phy$', full.names=T)
names(genetree_paths) <- sapply(genetree_paths, function(x) paste(strsplit(basename(x),'.', fixed=T)[[1]][[1]], sep='_') )

genetrees <- lapply(genetree_paths, function(x) ape::read.tree(x))

species_strings <- c(Callithrix_jacchus='ENSCJAG',Homo_sapiens='ENSG',Microcebus_murinus='ENSMICG',Macaca_mulatta='ENSMMUG',Papio_hamadryas='ENSPANG',Pongo_abelii='ENSPPYG')

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

species_tree <- ape::read.tree('data/primate_allometry/Phylo tree files/pruned_elevenspecies.phy')
species_tree$tip.label[species_tree$tip.label == 'Pongo_pygmaeus'] <- 'Pongo_abelii'

sample_data_filt <- sample_data[sample_data$Animal.ID %in% sub('_.*','',colnames(allcounts)),]
sample_data_filt$body_mass_log <- log(as.numeric(sample_data_filt$Body.Mass..g.))
body_mass_log_mean <- mean(sample_data_filt$body_mass_log)
body_mass_log_sd <- sd(sample_data_filt$body_mass_log)
sample_data_filt$body_mass_log_std <- (sample_data_filt$body_mass_log - body_mass_log_mean) / body_mass_log_sd

cat('Fitting models\n')
fits <- parallel::mclapply(rownames(allcounts), function(x) {
  cat(round(which(x==rownames(allcounts)) / nrow(allcounts),2) * 100, '%: ', x,'\n')
  dat <- data.frame(individual=sub('_.*','',colnames(allcounts)), treatment=factor(sub('.*_', '', colnames(allcounts)), levels=c('Null','LPS')), count = unname(t(allcounts[x,])))
  dat$size_factor <- apply(allcounts,2,function(x) sum(log(x[x>0])) / length(x))
  dat$species <- as.factor(sapply(dat$individual, function(z) sample_data_filt$genus_species[sample_data_filt$Animal.ID==z]))
  dat$body_mass_log_std <- sapply(dat$individual, function(z) sample_data_filt$body_mass_log_std[sample_data_filt$Animal.ID==z])

  if(sd(allcounts[x,])>0) {
    lc <- INLA::inla.make.lincomb(body_mass_log_std=1,'body_mass_log_std:treatmentLPS'=1)
    fit <- tryCatch(phyr::pglmm(count ~ offset(size_factor) + body_mass_log_std * treatment + (1 | individual) + (1 | species__) + (1 | species@treatment) + (1 | species__@treatment), family = "poisson", cov_ranef = list(species = species_tree), data=dat, bayes=TRUE, bayes_options=list(lincomb=lc)),
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

