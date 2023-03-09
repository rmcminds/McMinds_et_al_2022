
output_prefix <- path.expand('outputs/primates/')

# genes identified via lit review with specific hypotheses
focal_gene_names <- c("AKT1", "CASP4", "CASP8", "CCL3", "CCL4", "CD14", "CD69", "CD80", "CD86", "EIF2AK2", "PTK2", "FKBP5", "Gadd45B", "HSPA1A", "IFNG", "IL-10", "Il-12A", "Il-12B", "IL-17A", "IL-17B", "IL-17C", "IL-17F", "IL-1b", "IL-1R1", "IL-1R2", "IL-6", "IL-8", "IL1a", "IRAK4", "IRF3", "LY96", "MAP2K1", "MAP4K5", "MAPK8", "PARP14", "PARP9", "PBEF", "PEBP4", "PTGS2", "SOD1", "SOD2", "SORT1", "TLR2", "TLR4", "TNFalpha", "TRAF1")

load(file.path(output_prefix,'03a_orthologs_phyr_results_noINLA.RData'))

curgene <- unique(ensembl2ext$ensembl_gene_id[ensembl2ext$external_gene_name == 'TLR4' & grepl('ENSG[[:digit:]]', ensembl2ext$ensembl_gene_id)])

curgene <- unique(ensembl2ext$ensembl_gene_id[ensembl2ext$external_gene_name == 'IL1B' & grepl('ENSG[[:digit:]]', ensembl2ext$ensembl_gene_id)])

plot(fits[[curgene]]$dat$body_mass_log_sp_std, log(fits[[curgene]]$dat$count) - fits[[curgene]]$dat$norm_factor, col=fits[[curgene]]$dat$treatment, pch=(1:nlevels(fits[[curgene]]$dat$species))[fits[[curgene]]$dat$species])
plot(fits[[curgene]]$dat$body_mass_log_diff_std, log(fits[[curgene]]$dat$count) - fits[[curgene]]$dat$norm_factor, col=fits[[curgene]]$dat$treatment)
plot(body_mass_log_diff_sd * fits[[curgene]]$dat$body_mass_log_diff_std + body_mass_log_sd * fits[[curgene]]$dat$body_mass_log_sp_std, log(fits[[curgene]]$dat$count) - fits[[curgene]]$dat$norm_factor, col=fits[[curgene]]$dat$treatment)

hi <- t(sapply(unique(fits[[curgene]]$dat$individual), \(x) {
  datnow <- fits[[curgene]]$dat
  c(datnow[datnow$individual==x & datnow$treatment=='LPS','body_mass_log_sp_std'], 
    (log(datnow[datnow$individual==x & datnow$treatment=='LPS','count']) - datnow[datnow$individual==x & datnow$treatment=='LPS','norm_factor']) - 
      (log(datnow[datnow$individual==x & datnow$treatment=='Null','count']) - datnow[datnow$individual==x & datnow$treatment=='Null','norm_factor']))
}))

plot(hi, pch=(1:nlevels(as.factor(sample_data_filt$Family)))[as.factor(sample_data_filt$Family)[match(rownames(hi), sample_data_filt$Animal.ID)]])

plot(body_mass_log_sd * dat_ortho$body_mass_log_sp_std, log(dat_ortho$count) - dat_ortho$norm_factor, col=dat_ortho$treatment, main='Total immune gene expression', xlab='log body size', ylab='log normalized counts')
legend(x='topleft',legend=c('Null','LPS'),col=c('black','red'),lty=1)

plot(body_mass_log_sd * dat_ortho$body_mass_log_sp_std, body_mass_log_diff_sd * dat_ortho$body_mass_log_diff_std)

x1 <- (body_mass_log_diff_sd * fits[[curgene]]$dat$body_mass_log_diff_std)
x2 <- x1 * (fits[[curgene]]$dat$treatment == 'LPS')
y <- (log(fits[[curgene]]$dat$count) - fits[[curgene]]$dat$norm_factor)
curgenein <- lm(y~x1+x2)$coefficients[[1]]
curgeneco1 <- lm(y~x1+x2)$coefficients[[2]]
curgeneco2 <- lm(y~x1+x2)$coefficients[[3]]
adjustedcurgene <- y - (curgeneco1 * x1) - (curgeneco2 * x2)
plot(body_mass_log_sd * fits[[curgene]]$dat$body_mass_log_sp_std, adjustedcurgene, col=fits[[curgene]]$dat$treatment)
##

# filter out genes that didn't converge
fitsFilt <- fits[!is.na(sapply(fits,function(x) x$fit))]

natural_to_common <- log(exp(1),10)
mean_rescale <- natural_to_common * body_mass_log_mean
sd_rescale <- natural_to_common * body_mass_log_sd
rescaled_coefficients <- t(sapply(fitsFilt, function(fit) {
  ## convert both from natural to common log and from standardized log body size to actual log body size
  coefs_old <- fit$coefs * natural_to_common
  cis_old <- fit$coef_cis * natural_to_common
  int_diff <- coefs_old[[2]] * (mean_rescale / sd_rescale)
  lps_diff <- coefs_old[[4]] * (mean_rescale / sd_rescale)

  coefs_new <- c(constitutive_intercept.plotting          = coefs_old[[1]] - int_diff,
                 induced_intercept.plotting               = coefs_old[[3]] - lps_diff,
                 constitutive_mean.estimate_norm          = coefs_old[[1]],
                 constitutive_mean.l95_norm               = cis_old[1,1],
                 constitutive_mean.u95_norm               = cis_old[1,2],
                 constitutive_allometry.estimate          = coefs_old[[2]] / sd_rescale,
                 constitutive_allometry.l95               = cis_old[2,1] / sd_rescale,
                 constitutive_allometry.u95               = cis_old[2,2] / sd_rescale,
                 constitutive_allometry.estimate_norm     = coefs_old[[2]],
                 constitutive_allometry.l95_norm          = cis_old[2,1],
                 constitutive_allometry.u95_norm          = cis_old[2,2],
                 conserved_induction.estimate_norm        = coefs_old[[3]],
                 conserved_induction.l95_norm             = cis_old[3,1],
                 conserved_induction.u95_norm             = cis_old[3,2],
                 allometric_induction.estimate            = coefs_old[[4]] / sd_rescale,
                 allometric_induction.l95                 = cis_old[4,1] / sd_rescale,
                 allometric_induction.u95                 = cis_old[4,2] / sd_rescale,
                 allometric_induction.estimate_norm       = coefs_old[[4]],
                 allometric_induction.l95_norm            = cis_old[4,1],
                 allometric_induction.u95_norm            = cis_old[4,2],
                 induced_slope.estimate_lincomb           = coefs_old[[5]] / sd_rescale,
                 induced_slope.l95_lincomb                = cis_old[5,1] / sd_rescale,
                 induced_slope.u95_lincomb                = cis_old[5,2] / sd_rescale)

  return(coefs_new)

}))

modules <- read.table('data/primate_allometry/hawash/immune_modules.txt', sep='\t', quote='', header=T)

gene_names2module <- sapply(gene_names,
                            function(x) sapply(strsplit(x,'/'),
                                               function(y) modules$IIG_class[modules$HGNC_symbol %in% y]))
gene_names2module[sapply(gene_names2module,length) == 0] <- 'none'
gene_names2module <- simplify2array(gene_names2module)

gnm <- gene_names[gene_names2module != 'none']
missing <- modules[!modules$HGNC_symbol %in% unlist(strsplit(gnm, '/')),]

genes_modules <- data.frame(module = gene_names2module[names(gene_names) %in% names(fitsFilt)],
                            gene   = gene_names[names(gene_names) %in% names(fitsFilt)],
                            up     = as.numeric(sapply(fitsFilt, function(x) x$body_mass_LPS_significant & x$coefs[[4]] > 0)),
                            down   = as.numeric(sapply(fitsFilt, function(x) x$body_mass_LPS_significant & x$coefs[[4]] < 0)))
genes_modules <- genes_modules[rownames(genes_modules) %in% names(gene_names)[gene_names2module != 'none'],]

save.image(file.path(output_prefix,'03b_orthologs_results.RData'))

subtree <- ape::drop.tip(species_tree, species_tree$tip.label[!species_tree$tip.label %in% names(species_strings)])
subdat <- sample_data[sample_data$Animal.ID %in% sub('_.*','',unlist(sapply(filtered,colnames))),]
subdat$Body.Mass..g. <- as.numeric(subdat$Body.Mass..g.)
subdat$genus_species <- factor(subdat$genus_species, levels=subtree$tip.label)

## plot species phylogeny with individual body sizes
pdf(file.path(output_prefix,paste0('03b_orthologs_phylomass.pdf')), width=3,height=3)
treedat <- data.frame(id=subdat$genus_species, size=log(subdat$Body.Mass..g.,10))
plo <- ggtree::ggtree(subtree) + ggtree::geom_tiplab()
plo2 <- ggtree::facet_plot(plo, panel='dot', data=treedat, geom=ggtree::geom_point, ggtree::aes(x=size))
plo2 + ggtree::theme_tree2()
dev.off()

## identify in my orthologs matches to the focal genes
focal_genes <- sapply(focal_gene_names, function(x) {
  searchterm <- paste0('(^|/)',x,'($|/)|(^|/)',gsub('-','',x),'($|/)')
  y <- names(gene_names)[grep(searchterm,gene_names,ignore.case=TRUE)]
  if(length(y) == 0) y <- NA
  return(y)
})

focal_genes_missing <- names(focal_genes[is.na(focal_genes)])
focal_genes <- focal_genes[!is.na(focal_genes)]
focal_genes_nonconvergence <- focal_genes[focal_genes %in% names(fits[is.na(sapply(fits,function(x) x$fit))])]
focal_genes <- focal_genes[!focal_genes %in% focal_genes_nonconvergence]

focal_summaries <- as.data.frame(rescaled_coefficients[rownames(rescaled_coefficients) %in% focal_genes,])
focal_summaries$module <- gene_names2module[match(rownames(focal_summaries),names(gene_names2module))]
rownames(focal_summaries) <- names(focal_genes)[match(rownames(focal_summaries), focal_genes)]
write.table(focal_summaries, file=file.path(output_prefix,'03b_orthologs_focal_summaries.txt'), sep='\t', quote=FALSE)

## plot semi-raw relative abundance for the focal genes
rel2abs <- 0.25
gene <- 1
for(x in 1:ceiling(length(focal_genes)/7)) {
  pdf(file.path(output_prefix,paste0('03b_orthologs_raw_focal_genes_',x,'.pdf')), width=16, height=8)
  par(mfrow=c(2,4))
  for(orth in focal_genes[gene:min((gene+6),length(focal_genes))]) {

    newdat <- fits[[orth]]$data
    newdat$logrel <- log(newdat$count+0.5) - newdat$size_factor
    newdat$logrel10 <- newdat$logrel * natural_to_common
    newdat$body_mass_log <- newdat$body_mass_log_std * body_mass_log_sd + body_mass_log_mean
    newdat$body_mass_log10 <- newdat$body_mass_log * natural_to_common

    coefs <- rescaled_coefficients[orth,]

    iso_intercept <- -mean_rescale * rel2abs + coefs[['constitutive_intercept.plotting']] + mean_rescale * coefs[['constitutive_allometry.estimate']] + 0.5 * (coefs[['induced_intercept.plotting']] + mean_rescale * coefs[['allometric_induction.estimate']])

    plot(newdat[,c('body_mass_log10','logrel10')], col=c('blue','red')[newdat$treatment], main=gene_names[[orth]], xlab='Log10 g body mass', ylab='Log10 relative abundance', pch=(1:nlevels(newdat$species))[newdat$species], col.main = c('black','red')[fits[[orth]]$body_mass_LPS_significant + 1])
    abline(coef=c(iso_intercept,rel2abs), lty=2)
    abline(coef=coefs[c('constitutive_intercept.plotting','constitutive_allometry.estimate')], col='blue')
    abline(coef=coefs[c('constitutive_intercept.plotting','constitutive_allometry.estimate')] + coefs[c('induced_intercept.plotting','allometric_induction.estimate')], col='red')

  }
  plot.new()
  legend(x='bottomleft',legend=c('Null','LPS'),col=c('blue','red'),lty=1)
  legend(x='topleft',legend=levels(newdat$species), pch=(1:nlevels(newdat$species)), lty = rep(NULL,nlevels(newdat$species)))

  dev.off()

  gene <- gene + 7
}

## re-plot scatter plots for more-focused set of genes TLR4, IL1B, and CD69 for figure 1
pdf(file.path(output_prefix,paste0('Figure_1_raw.pdf')), width=4, height=16)
par(mfrow=c(4,1))
for(gene in c('TLR4', 'IL1B', 'CD69')) {

  orth <- names(gene_names)[grep(gene,gene_names)]
  newdat <- fits[[orth]]$data
  newdat$logrel <- log(newdat$count+0.5) - newdat$size_factor
  newdat$logrel10 <- newdat$logrel * natural_to_common
  newdat$body_mass_log <- newdat$body_mass_log_std * body_mass_log_sd + body_mass_log_mean
  newdat$body_mass_log10 <- newdat$body_mass_log * natural_to_common

  coefs <- rescaled_coefficients[orth,]

  iso_intercept <- -mean_rescale * rel2abs + coefs[['constitutive_intercept.plotting']] + mean_rescale * coefs[['constitutive_allometry.estimate']] + 0.5 * (coefs[['induced_intercept.plotting']] + mean_rescale * coefs[['allometric_induction.estimate']])

  plot(newdat[,c('body_mass_log10','logrel10')], col=c('blue','red')[newdat$treatment], main=gene_names[[orth]], xlab='Log10 g body mass', ylab='Log10 relative abundance', pch=(1:nlevels(newdat$species))[newdat$species], col.main = c('black','red')[fits[[orth]]$body_mass_LPS_significant + 1])
  abline(coef=c(iso_intercept,rel2abs), lty=2)
  abline(coef=coefs[c('constitutive_intercept.plotting','constitutive_allometry.estimate')], col='blue')
  abline(coef=coefs[c('constitutive_intercept.plotting','constitutive_allometry.estimate')] + coefs[c('induced_intercept.plotting','allometric_induction.estimate')], col='red')

}
plot.new()
legend(x='bottomleft',legend=c('Null','LPS'),col=c('blue','red'),lty=1)
legend(x='topleft',legend=levels(newdat$species), pch=(1:nlevels(newdat$species)), lty = rep(NULL,nlevels(newdat$species)))

dev.off()


## wrangle coefficients for plotting and meta-stats
newcoefs <- reshape2::melt(rescaled_coefficients, id.vars=row.names)
colnames(newcoefs) <- c('gene','variable','value')
newcoefs$quantity <- sapply(as.character(newcoefs$variable),function(x) strsplit(x,'\\.')[[1]][[2]])
newcoefs$variable <- sapply(as.character(newcoefs$variable),function(x) strsplit(x,'\\.')[[1]][[1]])

strong_hyper_const <- newcoefs[newcoefs$quantity == 'l95' & newcoefs$value > 0.25 & newcoefs$variable == 'constitutive_allometry',]
strong_hyper_ind <- newcoefs[newcoefs$quantity == 'l95' & newcoefs$value > 0.25 & newcoefs$variable == 'allometric_induction',]

newcoefs <- merge(newcoefs,genes_modules, by.x='gene', by.y=0, all=TRUE)
newcoefs$module[is.na(newcoefs$module)] <- 'non-immune'
module_means <- sapply(levels(newcoefs$variable), function(effname) sapply(unique(newcoefs$module), function(x) {
  mean(newcoefs[newcoefs$quantity=='estimate_norm' & newcoefs$module==x & newcoefs$variable==effname,'value'],na.rm=TRUE)
}))
newcoefs$modules_merged <- newcoefs$module
newcoefs$modules_merged[!newcoefs$modules_merged %in% c('non-immune','sensor','effector')] <- 'bureaucracy'
modules_merged_means <- sapply(levels(newcoefs$variable), function(effname) sapply(unique(newcoefs$modules_merged), function(x) {
  mean(newcoefs[newcoefs$quantity=='estimate_norm' & newcoefs$modules_merged==x & newcoefs$variable==effname,'value'],na.rm=TRUE)
}))
newcoefs$modules_merged <- as.factor(newcoefs$modules_merged)
species_counts <- sapply(filtered,function(x) rowSums(x))
rownames(species_counts) <- rownames(allcounts)
gre1 <- apply(species_counts>1,1,function(x)sum(x)>2) ### is this objective??
subs <- names(gre1)[gre1]
newcoefs <- newcoefs[newcoefs$variable != 'constitutive_mean' & newcoefs$gene %in% subs,]

newcoefs$module <- relevel(factor(newcoefs$module), 'non-immune')
newcoefs$modules_merged <- factor(newcoefs$modules_merged, levels=c('non-immune','sensor', 'bureaucracy', 'effector'))
newcoefs$variable <- factor(newcoefs$variable, levels=c("constitutive_intercept", "conserved_induction","constitutive_allometry","induced_intercept","allometric_induction","induced_slope"))
###

## do genes in different modules have consistent responses to LPS or body size
sink(file.path(output_prefix,'03b_orthologs_module_lm_alleffects_allmodules.txt'))
print(summary(lm(value~0+modules_merged:variable,data=droplevels(newcoefs[newcoefs$quantity=='estimate_norm',]))))
sink()
pdf(file.path(output_prefix,'03b_orthologs_module_beanplots_alleffects_allmodules.pdf'))
beanplot::beanplot(value~modules_merged:variable, data=droplevels(newcoefs[newcoefs$quantity=='estimate_norm',]), las=2, col=rep(list("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A"),3), main='Per-gene linear model coefficients summarized per module', ylim=c(-1,1))
abline(h=0)
dev.off()

## do genes in different modules have consistent responses to LPS or body size (include the linear combination quantity)
sink(file.path(output_prefix,'03b_orthologs_module_lm_alleffects_allmodules_wlincomb.txt'))
print(summary(lm(value~0+modules_merged:variable,data=droplevels(newcoefs[newcoefs$quantity %in% c('estimate_norm','estimate_lincomb'),]))))
sink()
pdf(file.path(output_prefix,'03b_orthologs_module_beanplots_alleffects_allmodules_wlincomb.pdf'))
beanplot::beanplot(value~modules_merged:variable, data=droplevels(newcoefs[newcoefs$quantity %in% c('estimate_norm','estimate_lincomb'),]), las=2, col=rep(list("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A"),3), main='Per-gene linear model coefficients summarized per module', what=c(1,1,1,0), ylim=c(-1,1))
abline(h=0)
dev.off()

## are allometries significant when converted to predicted absolute concentrations
absolute_coefs <- droplevels(newcoefs[newcoefs$quantity %in% c('estimate','estimate_lincomb') & newcoefs$variable != 'allometric_induction',])
absolute_coefs$value <- absolute_coefs$value - 0.25
sink(file.path(output_prefix,'03b_orthologs_module_lm_allmodules_absolute_estimates.txt'))
print(summary(lm(value~0+modules_merged:variable,data=absolute_coefs)))
sink()


## do immune-annotated genes in different modules have consistent responses to LPS or body size
tempdat <- droplevels(newcoefs[newcoefs$module != 'non-immune' & newcoefs$quantity=='estimate_norm',])
pdf(file.path(output_prefix,'03b_orthologs_module_beanplots_alleffects_onlyImmune.pdf'))
beanplot::beanplot(value~modules_merged:variable,data=tempdat, las=2, col=rep(list("#E41A1C", "#377EB8", "#4DAF4A"),3), main='All immune genes')
abline(h=0)
legend("bottomright", bty="n", levels(tempdat$modules_merged), fill = c("#E41A1C", "#377EB8", "#4DAF4A"))
dev.off()

## do genes that respond positively or negatively to LPS have consistent responses to body size or the interaction
upinlps <- newcoefs[newcoefs$variable == 'conserved_induction' & newcoefs$quantity == 'estimate_norm' & newcoefs$value > 0,'gene']
downinlps <- newcoefs[newcoefs$variable == 'conserved_induction' & newcoefs$quantity == 'estimate_norm' & newcoefs$value < 0,'gene']
newcoefs$LPSdirection <- sapply(newcoefs$gene, function(x) if(x %in% upinlps) {return('up')} else if(x %in% downinlps) {return('down')} else {return('nonsig')})

sink(file.path(output_prefix,'03b_orthologs_module_lm_alleffects_allmodules_per_LPS_response.txt'))
print(summary(lm(value~0+modules_merged:LPSdirection:variable,data=droplevels(newcoefs[newcoefs$variable != 'conserved_induction' & newcoefs$quantity=='estimate_norm' & newcoefs$LPSdirection!='nonsig',]))))
sink()

tempdat <- droplevels(newcoefs[newcoefs$gene %in% c(upinlps,downinlps) & newcoefs$variable != 'conserved_induction' & newcoefs$quantity=='estimate_norm',])
pdf(file.path(output_prefix,'03b_orthologs_module_beanplots_LPS_response_allometricEffects_allmodules.pdf'))
beanplot::beanplot(value~modules_merged:LPSdirection:variable,data=tempdat, las=2, col=rep(list("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A"),2), main='Sig-LPS gene allometry')
abline(h=0)
legend("bottomright", bty="n", levels(tempdat$modules_merged), fill = c("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A"))
dev.off()

## make plots that can be combined to show main LPS effects and the LPS-specific allometries (figure 3)
pdf(file.path(output_prefix,'03b_orthologs_module_beanplots_LPS_response_allometricEffects_allmodules_nonsymm.pdf'))
beanplot::beanplot(value~LPSdirection:modules_merged:variable,data=tempdat, las=2, col=rep(list("#C28DC8", "#733A78", "#F39192", "#B31417", "#89B7DC", "#2A618D", "#BFE3BF", "#3A8939"),2), main='Gene allometry by mean LPS response', side='both', what=c(1,1,1,0), ylim=c(-1,1))
abline(h=0)
legend("bottomright", bty="n", levels(tempdat$modules_merged), fill = c("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A"))
dev.off()

tempdat <- droplevels(newcoefs[newcoefs$variable != 'allometric_induction' & newcoefs$quantity=='estimate_norm',])
pdf(file.path(output_prefix,'03b_orthologs_module_beanplots_alleffects_allmodules_tocombine_w_LPS_nonsymm.pdf'))
beanplot::beanplot(value~modules_merged:variable, data=tempdat, las=2, col=list("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A"), main='alleffects', what=c(1,1,1,0), ylim=c(-1,1))
abline(h=0)
dev.off()
##

##do same as above but include linear combination quantity

sink(file.path(output_prefix,'03b_orthologs_module_lm_alleffects_allmodules_per_LPS_response_wlincomb.txt'))
print(summary(lm(value~0+modules_merged:LPSdirection:variable,data=droplevels(newcoefs[newcoefs$variable != 'conserved_induction' & newcoefs$quantity %in% c('estimate_norm','estimate_lincomb') & newcoefs$LPSdirection!='nonsig',]))))
sink()

tempdat <- droplevels(newcoefs[newcoefs$gene %in% c(upinlps,downinlps) & newcoefs$variable != 'conserved_induction' & newcoefs$quantity %in% c('estimate_norm','estimate_lincomb'),])
pdf(file.path(output_prefix,'03b_orthologs_module_beanplots_LPS_response_allometricEffects_allmodules_wlincomb.pdf'))
beanplot::beanplot(value~modules_merged:LPSdirection:variable,data=tempdat, las=2, col=rep(list("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A"),2), main='Sig-LPS gene allometry')
abline(h=0)
legend("bottomright", bty="n", levels(tempdat$modules_merged), fill = c("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A"))
dev.off()

## make plots that can be combined to show main LPS effects and the LPS-specific allometries (figure 3)
pdf(file.path(output_prefix,'03b_orthologs_module_beanplots_LPS_response_allometricEffects_allmodules_nonsymm_wlincomb.pdf'))
beanplot::beanplot(value~LPSdirection:modules_merged:variable,data=tempdat, las=2, col=rep(list("#C28DC8", "#733A78", "#F39192", "#B31417", "#89B7DC", "#2A618D", "#BFE3BF", "#3A8939"),2), main='Gene allometry by mean LPS response', side='both', what=c(1,1,1,0), ylim=c(-1,1))
abline(h=0)
legend("bottomright", bty="n", levels(tempdat$modules_merged), fill = c("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A"))
dev.off()

tempdat <- droplevels(newcoefs[newcoefs$variable != 'allometric_induction' & newcoefs$quantity %in% c('estimate_norm','estimate_lincomb'),])
pdf(file.path(output_prefix,'03b_orthologs_module_beanplots_alleffects_allmodules_tocombine_w_LPS_nonsymm_wlincomb.pdf'))
beanplot::beanplot(value~modules_merged:variable, data=tempdat, las=2, col=list("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A"), main='alleffects', what=c(1,1,1,0), ylim=c(-1,1))
abline(h=0)
dev.off()
##


## do immune-annotated that respond positively or negatively to LPS have consistent responses to body size or the interaction
tempdat <- droplevels(newcoefs[newcoefs$gene %in% c(upinlps,downinlps) & newcoefs$variable != 'conserved_induction' & newcoefs$module != 'non-immune' & newcoefs$quantity=='estimate_norm',])
pdf(file.path(output_prefix,'03b_orthologs_module_beanplots_LPS_response_allometricEffects_onlyImmune.pdf'))
beanplot::beanplot(value~modules_merged:LPSdirection:variable, data=tempdat, las=2, col=rep(list("#E41A1C", "#377EB8", "#4DAF4A"),2), main='Sig-LPS gene allometry')
abline(h=0)
legend("bottomright", bty="n", levels(tempdat$modules_merged), fill = c("#E41A1C", "#377EB8", "#4DAF4A"))
dev.off()

## do immune-annotated genes that respond positively to LPS have consistent responses to body size or the interaction
tempdat <- droplevels(newcoefs[newcoefs$module != 'non-immune' & newcoefs$quantity=='estimate_norm' & newcoefs_onlyimmune$gene %in% upinlps & newcoefs_onlyimmune$variable != 'conserved_induction',])
pdf(file.path(output_prefix,'03b_orthologs_module_beanplots_onlyLPSgenes_allometricEffects.pdf'))
beanplot::beanplot(value~modules_merged:variable,data=tempdat, las=2, col=rep(list("#E41A1C", "#377EB8", "#4DAF4A"),2), main='Positive-LPS gene allometry')
abline(h=0)
legend("bottomright", bty="n", levels(tempdat$modules_merged), fill = c("#E41A1C", "#377EB8", "#4DAF4A"))
dev.off()

## do immune-annotated genes that respond positively or negatively to LPS have consistent responses to body size or the interaction
tempdat <- droplevels(newcoefs[newcoefs$module != 'non-immune' & newcoefs$quantity=='estimate_norm' & newcoefs$variable != 'conserved_induction',])
pdf(file.path(output_prefix,'03b_orthologs_module_beanplots_onlyLPSgenesPOSNEG_allometricEffects.pdf'))
beanplot::beanplot(value~modules_merged:LPSdirection:variable,data=tempdat, las=2, col=rep(list("#E41A1C", "#377EB8", "#4DAF4A"),2), main='Sig-LPS gene allometry')
abline(h=0)
legend("bottomright", bty="n", levels(tempdat$modules_merged), fill = c("#E41A1C", "#377EB8", "#4DAF4A"))
dev.off()

## do immune-annotated genes that respond /significantly/ positively or negatively to LPS have consistent responses to body size or the interaction
tempdat <- droplevels(newcoefs[newcoefs$module != 'non-immune' & newcoefs$quantity=='estimate_norm' & newcoefs$variable != 'conserved_induction' & newcoefs$LPSdirection != 'nonsig',])
pdf(file.path(output_prefix,'03b_orthologs_module_beanplots_onlyLPSgenesPOSNEG_nononsig_allometricEffects.pdf'))
beanplot::beanplot(value~modules_merged:LPSdirection:variable,data=tempdat, las=2, col=rep(list("#E41A1C", "#377EB8", "#4DAF4A"),2), main='Sig-LPS gene allometry')
abline(h=0)
legend("bottomright", bty="n", levels(tempdat$modules_merged), fill = c("#E41A1C", "#377EB8", "#4DAF4A"))
dev.off()

newcoefs_merged_sig <- list()
for(eff in levels(newcoefs$variable)) {
  genes2keep4 <- sapply(unique(newcoefs$gene), function(x) prod(sign(newcoefs[newcoefs$quantity %in% c('l95_norm','u95_norm') & newcoefs$gene == x & newcoefs$variable==eff,'value'])) > 0)
  newcoefs_merged_sig[[eff]] <- droplevels(newcoefs[newcoefs$gene %in% genes2keep4 & newcoefs$variable == eff,])
}


