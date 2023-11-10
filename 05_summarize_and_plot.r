
output_prefix <- path.expand('outputs/primates_20230314_mixed/')

load(file.path(output_prefix, '04_phyr_results.RData'))

output_prefix <- path.expand('outputs/primates_20230314_mixed/05_summaries_and_plots')

dir.create(output_prefix)

ln2log10 <- log(2) / log(2,base=10)

## define some plotting functions
# plot the log of the normalized counts as a function of species' mean body size
plot_species_mean_raw <- function(dat, restab, plot_legends = TRUE, plot_title = 'Total immune gene expression', font.main=2, plot_axes=TRUE) {
  
  lo <- levels(as.factor(dat$species))
  matchlo <- match(names(sort(body_mass_log_sp_means)), lo)
  
  if(plot_axes) {
    plot(exp(body_mass_log_center + body_mass_log_sd * dat$body_mass_log_sp_std) / 1000, 
         log(dat$count + 0.5) - dat$norm_factor, 
         col  = dat$treatment, 
         main = plot_title,
         xlab = 'Species\' mean body mass (kg)', 
         ylab = 'Natural log of normalized counts', 
         pch  = (1:9)[dat$species],
         log  = 'x',
         asp = 1, 
         font.main = font.main)
  } else {
    plot(exp(body_mass_log_center + body_mass_log_sd * dat$body_mass_log_sp_std) / 1000, 
         log(dat$count + 0.5) - dat$norm_factor, 
         col  = dat$treatment, 
         xaxt = 'n',
         yaxt = 'n',
         ann  = FALSE,
         pch  = (1:9)[dat$species],
         log  = 'x',
         cex = 0.5,
         asp = 1, 
         font.main = font.main)    
  }
  if(plot_legends) {
    legend(x         = 'topleft', 
           legend    = sub('_', ' ', lo[matchlo[1:5]]), 
           pch       = (1:9)[matchlo[1:5]],
           text.font = 3)
    legend(x         = 'bottomright',
           legend    = sub('_', ' ', lo[matchlo[6:9]]), 
           pch       = (1:9)[matchlo[6:9]],
           text.font = 3)
    legend(x      = 'bottomleft',
           title  = 'Slopes',
           legend = c(paste0('LPS(+): β₂+β₃ = ', format(round(restab['evo_allometry_LPS','mean'],2),nsmall=2),          ' (', format(round(restab['evo_allometry_LPS','0.025quant'],2),nsmall=2),          '–', format(round(restab['evo_allometry_LPS','0.975quant'],2),nsmall=2),          ')'),
                      paste0('Δ–LPS: β₃ = ',     format(round(restab['evo_allometry_LPS_response','mean'],2),nsmall=2), ' (', format(round(restab['evo_allometry_LPS_response','0.025quant'],2),nsmall=2), '–', format(round(restab['evo_allometry_LPS_response','0.975quant'],2),nsmall=2), ')'),
                      paste0('LPS(–): β₂ = ',    format(round(restab['evo_allometry_baseline','mean'],2),nsmall=2),     ' (', format(round(restab['evo_allometry_baseline','0.025quant'],2),nsmall=2),     '–', format(round(restab['evo_allometry_baseline','0.975quant'],2),nsmall=2),     ')')),
           col    = c('red','black','white'),
           lty    = 1)
  }
  ## add the fit line from the model for Null samples
  abline(a = restab['(Intercept)','mean'] - restab['body_mass_log_sp_std','mean'] * (body_mass_log_center - log(1000)) / body_mass_log_sd,
         b = ln2log10 * restab['evo_allometry_baseline','mean'],
         col = scales::alpha('black',0.8))
  ## and the corresponding credible intervals
  matlines(exp(body_mass_log_center + body_mass_log_sd * seq_bmsp) / 1000, restab[7:56, c('0.025quant','0.975quant')], lty=2, col=scales::alpha('black',0.3))

  ## add the fit line from the model for LPS samples
  abline(a = restab['(Intercept)', 'mean'] + restab['treatmentLPS','mean'] - (restab['body_mass_log_sp_std', 'mean'] + restab['body_mass_log_sp_std:treatmentLPS','mean']) * (body_mass_log_center - log(1000)) / body_mass_log_sd,
         b = ln2log10 * restab['evo_allometry_LPS','mean'],
         col = scales::alpha('red',0.8))
  ## and the corresponding credible intervals
  matlines(exp(body_mass_log_center + body_mass_log_sd * seq_bmsp) / 1000, restab[57:106, c('0.025quant','0.975quant')], lty=2, col=scales::alpha('red',0.3))

}
#

# plot the log of the normalized counts as a function of individual body sizes normalized by their species' mean body size
plot_species_dev_raw <- function(dat, restab, plot_legends = TRUE, plot_title = 'Total immune gene expression', font.main = 2) {
  
  lo <- levels(as.factor(dat$species))
  matchlo <- match(names(sort(body_mass_log_sp_means)), lo)
  
  ymin <- min(c(log(dat$count + 0.5) - dat$norm_factor, restab['(Intercept)','mean'], restab['(Intercept)', 'mean'] + restab['treatmentLPS','mean']))
  ymax <- max(c(log(dat$count + 0.5) - dat$norm_factor, restab['(Intercept)','mean'], restab['(Intercept)', 'mean'] + restab['treatmentLPS','mean']))
  
  plot(exp(body_mass_log_diff_sd * dat$body_mass_log_diff_std), 
       log(dat$count + 0.5) - dat$norm_factor, 
       col = dat$treatment, 
       pch = (1:9)[dat$species], 
       main = plot_title,
       xlab = 'log individual mass normalized by species mean body mass', 
       ylab = 'log normalized counts',
       ylim = c(ymin, ymax),
       log = 'x',
       asp = 1, 
       font.main = font.main)
  if(plot_legends) {
    legend(x      = 'topleft', 
           legend = sub('_', ' ', lo[matchlo[1:5]]), 
           pch    = (1:9)[matchlo[1:5]],
           text.font = 3)
    legend(x      = 'bottomright', 
           legend = sub('_', ' ', lo[matchlo[6:9]]), 
           pch    = (1:9)[matchlo[6:9]],
           text.font = 3)
    legend(x      = 'bottomleft',
           title  = 'Slopes',
           legend = c(paste0('LPS(+): β₂+β₃ = ', format(round(restab['intra_allometry_LPS','mean'],2),nsmall=2),          ' (', format(round(restab['intra_allometry_LPS','0.025quant'],2),nsmall=2),          '–', format(round(restab['intra_allometry_LPS','0.975quant'],2),nsmall=2),          ')'),
                      paste0('Δ–LPS: β₃ = ',     format(round(restab['intra_allometry_LPS_response','mean'],2),nsmall=2), ' (', format(round(restab['intra_allometry_LPS_response','0.025quant'],2),nsmall=2), '–', format(round(restab['intra_allometry_LPS_response','0.975quant'],2),nsmall=2), ')'),
                      paste0('LPS(–): β₂ = ',    format(round(restab['intra_allometry_baseline','mean'],2),nsmall=2),     ' (', format(round(restab['intra_allometry_baseline','0.025quant'],2),nsmall=2),     '–', format(round(restab['intra_allometry_baseline','0.975quant'],2),nsmall=2),     ')')),
           col    = c('red','black','white'),
           lty    = 1)
  }
  ## add the fit line from the model for Null samples
  abline(a = restab['(Intercept)','mean'],
         b = ln2log10 * restab['intra_allometry_baseline','mean'],
         col = scales::alpha('black',0.8))
  ## and the corresponding credible intervals
  matlines(exp(body_mass_log_diff_sd * seq_bmdev), restab[157:206, c('0.025quant','0.975quant')], lty=2, col=scales::alpha('black',0.3))

  ## add the fit line from the model for LPS samples
  abline(a = restab['(Intercept)', 'mean'] + restab['treatmentLPS','mean'],
         b = ln2log10 * restab['intra_allometry_LPS','mean'],
         col = scales::alpha('red',0.8))
  ## and the corresponding credible intervals
  matlines(exp(body_mass_log_diff_sd * seq_bmdev), restab[207:256, c('0.025quant','0.975quant')], lty=2, col=scales::alpha('red',0.3))

}
#

# plot the difference between LPS and Null as a function of species' mean sizes
plot_species_mean_diff <- function(dat, restab, plot_legends = TRUE, plot_title = 'Total immune gene expression difference between LPS and Null', font.main = 2) {
  
  lo <- levels(as.factor(dat$species))
  matchlo <- match(names(sort(body_mass_log_sp_means)), lo)
  
  temp <- t(sapply(unique(dat$individual), \(x) {
    c(dat[dat$individual==x & dat$treatment=='LPS','body_mass_log_sp_std'], 
      (log(dat[dat$individual==x & dat$treatment=='LPS','count'] + 0.5) - dat[dat$individual==x & dat$treatment=='LPS','norm_factor']) - 
        (log(dat[dat$individual==x & dat$treatment=='Null','count'] + 0.5) - dat[dat$individual==x & dat$treatment=='Null','norm_factor']))
  }))
  
  plot(exp(body_mass_log_center + body_mass_log_sd * temp[,1]) / 1000, 
       temp[,2],
       pch  = (1:9)[as.factor(sample_data_filt$genus_species)[match(rownames(temp), sample_data_filt$Animal.ID)]], 
       main = plot_title, 
       xlab = 'Species\' mean body mass (kg)', 
       ylab = 'difference in natural log of normalized counts between LPS and Null',
       log  = 'x',
       asp = 1, 
       font.main = font.main)
  if(plot_legends) {
    legend(x      = 'topleft', 
           legend = sub('_', ' ', lo[matchlo]), 
           pch    = (1:9)[matchlo],
           text.font = 3)
  }
  abline(a = restab['treatmentLPS','mean'] - restab['body_mass_log_sp_std:treatmentLPS','mean'] * (body_mass_log_center - log(1000)) / body_mass_log_sd,
         b = ln2log10 * restab['evo_allometry_LPS_response','mean'])
  matlines(exp(body_mass_log_center + body_mass_log_sd * seq_bmsp) / 1000, restab[107:156, c('0.025quant','0.975quant')], lty=2, col = scales::alpha('black',0.3))
  
}
#

# plot the difference between LPS and Null as a function of individual body sizes normalized by their species' mean body size
plot_species_dev_diff <- function(dat, restab, plot_legends = TRUE, plot_title = 'Total immune gene expression difference between LPS and Null', font.main = 2) {
  
  lo <- levels(as.factor(dat$species))
  matchlo <- match(names(sort(body_mass_log_sp_means)), lo)
  
  temp <- t(sapply(unique(dat$individual), \(x) {
    c(dat[dat$individual==x & dat$treatment=='LPS','body_mass_log_diff_std'], 
      (log(dat[dat$individual==x & dat$treatment=='LPS','count'] + 0.5) - dat[dat$individual==x & dat$treatment=='LPS','norm_factor']) - 
        (log(dat[dat$individual==x & dat$treatment=='Null','count'] + 0.5) - dat[dat$individual==x & dat$treatment=='Null','norm_factor']))
  }))
  
  plot(exp(body_mass_log_diff_sd * temp[,1]), 
       temp[,2],
       pch  = (1:9)[as.factor(sample_data_filt$genus_species)[match(rownames(temp), sample_data_filt$Animal.ID)]], 
       main = plot_title, 
       xlab = 'log individual mass normalized by species mean body mass', 
       ylab = 'difference in natural log of normalized counts between LPS and Null',
       log  = 'x',
       asp = 1, 
       font.main = font.main)
  if(plot_legends) {
    legend(x      = 'topleft', 
           legend = sub('_', ' ', lo[matchlo]), 
           pch    = (1:9)[matchlo],
           text.font = 3)
  }
  abline(a = restab['treatmentLPS','mean'],
         b = ln2log10 * restab['intra_allometry_LPS_response','mean'])
  matlines(exp(body_mass_log_diff_sd * seq_bmdev), restab[257:306, c('0.025quant','0.975quant')], lty=2, col=scales::alpha('black',0.3))

}
#
##

## cluster genes by their responses

clustertab <- do.call(rbind, lapply(names(fits), \(x) {
  
  data.frame(row.names      = x, 
             HGNC           = unique(ensembl2ext[ensembl2ext$ensembl_gene_id == x, 'external_gene_name']),
             Deschamps      = unique(ensembl2ext[ensembl2ext$ensembl_gene_id == x, 'module']),
             Null_mean      = fits[[x]]$restab[1, 'mean'],
             response_mean  = fits[[x]]$restab[2, 'mean'],
             LPS_mean       = fits[[x]]$restab[3, 'mean'],
             Null_025       = fits[[x]]$restab[1, '0.025quant'],
             response_025   = fits[[x]]$restab[2, '0.025quant'],
             LPS_025        = fits[[x]]$restab[3, '0.025quant'],
             Null_975       = fits[[x]]$restab[1, '0.975quant'],
             response_975   = fits[[x]]$restab[2, '0.975quant'],
             LPS_975        = fits[[x]]$restab[3, '0.975quant'])
  
}))

clusts <- hclust(dist(clustertab[,3:8], method='euclidean'))
clustertab <- clustertab[clusts$order,]

write.table(clustertab, file = file.path(output_prefix,'05_gene_summaries.txt'), sep='\t', quote=FALSE, row.names=TRUE)

anysig <- apply(clustertab[,6:11], 1, \(x) any(c(x[1:3] > 0, x[4:6] < 0)))
numsig <- sum(anysig)

clustertab_immune <- clustertab[clustertab[,2] != 'no_module_annotation',]
clustertab_nonimmune <- clustertab[clustertab[,2] == 'no_module_annotation',]

immu_contingency <- matrix(c(sum(clustertab_immune$Null_025 > 0), sum(clustertab_immune$Null_025 <= 0 & clustertab_immune$Null_975 >= 0), sum(clustertab_immune$Null_975 < 0),
                             sum(clustertab_immune$response_025 > 0), sum(clustertab_immune$response_025 <= 0 & clustertab_immune$response_975 >= 0), sum(clustertab_immune$response_975 < 0),
                             sum(clustertab_immune$LPS_025 > 0), sum(clustertab_immune$LPS_025 <= 0 & clustertab_immune$LPS_975 >= 0), sum(clustertab_immune$LPS_975 < 0)),
                           nrow = 3,
                           dimnames = list(c('hyper','iso','hypo'),c('baseline','response','lps')))

noni_contingency <- matrix(c(sum(clustertab_nonimmune$Null_025 > 0), sum(clustertab_nonimmune$Null_025 <= 0 & clustertab_nonimmune$Null_975 >= 0), sum(clustertab_nonimmune$Null_975 < 0),
                             sum(clustertab_nonimmune$response_025 > 0), sum(clustertab_nonimmune$response_025 <= 0 & clustertab_nonimmune$response_975 >= 0), sum(clustertab_nonimmune$response_975 < 0),
                             sum(clustertab_nonimmune$LPS_025 > 0), sum(clustertab_nonimmune$LPS_025 <= 0 & clustertab_nonimmune$LPS_975 >= 0), sum(clustertab_nonimmune$LPS_975 < 0)),
                           nrow = 3,
                           dimnames = list(c('hyper','iso','hypo'),c('baseline','response','lps')))

hyper_sigtab_Null <- matrix(c(immu_contingency[1,1], sum(immu_contingency[2:3,1]), noni_contingency[1,1], sum(noni_contingency[2:3,1])), nrow=2, dimnames=list(c('hyper','not_hyper'), c('immune','not_immune')))
iso_sigtab_Null <- matrix(c(immu_contingency[2,1], sum(immu_contingency[c(1,3),1]), noni_contingency[2,1], sum(noni_contingency[c(1,3),1])), nrow=2, dimnames=list(c('iso','not_iso'), c('immune','not_immune')))
hypo_sigtab_Null <- matrix(c(immu_contingency[3,1], sum(immu_contingency[1:2,1]), noni_contingency[3,1], sum(noni_contingency[1:2,1])), nrow=2, dimnames=list(c('hypo','not_hypo'), c('immune','not_immune')))

hyper_sigtab_response <- matrix(c(immu_contingency[1,2], sum(immu_contingency[2:3,2]), noni_contingency[1,2], sum(noni_contingency[2:3,2])), nrow=2, dimnames=list(c('hyper','not_hyper'), c('immune','not_immune')))
iso_sigtab_response <- matrix(c(immu_contingency[2,2], sum(immu_contingency[c(1,3),2]), noni_contingency[2,2], sum(noni_contingency[c(1,3),2])), nrow=2, dimnames=list(c('iso','not_iso'), c('immune','not_immune')))
hypo_sigtab_response <- matrix(c(immu_contingency[3,2], sum(immu_contingency[1:2,2]), noni_contingency[3,2], sum(noni_contingency[1:2,2])), nrow=2, dimnames=list(c('hypo','not_hypo'), c('immune','not_immune')))

hyper_sigtab_LPS <- matrix(c(immu_contingency[1,3], sum(immu_contingency[2:3,3]), noni_contingency[1,3], sum(noni_contingency[2:3,3])), nrow=2, dimnames=list(c('hyper','not_hyper'), c('immune','not_immune')))
iso_sigtab_LPS <- matrix(c(immu_contingency[2,3], sum(immu_contingency[c(1,3),3]), noni_contingency[2,3], sum(noni_contingency[c(1,3),3])), nrow=2, dimnames=list(c('iso','not_iso'), c('immune','not_immune')))
hypo_sigtab_LPS <- matrix(c(immu_contingency[3,3], sum(immu_contingency[1:2,3]), noni_contingency[3,3], sum(noni_contingency[1:2,3])), nrow=2, dimnames=list(c('hypo','not_hypo'), c('immune','not_immune')))

t1 <- fisher.test(hyper_sigtab_Null) ## immune annotated genes are NOT differently likely to be hypermetric in Null
t2 <- fisher.test(iso_sigtab_Null) ## immune annotated genes are NOT differently likely to be isometric in Null
t3 <- fisher.test(hypo_sigtab_Null) ## immune annotated genes are NOT differently likely to be hypometric in Null

t4 <- fisher.test(hyper_sigtab_response) ## immune annotated genes ARE more likely to have hypermetric responses
t5 <- fisher.test(iso_sigtab_response) ## immune annotated genes are NOT differently likely to have isometric responses
t6 <- fisher.test(hypo_sigtab_response) ## immune annotated genes ARE less likely to have hypometric responses

t7 <- fisher.test(hyper_sigtab_LPS) ## immune annotated genes ARE more likely to be hypermetric in LPS
t8 <- fisher.test(iso_sigtab_LPS) ## immune annotated genes ARE less likely to be isometric in LPS
t9 <- fisher.test(hypo_sigtab_LPS) ## immune annotated genes are NOT differently likely to be hypometric in LPS

full_contingency <- data.frame(allometry=rep(c('hyper','iso','hypo'),3), beta=c(rep('baseline',3),rep('response',3),rep('LPS',3)), immune=c(immu_contingency), nonimmune=c(noni_contingency), pvalue=sprintf('%.4f', c(t1$p.value,t2$p.value,t3$p.value,t4$p.value,t5$p.value,t6$p.value,t7$p.value,t8$p.value,t9$p.value)))
write.table(full_contingency, file = file.path(output_prefix,'05_contingency.txt'), sep='\t', quote=FALSE, row.names=FALSE)

##

# genes identified via lit review with specific hypotheses
focal_gene_names <- c("AKT1", "AP1", "CASP4", "CASP8", "CCL2", "CCL3", "CCL4", "CCL5", 
                      "CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CD14", "CD36", "CD40", 
                      "CD40LG", "CD58", "CD69", "CD80", "CD86", "CTLA4", "CXCL8", "CXCR1", 
                      "CXCR2", "EIF2AK2", "FKBP5", "FOS", "FOSL1", "FOSL2", "Gadd45B", 
                      "HSPA1A", "ICAM1", "ICAM2", "ICAM3", "ICAM4", "IFNG", "IL10", 
                      "IL10RA", "IL10RB", "IL12A", "IL12B", "IL15RA", "IL17A", "IL17B", 
                      "IL17C", "IL17F", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL1R2", 
                      "IL1RAP", "IL1RL1", "IL1RL2", "IL21R", "IL23R", "IL27RA", "IL2RA", 
                      "IL2RB", "IL2RG", "IL4R", "IL5RA", "IL6", "IL6R", "IL7R", "IL9R", 
                      "IRAK1", "IRAK4", "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF7", 
                      "IRF8", "ITGAL", "ITGAM", "ITGB2", "JUN", "LIFR", "LILRA1", "LILRB1", 
                      "LILRB2", "LILRB3", "LILRB4", "LRP1", "LY96", "MAP2K1", "MAP4K5", 
                      "MAPK8", "MYD88", "NFAT5", "NFATC1", "NFATC2", "NFATC3", "NFATC4", 
                      "NFKB1", "NFKB2", "NFKBIA", "NFKBIB", "NFKBIE", "NFKBIZ", "PARP14", 
                      "PARP9", "PBEF", "PEBP4", "PTGS2", "PTK2", "PTPRC", "PTPRCAP", 
                      "REL", "RELA", "RELB", "SELL", "SIRPA", "SOD1", "SOD2", "SORT1", 
                      "STAT1", "STAT3", "STAT5", "STAT6", "TBK1", "TLR1", "TLR10", 
                      "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TLR9", 
                      "TNF", "TRAF1", "TRAF6", "TREM1", "VCAM1")

## identify in my orthologs matches to the focal genes
focal_genes <- sapply(focal_gene_names, function(x) {
  searchterm <- paste0('(^|/)', x, '($|/)|(^|/)', gsub('-','',x), '($|/)')
  matches <- grep(searchterm, ensembl2ext$external_gene_name, ignore.case=TRUE)
  if(length(matches) > 0) {
    y <- ensembl2ext[matches[[1]], 'ensembl_gene_id']
  } else {
    y <- NA
  }
  return(y)
})

focal_genes_missing <- names(focal_genes[is.na(focal_genes) | (!focal_genes %in% names(fits))])
focal_genes <- focal_genes[!names(focal_genes) %in% focal_genes_missing]
focal_genes_nonconvergence <- focal_genes[focal_genes %in% names(fits[is.na(sapply(fits,function(x) x$fit))])]
focal_genes <- focal_genes[!focal_genes %in% focal_genes_nonconvergence]

focal_summaries <- do.call(rbind, lapply(names(focal_genes), \(x) {
  cbind(x, focal_genes[[x]], unique(ensembl2ext[ensembl2ext$ensembl_gene_id == focal_genes[[x]],'module']), rownames(fits[[focal_genes[[x]]]]$restab), fits[[focal_genes[[x]]]]$restab[,c('mean','0.025quant','0.975quant')])
}))
colnames(focal_summaries) <- c('gene_name','ensembl_id','module','coefficient','mean','0.025quant','0.975quant')
rownames(focal_summaries) <- NULL

write.table(focal_summaries, file = file.path(output_prefix,'05_focal_summaries.txt'), sep='\t', quote=FALSE, row.names=FALSE)

## make summaries ready to copy in to manuscript, and convert standardized effect sizes to natural effect sizes
formatted_summaries <- t(sapply(unique(focal_summaries$gene_name), \(x) {
  
  patterns <- sapply(c('evo_allometry_baseline', 'evo_allometry_LPS_response', 'evo_allometry_LPS'), \(y) {
    if(sign(prod(focal_summaries$'0.025quant'[focal_summaries$gene_name==x & focal_summaries$coefficient %in% y], focal_summaries$'0.975quant'[focal_summaries$gene_name==x & focal_summaries$coefficient %in% y])) != 1) {
      return(0)
    } else {
      return(ifelse(focal_summaries$mean[focal_summaries$gene_name==x & focal_summaries$coefficient %in% y] > 0, '+', '-'))
    }
  })
  
  c(x,
    unique(focal_summaries$module[focal_summaries$gene_name==x & focal_summaries$coefficient == 'evo_allometry_baseline']),
    paste0(patterns, collapse=','),
    paste0(sprintf('%.2f', focal_summaries$mean[focal_summaries$gene_name==x & focal_summaries$coefficient == 'evo_allometry_baseline']), 
           '[](', 
           sprintf('%.2f', focal_summaries$'0.025quant'[focal_summaries$gene_name==x & focal_summaries$coefficient == 'evo_allometry_baseline']), 
           ', ', 
           sprintf('%.2f', focal_summaries$'0.975quant'[focal_summaries$gene_name==x & focal_summaries$coefficient == 'evo_allometry_baseline']), 
           ')'),
    paste0(sprintf('%.2f', focal_summaries$mean[focal_summaries$gene_name==x & focal_summaries$coefficient == 'evo_allometry_LPS_response']), 
           '[](', 
           sprintf('%.2f', focal_summaries$'0.025quant'[focal_summaries$gene_name==x & focal_summaries$coefficient == 'evo_allometry_LPS_response']), 
           ', ', 
           sprintf('%.2f', focal_summaries$'0.975quant'[focal_summaries$gene_name==x & focal_summaries$coefficient == 'evo_allometry_LPS_response']), 
           ')'),
    paste0(sprintf('%.2f', focal_summaries$mean[focal_summaries$gene_name==x & focal_summaries$coefficient == 'evo_allometry_LPS']), 
           '[](', 
           sprintf('%.2f', focal_summaries$'0.025quant'[focal_summaries$gene_name==x & focal_summaries$coefficient == 'evo_allometry_LPS']), 
           ', ', 
           sprintf('%.2f', focal_summaries$'0.975quant'[focal_summaries$gene_name==x & focal_summaries$coefficient == 'evo_allometry_LPS']), 
           ')'))
}))
colnames(formatted_summaries) <- c('Gene','Deschamps annotation','Pattern','Constitutive Allometry[](β₂)','Response Allometry[](β₃)','Sepsis Allometry[](β₂+β₃)')

formatted_summaries <- formatted_summaries[order(sapply(strsplit(formatted_summaries[,3], ''), \(x) gsub('-',1,gsub('+',-1,paste0(rev(x), collapse='')))), formatted_summaries[,2]),]

write.table(formatted_summaries, file = file.path(output_prefix,'05_formatted_summaries.txt'), sep='\t', quote=FALSE, row.names=FALSE)
## for final table, just replace all instances of '[]' with line breaks, and embolden cells with significance

## plot species phylogeny with individual body sizes
cairo_pdf(file.path(output_prefix,paste0('05_orthologs_phylomass.pdf')), width=4,height=3)
treedat1 <- data.frame(id=sample_data$genus_species, size=log(as.numeric(sample_data$Body.Mass..g.)/1000,2))
treedat1 <- treedat1[sample_data$Genus != 'Microcebus' | sample_data$Animal.ID == 'mmurPool',]
treedat2 <- data.frame(id=names(body_mass_log_sp_means), size=log(exp(body_mass_log_sp_means)/1000,2))

plo <- ggtree::ggtree(species_tree) + ggtree::geom_tiplab()
plo2 <- ggtree::facet_plot(plo, panel='dot', data=treedat1, geom=\(...) ggtree::geom_point(..., position=ggplot2::position_jitter(0,0.3)), pch=1, ggtree::aes(x=size)) + ggtree::geom_facet(data=treedat2, geom=ggtree::geom_point, pch=3, col='red', panel='dot', ggtree::aes(x=size))
plo2 + ggtree::theme_tree2() + ggplot2::scale_x_continuous(breaks=log(c(0.05, 1, 20, 125),2), labels=c('0.05', '1.00', '20.0', '125'))
dev.off()

## plot log normalized counts against species mean mass
cairo_pdf(file.path(output_prefix,'05_total_orthologous_immune.pdf'), width=8, height=8)
plot_species_mean_raw(dat_ortho, ortho_immune_restab)
dev.off()
cairo_pdf(file.path(output_prefix,'05_total_immune.pdf'), width=8, height=8)
plot_species_mean_raw(dat_og, og_immune_restab)
dev.off()

## plot on natural scale
dat_ortho$counts_norm <- exp(log(dat_ortho$count + 0.5) - dat_ortho$norm_factor)
dat_ortho_agg <- aggregate(dat_ortho, by = list(species = dat_ortho$species, treatment = dat_ortho$treatment), FUN = mean)
cairo_pdf(file.path(output_prefix,'05_total_orthologous_immune_naturalScale.pdf'), width=8, height=8)
plot(exp(body_mass_log_center + body_mass_log_sd * dat_ortho_agg$body_mass_log_sp_std - log(1000)), 
     dat_ortho_agg$counts_norm, 
     col  = dat_ortho_agg$treatment, 
     xlab = 'Species\' mean body mass (kg)', 
     ylab = 'Normalized counts', 
     pch  = (1:9)[dat_ortho_agg$species])
dev.off()

cairo_pdf(file.path(output_prefix,'05_total_orthologous_immune_naturalScaleX.pdf'), width=8, height=8)
plot(exp(body_mass_log_center + body_mass_log_sd * dat_ortho_agg$body_mass_log_sp_std - log(1000)), 
     log(dat_ortho_agg$counts_norm), 
     col  = dat_ortho_agg$treatment, 
     xlab = 'Species\' mean body mass (kg)', 
     ylab = 'Normalized counts', 
     pch  = (1:9)[dat_ortho_agg$species])
dev.off()


## plot diff in log normalized counts against species mean mass
cairo_pdf(file.path(output_prefix,'05_total_orthologous_immune_diff.pdf'), width=8, height=8)
plot_species_mean_diff(dat_ortho, ortho_immune_restab)
dev.off()
cairo_pdf(file.path(output_prefix,'05_total_immune_diff.pdf'), width=8, height=8)
plot_species_mean_diff(dat_og, og_immune_restab)
dev.off()

cairo_pdf(file.path(output_prefix,'05_total_orthologous_immune_intra.pdf'), width=8, height=8)
plot_species_dev_raw(dat_ortho, ortho_immune_restab)
dev.off()
cairo_pdf(file.path(output_prefix,'05_total_immune_intra.pdf'), width=8, height=8)
plot_species_dev_raw(dat_og, og_immune_restab)
dev.off()

cairo_pdf(file.path(output_prefix,'05_total_orthologous_immune_intra_diff.pdf'), width=8, height=8)
plot_species_dev_diff(dat_ortho, ortho_immune_restab)
dev.off()
cairo_pdf(file.path(output_prefix,'05_total_immune_intra_diff.pdf'), width=8, height=8)
plot_species_dev_diff(dat_og, og_immune_restab)
dev.off()

cairo_pdf(file.path(output_prefix,'05_total_orthologous_non-immune.pdf'), width=8, height=8)
plot_species_mean_raw(dat_ortho_nonimmune, ortho_nonimmune_restab, plot_title = 'Total non-immune gene expression')
dev.off()

cairo_pdf(file.path(output_prefix,'05_total_orthologous_non-immune_inset.pdf'), width=3.5, height=4)
plot_species_mean_raw(dat_ortho_nonimmune, ortho_nonimmune_restab, plot_legends = FALSE, plot_axes=FALSE)
dev.off()

gene <- 1
for(x in 1:ceiling(length(focal_genes)/7)) {
  cairo_pdf(file.path(output_prefix,paste0('05_raw_focal_genes_',x,'.pdf')), width=16, height=8)
  par(mfrow=c(2,4))
  for(g in names(focal_genes)[gene:min((gene+6),length(focal_genes))]) {

    curgene <- focal_genes[[g]]
      
    ## plot the gene abundance as function of species mean size
    plot_species_mean_raw(fits[[curgene]]$dat, fits[[curgene]]$restab, FALSE, g)

  }
  plot.new()
  lo <- levels(as.factor(fits[[1]]$dat$species))
  matchlo <- match(names(sort(body_mass_log_sp_means)), lo)
  legend(x = 'bottomleft', legend = c('Null','LPS'), col = c('black','red'), lty = 1)
  legend(x = 'topleft', legend = lo[matchlo], pch = (1:9)[matchlo])
  dev.off()
  
  cairo_pdf(file.path(output_prefix,paste0('05_raw_focal_genes_dev_',x,'.pdf')), width=16, height=8)
  par(mfrow=c(2,4))
  for(g in names(focal_genes)[gene:min((gene+6),length(focal_genes))]) {

    curgene <- focal_genes[[g]]
      
    ## plot the gene abundance as function of species mean size
    plot_species_dev_raw(fits[[curgene]]$dat, fits[[curgene]]$restab, FALSE, g)

  }
  plot.new()
  lo <- levels(as.factor(fits[[1]]$dat$species))
  matchlo <- match(names(sort(body_mass_log_sp_means)), lo)
  legend(x = 'bottomleft', legend = c('Null','LPS'), col = c('black','red'), lty = 1)
  legend(x = 'topleft', legend = lo[matchlo], pch = (1:9)[matchlo])
  dev.off()

  gene <- gene + 7
}
##

## figure 3

cairo_pdf(file.path(output_prefix, '05_figure_3_raw.pdf'), width=9, height=6.5)
par(mfrow=c(2,3), pty='s')

## plot the gene abundance as function of species mean size
plot_species_mean_raw(fits[[focal_genes[['TLR4']]]]$dat, fits[[focal_genes[['TLR4']]]]$restab, FALSE, 'TLR4', font.main=4)
legend('topleft', "D.", bty='n', text.font=2, x.intersp=0)

plot.new()
lo <- levels(as.factor(fits[[1]]$dat$species))
matchlo <- match(names(sort(body_mass_log_sp_means)), lo)
legend(x = 'bottomleft', legend = sub('_', ' ', lo[matchlo]), pch = (1:9)[matchlo], text.font = 3)
legend(x = 'bottomright', legend = c('Null','LPS'), col = c('black','red'), lty = 1)

plot_species_mean_raw(fits[[focal_genes[['CTLA4']]]]$dat, fits[[focal_genes[['CTLA4']]]]$restab, FALSE, 'CTLA4', font.main=4)
legend('topleft', "E.", bty='n', text.font=2, x.intersp=0)
for(g in c('IFNG','IL6','IL1B')) {

  curgene <- focal_genes[[g]]
    
  ## plot the gene abundance as function of species mean size
  plot_species_mean_raw(fits[[curgene]]$dat, fits[[curgene]]$restab, FALSE, g, font.main=4)
  legend('topleft', c(IFNG='F.', IL6='G.', IL1B='H.')[[g]], bty='n', text.font=2, x.intersp=0)

}
dev.off()

##


## make various plots for each focal gene
for(g in names(focal_genes)) {
  
  curgene <- focal_genes[[g]]
  
  lo <- levels(as.factor(fits[[curgene]]$dat$species))
  
  ## plot the gene abundance as function of species mean size
  cairo_pdf(file.path(output_prefix, paste0('05_',g,'_species_means.pdf')), width=8, height=8)
  plot_species_mean_raw(fits[[curgene]]$dat, fits[[curgene]]$restab, FALSE, g)
  dev.off()
  
  ## plot gene abundance as function of difference between individual size and species mean size
  cairo_pdf(file.path(output_prefix, paste0('05_',g,'_individual_norm_by_species_means.pdf')), width=8, height=8)
  plot_species_dev_raw(fits[[curgene]]$dat, fits[[curgene]]$restab, FALSE, g)
  dev.off()
  
  ## plot gene abundance as function of individual size
  cairo_pdf(file.path(output_prefix, paste0('05_',g,'_individual_mass.pdf')), width=8, height=8)
  plot(body_mass_log_diff_sd * fits[[curgene]]$dat$body_mass_log_diff_std + body_mass_log_sd * fits[[curgene]]$dat$body_mass_log_sp_std, 
       log(fits[[curgene]]$dat$count) - fits[[curgene]]$dat$norm_factor, 
       col  = fits[[curgene]]$dat$treatment, 
       pch  = (1:nlevels(fits[[curgene]]$dat$species))[fits[[curgene]]$dat$species], 
       main = paste0(g, ' gene expression'), 
       xlab = 'log individual body size', 
       ylab = 'log normalized counts')
  legend(x      = 'topleft', 
         legend = lo[match(names(sort(body_mass_log_sp_means)), lo)], 
         pch    = (1:length(lo))[match(names(sort(body_mass_log_sp_means)), lo)])
  dev.off()
    
  ## plot the difference between LPS and Null as a function of species' mean sizes
  cairo_pdf(file.path(output_prefix, paste0('05_',g,'_species_means_LPS_vs_Null.pdf')), width=8, height=8)
  plot_species_mean_diff(fits[[curgene]]$dat, fits[[curgene]]$restab, FALSE, g)
  dev.off()
  ##
  
  ## plot the difference between LPS and Null as a function of the difference between the individual size and the species' mean size
  cairo_pdf(file.path(output_prefix, paste0('05_',g,'_individual_norm_by_species_means_LPS_vs_Null.pdf')), width=8, height=8)
  plot_species_dev_diff(fits[[curgene]]$dat, fits[[curgene]]$restab, FALSE, g)
  dev.off()
  ##
  
}
