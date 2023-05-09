
nthreads <- 7
raw_data_prefix <- path.expand('raw_data/20221215_primate_allometry/')
output_prefix <- path.expand('outputs/primates_20230314_mixed/')

ensembl_species <- c('Callithrix_jacchus', 'Homo_sapiens', 'Macaca_mulatta', 'Microcebus_murinus', 'Papio_anubis', 'Pongo_abelii')
ncbi_species <- c('Daubentonia_madagascariensis', 'Lemur_catta', 'Sapajus_apella')

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

martlist <- list(Callithrix_jacchus='cjacchus_gene_ensembl', Homo_sapiens='hsapiens_gene_ensembl', Macaca_mulatta='mmulatta_gene_ensembl', Microcebus_murinus='mmurinus_gene_ensembl', Papio_anubis='panubis_gene_ensembl', Pongo_abelii='pabelii_gene_ensembl')

orthologs <- read.table(Sys.glob(file.path(output_prefix,'02_find_orthologs','of_out','*','Orthogroups','Orthogroups.tsv')), header=TRUE, row.names = 1, sep='\t')
orthologs_1_1 <- orthologs[!apply(orthologs, 1, \(x) any(grepl(',',x) | (nchar(x) == 0))),]
orthologs_1_1[,colnames(orthologs_1_1) %in% ncbi_species] <- apply(orthologs_1_1[,colnames(orthologs_1_1) %in% ncbi_species], 2, \(x) sapply(x, \(y) paste(strsplit(y,'.',fixed=TRUE)[[1]][1:2],collapse='.')))
orthologs_1_1[,ensembl_species] <- sapply(ensembl_species, \(x) {
  for(i in 1:numtries) {
    try({
    ee <- biomaRt::getBM(attributes = c('ensembl_peptide_id_version', 'ensembl_gene_id'), 
                         filters = 'ensembl_peptide_id_version', 
                         values = orthologs_1_1[,x], 
                         mart = biomaRt::useMart('ENSEMBL_MART_ENSEMBL', martlist[[x]], version='Ensembl Genes 109'), 
                         uniqueRows = TRUE)
    break
    }, silent = FALSE)
  }
  rownames(ee) <- ee[,'ensembl_peptide_id_version']
  return(ee[orthologs_1_1[,x], 'ensembl_gene_id'])
})
orthologs_1_1 <- orthologs_1_1[!orthologs_1_1[,'Homo_sapiens'] %in% human_globin_gene_ids,]

for(i in 1:numtries) {
  try({
    ensembl2ext <- biomaRt::getBM(attributes = c('ensembl_peptide_id_version', 'ensembl_gene_id', 'external_gene_name', 'go_id'), 
                                  filters = 'ensembl_peptide_id_version', 
                                  values = unlist(strsplit(orthologs[,'Homo_sapiens'],', ')), 
                                  mart = biomaRt::useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', version='Ensembl Genes 109'), 
                                  uniqueRows = TRUE)
  break
  }, silent = FALSE)
}

## feels weird to round decimals for poisson error, but since data were counts at one point, zeros are possible, so doesn't make sense to just log-transform and use gaussian error; and error probably still scales like poisson such that low counts are less meaningful (and the decimal values rounded off would contribute negligible information) (cite Z1000 tximport paper)
raw_txi <- c(lapply(ncbi_species, \(x) {
  cat(paste0('Importing ', x, ' counts\n'))
  cfiles <- countfiles[grep(x, names(countfiles), ignore.case = TRUE)]
  tnames <- unique(do.call(rbind, lapply(cfiles, read.table, header=TRUE, sep='\t'))$Name)
  tx2gene <- data.frame(transcript=tnames, gene=sapply(tnames, \(y) paste(strsplit(y,'.',fixed=TRUE)[[1]][1:2],collapse='.')))
  txi <- tximport::tximport(cfiles, type = 'salmon', tx2gene = tx2gene)
  return(txi)
}), lapply(ensembl_species, \(x) {
      cat(paste0('Retrieving ensembl genes for ', x, ' \n'))
      cfiles <- countfiles[grep(x, names(countfiles), ignore.case = TRUE)]
      tnames <- unique(do.call(rbind, lapply(cfiles, read.table, header=TRUE, sep='\t'))$Name)
      for(i in 1:numtries) {
        try({
          ee <- biomaRt::getBM(attributes = c('ensembl_transcript_id_version','ensembl_gene_id'), 
                               filters    = 'ensembl_transcript_id_version', 
                               values     = tnames, 
                               mart       = biomaRt::useMart('ENSEMBL_MART_ENSEMBL', martlist[[x]], version='Ensembl Genes 109'), 
                               uniqueRows = TRUE)
          break
        }, silent = FALSE)
      }
      tx2gene <- ee[,c('ensembl_transcript_id_version', 'ensembl_gene_id')]
      txi <- tximport::tximport(cfiles, type = 'salmon', tx2gene = tx2gene)
      return(txi)
}))
names(raw_txi) <- c(ncbi_species, ensembl_species)
##

## import reference data for species average body sizes
body_size_ref <- read.csv(file.path(raw_data_prefix, 'gyz043_suppl_Supplement_Data.csv'))
body_size_ref$genus_species <- paste(body_size_ref$genus, body_size_ref$species, sep='_')
body_size_ref$genus_species[body_size_ref$genus_species == 'Cebus_apella'] <- 'Sapajus_apella' 
##

## calculate differences of individual sizes from species means and log-transform
sample_data_filt$body_mass_log <- log(as.numeric(sample_data_filt$Body.Mass..g.))

body_mass_log_sp_means <- sapply(unique(sample_data_filt$genus_species), \(x) log(as.numeric(body_size_ref$Mean_body_mass_g[body_size_ref$genus_species == x])))

sample_data_filt$body_mass_log_sp <- body_mass_log_sp_means[sample_data_filt$genus_species]
sample_data_filt$body_mass_log_diff <- sample_data_filt$body_mass_log - sample_data_filt$body_mass_log_sp
##

## calculate ancestral body size (somewhat arbitrary for our purposes but could help interpret 'main effects')
body_mass_log_center <- ape::ace(body_mass_log_sp_means, species_tree_norm)$ace[[1]]
##

## standardize masses for model input
body_mass_log_sd <- sd(sample_data_filt$body_mass_log_sp)
body_mass_log_diff_sd <- sd(sample_data_filt$body_mass_log_diff)

sample_data_filt$body_mass_log_sp_std <- (sample_data_filt$body_mass_log_sp - body_mass_log_center) / body_mass_log_sd
sample_data_filt$body_mass_log_diff_std <- sample_data_filt$body_mass_log_diff / body_mass_log_diff_sd
##

ortholog_mat <- orthologs_1_1[,names(raw_txi)]
ortholog_mat <- ortholog_mat[ortholog_mat[,1] %in% rownames(raw_txi[[1]]$abundance), ]
txi_ortho <- list(abundance = raw_txi[[1]]$abundance[ortholog_mat[,1],], 
                  counts    = raw_txi[[1]]$counts[ortholog_mat[,1],],
                  length    = raw_txi[[1]]$length[ortholog_mat[,1],],
                  countsFromAbundance = raw_txi[[1]]$countsFromAbundance)
for(i in 2:length(raw_txi)) {
  ortholog_mat <- ortholog_mat[ortholog_mat[,i] %in% rownames(raw_txi[[i]]$abundance), ]

  txi_ortho$abundance <- cbind(txi_ortho$abundance[ortholog_mat[,1],], raw_txi[[i]]$abundance[ortholog_mat[,i],])
  txi_ortho$counts <- cbind(txi_ortho$counts[ortholog_mat[,1],], raw_txi[[i]]$counts[ortholog_mat[,i],])
  txi_ortho$length <- cbind(txi_ortho$length[ortholog_mat[,1],], raw_txi[[i]]$length[ortholog_mat[,i],])
} 
rownames(txi_ortho$abundance) <- rownames(txi_ortho$counts) <- rownames(txi_ortho$length) <- ortholog_mat[match(rownames(txi_ortho$length), ortholog_mat[,1]), 'Homo_sapiens']

## get per-gene-and-sample normalization factors for 1:1 orthologs, and per-sample size factors based on orthologs

des <- DESeq2:::DESeqDataSetFromTximport(txi_ortho, data.frame(Int=rep(1,ncol(txi_ortho$counts))), ~1)
des <- DESeq2::estimateSizeFactors(des)

sfs <- apply(SummarizedExperiment::assay(des, 'normalizationFactors'), 2, \(x) exp(mean(log(x))))

##

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
dat_ortho <- dat_ortho[dat_ortho$individual %in% sample_data_filt$Animal.ID,]
dat_ortho$body_mass_log_sp_std <- sapply(dat_ortho$individual, function(z) sample_data_filt$body_mass_log_sp_std[sample_data_filt$Animal.ID==z])
dat_ortho$body_mass_log_diff_std <- sapply(dat_ortho$individual, function(z) sample_data_filt$body_mass_log_diff_std[sample_data_filt$Animal.ID==z])

## define phyr formula
phy_formula <- count ~ offset(norm_factor) + body_mass_log_sp_std * treatment + body_mass_log_diff_std + body_mass_log_diff_std:treatment + (1 | individual) + (1 | species__) + (1 | species@treatment) + (1 | species__@treatment)

## define linear combinations of effects of interest
lc1 <- INLA::inla.make.lincomb('body_mass_log_sp_std'                = 1/body_mass_log_sd) 
names(lc1) = "evo_allometry_baseline"

lc2 <- INLA::inla.make.lincomb('body_mass_log_sp_std:treatmentLPS'   = 1/body_mass_log_sd) 
names(lc2) = "evo_allometry_LPS_response"

lc3 <- INLA::inla.make.lincomb('body_mass_log_sp_std'                = 1/body_mass_log_sd, 
                               'body_mass_log_sp_std:treatmentLPS'   = 1/body_mass_log_sd) 
names(lc3) = "evo_allometry_LPS"

lc4 <- INLA::inla.make.lincomb('body_mass_log_diff_std'              = 1/body_mass_log_diff_sd)
names(lc4) = "intra_allometry_baseline"

lc5 <- INLA::inla.make.lincomb('treatmentLPS:body_mass_log_diff_std' = 1/body_mass_log_diff_sd)
names(lc5) = "intra_allometry_LPS_response"

lc6 <- INLA::inla.make.lincomb('body_mass_log_diff_std'              = 1/body_mass_log_diff_sd, 
                               'treatmentLPS:body_mass_log_diff_std' = 1/body_mass_log_diff_sd)
names(lc6) = "intra_allometry_LPS"


cat('Modeling immune orthologs\n')
ortho_immune_res <- phyr::pglmm(formula       = phy_formula, 
                                family        = "poisson", 
                                cov_ranef     = list(species = species_tree_norm), 
                                data          = dat_ortho, 
                                bayes         = TRUE, 
                                bayes_options = list(lincomb = c(lc1,lc2,lc3,lc4,lc5,lc6)))
ortho_immune_restab <- ortho_immune_res$inla.model$summary.lincomb.derived
##

## do same but for non-immune, to show it's not a general artifact of the summarization
geneGroups <- rownames(txi_ortho$abundance)
geneGroups[!geneGroups %in% immuneGenes] <- 'nonimmune'
txi_ortho_nonimmune_pool <- txi_ortho
txi_ortho_nonimmune_pool$abundance <- rowsum(txi_ortho$abundance, geneGroups)
txi_ortho_nonimmune_pool$counts <- rowsum(txi_ortho$counts, geneGroups)
txi_ortho_nonimmune_pool$length <- rowsum(txi_ortho$abundance * txi_ortho$length, geneGroups) / txi_ortho_nonimmune_pool$abundance ## weighted arithmetic mean length
txi_ortho_nonimmune_pool$length[rownames(txi_ortho_nonimmune_pool$length) != 'nonimmune',] <- txi_ortho$length[rownames(txi_ortho_nonimmune_pool$length)[rownames(txi_ortho_nonimmune_pool$length) != 'nonimmune'],]

## get normalization factors
des_ortho_nonimmune_pool <- DESeq2:::DESeqDataSetFromTximport(txi_ortho_nonimmune_pool, data.frame(Int=rep(1,ncol(txi_ortho_nonimmune_pool$counts))), ~1)
des_ortho_nonimmune_pool <- DESeq2::estimateSizeFactors(des_ortho_nonimmune_pool)

## model total differential expression of non-immune genes
dat_ortho_nonimmune <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des_ortho_nonimmune_pool)), '_'), \(x) x[[3]]), 
                                  species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des_ortho_nonimmune_pool)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                                  treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des_ortho_nonimmune_pool)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                                  norm_factor = log(DESeq2::normalizationFactors(des_ortho_nonimmune_pool)['nonimmune',]),
                                  count       = DESeq2::counts(des_ortho_nonimmune_pool)['nonimmune',])
dat_ortho_nonimmune <- dat_ortho_nonimmune[dat_ortho_nonimmune$individual %in% sample_data_filt$Animal.ID,]
dat_ortho_nonimmune$body_mass_log_sp_std <- sapply(dat_ortho_nonimmune$individual, function(z) sample_data_filt$body_mass_log_sp_std[sample_data_filt$Animal.ID==z])
dat_ortho_nonimmune$body_mass_log_diff_std <- sapply(dat_ortho_nonimmune$individual, function(z) sample_data_filt$body_mass_log_diff_std[sample_data_filt$Animal.ID==z])

ortho_nonimmune_res <- phyr::pglmm(formula       = phy_formula, 
                                   family        = "poisson", 
                                   cov_ranef     = list(species = species_tree_norm), 
                                   data          = dat_ortho_nonimmune, 
                                   bayes         = TRUE, 
                                   bayes_options = list(lincomb = c(lc1,lc2,lc3,lc4,lc5,lc6)))
ortho_nonimmune_restab <- ortho_nonimmune_res$inla.model$summary.lincomb.derived
##

## import counts for entire orthogroups, not just 1:1 orthologs
txi_all_tmp <- list()
for(x in ncbi_species) {
  
  cfiles <- countfiles[grep(x, names(countfiles), ignore.case = TRUE)]
  tnames <- unique(do.call(rbind, lapply(cfiles, read.table, header=TRUE, sep='\t'))$Name)
  tmp <- data.frame(transcript=tnames, gene.species=sapply(tnames, \(y) paste(strsplit(y,'.',fixed=TRUE)[[1]][1:2],collapse='.')))
  tmp$gene <- sapply(tmp$gene.species, \(y) {
    tmpidx <- grep(paste0(y,'.'), orthologs[,x],fixed=TRUE)
    if(length(tmpidx) == 1) {
      return(rownames(orthologs)[[tmpidx]])
    } else {
      return(y)
    }
  })
  tx2gene <- data.frame(transcript = tmp$transcript, gene=tmp$gene)

  txi_all_tmp[[x]] <- tximport::tximport(cfiles, type = 'salmon', tx2gene = tx2gene)
  
}

for(x in ensembl_species) {
  
  cfiles <- countfiles[grep(x, names(countfiles), ignore.case = TRUE)]
  tnames <- unique(do.call(rbind, lapply(cfiles, read.table, header=TRUE, sep='\t'))$Name)
  
  for(i in 1:numtries) {
    try({
      ee_t <- biomaRt::getBM(attributes = c('ensembl_transcript_id_version','ensembl_gene_id'), 
                             filters    = 'ensembl_transcript_id_version', 
                             values     = tnames, 
                             mart       = biomaRt::useMart('ENSEMBL_MART_ENSEMBL', martlist[[x]], version='Ensembl Genes 109'), 
                             uniqueRows = TRUE)
      break
    }, silent = FALSE)
  }
  
  for(i in 1:numtries) {
    try({
      ee_p <- biomaRt::getBM(attributes = c('ensembl_peptide_id_version', 'ensembl_gene_id'), 
                             filters = 'ensembl_peptide_id_version', 
                             values = unlist(strsplit(orthologs[,x],', ')), 
                             mart = biomaRt::useMart('ENSEMBL_MART_ENSEMBL', martlist[[x]], version='Ensembl Genes 109'), 
                             uniqueRows = TRUE)
    break
    }, silent = FALSE)
  }
  
  ee_p$og <- sapply(ee_p$ensembl_peptide_id_version, \(y) {
    tmpidx <- grep(paste0(y,'\\b'), orthologs[,x])
    if(length(tmpidx) == 1) {
      return(rownames(orthologs)[[tmpidx]])
    } else {
      return(y)
    }
  })
  
  ee <- merge(ee_t,ee_p, by = 'ensembl_gene_id', all.x = TRUE)
  
  tx2gene <- data.frame(transcript = ee$ensembl_transcript_id_version, gene = ee$og)
  
  tx2gene$gene[is.na(tx2gene$gene)] <- tx2gene$transcript[is.na(tx2gene$gene)]
  
  txi_all_tmp[[x]] <- tximport::tximport(cfiles, type = 'salmon', tx2gene = tx2gene)

}

txi_all <- txi_all_tmp[[1]]
for(s in 2:length(txi_all_tmp)) {
  
  txi_all[[1]] <- merge(txi_all[[1]], txi_all_tmp[[s]][[1]], by='row.names')
  rownames(txi_all[[1]]) <- txi_all[[1]][,1]
  txi_all[[1]] <- as.matrix(txi_all[[1]][,-1])
  txi_all[[2]] <- merge(txi_all[[2]], txi_all_tmp[[s]][[2]], by='row.names')
  rownames(txi_all[[2]]) <- txi_all[[2]][,1]
  txi_all[[2]] <- as.matrix(txi_all[[2]][,-1])
  txi_all[[3]] <- merge(txi_all[[3]], txi_all_tmp[[s]][[3]], by='row.names')
  rownames(txi_all[[3]]) <- txi_all[[3]][,1]
  txi_all[[3]] <- as.matrix(txi_all[[3]][,-1])
  
}

## create txi_all_immune_pool analogous to above

immunePeps <- unique(ensembl2ext$ensembl_peptide_id_version[ensembl2ext$module != 'no_module_annotation'])
immune_OGs <- rownames(orthologs)[sapply(rownames(orthologs), \(x) any(unlist(strsplit(orthologs[x,'Homo_sapiens'],', ')) %in% immunePeps))]

## pool immune-annotated orthogroups
geneGroups <- rownames(txi_all$abundance)
geneGroups[geneGroups %in% immune_OGs] <- 'immune' ## change this part so that if any member of an orthogroup is immune, the whole thing is
txi_og_immune_pool <- txi_all
txi_og_immune_pool$abundance <- rowsum(txi_all$abundance, geneGroups)
txi_og_immune_pool$counts <- rowsum(txi_all$counts, geneGroups)
txi_og_immune_pool$length <- rowsum(txi_all$abundance * txi_all$length, geneGroups) / txi_og_immune_pool$abundance ## weighted arithmetic mean length
txi_og_immune_pool$length[rownames(txi_og_immune_pool$length) != 'immune',] <- txi_all$length[rownames(txi_og_immune_pool$length)[rownames(txi_og_immune_pool$length) != 'immune'],]
##

des_og <- DESeq2:::DESeqDataSetFromTximport(txi_og_immune_pool, data.frame(Int=rep(1,ncol(txi_og_immune_pool$counts))), ~1)

## model total differential expression of immune genes
dat_og <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des_og)), '_'), \(x) x[[3]]), 
                     species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des_og)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                     treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des_og)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                     norm_factor = log(sfs[colnames(DESeq2::counts(des_og))]),
                     count       = DESeq2::counts(des_og)['immune',])
dat_og <- dat_og[dat_og$individual %in% sample_data_filt$Animal.ID,]
dat_og$body_mass_log_sp_std <- sapply(dat_og$individual, function(z) sample_data_filt$body_mass_log_sp_std[sample_data_filt$Animal.ID==z])
dat_og$body_mass_log_diff_std <- sapply(dat_og$individual, function(z) sample_data_filt$body_mass_log_diff_std[sample_data_filt$Animal.ID==z])

cat('Modeling immune orthogroups\n')
og_immune_res <- phyr::pglmm(formula       = phy_formula, 
                             family        = "poisson", 
                             cov_ranef     = list(species = species_tree_norm), 
                             data          = dat_og, 
                             bayes         = TRUE, 
                             bayes_options = list(lincomb = c(lc1,lc2,lc3,lc4,lc5,lc6)))
og_immune_restab <- og_immune_res$inla.model$summary.lincomb.derived
##

##

##
sensors <- unique(ensembl2ext$ensembl_gene_id[ensembl2ext$module == 'sensor'])
effectors <- unique(ensembl2ext$ensembl_gene_id[ensembl2ext$module == 'effector'])
bureaucracy <- unique(ensembl2ext$ensembl_gene_id[!ensembl2ext$module %in% c('sensor','effector','no_module_annotation')])

## pool immune-annotated 1:1 orthologous genes
geneGroups <- rownames(txi_ortho$abundance)
geneGroups[geneGroups %in% sensors] <- 'sensor'
geneGroups[geneGroups %in% effectors] <- 'effector'
geneGroups[geneGroups %in% bureaucracy] <- 'bureaucracy'

txi_modules_pool <- txi_ortho
txi_modules_pool$abundance <- rowsum(txi_ortho$abundance, geneGroups)
txi_modules_pool$counts <- rowsum(txi_ortho$counts, geneGroups)
txi_modules_pool$length <- rowsum(txi_ortho$abundance * txi_ortho$length, geneGroups) / txi_modules_pool$abundance ## weighted arithmetic mean length
txi_modules_pool$length[!rownames(txi_modules_pool$length) %in% c('sensor','effector','bureaucracy'),] <- txi_ortho$length[rownames(txi_modules_pool$length)[!rownames(txi_modules_pool$length) %in% c('sensor','effector','bureaucracy')],]

## get normalization factors
des_modules_pool <- DESeq2:::DESeqDataSetFromTximport(txi_modules_pool, data.frame(Int=rep(1,ncol(txi_modules_pool$counts))), ~1)
des_modules_pool <- DESeq2::estimateSizeFactors(des_modules_pool)

module_fits <- lapply(c('sensor','effector','bureaucracy'), function(x) {
## model total differential expression of immune genes
  dat <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des_modules_pool)), '_'), \(x) x[[3]]), 
                    species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des_modules_pool)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                    treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des_modules_pool)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                    norm_factor = log(DESeq2::normalizationFactors(des_modules_pool)[x,]),
                    count       = DESeq2::counts(des_modules_pool)[x,])
  dat <- dat[dat$individual %in% sample_data_filt$Animal.ID,]
  dat$body_mass_log_sp_std <- sapply(dat$individual, function(z) sample_data_filt$body_mass_log_sp_std[sample_data_filt$Animal.ID==z])
  dat$body_mass_log_diff_std <- sapply(dat$individual, function(z) sample_data_filt$body_mass_log_diff_std[sample_data_filt$Animal.ID==z])

  res <- phyr::pglmm(formula       = phy_formula, 
                     family        = "poisson", 
                     cov_ranef     = list(species = species_tree_norm), 
                     data          = dat, 
                     bayes         = TRUE, 
                     bayes_options = list(lincomb = c(lc1,lc2,lc3,lc4,lc5,lc6)))
  restab <- res$inla.model$summary.lincomb.derived
  
  return(list(dat=dat, res=res, restab=restab))

})
names(module_fits) <- c('sensor','effector','bureaucracy')
##

## fit per-gene models
cat('Fitting models\n')
numtriesfit <- 15
fits <- parallel::mclapply(rownames(des), function(x) {
  cat(round(which(x == rownames(des)) / nrow(des),2) * 100, '%: ', x,'\n')
  dat <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des)), '_'), \(x) x[[3]]), 
                    species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                    treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                    norm_factor = log(DESeq2::normalizationFactors(des)[x,]),
                    count       = DESeq2::counts(des)[x,])
  dat <- dat[dat$individual %in% sample_data_filt$Animal.ID,]
  dat$body_mass_log_sp_std <- sapply(dat$individual, function(z) sample_data_filt$body_mass_log_sp_std[sample_data_filt$Animal.ID==z])
  dat$body_mass_log_diff_std <- sapply(dat$individual, function(z) sample_data_filt$body_mass_log_diff_std[sample_data_filt$Animal.ID==z])
  
  if(sd(DESeq2::counts(des)[x,]) > 0) {
    
    fit <- NA
    for(i in 1:numtriesfit) {
      try({
        
        fit <- phyr::pglmm(formula       = phy_formula, 
                           family        = "poisson", 
                           cov_ranef     = list(species = species_tree_norm), 
                           data          = dat, 
                           bayes         = TRUE, 
                           bayes_options = list(lincomb = c(lc1,lc2,lc3,lc4,lc5,lc6)))
        
        if(!'(Intercept)' %in% rownames(fit$inla.model$summary.lincomb.derived)) {
          stop(paste0('Gene ', rownames(des)[[x]], ' has incomplete summary on try ', i))
        }
        
        break
        
      }, silent = FALSE)
    }
    
    if(!all(is.na(fit))) {
      restab <- fit$inla.model$summary.lincomb.derived
      fit$inla.model <- NULL
    } else {
      restab <- NA
    }
    
  } else {
    
    fit <- NA
    restab <- NA
    
  }
  
  return(list(data=dat, fit=fit, restab=restab))
  
}, mc.cores = nthreads)
names(fits) <- rownames(des)

save.image(file.path(output_prefix, '04_phyr_results.RData'))

