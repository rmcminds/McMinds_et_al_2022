
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

## identify database names for all species
martlist <- list(Callithrix_jacchus='cjacchus_gene_ensembl', Homo_sapiens='hsapiens_gene_ensembl', Macaca_mulatta='mmulatta_gene_ensembl', Microcebus_murinus='mmurinus_gene_ensembl', Papio_anubis='panubis_gene_ensembl', Pongo_abelii='pabelii_gene_ensembl')

## create globin gene list for later filtering (these were targeted for removal during library prep so are likely unreliable, and we don't want them to influence e.g. size factors)
# https://doi.org/10.1038/s41598-020-62801-6 Table S1
human_globin_gene_ids <- c('ENSG00000206172', 'ENSG00000188536', 'ENSG00000244734', 'ENSG00000229988', 'ENSG00000223609', 'ENSG00000213931', 'ENSG00000213934', 'ENSG00000196565', 'ENSG00000206177', 'ENSG00000086506', 'ENSG00000130656', 'ENSG00000206178') 

## identify count files
countfiles <- Sys.glob(file.path(output_prefix, '/03_generate_counts/*/quant.sf'))
samples <- sapply(basename(dirname(countfiles)), \(x) paste(strsplit(x,'_')[[1]][3:4], collapse='_'))
names(countfiles) <- names(samples)
sample_data_9s <- sample_data[sample_data$Animal.ID %in% sub('_.*', '', samples),]

## import orthologs
orthologs <- read.table(Sys.glob(file.path(output_prefix,'02_find_orthologs','of_out','*','Phylogenetic_Hierarchical_Orthogroups','N0.tsv')), header=TRUE, row.names = 1, sep='\t')

## import counts
## feels weird to round decimals for poisson error, but since data were counts at one point, zeros are possible, so doesn't make sense to just log-transform and use gaussian error; and error probably still scales like poisson such that low counts are less meaningful (and the decimal values rounded off would contribute negligible information) (cite Z1000 tximport paper)
cat('Importing counts\n')
txi_genes_ncbi_tmp <- list()
txi_orthogroups_ncbi_tmp <- list()
for(x in ncbi_species) {
  
  # find all count files that belong to the current species
  cfiles <- countfiles[grep(x, names(countfiles), ignore.case = TRUE)]
  
  # collect all transcript names
  tnames <- unique(do.call(rbind, lapply(cfiles, read.table, header=TRUE, sep='\t'))$Name)
  
  # create counts per gene
  tx2gene <- data.frame(transcript=tnames, gene=sapply(tnames, \(y) paste(strsplit(y,'.',fixed=TRUE)[[1]][1:2],collapse='.')))
  txi_genes_ncbi_tmp[[x]] <- tximport::tximport(cfiles, type = 'salmon', tx2gene = tx2gene)
  
  # create counts per orthogroup
  # create intermediate tx2gene mapping table by extracting the specific gene names, which are a subset of transcript names
  tmp <- data.frame(transcript   = tnames, 
                    gene.species = sapply(tnames, \(y) paste(strsplit(y, '.', fixed = TRUE)[[1]][1:2], collapse = '.')))
  
  
  tmp$gene <- unlist(parallel::mclapply(tmp$gene.species, \(y) {
    ## if the 'gene.species' value exists in 'orthologs' with a period appended to it, then orthofinder has placed it in an orthogroup
    tmpidx <- grep(paste0(y, '.'), orthologs[,x], fixed = TRUE)
    if(length(tmpidx) == 1) {
      return(rownames(orthologs)[[tmpidx]])
    } else {
      return(y)
    }
  }, mc.cores = nthreads))
  tx2og <- data.frame(transcript = tmp$transcript, gene=tmp$gene)

  txi_orthogroups_ncbi_tmp[[x]] <- tximport::tximport(cfiles, type = 'salmon', tx2gene = tx2og)
  
}

txi_genes_ensembl_tmp <- list()
txi_orthogroups_ensembl_tmp <- list()
numtries <- 5
for(x in ensembl_species) {
  
  cfiles <- countfiles[grep(x, names(countfiles), ignore.case = TRUE)]
  tnames <- unique(do.call(rbind, lapply(cfiles, read.table, header=TRUE, sep='\t'))$Name)
  
  for(i in 1:numtries) {
    try({
      ee_t <- biomaRt::getBM(attributes = c('ensembl_transcript_id_version','ensembl_gene_id'), 
                             filters    = 'ensembl_transcript_id_version', 
                             values     = tnames, 
                             mart       = biomaRt::useMart('ENSEMBL_MART_ENSEMBL', martlist[[x]], host='https://feb2023.archive.ensembl.org'), 
                             uniqueRows = TRUE)
      break
    }, silent = FALSE)
  }
  
  tx2gene <- ee_t[,c('ensembl_transcript_id_version', 'ensembl_gene_id')]
  txi_genes_ensembl_tmp[[x]] <- tximport::tximport(cfiles, type = 'salmon', tx2gene = tx2gene)
  
  for(i in 1:numtries) {
    try({
      ee_p <- biomaRt::getBM(attributes = c('ensembl_peptide_id_version', 'ensembl_gene_id'), 
                             filters = 'ensembl_peptide_id_version', 
                             values = unlist(strsplit(orthologs[,x],', ')), 
                             mart = biomaRt::useMart('ENSEMBL_MART_ENSEMBL', martlist[[x]], host='https://feb2023.archive.ensembl.org'), 
                             uniqueRows = TRUE)
    break
    }, silent = FALSE)
  }
  
  ee_p$og <- unlist(parallel::mclapply(ee_p$ensembl_peptide_id_version, \(y) {
    tmpidx <- grep(paste0(y,'\\b'), orthologs[,x])
    if(length(tmpidx) == 1) {
      return(rownames(orthologs)[[tmpidx]])
    } else {
      return(y)
    }
  }, mc.cores = nthreads))
  
  ee <- merge(ee_t,ee_p, by = 'ensembl_gene_id', all.x = TRUE)
  
  # 'gene' is the entire orthogroup, so tximport considers all gene duplicates as if they were different isoforms
  tx2og <- data.frame(transcript = ee$ensembl_transcript_id_version, gene = ee$og)
  tx2og$gene[is.na(tx2og$gene)] <- tx2og$transcript[is.na(tx2og$gene)]
  
  txi_orthogroups_ensembl_tmp[[x]] <- tximport::tximport(cfiles, type = 'salmon', tx2gene = tx2og)

}
##

cat('Retrieving ensembl data\n')
for(i in 1:numtries) {
  try({
    ensembl2ext <- biomaRt::getBM(attributes = c('ensembl_peptide_id_version', 'ensembl_gene_id', 'external_gene_name', 'go_id'), 
                                  filters = 'ensembl_peptide_id_version', 
                                  values = unlist(strsplit(orthologs[,'Homo_sapiens'],', ')), 
                                  mart = biomaRt::useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', host='https://feb2023.archive.ensembl.org'), 
                                  uniqueRows = TRUE)
    break
  }, silent = FALSE)
}

human_globin_peps <- unique(ensembl2ext$ensembl_peptide_id_version[ensembl2ext$ensembl_gene_id %in% human_globin_gene_ids])
human_globin_OGs <- rownames(orthologs)[sapply(rownames(orthologs), \(x) any(unlist(strsplit(orthologs[x,'Homo_sapiens'],', ')) %in% human_globin_peps))]

## import immune annotations
## all ensembl identifiers corresponding to immune annotations
modules <- read.table(file.path(raw_data_prefix, 'immune_modules.txt'), sep='\t', quote='', header=T)
ensembl2ext$module <- modules$IIG_class[match(ensembl2ext$external_gene_name, modules$HGNC_symbol)]
ensembl2ext$module[is.na(ensembl2ext$module)] <- 'no_module_annotation'

immuneGenes <- unique(ensembl2ext$ensembl_gene_id[ensembl2ext$module != 'no_module_annotation'])
immunePeps <- unique(ensembl2ext$ensembl_peptide_id_version[ensembl2ext$module != 'no_module_annotation'])
immune_OGs <- rownames(orthologs)[sapply(rownames(orthologs), \(x) any(unlist(strsplit(orthologs[x,'Homo_sapiens'],', ')) %in% immunePeps))]

## import reference data for species average body sizes
body_size_ref <- read.csv(file.path(raw_data_prefix, 'gyz043_suppl_Supplement_Data.csv'))
body_size_ref$genus_species <- paste(body_size_ref$genus, body_size_ref$species, sep='_')
body_size_ref$genus_species[body_size_ref$genus_species == 'Cebus_apella'] <- 'Sapajus_apella' 
##

####### begin 9-species analyses

## calculate differences of individual sizes from species means and log-transform
sample_data_9s$body_mass_log <- log(as.numeric(sample_data_9s$Body.Mass..g.))

body_mass_log_sp_means <- sapply(unique(sample_data_9s$genus_species), \(x) log(as.numeric(body_size_ref$Mean_body_mass_g[body_size_ref$genus_species == x])))

sample_data_9s$body_mass_log_sp <- body_mass_log_sp_means[sample_data_9s$genus_species]
sample_data_9s$body_mass_log_diff <- sample_data_9s$body_mass_log - sample_data_9s$body_mass_log_sp
##

## calculate ancestral body size (somewhat arbitrary for our purposes but could help interpret 'main effects')
body_mass_log_center <- ape::ace(body_mass_log_sp_means, species_tree_norm)$ace[[1]]
##

## standardize masses for model input
body_mass_log_sd <- sd(sample_data_9s$body_mass_log_sp)
body_mass_log_diff_sd <- sd(sample_data_9s$body_mass_log_diff)

sample_data_9s$body_mass_log_sp_std <- (sample_data_9s$body_mass_log_sp - body_mass_log_center) / body_mass_log_sd
sample_data_9s$body_mass_log_diff_std <- sample_data_9s$body_mass_log_diff / body_mass_log_diff_sd
##

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

## more linear combinations for plotting credible intervals
range_bmsp <- range(sample_data_9s$body_mass_log_sp_std)
seq_bmsp <- seq(range_bmsp[[1]], range_bmsp[[2]], length.out=50)

ci1 <- INLA::inla.make.lincombs('body_mass_log_sp_std' = seq_bmsp,
                                '(Intercept)'          = rep(1,50))
names(ci1) <- paste0('evo_allometry_baseline_prediction_', 1:50)

ci2 <- INLA::inla.make.lincombs('body_mass_log_sp_std' = seq_bmsp,
                                '(Intercept)'          = rep(1,50),
                                'treatmentLPS'         = rep(1,50),
                                'body_mass_log_sp_std:treatmentLPS' = seq_bmsp)
names(ci2) <- paste0('evo_allometry_LPS_prediction_', 1:50)

ci3 <- INLA::inla.make.lincombs('treatmentLPS'         = rep(1,50),
                                'body_mass_log_sp_std:treatmentLPS' = seq_bmsp)
names(ci3) <- paste0('evo_allometry_LPS_response_prediction_', 1:50)

range_bmdev <- range(sample_data_9s$body_mass_log_diff_std)
seq_bmdev <- seq(range_bmdev[[1]], range_bmdev[[2]], length.out=50)

ci4 <- INLA::inla.make.lincombs('body_mass_log_diff_std' = seq_bmdev,
                                '(Intercept)'          = rep(1,50))
names(ci4) <- paste0('intra_allometry_baseline_prediction_', 1:50)

ci5 <- INLA::inla.make.lincombs('body_mass_log_diff_std' = seq_bmdev,
                                '(Intercept)'          = rep(1,50),
                                'treatmentLPS'         = rep(1,50),
                                'treatmentLPS:body_mass_log_diff_std' = seq_bmdev)
names(ci5) <- paste0('intra_allometry_LPS_prediction_', 1:50)

ci6 <- INLA::inla.make.lincombs('treatmentLPS'         = rep(1,50),
                                'treatmentLPS:body_mass_log_diff_std' = seq_bmdev)
names(ci6) <- paste0('intra_allometry_LPS_response_prediction_', 1:50)

all_lincombs <- c(lc1,lc2,lc3,lc4,lc5,lc6,ci1,ci2,ci3,ci4,ci5,ci6)


## identify 1:1 orthologs for all 9 species
# any row (orthogroup) that has a comma has duplications in at least one species. Any row with a column with no characters means the orthogroup doesn't exist in at least one species
orthologs_9s_11 <- orthologs[!apply(orthologs[,-c(1,2)], 1, \(x) any(grepl(',',x) | (nchar(x) == 0))),] 
# edit names to remove peptide ID suffixes for ncbi species (making them just gene IDs)
orthologs_9s_11[,colnames(orthologs_9s_11) %in% ncbi_species] <- apply(orthologs_9s_11[,colnames(orthologs_9s_11) %in% ncbi_species], 2, \(x) sapply(x, \(y) paste(strsplit(y,'.',fixed=TRUE)[[1]][1:2],collapse='.')))
# replace peptide ids with gene ids for ensembl species
orthologs_9s_11[,ensembl_species] <- sapply(ensembl_species, \(x) {
  for(i in 1:numtries) {
    try({
      ee <- biomaRt::getBM(attributes = c('ensembl_peptide_id_version', 'ensembl_gene_id'), 
                           filters = 'ensembl_peptide_id_version', 
                           values = orthologs_9s_11[,x], 
                           mart = biomaRt::useMart('ENSEMBL_MART_ENSEMBL', martlist[[x]], host='https://feb2023.archive.ensembl.org'), 
                           uniqueRows = TRUE)
      break
    }, silent = FALSE)
  }
  rownames(ee) <- ee[,'ensembl_peptide_id_version']
  return(ee[orthologs_9s_11[,x], 'ensembl_gene_id'])
})
# remove globins
orthologs_9s_11 <- orthologs_9s_11[!orthologs_9s_11[,'Homo_sapiens'] %in% human_globin_gene_ids,] 

## prepare count data for 1:1 ortholog analysis
txi_9s_genes_tmp <- c(txi_genes_ncbi_tmp, txi_genes_ensembl_tmp)
names(txi_9s_genes_tmp) <- c(ncbi_species, ensembl_species)

ortholog_9s_11_mat <- orthologs_9s_11[,names(txi_9s_genes_tmp)]
ortholog_9s_11_mat <- ortholog_9s_11_mat[ortholog_9s_11_mat[,1] %in% rownames(txi_9s_genes_tmp[[1]]$abundance), ]
txi_9s_11_genes <- list(abundance = txi_9s_genes_tmp[[1]]$abundance[ortholog_9s_11_mat[,1],], 
                        counts    = txi_9s_genes_tmp[[1]]$counts[ortholog_9s_11_mat[,1],],
                        length    = txi_9s_genes_tmp[[1]]$length[ortholog_9s_11_mat[,1],],
                        countsFromAbundance = txi_9s_genes_tmp[[1]]$countsFromAbundance)
for(i in 2:length(txi_9s_genes_tmp)) {
  ortholog_9s_11_mat <- ortholog_9s_11_mat[ortholog_9s_11_mat[,i] %in% rownames(txi_9s_genes_tmp[[i]]$abundance), ]

  txi_9s_11_genes$abundance <- cbind(txi_9s_11_genes$abundance[ortholog_9s_11_mat[,1],], txi_9s_genes_tmp[[i]]$abundance[ortholog_9s_11_mat[,i],])
  txi_9s_11_genes$counts <- cbind(txi_9s_11_genes$counts[ortholog_9s_11_mat[,1],], txi_9s_genes_tmp[[i]]$counts[ortholog_9s_11_mat[,i],])
  txi_9s_11_genes$length <- cbind(txi_9s_11_genes$length[ortholog_9s_11_mat[,1],], txi_9s_genes_tmp[[i]]$length[ortholog_9s_11_mat[,i],])
} 
rownames(txi_9s_11_genes$abundance) <- rownames(txi_9s_11_genes$counts) <- rownames(txi_9s_11_genes$length) <- ortholog_9s_11_mat[match(rownames(txi_9s_11_genes$length), ortholog_9s_11_mat[,1]), 'Homo_sapiens']

## get per-gene-and-sample normalization factors for 1:1 orthologs, and per-sample size factors based on orthologs
des_9s_11_genes <- DESeq2:::DESeqDataSetFromTximport(txi_9s_11_genes, data.frame(Int=rep(1,ncol(txi_9s_11_genes$counts))), ~1)
des_9s_11_genes <- DESeq2::estimateSizeFactors(des_9s_11_genes)

sfs_9s_11_genes <- apply(SummarizedExperiment::assay(des_9s_11_genes, 'normalizationFactors'), 2, \(x) exp(mean(log(x))))
##

## pool immune-annotated 1:1 orthologous genes
geneGroups_9s_11_genes_immune <- rownames(txi_9s_11_genes$abundance)
geneGroups_9s_11_genes_immune[geneGroups_9s_11_genes_immune %in% immuneGenes] <- 'immune'
txi_9s_11_genes_immunePool <- txi_9s_11_genes
txi_9s_11_genes_immunePool$abundance <- rowsum(txi_9s_11_genes$abundance, geneGroups_9s_11_genes_immune)
txi_9s_11_genes_immunePool$counts <- rowsum(txi_9s_11_genes$counts, geneGroups_9s_11_genes_immune)
txi_9s_11_genes_immunePool$length <- rowsum(txi_9s_11_genes$abundance * txi_9s_11_genes$length, geneGroups_9s_11_genes_immune) / txi_9s_11_genes_immunePool$abundance ## weighted arithmetic mean length
txi_9s_11_genes_immunePool$length[rownames(txi_9s_11_genes_immunePool$length) != 'immune',] <- txi_9s_11_genes$length[rownames(txi_9s_11_genes_immunePool$length)[rownames(txi_9s_11_genes_immunePool$length) != 'immune'],]

## get normalization factors
des_9s_11_genes_immunePool <- DESeq2:::DESeqDataSetFromTximport(txi_9s_11_genes_immunePool, data.frame(Int=rep(1,ncol(txi_9s_11_genes_immunePool$counts))), ~1)
des_9s_11_genes_immunePool <- DESeq2::estimateSizeFactors(des_9s_11_genes_immunePool)

## model total differential expression of immune genes
dat_9s_11_genes_immunePool <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des_9s_11_genes_immunePool)), '_'), \(x) x[[3]]), 
                                         species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des_9s_11_genes_immunePool)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                                         treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des_9s_11_genes_immunePool)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                                         norm_factor = log(DESeq2::normalizationFactors(des_9s_11_genes_immunePool)['immune',]),
                                         count       = DESeq2::counts(des_9s_11_genes_immunePool)['immune',])
dat_9s_11_genes_immunePool <- dat_9s_11_genes_immunePool[dat_9s_11_genes_immunePool$individual %in% sample_data_9s$Animal.ID,]
dat_9s_11_genes_immunePool$body_mass_log_sp_std <- sapply(dat_9s_11_genes_immunePool$individual, function(z) sample_data_9s$body_mass_log_sp_std[sample_data_9s$Animal.ID==z])
dat_9s_11_genes_immunePool$body_mass_log_diff_std <- sapply(dat_9s_11_genes_immunePool$individual, function(z) sample_data_9s$body_mass_log_diff_std[sample_data_9s$Animal.ID==z])

cat('Modeling immune orthologs\n')
res_9s_11_genes_immunePool <- phyr::pglmm(formula       = phy_formula, 
                                          family        = "poisson", 
                                          cov_ranef     = list(species = species_tree_norm), 
                                          data          = dat_9s_11_genes_immunePool, 
                                          bayes         = TRUE, 
                                          bayes_options = list(lincomb = all_lincombs, 
                                                               control.compute = list(config = TRUE)))
restab_9s_11_genes_immunePool <- res_9s_11_genes_immunePool$inla.model$summary.lincomb.derived
##

## do same but for non-immune, to show it's not a general artifact of the summarization
geneGroups_9s_11_genes_nonimmune <- rownames(txi_9s_11_genes$abundance)
geneGroups_9s_11_genes_nonimmune[!geneGroups_9s_11_genes_nonimmune %in% immuneGenes] <- 'nonimmune'
txi_9s_11_genes_nonimmunePool <- txi_9s_11_genes
txi_9s_11_genes_nonimmunePool$abundance <- rowsum(txi_9s_11_genes$abundance, geneGroups_9s_11_genes_nonimmune)
txi_9s_11_genes_nonimmunePool$counts <- rowsum(txi_9s_11_genes$counts, geneGroups_9s_11_genes_nonimmune)
txi_9s_11_genes_nonimmunePool$length <- rowsum(txi_9s_11_genes$abundance * txi_9s_11_genes$length, geneGroups_9s_11_genes_nonimmune) / txi_9s_11_genes_nonimmunePool$abundance ## weighted arithmetic mean length
txi_9s_11_genes_nonimmunePool$length[rownames(txi_9s_11_genes_nonimmunePool$length) != 'nonimmune',] <- txi_9s_11_genes$length[rownames(txi_9s_11_genes_nonimmunePool$length)[rownames(txi_9s_11_genes_nonimmunePool$length) != 'nonimmune'],]

## get normalization factors
des_9s_11_genes_nonimmunePool <- DESeq2:::DESeqDataSetFromTximport(txi_9s_11_genes_nonimmunePool, data.frame(Int=rep(1,ncol(txi_9s_11_genes_nonimmunePool$counts))), ~1)
des_9s_11_genes_nonimmunePool <- DESeq2::estimateSizeFactors(des_9s_11_genes_nonimmunePool)

## model total differential expression of non-immune genes
dat_9s_11_genes_nonimmunePool <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des_9s_11_genes_nonimmunePool)), '_'), \(x) x[[3]]), 
                                            species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des_9s_11_genes_nonimmunePool)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                                            treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des_9s_11_genes_nonimmunePool)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                                            norm_factor = log(DESeq2::normalizationFactors(des_9s_11_genes_nonimmunePool)['nonimmune',]),
                                            count       = DESeq2::counts(des_9s_11_genes_nonimmunePool)['nonimmune',])
dat_9s_11_genes_nonimmunePool <- dat_9s_11_genes_nonimmunePool[dat_9s_11_genes_nonimmunePool$individual %in% sample_data_9s$Animal.ID,]
dat_9s_11_genes_nonimmunePool$body_mass_log_sp_std <- sapply(dat_9s_11_genes_nonimmunePool$individual, function(z) sample_data_9s$body_mass_log_sp_std[sample_data_9s$Animal.ID==z])
dat_9s_11_genes_nonimmunePool$body_mass_log_diff_std <- sapply(dat_9s_11_genes_nonimmunePool$individual, function(z) sample_data_9s$body_mass_log_diff_std[sample_data_9s$Animal.ID==z])

res_9s_11_genes_nonimmunePool <- phyr::pglmm(formula       = phy_formula, 
                                             family        = "poisson", 
                                             cov_ranef     = list(species = species_tree_norm), 
                                             data          = dat_9s_11_genes_nonimmunePool, 
                                             bayes         = TRUE, 
                                             bayes_options = list(lincomb = all_lincombs, 
                                                                  control.compute = list(config = TRUE)))
restab_9s_11_genes_nonimmunePool <- res_9s_11_genes_nonimmunePool$inla.model$summary.lincomb.derived
##

## prepare count data for orthogroup analysis
txi_9s_orthogroups_tmp <- c(txi_orthogroups_ncbi_tmp, txi_orthogroups_ensembl_tmp)
txi_9s_orthogroups <- txi_9s_orthogroups_tmp[[1]]
for(s in 2:length(txi_9s_orthogroups_tmp)) {
  
  txi_9s_orthogroups[[1]] <- merge(txi_9s_orthogroups[[1]], txi_9s_orthogroups_tmp[[s]][[1]], by='row.names')
  rownames(txi_9s_orthogroups[[1]]) <- txi_9s_orthogroups[[1]][,1]
  txi_9s_orthogroups[[1]] <- as.matrix(txi_9s_orthogroups[[1]][,-1])
  
  txi_9s_orthogroups[[2]] <- merge(txi_9s_orthogroups[[2]], txi_9s_orthogroups_tmp[[s]][[2]], by='row.names')
  rownames(txi_9s_orthogroups[[2]]) <- txi_9s_orthogroups[[2]][,1]
  txi_9s_orthogroups[[2]] <- as.matrix(txi_9s_orthogroups[[2]][,-1])
  
  txi_9s_orthogroups[[3]] <- merge(txi_9s_orthogroups[[3]], txi_9s_orthogroups_tmp[[s]][[3]], by='row.names')
  rownames(txi_9s_orthogroups[[3]]) <- txi_9s_orthogroups[[3]][,1]
  txi_9s_orthogroups[[3]] <- as.matrix(txi_9s_orthogroups[[3]][,-1])
  
}

## remove globins
txi_9s_orthogroups$abundance <- txi_9s_orthogroups$abundance[!rownames(txi_9s_orthogroups$abundance) %in% human_globin_OGs,]
txi_9s_orthogroups$counts <- txi_9s_orthogroups$counts[!rownames(txi_9s_orthogroups$counts) %in% human_globin_OGs,]
txi_9s_orthogroups$length <- txi_9s_orthogroups$length[!rownames(txi_9s_orthogroups$length) %in% human_globin_OGs,]

des_9s_orthogroups <- DESeq2:::DESeqDataSetFromTximport(txi_9s_orthogroups, data.frame(Int=rep(1,ncol(txi_9s_orthogroups$counts))), ~1)

## create txi_9s_orthogroups_immunePool analogous to above

## pool immune-annotated orthogroups
geneGroups_9s_11_orthogroups_immune <- rownames(txi_9s_orthogroups$abundance)
geneGroups_9s_11_orthogroups_immune[geneGroups_9s_11_orthogroups_immune %in% immune_OGs] <- 'immune' ## if any member of an orthogroup is immune, the whole thing is
txi_9s_orthogroups_immunePool <- txi_9s_orthogroups
txi_9s_orthogroups_immunePool$abundance <- rowsum(txi_9s_orthogroups$abundance, geneGroups_9s_11_orthogroups_immune)
txi_9s_orthogroups_immunePool$counts <- rowsum(txi_9s_orthogroups$counts, geneGroups_9s_11_orthogroups_immune)
txi_9s_orthogroups_immunePool$length <- rowsum(txi_9s_orthogroups$abundance * txi_9s_orthogroups$length, geneGroups_9s_11_orthogroups_immune) / txi_9s_orthogroups_immunePool$abundance ## weighted arithmetic mean length
txi_9s_orthogroups_immunePool$length[rownames(txi_9s_orthogroups_immunePool$length) != 'immune',] <- txi_9s_orthogroups$length[rownames(txi_9s_orthogroups_immunePool$length)[rownames(txi_9s_orthogroups_immunePool$length) != 'immune'],]
##

des_9s_orthogroups_immunePool <- DESeq2:::DESeqDataSetFromTximport(txi_9s_orthogroups_immunePool, data.frame(Int=rep(1,ncol(txi_9s_orthogroups_immunePool$counts))), ~1)

## model total differential expression of immune genes
dat_9s_orthogroups_immunePool <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des_9s_orthogroups_immunePool)), '_'), \(x) x[[3]]), 
                                            species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des_9s_orthogroups_immunePool)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                                            treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des_9s_orthogroups_immunePool)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                                            norm_factor = log(sfs_9s_11_genes[colnames(DESeq2::counts(des_9s_orthogroups_immunePool))]),
                                            count       = DESeq2::counts(des_9s_orthogroups_immunePool)['immune',])
dat_9s_orthogroups_immunePool <- dat_9s_orthogroups_immunePool[dat_9s_orthogroups_immunePool$individual %in% sample_data_9s$Animal.ID,]
dat_9s_orthogroups_immunePool$body_mass_log_sp_std <- sapply(dat_9s_orthogroups_immunePool$individual, function(z) sample_data_9s$body_mass_log_sp_std[sample_data_9s$Animal.ID==z])
dat_9s_orthogroups_immunePool$body_mass_log_diff_std <- sapply(dat_9s_orthogroups_immunePool$individual, function(z) sample_data_9s$body_mass_log_diff_std[sample_data_9s$Animal.ID==z])

cat('Modeling immune orthogroups\n')
res_9s_orthogroups_immunePool <- phyr::pglmm(formula       = phy_formula, 
                                             family        = "poisson", 
                                             cov_ranef     = list(species = species_tree_norm), 
                                             data          = dat_9s_orthogroups_immunePool, 
                                             bayes         = TRUE, 
                                             bayes_options = list(lincomb = all_lincombs, 
                                                                  control.compute = list(config = TRUE)))
restab_9s_orthogroups_immunePool <- res_9s_orthogroups_immunePool$inla.model$summary.lincomb.derived
##

## do same but for non-immune, to show it's not a general artifact of the summarization
geneGroups_9s_11_orthogroups_nonimmune <- rownames(txi_9s_orthogroups$abundance)
geneGroups_9s_11_orthogroups_nonimmune[!geneGroups_9s_11_orthogroups_nonimmune %in% immune_OGs] <- 'nonimmune' 
txi_9s_orthogroups_nonimmunePool <- txi_9s_orthogroups
txi_9s_orthogroups_nonimmunePool$abundance <- rowsum(txi_9s_orthogroups$abundance, geneGroups_9s_11_orthogroups_nonimmune)
txi_9s_orthogroups_nonimmunePool$counts <- rowsum(txi_9s_orthogroups$counts, geneGroups_9s_11_orthogroups_nonimmune)
txi_9s_orthogroups_nonimmunePool$length <- rowsum(txi_9s_orthogroups$abundance * txi_9s_orthogroups$length, geneGroups_9s_11_orthogroups_nonimmune) / txi_9s_orthogroups_nonimmunePool$abundance ## weighted arithmetic mean length
txi_9s_orthogroups_nonimmunePool$length[rownames(txi_9s_orthogroups_nonimmunePool$length) != 'nonimmune',] <- txi_9s_orthogroups$length[rownames(txi_9s_orthogroups_nonimmunePool$length)[rownames(txi_9s_orthogroups_nonimmunePool$length) != 'nonimmune'],]
##

des_9s_orthogroups_nonimmunePool <- DESeq2:::DESeqDataSetFromTximport(txi_9s_orthogroups_nonimmunePool, data.frame(Int=rep(1,ncol(txi_9s_orthogroups_nonimmunePool$counts))), ~1)

## model total differential expression of immune genes
dat_9s_orthogroups_nonimmunePool <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des_9s_orthogroups_nonimmunePool)), '_'), \(x) x[[3]]), 
                                               species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des_9s_orthogroups_nonimmunePool)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                                               treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des_9s_orthogroups_nonimmunePool)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                                               norm_factor = log(sfs_9s_11_genes[colnames(DESeq2::counts(des_9s_orthogroups_nonimmunePool))]),
                                               count       = DESeq2::counts(des_9s_orthogroups_nonimmunePool)['nonimmune',])
dat_9s_orthogroups_nonimmunePool <- dat_9s_orthogroups_nonimmunePool[dat_9s_orthogroups_nonimmunePool$individual %in% sample_data_9s$Animal.ID,]
dat_9s_orthogroups_nonimmunePool$body_mass_log_sp_std <- sapply(dat_9s_orthogroups_nonimmunePool$individual, function(z) sample_data_9s$body_mass_log_sp_std[sample_data_9s$Animal.ID==z])
dat_9s_orthogroups_nonimmunePool$body_mass_log_diff_std <- sapply(dat_9s_orthogroups_nonimmunePool$individual, function(z) sample_data_9s$body_mass_log_diff_std[sample_data_9s$Animal.ID==z])

cat('Modeling immune orthogroups\n')
res_9s_orthogroups_nonimmunePool <- phyr::pglmm(formula       = phy_formula, 
                                                family        = "poisson", 
                                                cov_ranef     = list(species = species_tree_norm), 
                                                data          = dat_9s_orthogroups_nonimmunePool, 
                                                bayes         = TRUE, 
                                                bayes_options = list(lincomb = all_lincombs, 
                                                                     control.compute = list(config = TRUE)))
restab_9s_orthogroups_nonimmunePool <- res_9s_orthogroups_nonimmunePool$inla.model$summary.lincomb.derived
##

save.image(file.path(output_prefix, '04_phyr_results_20231112_9s_pre-genes.RData'))

## fit per-gene models
numtriesfit <- 15

cat('Fitting 1:1 ortholog gene models\n')
fits_9s_11_genes <- parallel::mclapply(rownames(des_9s_11_genes), \(x) {
  cat(round(which(x == rownames(des_9s_11_genes)) / nrow(des_9s_11_genes),2) * 100, '%: ', x,'\n')
  dat <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des_9s_11_genes)), '_'), \(x) x[[3]]), 
                    species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des_9s_11_genes)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                    treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des_9s_11_genes)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                    norm_factor = log(DESeq2::normalizationFactors(des_9s_11_genes)[x,]),
                    count       = DESeq2::counts(des_9s_11_genes)[x,])
  dat <- dat[dat$individual %in% sample_data_9s$Animal.ID,]
  dat$body_mass_log_sp_std <- sapply(dat$individual, \(z) sample_data_9s$body_mass_log_sp_std[sample_data_9s$Animal.ID==z])
  dat$body_mass_log_diff_std <- sapply(dat$individual, \(z) sample_data_9s$body_mass_log_diff_std[sample_data_9s$Animal.ID==z])
  
  if(sd(DESeq2::counts(des_9s_11_genes)[x,]) > 0) {
    
    fit <- NA
    for(i in 1:numtriesfit) {
      try({
        
        fit <- phyr::pglmm(formula       = phy_formula, 
                           family        = "poisson", 
                           cov_ranef     = list(species = species_tree_norm), 
                           data          = dat, 
                           bayes         = TRUE, 
                           bayes_options = list(lincomb = all_lincombs))
        
        if(!'(Intercept)' %in% rownames(fit$inla.model$summary.lincomb.derived)) {
          stop(paste0('Gene ', rownames(des_9s_11_genes)[[x]], ' has incomplete summary on try ', i))
        }
        
        break
        
      }, silent = FALSE)
    }
    
    if(!all(is.na(fit))) {
      restab <- fit$inla.model$summary.lincomb.derived
      fit$inla.model <- NULL ## this is simply too large to keep for all genes
    } else {
      restab <- NA
    }
    
  } else {
    
    fit <- NA
    restab <- NA
    
  }
  
  return(list(data=dat, fit=fit, restab=restab))
  
}, mc.cores = nthreads)
names(fits_9s_11_genes) <- rownames(des_9s_11_genes)

cat('Fitting all orthogroup models\n')
## almost identical to 1:1 procedure, except using size factors instead of normalization factors due to potential issues with messy genome assemblies
fits_9s_orthogroups <- parallel::mclapply(rownames(des_9s_orthogroups), \(x) {
  cat(round(which(x == rownames(des_9s_orthogroups)) / nrow(des_9s_orthogroups),2) * 100, '%: ', x,'\n')
  dat <- data.frame(individual  = sapply(strsplit(colnames(DESeq2::counts(des_9s_orthogroups)), '_'), \(x) x[[3]]), 
                    species     = as.factor(sapply(strsplit(colnames(DESeq2::counts(des_9s_orthogroups)), '_'), \(x) paste(tools::toTitleCase(x[1]), x[2], sep='_'))),
                    treatment   = factor(sapply(strsplit(colnames(DESeq2::counts(des_9s_orthogroups)), '_'), \(x) x[[4]]), levels=c('Null','LPS')),
                    norm_factor = log(sfs_9s_11_genes[colnames(DESeq2::counts(des_9s_orthogroups_immunePool))]),
                    count       = DESeq2::counts(des_9s_orthogroups)[x,])
  dat <- dat[dat$individual %in% sample_data_9s$Animal.ID,]
  dat$body_mass_log_sp_std <- sapply(dat$individual, \(z) sample_data_9s$body_mass_log_sp_std[sample_data_9s$Animal.ID==z])
  dat$body_mass_log_diff_std <- sapply(dat$individual, \(z) sample_data_9s$body_mass_log_diff_std[sample_data_9s$Animal.ID==z])
  
  if(sd(DESeq2::counts(des_9s_orthogroups)[x,]) > 0) {
    
    fit <- NA
    for(i in 1:numtriesfit) {
      try({
        
        fit <- phyr::pglmm(formula       = phy_formula, 
                           family        = "poisson", 
                           cov_ranef     = list(species = species_tree_norm), 
                           data          = dat, 
                           bayes         = TRUE, 
                           bayes_options = list(lincomb = all_lincombs))
        
        if(!'(Intercept)' %in% rownames(fit$inla.model$summary.lincomb.derived)) {
          stop(paste0('Gene ', rownames(des_9s_orthogroups)[[x]], ' has incomplete summary on try ', i))
        }
        
        break
        
      }, silent = FALSE)
    }
    
    if(!all(is.na(fit))) {
      restab <- fit$inla.model$summary.lincomb.derived
      fit$inla.model <- NULL ## this is simply too large to keep for all genes
    } else {
      restab <- NA
    }
    
  } else {
    
    fit <- NA
    restab <- NA
    
  }
  
  return(list(data=dat, fit=fit, restab=restab))
  
}, mc.cores = nthreads)
names(fits_9s_orthogroups) <- rownames(des_9s_orthogroups)


save.image(file.path(output_prefix, '04_phyr_results_20231112_9s.RData'))

