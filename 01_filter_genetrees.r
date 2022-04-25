thresh <- 4 #number of our species that the tree needs to contain in order to be kept
nthreads <- 20

samples <- list.files(path='inputs/reference_counts_cleaned', full.names=T, recursive=F)
names(samples) <- sapply(samples, function(x) {

  bn <- basename(x)
  paste(strsplit(bn, '_')[[1]][c(1,3)], collapse='_')

})

sample_data <- read.table('inputs/sample_data.txt', sep='\t', header=T)
## summarize pooled Microcebus murinus samples into one pseudo-sample
sample_data <- rbind(sample_data,apply(sample_data[sample_data$Genus=='Microcebus',], 2, function(x) {temp <- unique(x); if(length(temp)==1) return(temp) else return(mean(as.numeric(temp)))}))
sample_data$Animal.ID[nrow(sample_data)] <- 'Mouse'
##

sample_data$genus_species <- apply(sample_data[,c('Genus','Species')], 1, paste, collapse='_')
indivs_per_species <- lapply(unique(sample_data$genus_species), function(x) sample_data$Animal.ID[sample_data$genus_species == x])
names(indivs_per_species) <- unique(sample_data$genus_species)

## manually removing species that we know were mapped to inadequate genome assemblies
indivs_per_species <- indivs_per_species[!names(indivs_per_species) %in% c('Daubentonia_madagascariensis','Lemur_catta','Sapajus_apella')]
##

raw <- list()
for(species in names(indivs_per_species)) {
  for(indiv in indivs_per_species[[species]]) {
    for(nl in c('LPS','Null')) {
      nl2 <- paste(indiv, nl, sep='_')
      if(nl2 %in% names(samples)) {
        temp <- read.table(Sys.glob(file.path('inputs/reference_counts_cleaned', paste0(indiv, '*', nl, '*_counts.txt'))), header=T, sep='\t', row.names=1)[, 6, drop=F]
        colnames(temp) <- nl2
        if(species %in% names(raw)) {
          raw[[species]] <- merge(raw[[species]], temp, by=0)
          rownames(raw[[species]]) <- raw[[species]][,1]
          raw[[species]] <- raw[[species]][,-1]
        }
        else {
          raw[[species]] <- temp
        }
      }
    }
  }
}

allgenes <- lapply(raw, function(x) rownames(x))

genetrees <- list.files(path='outputs/primates/organized_ensembl_trees', pattern='*.xml$', full.names=T)
names(genetrees) <- sapply(genetrees, function(x) {

  bn <- basename(x)
  paste(strsplit(bn, '_')[[1]][[1]], strsplit(bn,'.', fixed=T)[[1]][[3]], sep='_')

})

dir.create('outputs/primates/filtered_ensembl_trees', showWarnings = F)

treenames <- parallel::mclapply(names(genetrees), function(tree) {
  hi <- rphyloxml::read_phyloxml(genetrees[[tree]], options='HUGE')
  thisxml <- xml2::read_xml(genetrees[[tree]], options='HUGE')
  hi[[1]]$edges[,'branch_length'] <- xml2::xml_double(xml2::xml_find_all(thisxml, "//@branch_length"))
  newtree <- rphyloxml::phyloxml2phylo(hi)[[1]]

  matchtips <- lapply(allgenes, function(x) intersect(x, newtree$tip.label))
  present_in <- sum(sapply(matchtips, length) > 0)
  cat(round(which(tree == names(genetrees)) / length(genetrees), 2) * 100, '%: ', tree, ' present in ', present_in, '\n', sep='')
  keep <- present_in >= thresh
  if(keep) {
    newtree <- ape::keep.tip(newtree,unlist(matchtips))
    ape::write.tree(newtree,paste0('outputs/primates/filtered_ensembl_trees/',tree,'.phy'))
  }

  treename <- tryCatch(xml2::as_list(thisxml)$phyloxml$phylogeny$property[[1]], error = function(e) NA)

  return(treename)
}, mc.cores=nthreads)
names(treenames) <- names(genetrees)

save.image('outputs/primates/filter_genetrees.RData')
