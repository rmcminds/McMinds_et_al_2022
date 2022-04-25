
load('outputs/primates/filter_genetrees.RData')

output_prefix <- path.expand('outputs/primates/')

genetree_paths <- list.files(path='outputs/primates/filtered_ensembl_trees', pattern='*.phy$', full.names=T)
names(genetree_paths) <- sapply(genetree_paths, function(x) paste(strsplit(basename(x),'.', fixed=T)[[1]][[1]], sep='_') )

genetrees <- lapply(genetree_paths, function(x) ape::read.tree(x))

genefamilies <- lapply(genetrees, function(x) x$tip.label)

martlist <- list(ENSCJAG='cjacchus_gene_ensembl', ENSG='hsapiens_gene_ensembl', ENSMICG='mmurinus_gene_ensembl', ENSMMUG='mmulatta_gene_ensembl', ENSPANG='panubis_gene_ensembl', ENSPPYG='pabelii_gene_ensembl')
marts <- lapply(martlist, function(y) biomaRt::useMart('ensembl', y))

genefamilies_ul <- unlist(genefamilies)

allspecs <- gsub('[0-9]','',genefamilies_ul)
gene_data <- lapply(unique(allspecs), function(y) biomaRt::getBM(attributes=c('ensembl_gene_id','external_gene_name','go_id'), filters='ensembl_gene_id', values=genefamilies_ul[allspecs==y], mart=marts[[y]], uniqueRows=TRUE))

save.image(file.path(output_prefix,'ensembl2GO.Rdata'))
