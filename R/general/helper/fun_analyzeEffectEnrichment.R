
require(clusterProfiler)
require(org.Hs.eg.db)

if(me == "hru"){

  hallmark   <- read.gmt("/Users/hru/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt")
  tf         <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c3.tft.v7.0.symbols.gmt")
  go         <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c5.bp.v7.0.symbols.gmt")
  immunology <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c7.all.v7.0.symbols.gmt")


}


if(me == "janihuuh"){

  hallmark   <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt")
  tf         <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c3.tft.v7.0.symbols.gmt")
  go         <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c5.bp.v7.0.symbols.gmt")
  immunology <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c7.all.v7.0.symbols.gmt")

}




## == No need to modify

plotHypergeometric <- function(genes_df, universe_df, term_df){

  require(clusterProfiler)

  if(nrow(genes_df) == 0) return(NULL)

  # Enrichment in hypergeometric test
  # in: de_df with BM-annotations, direction (Up/Down) and universe to count the enrichment
  # from, e.g. all genes available (~60 000) or all genes expressed in the data (usually ~20 000)

  # out: df with enrichment results

  enrich <- enricher(genes_df$gene, universe = universe_df, TERM2GENE = term_df)

  if(table(enrich@result$p.adjust < 0.05) %>% length() > 1){
    heatplot(enrich)
  }

  else(NULL)

}

getHypergeometric <- function(genes_df, universe_df, term_df){

  require(clusterProfiler)

  if(nrow(genes_df) == 0) return(NULL)

  # Enrichment in hypergeometric test
  # in: de_df with BM-annotations, direction (Up/Down) and universe to count the enrichment
  # from, e.g. all genes available (~60 000) or all genes expressed in the data (usually ~20 000)

  # out: df with enrichment results

  enrich <- enricher(genes_df$gene, universe = universe_df, TERM2GENE = term_df)
  enrich <- do.call(rbind, enrich@result) %>% t %>% as.data.frame()
  enrich[,c(5:7, 9)] <- sapply(enrich[,c(5:7, 9)], function(x) {as.numeric(as.character(x))})
  return(enrich)

}
