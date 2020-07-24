
################## Study the enrichment of antigen-specific clusters to responders

filterClusters      <- function(df){

  df_raw <- df %>% pull(index) %>% unique() %>% length()

  # df <- df %>% filter(number_subject > 2 & vb_score < 0.05 & Fisher_score < 0.05 & number_unique_cdr3 > 4)
  df <- df %>% filter(number_subject > 2 & number_unique_cdr3 > 9)

  df_new<- df %>% pull(index) %>% unique() %>% length()

  message(paste("Removed", df_raw - df_new, "of", df_raw, paste0("(", round(c(df_raw - df_new) / df_raw, 3), ")"), "clusters"))

  return(df)

}

getPatientInfoGliph <- function(df){

  ## Get patient info from Gliph-files provided by Liang

  extractGliphName <- function(str1){

    sub("\\:.*", "", str1)

  }


  df <- df %>% mutate(name = extractGliphName(df$Sample))
  meta_clinical_temp <- meta_clinical %>% dplyr::rename(name = `Study ID`)
  df <- df %>% left_join(meta_clinical_temp, by = "name")

  return(df)

}

gliphEnrichment     <- function(gliph_df, alternative = "greater"){

  getEnrichment <- function(df, cluster, alternative = "greater", r_clonotypes_total = r_clonotypes_total, n_clonotypes_total = n_clonotypes_total){

    message(cluster)

    df_temp              <- df %>% filter(index == cluster)

    n_clonotypes_cluster <- df_temp %>% filter(overall == "N") %>% pull(n)
    r_clonotypes_cluster <- df_temp %>% filter(overall == "R") %>% pull(n)

    if(length(n_clonotypes_cluster) == 0){n_clonotypes_cluster <- 0}
    if(length(r_clonotypes_cluster) == 0){r_clonotypes_cluster <- 0}

    fisher_df <- matrix(c(r_clonotypes_cluster, r_clonotypes_total - r_clonotypes_cluster, n_clonotypes_cluster, n_clonotypes_total - n_clonotypes_cluster), ncol = 2) %>%
      fisher.test(alternative = alternative) %>% broom::tidy() %>%
      mutate(r_clonotypes_cluster = r_clonotypes_cluster, n_clonotypes_cluster = n_clonotypes_cluster, r_clonotypes_total = r_clonotypes_total - r_clonotypes_cluster, n_clonotypes_total = n_clonotypes_total - n_clonotypes_cluster, index = cluster)

    return(fisher_df)

  }

  ## Do enrichment testing for each cluster. Which cluster contains more clonotypes from R vs N patients?
  universum          <- gliph_df %>% group_by(overall) %>% summarise(n = n())
  n_clonotypes_total <- universum %>% filter(overall == "N") %>% pull(n)
  r_clonotypes_total <- universum %>% filter(overall == "R") %>% pull(n)

  df        <- gliph_df %>% group_by(index, overall) %>% summarise(n = n()) %>% mutate(prop = n/sum(n))
  clusters  <- df$index %>% unique
  fisher_df <- lapply(clusters, getEnrichment, df = df, alternative = alternative, r_clonotypes_total = r_clonotypes_total, n_clonotypes_total = n_clonotypes_total) %>% rbind_list %>% mutate(p.adj = p.adjust(p.value, method = "BH"))

  return(fisher_df)

}

modifyFisherOutput  <- function(df){

  df %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% filter(p.value < 0.05)

}

extractGliphName    <- function(str1){
  strsplit(str1, "[:]")[[1]][1]
}

getGliphFreq        <- function(df, dataset_name){

  ## Look how the epitope-specific clusters change between the time points.
  ## For this information, we need to open the original files and fetch the frequency info
  ## of each clonotype in that time point

  ## input:
  # df = gliph output df, preferably filtered version

  df$clonotype   <- paste0(df$V, "_", df$TcRb)
  gliph_patients <- unique(df$Sample)

  gliph_list <- NULL
  i <- 1

  for(gliph_patient in gliph_patients){

    # gliph_patient = gliph_patients[1]
    message(gliph_patient)

    patient       <- extractGliphName(gliph_patient)
    gliph_pt_df   <- df %>% filter(Sample %in% gliph_patient)
    patient_files <- list.files(paste0("data/unselected_TCRseq/", dataset_name), pattern = patient, full.names = T)

    pre_file     <- grep("0m", patient_files, value = T)
    post_3m_file <- grep("3m", patient_files, value = T)

    post_1m_file <- grep("1m", patient_files, value = T)
    post_6m_file <- grep("6m", patient_files, value = T)

    ## Make sure to take the PB samples from Helsinki cohort
    if(dataset_name == "Helsinki"){
      pre_file         <- grep("PB", pre_file, value = T)
      post_3m_file     <- grep("PB", post_3m_file, value = T)
      post_1m_file     <- grep("PB", post_1m_file, value = T)
      post_6m_file     <- grep("PB", post_6m_file, value = T)
    }


    if(length(pre_file) > 0){

      pre_df <- fread(pre_file) %>% mutate(clonotype = paste0(v, "_", cdr3aa)) %>% dplyr::select(clonotype, freq, cdr3aa)
      pre_df <- aggregate(freq ~ clonotype, pre_df, sum)
      gliph_pt_df <- merge(gliph_pt_df, pre_df, by = "clonotype", all.x = T) %>% mutate(freq = ifelse(is.na(freq), 0, freq)) %>% dplyr::rename(pre_freq = "freq")


    }

    if(length(post_3m_file) > 0){

      post_df <- fread(post_3m_file) %>% mutate(clonotype = paste0(v, "_", cdr3aa)) %>% dplyr::select(clonotype, freq, cdr3aa)
      post_df <- aggregate(freq ~ clonotype, post_df, sum)
      gliph_pt_df <- merge(gliph_pt_df, post_df, by = "clonotype", all.x = T) %>% mutate(freq = ifelse(is.na(freq), 0, freq)) %>% dplyr::rename(post_3m_freq = "freq")

    }

    if(length(post_1m_file) > 0){

      post_df <- fread(post_file) %>% mutate(clonotype = paste0(v, "_", cdr3aa)) %>% mutate(freq = replace(freq,is.na(freq),0)) %>% dplyr::select(clonotype, freq)
      post_df <- aggregate(freq ~ clonotype, post_df, sum)
      gliph_pt_df <- merge(gliph_pt_df, post_df, by = "clonotype", all.x = T) %>% mutate(freq = ifelse(is.na(freq), 0, freq)) %>% dplyr::rename(post_1m_freq = "freq")

    }

    if(length(post_6m_file) > 0){

      post_df <- fread(post_file) %>% mutate(clonotype = paste0(v, "_", cdr3aa)) %>% mutate(freq = replace(freq,is.na(freq),0)) %>% dplyr::select(clonotype, freq)
      post_df <- aggregate(freq ~ clonotype, post_df, sum)
      gliph_pt_df <- merge(gliph_pt_df, post_df, by = "clonotype", all.x = T) %>% mutate(freq = ifelse(is.na(freq), 0, freq)) %>% dplyr::rename(post_6m_freq = "freq")

    }

    gliph_list[[i]] <- gliph_pt_df
    i <- i + 1

  }

  df_with_freqs <- gliph_list %>% rbindlist(fill = T) %>% arrange(desc(index))

  return(df_with_freqs)

}




######################### Study the differences in proliferation of antigen-specific clusters after therapy initiation


getFreqForGliph <- function(df, project_name){

  ## Look how the epitope-specific clusters change between the time points.
  ## For this information, we need to open the original files and fetch the frequency info
  ## of each clonotype in that time point

  ## input:
  # df = gliph output df, preferably filtered version

  ## In order to get the right freq, we need to open the data
  files   <- list.files(paste0("data/unselected_TCRseq/", project_name), full.names = T)
  df_list <- NULL
  i <- 1

  # name_temp = unique(df$name)[1]

  for(name_temp in unique(df$name)){

    message(name_temp)
    df_name    <- filter(df, name == name_temp) %>% mutate(clonotype = paste0(V, TcRb))

    files_name <- grep(name_temp, files, value = T)
    file_0m    <- grep("0m", files_name, value = T)
    file_3m    <- grep("3m", files_name, value = T)

    ## Make sure to take the PB samples from Helsinki cohort
    if(dataset_name == "Helsinki"){
      file_0m     <- grep("PB", file_0m, value = T)
      file_3m     <- grep("PB", file_3m, value = T)
    }


    ## Specific for this function, reduces repetition in the if-else statements
    openFile <- function(file, df_name){

      df_temp <- fread(file) %>% mutate(clonotype = paste0(v,cdr3aa)) %>% dplyr::select(clonotype, freq)
      df_temp <- aggregate(freq ~clonotype, df_temp, sum)

      df_name <- merge(df_name, df_temp, all.x = T, by = "clonotype") %>% mutate(freq = replace(freq, is.na(freq), 0))

    }

    if(length(file_0m) > 0){df_name <- openFile(file = file_0m, df_name) %>% dplyr::rename("freq_0m" = freq)}
    if(length(file_3m) > 0){df_name <- openFile(file = file_3m, df_name) %>% dplyr::rename("freq_3m" = freq)}

    df_list[[i]] <- df_name
    i <- i + 1

  }

  df_tot <- df_list %>% rbindlist(fill = TRUE) %>% group_by(pattern)
  return(df_tot)

}

getProliferationDEGliph <- function(df){

  # df = tumeh_freq

  ## Test if there's a difference between the proliferation between R and N
  # df <- data.frame(aggregate(freq_0m ~ name + overall + pattern, df, sum, drop = F), freq_3m = aggregate(freq_3m ~ name + overall + pattern, df, sum, drop = F)[,4])

  df_0m <- data.frame(aggregate(freq_0m ~ name + overall + pattern, df, sum, drop = T)) %>% mutate(temp = paste0(name, pattern))
  df_3m <- data.frame(aggregate(freq_3m ~ name + overall + pattern, df, sum, drop = T)) %>% mutate(temp = paste0(name, pattern))

  df <- merge(df_0m, select(df_3m, temp, freq_3m), by = "temp") %>% dplyr::select(-temp)
  df <- df %>% mutate(log2fc = log2(freq_3m / freq_0m))
  df <- df %>% dplyr::select(overall, log2fc, name, pattern) %>% melt %>% filter(is.finite(value))

  doTest <- function(df, pattern_temp){

    df_temp <- df %>% filter(pattern == pattern_temp)

    # df_temp <- data.frame(aggregate(freq_0m ~ name + overall, df_temp, sum), freq_3m = aggregate(freq_3m ~ name + overall, df_temp, sum)[,3])
    # df_temp <- df_temp %>% mutate(log2fc = log2(freq_3m / freq_0m))
    # df_temp <- df_temp %>% dplyr::select(overall, log2fc, name) %>% melt %>% filter(is.finite(value))

    if(length(unique(df_temp$overall)) > 1){

      median  <- df_temp %>% group_by(overall) %>% summarise(med = median(value)) %>% pull(med)
      test_df <- wilcox.test(value ~ overall, data = df_temp) %>% broom::glance() %>% mutate(pattern = pattern_temp, median_n = median[1], median_r = median[2])

      return(test_df)

    }

    # ggplot(df_temp, aes(overall,value)) + geom_boxplot() + geom_jitter() + ggsignif::geom_signif(comparisons = list(c("N", "R")))

  }

  test_df <- lapply(unique(df$pattern), function(x) doTest(df = df, pattern = x)) %>% rbindlist(fill=TRUE) %>% mutate(up = ifelse(median_n > median_r, "n", "r"))
  return(test_df)

}

plotProliferationGliph <- function(df, pattern_temp){


  ## Test if there's a difference between the proliferation between R and N
  df <- df %>% filter(pattern %in% pattern_temp)

  # df <- data.frame(aggregate(freq_0m ~ name + overall + pattern, df, sum), freq_3m = aggregate(freq_3m ~ name + overall + pattern, df, sum)[,4])
  df_0m <- data.frame(aggregate(freq_0m ~ name + overall + pattern, df, sum, drop = T)) %>% mutate(temp = paste0(name, pattern))
  df_3m <- data.frame(aggregate(freq_3m ~ name + overall + pattern, df, sum, drop = T)) %>% mutate(temp = paste0(name, pattern))

  df <- merge(df_0m, select(df_3m, temp, freq_3m), by = "temp") %>% dplyr::select(-temp)
  df <- df %>% mutate(log2fc = log2(freq_3m / freq_0m))
  df <- df %>% dplyr::select(overall, log2fc, name, pattern) %>% melt %>% filter(is.finite(value))

  ggplot(df, aes(overall,value, fill = overall)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(~pattern) +
    response_fill + geom_hline(yintercept = 0, linetype = "dotted") + facets_nice + theme(legend.position = "none") + labs(y = "log2fc")

}




######################### Basic gliph functions

readGliphFile <- function(filename){

  gliph_df   <- fread(filename, fill = T, stringsAsFactors = F)

  gliph_df$vb_score <- as.numeric(gliph_df$vb_score)
  gliph_df$Fisher_score <- as.numeric(gliph_df$Fisher_score)
  gliph_df$final_score <- as.numeric(gliph_df$final_score)
  gliph_df$number_unique_cdr3 <- as.numeric(gliph_df$number_unique_cdr3)

  gliph_df$number_subject <- as.numeric(gliph_df$number_subject)
  gliph_df$length_score <- as.numeric(gliph_df$length_score)
  gliph_df$cluster_size_score <- as.numeric(gliph_df$cluster_size_score)

  ## Remove everything after last numeric index
  max_index  <- max(as.numeric(gliph_df$index), na.rm = T)
  max_indexs <- which(gliph_df$index == max_index)
  gliph_df   <- gliph_df[1:max_indexs[length(max_indexs)], ]
  return(gliph_df)

}

plotGLIPH     <- function(filename, color = 'orange', filter_n = NA){

  require(igraph)

  ## Create edgelist from gliph output
  makeEL = function(df){

    if(nrow(df) > 1){
      el = df %>% as.data.frame() %>% pull(TcRb) %>% combn(m=2) %>% t %>% as.data.frame
      return(el)
    }

  }


  el = fread(filename, fill = T, stringsAsFactors = F) %>% dplyr::select(TcRb, pattern) %>% filter(pattern != '' | TcRb != '') %>% mutate(pattern = as.factor(pattern))
  if(!is.na(filter_n)){el = fread(filename, fill = T, stringsAsFactors = F) %>% mutate(pattern = as.factor(pattern)) %>% filter(index < filter_n) %>% filter(pattern != '' | TcRb != '')}


  el = split(el, f = el$pattern)
  el = lapply(el, makeEL) %>% rbindlist() %>% filter(as.character(V1) != as.character(V2))
  el = el[!duplicated(el), ]

  ##
  el2 <- el[sample(1:nrow(el), 1e3), ]

  ## From edgelist, create graph
  # g        = graph.data.frame(el, directed = F)
  g        = graph.data.frame(el2, directed = F)
  am       = get.adjacency(g,sparse=T)
  g_am     = graph_from_adjacency_matrix(am, mode = "undirected")
  g_am     = simplify(g_am, remove.multiple = F, remove.loops = T)

  ## Plot
  V(g_am)$size <- 3
  V(g_am)$frame.color <- "white"
  V(g_am)$color <- color
  V(g_am)$label <- ""
  E(g_am)$arrow.mode <- 0

  layout = layout_with_graphopt(g_am)
  plot.igraph(g_am, layout = layout)

}


######################### EPitope-specific (i.e. TCRGP-like) functions

plotFoldchange        <- function(foldchange_file){

  ## For epitope specific analyses
  limit = foldchange_file$foldchange[is.finite(foldchange_file$foldchange)] %>% abs %>% max %>% round(0)

  foldchange_file %>%
    ggplot(aes(overall,foldchange, fill = overall)) +
    geom_boxplot(outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = 3, color = 'black') +
    geom_jitter(size = 0.5) +
    facet_wrap(~species.x, ncol = 4) + facets_nice +
    ggsignif::geom_signif(comparisons = list(c('N', 'R'))) +
    response_fill + facets_nice + theme(legend.position = 'none') +
    labs(x = "", y = "log2(fc)") +
    ylim(-limit, limit)

}

plotStaticEpitopeSpec <- function(gliph_df){

  gliph_df = aggregate(freq ~ species + Sample + overall + io_stat, gliph_df, sum) %>% filter(species %in% interesting_species_gliph)

  gliph_df %>%
    ggplot(aes(overall,freq,fill=overall)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) +
    facet_wrap(~species, scales = "free", ncol = 4) + facets_nice +
    ggsignif::geom_signif(comparisons = list(c("N", "R"))) +
    response_fill +
    labs(x = "") + theme(legend.position = "none")

}

modifyForFoldchange   <- function(gliph_df, timepoint1, timepoint2){

  ## Analyze how much melanoma-specific clonotypes at therapy

  gliph_baseline  <- aggregate(freq ~ species + Sample + overall + io_stat, filter(gliph_df, timepoint == timepoint1), sum) %>% mutate(merge_id = paste0(Sample, species))
  gliph_followup  <- aggregate(freq ~ species + Sample, filter(gliph_df, timepoint == timepoint2), sum) %>% mutate(merge_id = paste0(Sample, species))

  gliph_tot       <- merge(gliph_baseline, gliph_followup, by = "merge_id") %>% mutate(foldchange = log2(freq.y / freq.x)) %>% filter(species.x %in% interesting_species_gliph)

  return(gliph_tot)

}

combineGliphResults   <- function(file1,file2){

  gliph_r <- fread(file1) %>% filter(target != "multi")
  gliph_n <- fread(file2) %>% filter(target != "multi")
  gliph   <- rbind(gliph_r, gliph_n)
  gliph$species[gliph$species == "HomoSapiens"] = "MELANOMA"

  return(gliph)

}
