
## sorted-files
pr1_tcrgp_raw <- lapply(list.files("tcrb_data/sorted/pr1_tolerable/", full.names = T), function(x) fread(x) %>% mutate(filename = extractFileName(x))) %>% rbindlist() %>% filter(count > 1) %>%
  group_by(filename) %>% mutate(freq = freq / sum(freq)) %>% mutate(subject = gsub(filename, pattern = ".txt", replacement = ""), key = paste(subject, cdr3aa))

pr1_tcrgp_raw %>% ggplot(aes(area = freq, fill = cdr3aa)) + treemapify::geom_treemap() + facet_wrap(~subject) + theme(legend.position = "none") + scale_fill_manual(values = getPalette3(772)) + facets_nice
ggsave("results/tcrgp/tree_pr1_raw.png", width = 8, height = 5.5)

## tcrgp file
pr1_tcrgp <- fread("results/tcrgp/input/pr1_tcrgp.txt") %>% mutate(key = paste(subject, cdr3b))
pr1_tcrgp %>% left_join(pr1_tcrgp_raw, by = "key")
pr1_tcrgp %>% ggplot(aes(area = freq, fill = cdr3aa)) + treemapify::geom_treemap() + facet_wrap(~subject) + theme(legend.position = "none") + scale_fill_manual(values = getPalette3(772))


## TCRGP input file
bfore_files    <- lapply(list.files("tcrb_data/unsorted/bfore/", full.names = T), FUN = function(x) fread(x) %>% mutate(filename = extractFileName(x))) %>% rbindlist()
dasa_files     <- lapply(list.files("tcrb_data/unsorted/dasastop/", full.names = T), FUN = function(x) fread(x) %>% mutate(filename = extractFileName(x))) %>% rbindlist()
cml_stop_files <- lapply(list.files("tcrb_data/unsorted/cml_sc/", full.names = T), FUN = function(x) fread(x) %>% mutate(filename = extractFileName(x))) %>% rbindlist()
cml_files      <- lapply(list.files("tcrb_data/unsorted/cml/", full.names = T), FUN = function(x) fread(x) %>% mutate(filename = extractFileName(x))) %>% rbindlist()
cml_tcrb_files <- rbind(bfore_files, dasa_files, cml_stop_files, cml_files)
fwrite(cml_tcrb_files, "/Users/janihuuh/Dropbox/cml_stop/results/tcrgp/cml_input.txt", sep = "\t", quote = F, row.names = F)
# data.frame(filename = list.files("tcrb_data/unsorted/dasastop/", full.names = T) %>% extractFileName() %>% gsub(pattern = ".txt", replacement = "")) %>% fwrite("data/clinical/dasastop_samples.txt", quote = F, row.names = F, sep = "\t")
cml_tcrb_files <- fread("/Users/janihuuh/Dropbox/cml_stop/results/tcrgp/cml_input.txt")



## Find antigen-specific TCRs, select 0.95 as a threshold
ag_specific           <- lapply(list.files("results/tcrgp/raw/tcrb/", full.names = T), function(x) fread(x)) %>% bind_cols()
colnames(ag_specific) <- do.call(lapply(list.files("results/tcrgp/raw/tcrb/", full.names = T), extractFileName) , what = "c") %>% gsub(pattern = ".csv", replacement =  "") %>% gsub(pattern = "cml_", replacement =  "")
ag_specific[ag_specific < 0.95] <- 0
ag_specific_df   <- cml_tcrb_files %>% bind_cols(ag_specific)
ag_specific_filt <- ag_specific_df # %>% filter(count > 1)

preds_df <- ag_specific_filt %>% dplyr::select(ELAGIGILTV_cdr3b_comb:YVLDHLIVV_cdr3b)
target <- apply(preds_df, 1, function(x) ifelse(max(x) == 0, 16, which.max(x)))
ag_specific_filt$target <- target
ag_specific_filt$target <- c(colnames(preds_df), "none")[ag_specific_filt$target]
ag_specific_filt$target[rowSums(preds_df) > 1] <- "multi"
ag_specific_filt$target[is.na(ag_specific_filt$target)] <- "none"

## Add meta info
extractTimepoint <- function(str1){
  strsplit(str1, "[_]")[[1]][2]
  # sub("\\_.*", "", str1)
}

ag_specific_filt$filename  <- gsub(pattern = ".txt", replacement = "", ag_specific_filt$filename)
ag_specific_filt$timepoint <- do.call(lapply(ag_specific_filt$filename, extractTimepoint), what = "c")
ag_specific_filt$type      <- ifelse(grepl("BM", ag_specific_filt$filename), "BM", "PB")
ag_specific_filt$name      <- extractName(ag_specific_filt$filename)

ag_specific_filt$timepoint[ag_specific_filt$timepoint == "d0"]       <- "0m"
ag_specific_filt$timepoint[ag_specific_filt$timepoint == "baseline"] <- "0m"
ag_specific_filt$timepoint[!ag_specific_filt$timepoint %in% c("0m", "3m", "6m", "12m", "relapse")] <- "none"
ag_specific_filt$timepoint <- factor(as.character(ag_specific_filt$timepoint), levels = c("0m", "3m", "6m", "12m", "relapse"))

## Add clin info
clinical_data              <- fread("data/clinical/tcrb_data_clinical.txt") %>% mutate(filename = gsub(filename, pattern = ".txt", replacement = ""))
ag_specific_filt           <- ag_specific_filt %>% left_join(clinical_data, by = "filename")
ag_specific_filt$timepoint <- ag_specific_filt$timepoint.y
ag_specific_filt$type      <- ag_specific_filt$type.y
ag_specific_filt <- ag_specific_filt %>% filter(!is.na(tki))
ag_specific_filt$timepoint <- factor(as.character(ag_specific_filt$timepoint), levels = c("dg", "6m", "12m", "baseline", "tf_1m", "tf_6m", "relapse", "2nd_baseline", "unknown"))

ag_specific_filt %>%
  group_by(filename,type,target) %>% summarise(tot = sum(freq), n = n()) %>%
  ggplot(aes(n,tot,fill=type)) + geom_point(shape=21, size = 3) +
  scale_x_log10() + scale_y_log10() +
  geom_smooth(method = "lm", color = "darkred", fill = "gray90") + ggpubr::stat_cor() + labs(x = "nClonotypes", y = "freq clonotypes") + scale_fill_manual(values = getPalette3(4)) + facet_wrap(~target, scales = "free") + facets_nice
ggsave("results/tcrgp/scatter_total.png", width = 9, height = 8)

ag_specific_filt %>% group_by(filename,type,target) %>% summarise(tot = sum(freq), n = n()) %>%
  filter(target != "none") %>%
  left_join(ag_specific_filt %>% group_by(filename,type,target) %>% summarise(tot = sum(freq), n = n()) %>% group_by(target) %>% summarise(median = median(tot)), by = "target") %>%
  ggplot(aes(reorder(target, median),tot,fill=type)) + geom_violin(draw_quantiles = 0.5) + facet_wrap(~type, scales = "free_y") +
  ggpubr::stat_compare_means() +
  facets_nice + scale_y_log10() + theme_classic(base_size = 12) + theme(legend.position = "none") + labs(x = "", y = "freq") +
  ggpubr::rotate_x_text(45) + facets_nice + scale_fill_manual(values = c("ivory", "salmon"))
ggsave("results/tcrgp/box_tcrgp_predictions.pdf", width = 7, height = 3.5)

ag_specific_filt %>%
  filter(timepoint == "dg") %>% filter(target == "pr1_cdr3b") %>%
  group_by(filename,type) %>% summarise(tot = sum(freq), n = n(), gini = ineq::Gini(freq)) %>%
  ggplot(aes(n,tot,fill=type,size=gini)) + geom_point(shape=21) +
  theme_classic(base_size = 12) +
  scale_x_log10() + scale_y_log10() + labs(x = "n anti-PR1 clonotypes", y = "freq anti-PR1 clonotypes") + stat_ellipse(aes(color=type),size=1) + scale_fill_manual(values = getPalette3(4)) + scale_color_manual(values = getPalette3(4))
ggsave("results/tcrgp/scatter_pb_bm.pdf", width = 5.5, height = 4)

ag_specific_filt %>%
  filter(type == "BM") %>% filter(target == "pr1_cdr3b") %>%
  mutate(timepoint = as.character(timepoint)) %>%
  mutate(timepoint = ifelse(timepoint == "unknown", "dg", timepoint)) %>%
  mutate(timepoint = factor(timepoint, levels = c("dg", "6m"))) %>%
  group_by(filename,timepoint) %>% summarise(tot = sum(freq), n = n(), gini = ineq::Gini(freq)) %>%
  ggplot(aes(n,tot,fill=timepoint,size=gini)) + geom_point(shape=21) +
  theme_classic(base_size = 12) +
  scale_x_log10() + scale_y_log10() + labs(x = "n anti-PR1 clonotypes", y = "freq anti-PR1 clonotypes") + stat_ellipse(aes(color=timepoint),size=1) + scale_fill_manual(values = getPalette(6)) + scale_color_manual(values = getPalette(6))
ggsave("results/tcrgp/scatter_bm_timepoint.pdf", width = 5.5, height = 4)

ag_specific_filt %>%
  filter(type != "BM") %>% filter(target == "pr1_cdr3b") %>%
  group_by(filename,timepoint) %>% summarise(tot = sum(freq), n = n(), gini = ineq::Gini(freq)) %>%
  ggplot(aes(n,tot,fill=timepoint,size=gini)) + geom_point(shape=21) +
  theme_classic(base_size = 12) +
  scale_x_log10() + scale_y_log10() + labs(x = "n anti-PR1 clonotypes", y = "freq anti-PR1 clonotypes")  + scale_fill_manual(values = getPalette4(8)) + scale_color_manual(values = getPalette4(8))
ggsave("results/tcrgp/scatter_bm_timepoint.pdf", width = 5.5, height = 4)


ag_specific_filt %>%
  filter(target == "pr1_cdr3b" & timepoint != "none") %>%
  mutate(timepoint = as.character(timepoint)) %>%
  mutate(timepoint = ifelse(timepoint == "unknown", "dg", timepoint)) %>%
  mutate(timepoint = factor(timepoint, levels = c("dg", "6m"))) %>%
  filter(type == "BM") %>%
  group_by(name,timepoint,type,target,design,tki) %>% summarise(tot = sum(freq), n = n()) %>%
  ggplot(aes(timepoint,tot,fill=timepoint)) + geom_path() + scale_y_log10() + theme_classic(base_size = 12) +
  ggpubr::rotate_x_text(angle = 0) + geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(label="p.format") + geom_jitter(size=0.5) + scale_fill_manual(values = getPalette4(4)) + theme(legend.position = "none") + labs(x = "", y = "freq of anti-PR1 clonotypes") # + geom_path(aes(group=name),alpha=0.5,color="lightgrey")
ggsave("results/tcrgp/box_pr1_freq_bm.pdf", width = 5.5, height = 4)

ag_specific_filt %>%
  filter(target == "pr1_cdr3b" & timepoint != "none") %>%
  mutate(timepoint = as.character(timepoint)) %>%
  mutate(timepoint = ifelse(timepoint == "unknown", "dg", timepoint)) %>%
  mutate(timepoint = factor(timepoint, levels = c("dg", "6m"))) %>%
  filter(type == "BM") %>%
  group_by(name,timepoint,type,target,design,tki) %>% summarise(tot = sum(freq), n = n()) %>%
  ggplot(aes(timepoint,tot,fill=timepoint)) + geom_path() + scale_y_log10() + theme_classic(base_size = 12) +
  ggpubr::rotate_x_text(angle = 0) + geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(label="p.format") + geom_jitter(size=0.5) + scale_fill_manual(values = getPalette4(4)) + facet_wrap(~tki) + facets_nice + theme(legend.position = "none") + labs(x = "", y = "freq of anti-PR1 clonotypes") # + geom_path(aes(group=name),alpha=0.5,color="lightgrey")
ggsave("results/tcrgp/box_pr1_freq_bm_tki.pdf", width = 9, height = 8)

ag_specific_filt %>%
  filter(target == "pr1_cdr3b" & timepoint != "none") %>%
  filter(type == "PB") %>% filter(timepoint != "2nd_baseline") %>% filter(!tki %in% c("bosutinib", "nilotinib")) %>% mutate(tki = factor(as.character(tki), levels = c("imatinib", "dasatinib"))) %>%
  group_by(name,timepoint,type,target,design,tki) %>% summarise(tot = sum(freq), n = n()) %>%
  ggplot(aes(timepoint,tot,fill=timepoint)) + geom_path() + scale_y_log10() + theme_classic(base_size = 12) + geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.5) + scale_fill_manual(values = getPalette4(8))  +
  ggpubr::stat_compare_means(label="p.format") + facets_nice + labs(x = "", y = "freq of anti-PR1 clonotypes") + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(angle = 45)
ggsave("results/tcrgp/box_pr1_freq_pb.pdf", width = 5, height = 4)

ag_specific_filt %>%
  filter(target == "pr1_cdr3b" & timepoint != "none") %>%
  filter(type == "PB") %>% filter(timepoint != "2nd_baseline") %>% filter(!tki %in% c("bosutinib", "nilotinib")) %>% mutate(tki = factor(as.character(tki), levels = c("imatinib", "dasatinib"))) %>%
  group_by(name,timepoint,type,target,design,tki) %>% summarise(tot = sum(freq), n = n()) %>%
  ggplot(aes(timepoint,tot,fill=timepoint)) + geom_path() + scale_y_log10() + theme_classic(base_size = 12) + geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.5) + scale_fill_manual(values = getPalette4(8)) + facet_wrap(~tki, scales = "free") +
  ggpubr::stat_compare_means(label="p.format") + facets_nice + labs(x = "", y = "freq of anti-PR1 clonotypes") + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(angle = 45)
ggsave("results/tcrgp/box_pr1_freq_pb_tki.pdf", width = 5.5, height = 4)

ag_specific_filt %>% ungroup() %>%
  filter(target == "pr1_cdr3b" & timepoint != "none") %>%
  filter(type == "PB") %>% filter(timepoint != "2nd_baseline") %>% filter(!tki %in% c("bosutinib", "nilotinib", "imatinib")) %>% mutate(tki = factor(as.character(tki), levels = c("dasatinib"))) %>%
  group_by(name,timepoint,type,target,design,tki) %>% summarise(tot = sum(freq), n = n()) %>%
  ggplot(aes(timepoint,tot,fill=timepoint)) + geom_path() + scale_y_log10() + theme_classic(base_size = 12) + geom_boxplot(outlier.shape = NA) + geom_jitter(size=0.5) + scale_fill_manual(values = getPalette4(8)) + facet_wrap(~tki, scales = "free") +
  ggpubr::stat_compare_means(label="p.format") + facets_nice + labs(x = "", y = "freq of anti-PR1 clonotypes") + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(angle = 45)


ag_specific_filt %>%
  filter(target == "pr1_cdr3b" & timepoint != "none") %>%
  mutate(timepoint = ifelse(timepoint %in% c("0m", "relapse"), as.character(timepoint), "follow-up")) %>%
  group_by(name,timepoint,type,target) %>% summarise(tot = sum(freq), n = n(), clonality = ineq::Gini(freq/sum(freq))) %>%
  ggplot(aes(n,tot,size = log10(clonality), color=timepoint)) + geom_point() + geom_path(aes(group=name),size=1) + ggpubr::stat_compare_means() + stat_ellipse()

df <- ag_specific_filt %>%
  filter(target == "pr1_cdr3b") %>%
  mutate(timepoint = ifelse(timepoint %in% c("0m", "relapse"), as.character(timepoint), "follow-up")) %>%
  filter(timepoint %in% c("0m", "follow-up")) %>%
  group_by(name,timepoint) %>% summarise(tot = sum(freq), n = n(), clonality = ineq::Gini(freq/sum(freq)))

names.to.keep <- df %>% group_by(name) %>% summarise(n = n()) %>% filter(n > 1) %>% pull(name)

df %>% filter(name %in% names.to.keep) %>%
  ggplot(aes(timepoint,clonality,fill=timepoint)) + geom_boxplot() + ggpubr::stat_compare_means(paired = T) + geom_path(aes(group=name),alpha=0.2) +
  scale_fill_manual(values = getPalette3(4)) + add_guide + labs(x = "", y = "clonality (Gini) \nof PR1-specific clonotypes")





