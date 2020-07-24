
## Analyze diversity on resampled to 10k reads
div_10k_healthy <- fread("/Users/janihuuh/Dropbox/cml_stop/tcrb_data/results/diversity/10k_sampled_w_healthy.divresity.aa.resampled.sm.txt")

div_10k_healthy <- div_10k_healthy %>%
  mutate(disease = ifelse(grepl("HIP", sample_id),  "healthy", "cml")) %>%
  mutate(disease = ifelse(grepl("Keck", sample_id),  "healthy", disease)) %>%
  mutate(sample_id = gsub("\\-", "\\_", sample_id)) %>%
  mutate(name = extractName(sample_id)) %>%
  mutate(clonality = 1 - normalizedShannonWienerIndex_mean) %>%
  dplyr::rename(cohort = metadata_blank,
                drug = `drug (1 imatinib), 2 (dasatinib), 3 (nilotinib)`,
                timepoint = `baseline (0), 1 month (1), 3 months (2), 6 months (3), relapse (4), baseline dastop (5)`)

div_10k_healthy$drug      <- plyr::revalue(as.factor(div_10k_healthy$drug), c("1" = "imatinib", "2" = "dasatinib", "3" = "nilotinib"))
div_10k_healthy$timepoint <- plyr::revalue(as.factor(div_10k_healthy$timepoint), c("0" = "baseline", "1" = "1mo", "2" = "3mo", "3" = "6mo", "4" = "relapse", "5" = "baseline dastop"))

div_columns <- c( "name", "sample_id", "disease", "reads", "diversity", "clonality", "cohort", "drug", "timepoint",
                  grep("mean", colnames(div_10k_healthy), value = T))

div_10k_healthy <- div_10k_healthy %>% select(div_columns) %>% select(-chaoE_mean)

div_10k_healthy %>% melt(id = c("name", "sample_id", "disease")) %>%
  ggplot(aes(disease,value,fill=disease)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) + facet_wrap(~variable, scales = "free", ncol = 5) + theme_bw() + ggsignif::geom_signif(comparisons = list(c("cml", "healthy"))) +
  scale_fill_manual(values = c("salmon", "lightgrey")) + theme(legend.position = "none")
ggsave("results/diversity/plots/box_total_10k_sampled.png", width = 12, height = 8)

div_10k_healthy %>%
  filter(!duplicated(name)) %>%
  melt(id = c("name", "sample_id", "disease")) %>%
  ggplot(aes(disease,value,fill=disease)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) + facet_wrap(~variable, scales = "free", ncol = 5) + theme_bw() + ggsignif::geom_signif(comparisons = list(c("cml", "healthy"))) +
  scale_fill_manual(values = c("salmon", "lightgrey")) + theme(legend.position = "none")
ggsave("results/diversity/plots/box_total_10k_sampled_unique.png", width = 12, height = 8)



## Only euroski, and their course as function of time
euroski <- div_10k_healthy %>% filter(cohort == "euroski")

ggplot(data = euroski, aes(timepoint, clonality,group=name,color=drug, fill=drug)) + geom_point(shape = 21, size = 3) + geom_path() + theme_bw() + scale_fill_manual(values = getPalette(4)) + scale_color_manual(values = getPalette(4)) +
  ggrepel::geom_text_repel(data = subset(euroski, timepoint == "baseline"), aes(timepoint, clonality, label=name), vjust = 0, nudge_y = 0, nudge_x = -0.2)
ggsave("results/diversity/plots/line_euroski_clonality_10k.png", width = 6, height = 4)

ggplot(data = euroski, aes(timepoint, clonality,group=name,color=drug, fill=drug)) + geom_point(shape = 21, size = 3) + geom_path() + theme_bw() + scale_fill_manual(values = getPalette(4)) + scale_color_manual(values = getPalette(4)) +
  ggrepel::geom_text_repel(data = subset(euroski, timepoint == "baseline"), aes(timepoint, clonality, label=name), vjust = 0, nudge_y = 0, nudge_x = -0.2) +
  ggsignif::geom_signif(comparisons = list(c("baseline", "1mo")), test.args = c(paired = T))
ggsave("results/diversity/plots/line_euroski_clonality_10k.png", width = 6, height = 4)


ggplot(data = euroski, aes(timepoint, clonality,group=name,color=drug, fill=drug)) + geom_point(shape = 21, size = 3) + geom_path() + theme_bw() + scale_fill_manual(values = getPalette(4)) + scale_color_manual(values = getPalette(4)) +
  ggrepel::geom_text_repel(data = subset(euroski, timepoint == "baseline"), aes(timepoint, clonality, label=name), vjust = 0, nudge_y = 0, nudge_x = -0.2) +
  ggsignif::geom_signif(comparisons = list(c("baseline", "1mo"))) + facet_wrap(~drug) + theme(legend.position = "none")
ggsave("results/diversity/plots/line_euroski_clonality_10k_2.png", width = 8, height = 4)

ggplot(data = subset(euroski, drug == "dasatinib"), aes(timepoint, clonality,group=name,color=drug, fill=drug)) + geom_point(shape = 21, size = 3) + geom_path() + theme_bw() + scale_fill_manual(values = getPalette(4)) + scale_color_manual(values = getPalette(4)) +
  ggsignif::geom_signif(comparisons = list(c("baseline", "1mo")), test.args = c(paired = T)) + facet_wrap(~drug)
ggsave("results/diversity/plots/line_euroski_clonality_10k_dasa.png", width = 6, height = 4)



div_10k_healthy %>% filter(cohort == "euroski") %>% filter(timepoint == "baseline") %>%
  ggplot(aes(drug, clonality,fill=drug)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) +
  ggsignif::geom_signif(comparisons = list(c("imatinib", "dasatinib"), c("imatinib", "nilotinib"), c("nilotinib", "dasatinib")), step_increase = 0.1) + theme_bw() +
  theme(legend.position = "none") + scale_fill_manual(values = getPalette(4))
ggsave("results/diversity/plots/box_euroski_clonality_10k_baseline.png", width = 5, height = 4)

div_10k_healthy %>% filter(cohort == "euroski" & timepoint == "baseline" | cohort == ".") %>%
  ggplot(aes(disease, clonality,fill=disease)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) +
  ggsignif::geom_signif(comparisons = list(c("healthy", "cml")), step_increase = 0.1) + theme_bw() +
  theme(legend.position = "none") + scale_fill_manual(values = getPalette(4))
ggsave("results/diversity/plots/box_euroski_clonality_10k_baseline_healthy.png", width = 5, height = 4)


## Only dastop
div_10k_healthy %>% filter(cohort == "dastop2" | cohort == ".") %>%
  ggplot(aes(disease, clonality,fill=disease)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) +
  ggsignif::geom_signif(comparisons = list(c("healthy", "cml")), step_increase = 0.1) + theme_bw() +
  theme(legend.position = "none") + scale_fill_manual(values = getPalette(4))
ggsave("results/diversity/plots/box_dasastop_clonality_10k_baseline_healthy.png", width = 5, height = 4)



## dastop2 vs euroski
div_10k_healthy %>% filter(cohort == "euroski" & timepoint == "baseline" | cohort == "." | cohort == "dastop2") %>%
  ggplot(aes(disease, clonality,fill=disease)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) +
  ggsignif::geom_signif(comparisons = list(c("healthy", "cml")), step_increase = 0.1) + theme_bw() +
  theme(legend.position = "none") + scale_fill_manual(values = getPalette(4))
ggsave("results/diversity/plots/box_dasastop_euroski_clonality_10k_baseline_healthy.png", width = 5, height = 4)

div_10k_healthy %>% filter(cohort == "euroski" & timepoint == "baseline" | cohort == "." | cohort == "dastop2") %>%
  ggplot(aes(cohort, clonality,fill=disease)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) +
  ggsignif::geom_signif(comparisons = list(c("dastop2", "euroski")), step_increase = 0.1) + theme_bw() +
  theme(legend.position = "none") + scale_fill_manual(values = getPalette(4))
ggsave("results/diversity/plots/box_dasastop_vs_euroski_clonality_10k_baseline_healthy.png", width = 5, height = 4)








## Analyze diversity without resampling to 10k reads
div_healthy <- fread("results/diversity/10k_sampled_w_healthy.diversity.aa.exact.txt")
div_healthy <- div_healthy %>%
  mutate(disease = ifelse(grepl("HIP", sample_id),  "healthy", "cml")) %>%
  mutate(disease = ifelse(grepl("Keck", sample_id),  "healthy", disease)) %>%
  mutate(sample_id = gsub("\\-", "\\_", sample_id)) %>%
  mutate(name = extractName(sample_id)) %>%
  mutate(clonality = 1 - normalizedShannonWienerIndex_mean)

div_columns <- c( "name", "sample_id", "disease", "reads", "diversity", "clonality",
                  grep("mean", colnames(div_healthy), value = T))
div_healthy <- div_healthy %>% select(div_columns) %>% select(-chaoE_mean)

div_healthy %>% melt(id = c("name", "sample_id", "disease")) %>%
  ggplot(aes(disease,value,fill=disease)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) + facet_wrap(~variable, scales = "free", ncol = 5) + theme_bw() + ggsignif::geom_signif(comparisons = list(c("cml", "healthy"))) +
  scale_fill_manual(values = c("salmon", "lightgrey")) + theme(legend.position = "none")
ggsave("results/diversity/plots/box_total_not_sampled.png", width = 12, height = 8)

div_healthy %>%
  filter(!duplicated(name)) %>%
  melt(id = c("name", "sample_id", "disease")) %>%
  ggplot(aes(disease,value,fill=disease)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) + facet_wrap(~variable, scales = "free", ncol = 5) + theme_bw() + ggsignif::geom_signif(comparisons = list(c("cml", "healthy"))) +
  scale_fill_manual(values = c("salmon", "lightgrey")) + theme(legend.position = "none")
ggsave("results/diversity/plots/box_total_not_sampled_unique.png", width = 12, height = 8)

div_healthy %>%
  ggplot(aes(reads,clonality)) + geom_point() + geom_smooth(method = "lm") + scale_x_log10()


