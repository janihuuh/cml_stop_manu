
vdjToGliph <- function(vdj_df){

  ## Write gliph files to the vdj files
  # @ param
  # input: df from vdj

  df <- vdj_df %>% dplyr::select(cdr3aa, v, j, name, count) %>%
    dplyr::rename("CDR3b"   = "cdr3aa",
                  "TRBV"    = "v",
                  "TRBJ"    = "j",
                  "Patient" =  name,
                  "Counts"  = count)

  return(df)

}

## Run Pr1-specific
pr1_specific_data <- lapply(list.files("tcrb_data/sorted/pr1_tolerable/", full.names = T, pattern = "PR1"), FUN = function(x) x %>% fread() %>% filter(count > 0) %>% mutate(name = extractFileName(x) %>% gsub(pattern = ".txt", replacement = ""))) %>% rbindlist()
pr1_specific_data <- pr1_specific_data %>% as.data.frame %>% vdjToGliph()

df_test <- pr1_specific_data %>% mutate(CDR3a = NA) %>% dplyr::select(CDR3b, TRBV, TRBJ, CDR3a, Patient, Counts) %>% mutate(Patient = paste0(Patient, ":CML")) %>% dplyr::rename("subject:condition" = "Patient",	"count" = "Counts")
fwrite(df_test, "/Users/janihuuh/Dropbox/applications/gliph2.1/pr1.txt", sep = "\t", quote = F, row.names = F, col.names = F)

hla_test <- data.frame(`subject:condition` = unique(df_test$`subject:condition`)) %>% mutate(allele1 = NA, allele2 = NA)
fwrite(hla_test, "/Users/janihuuh/Dropbox/applications/gliph2.1/pr1_hla.txt", sep = "\t", quote = F, row.names = F, col.names = F)

## To tcrgp
pr1_tcrgp <- lapply(list.files("tcrb_data/sorted/", full.names = T, pattern = "PR1"), FUN = function(x) x %>% fread() %>% filter(count > 2) %>% mutate(name = extractFileName(x) %>% gsub(pattern = ".txt", replacement = ""))) %>% rbindlist() %>% filter(!grepl("\\*", cdr3aa))
pr1_tcrgp <- pr1_tcrgp %>% mutate(epitope = "pr1", va = NA, cdr3a = NA) %>% dplyr::rename(subject=name,cdr3b=cdr3aa,vb=v) %>% dplyr::select(epitope,subject,va,vb,cdr3a,cdr3b)# %>% mutate(vb = NA)

## satu-filtered
pr1_tcrgp <- lapply(list.files("tcrb_data/sorted/pr1_tolerable/", full.names = T, pattern = "PR1"), FUN = function(x) x %>% fread() %>% filter(count > 0) %>% mutate(name = extractFileName(x) %>% gsub(pattern = ".txt", replacement = ""))) %>% rbindlist() %>% filter(!grepl("\\*", cdr3aa))
pr1_tcrgp <- pr1_tcrgp %>% mutate(epitope = "pr1", va = NA, cdr3a = NA) %>% dplyr::rename(subject=name,cdr3b=cdr3aa,vb=v) %>% dplyr::select(epitope,subject,va,vb,cdr3a,cdr3b)# %>% mutate(vb = NA)

## Remove public sequences
pr1_tcrgp <- pr1_tcrgp[!duplicated(pr1_tcrgp$cdr3b), ]

## Remove too short or too long seqs
# pr1_tcrgp <- pr1_tcrgp[nchar(pr1_tcrgp$cdr3b) >= 10 & nchar(pr1_tcrgp$cdr3b) <= 22, ]

## Remove non-canonical tcrs
# pr1_tcrgp <- pr1_tcrgp[grep("^C", pr1_tcrgp$cdr3b, invert = F), ]
# pr1_tcrgp <- pr1_tcrgp[grep("F$", pr1_tcrgp$cdr3b, invert = F), ]

pr1_tcrgp <- pr1_tcrgp[pr1_tcrgp$cdr3b %in% pr1_gliph$TcRb , ]
fwrite(pr1_tcrgp, "results/tcrgp/input/pr1_tcrgp.txt", sep = "\t", quote = F, row.names = F)


## Overall
pr1_data <- lapply(list.files("tcrb_data/sorted/pr1/", full.names = T, pattern = "PR1"), FUN = function(x) x %>% fread() %>% filter(count > 2) %>% mutate(name = extractFileName(x) %>% gsub(pattern = ".txt", replacement = ""))) %>% rbindlist() %>%
  filter(!grepl("\\*", cdr3aa))
fwrite(pr1_data, "results/tcrgp/input/pr1_no_singletons.txt", sep = "\t", quote = F, row.names = F)

pr1_data %>% group_by(name) %>% summarise(n = n()) %>%
  ggplot(aes(reorder(name, n),n)) + geom_bar(stat = "identity", fill = "lightgrey") + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "nClonotypes")
ggsave("results/tcrgp/bar_pr1.pdf", width = 5, height = 4)

pr1_data %>% mutate(sample = extractName(name)) %>% group_by(sample,name) %>% summarise(n=n()) %>%
  filter(sample %in% c("1624", "788", "822", "824")) %>%
  ggplot(aes(name, n)) + geom_point() + facet_wrap(~sample, scales = "free") + ggpubr::rotate_x_text(angle = 0) + ggpubr::rotate_x_text(angle = 45)

pr1_data %>% mutate(date = gsub(pattern = name, replacement = "", name))



## Analyze gliph results
pr1_gliph <- readGliphFile("results/pr1/pr1_gliph_tolerable.csv")
pr1_gliph <- pr1_gliph %>% filter(vb_score < 0.05 & number_unique_cdr3 >= 3)


pr1_gliph %>% group_by(pattern) %>% summarise(n = n()) %>% filter(n>4) %>%
  ggplot(aes(reorder(pattern,n),n)) + geom_bar(stat="identity",fill="lightgrey") + coord_flip() + labs(x = "", y = "TCRs with motif")
ggsave("results/tcr/pr1/bar_gliph.pdf", width = 3, height = 4)

## Logoplot
require(ggseqlogo)

plotLogo <- function(seq){

  ## Determine the commonest sequence
  max_n <- table(nchar(seq)) %>% which.max %>% names %>% as.numeric()
  seq_cut <- seq[nchar(seq) == max_n]
  seq_cut <- substr(seq_cut, 4, nchar(seq_cut) - 3)

  # ggplot() + geom_logo(seq_cut,  method = 'prob') + theme_logo()
  ggplot() + geom_logo(seq_cut,  method = 'bits') + theme_logo()

}

plotLogoList <- function(seq_list){

  ## Determine the commonest sequence
  getSeq <- function(seq){

    max_n <- table(nchar(seq)) %>% which.max %>% names %>% as.numeric()
    seq_cut <- seq[nchar(seq) == max_n]
    seq_cut <- substr(seq_cut, 4, nchar(seq_cut) - 3)

  }

  seq_list_cut <- lapply(seq_list, getSeq)

  ggseqlogo(seq_list_cut,  method = 'bits', ncol=5) + theme_logo()

}


logolist <- split(pr1_gliph, f = pr1_gliph$index)[1:20]
names(logolist) <- lapply(logolist, FUN = function(x) unique(x$pattern))
seq_list <- lapply(logolist, FUN = function(x) x %>% pull(TcRb))
plotLogoList(seq_list)
ggsave("results/tcr/pr1/logoplots_top20.pdf", width = 12, height = 8)

logolist <- split(pr1_gliph, f = pr1_gliph$index)[1:10]
names(logolist) <- lapply(logolist, FUN = function(x) unique(x$pattern))
seq_list <- lapply(logolist, FUN = function(x) x %>% pull(TcRb))
plotLogoList(seq_list)
ggsave("results/tcr/pr1/logoplots_top10.pdf", width = 12, height = 4)

df <- pr1_gliph[!duplicated(pr1_gliph$pattern), ]
ggplot(data = df, aes(number_subject, number_unique_cdr3, label = pattern)) + geom_jitter() + ggrepel::geom_text_repel(data = subset(df, df$number_subject > 2 & df$number_unique_cdr3), aes(number_subject, number_unique_cdr3, label = pattern))
ggsave("results/tcr/pr1/scatter_subject_unique_cdr3.pdf", width = 5, height = 4)


p <- NULL
i <- 1

for(pattern_temp in pr1_gliph$pattern){

  temp <- pr1_gliph %>% filter(pattern == pattern_temp)
  p[[i]] <- combn(temp$TcRb, 2) %>% t() %>% as.data.frame()
  i <- i + 1

}

el  <- p %>% rbindlist()
el  <- el[!duplicated(el), ]
el2 <- el

g    <- graph.data.frame(el2, directed = F)
am   <- get.adjacency(g,sparse=T)
g_am <- graph_from_adjacency_matrix(am, mode = "undirected")
g_am <- simplify(g_am, remove.multiple = F, remove.loops = T)

## Plot
V(g_am)$size <- 3
V(g_am)$frame.color <- "white"
V(g_am)$color <- "dodgerblue"
V(g_am)$label <- ""
E(g_am)$arrow.mode <- 0

layout = layout_with_graphopt(g_am)

pdf("results/tcrgp/network_gliph.pdf", width = 5, height = 4)
plot.igraph(g_am, layout = layout)
dev.off()

df <- fread("data/clinical/scrnaseq_clinical2.txt")
dim(df)
