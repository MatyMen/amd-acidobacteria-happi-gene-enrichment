library(dplyr)
library(readr)
library(readxl)
library(parallel)
library(happi)

genera <- c("Acidicapsa", "Granulicella", "JAKAAZ01", "Telmatobacter")
meta <- read_excel("Acido_newniche_Matilde_Andrea1.xlsx") %>% transmute(genome = Accession, group = ifelse(Niche == "Acid Mine Drainage", "AMD", "nonAMD"), completeness = checkm_completeness)

# HAPPI on one cluster; retry with more random starts (1, 5, 20)
run_happi <- function(y, target, quality) {
  for (ns in c(1, 5, 20)) {
    res <- try(happi::happi(outcome = y, covariate = model.matrix(~ target), h0_param = 2,
                            quality_var = quality, max_iterations = 1000, change_threshold = 0.001,
                            epsilon = 0.1, method = "splines", firth = TRUE, spline_df = 3,
                            nstarts = ns, seed = 13), silent = TRUE)
    if (!inherits(res, "try-error")) break
  }
  if (inherits(res, "try-error")) return(c(NA_real_, NA_real_))
  c(as.numeric(tail(res$loglik$pvalue[!is.na(res$loglik$pvalue)], 1)),            #pval
    as.numeric(tail(res$beta[!is.na(res$beta[, 1]), , drop = FALSE], 1)[1, 2]))   #beta of (AMD vs nonAMD)
}
cl <- makeCluster(max(1, detectCores() - 1))
clusterEvalQ(cl, library(happi)); clusterExport(cl, "run_happi")
dir.create("results", showWarnings = FALSE)
all_results <- list()

for (g in genera) {
  sel <- read_tsv(paste0("enrichment_results_", g, ".txt"), show_col_types = FALSE) %>% filter(adjusted_q_value < 0.01) %>%
    mutate(class = case_when(p_AMD > 0.6 & p_nonAMD < 0.1 ~ "AMD",p_nonAMD > 0.6 & p_AMD < 0.1 ~ "nonAMD")) %>% filter(!is.na(class))
  rtab <- list.files(pattern = paste0("^gene_presence_absence_", g, "(_[0-9]+)?\\.Rtab$"))[1]
  pa <- read_tsv(rtab, col_types = cols(Gene = col_character(), .default = col_integer()))
  pa <- pa[pa$Gene %in% sel$key, ]
  m <- t(as.matrix(pa[, -1])); colnames(m) <- pa$Gene

  #completeness of each genome
  info <- meta[match(rownames(m), meta$genome),]
  stopifnot(!anyNA(info$group))
  target <- as.integer(info$group == "AMD"); quality <- info$completeness
  cat(g,":", nrow(m), "genomes,", ncol(m), "clusters\n")

  #happi on every cluster
  out <- parLapply(cl, as.data.frame(m), run_happi, target, quality)
  res <- data.frame(genus = g, key = colnames(m),pvalue = sapply(out, `[`, 1), beta_AMD = sapply(out, `[`, 2))
  res$qvalue <- p.adjust(res$pvalue, method = "BH")
  res <- res %>% left_join(sel %>% select(key, class, enrichment_score, p_AMD, p_nonAMD), by = "key") %>% arrange(pvalue)
  write_tsv(res, paste0("results/happi_results_", g, ".tsv"))
  all_results[[g]] <- res
}
stopCluster(cl)
write_tsv(bind_rows(all_results), "results/happi_results_ALL.tsv")
