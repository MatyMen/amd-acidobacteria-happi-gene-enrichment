Purpose: test whether gene clusters selected from the anvi'o functional enrichment are really associated with AMD or non-AMD genomes, correcting for genome completeness (HAPPI).

Steps (run separately for each genus):
1. Select clusters from enrichment_results_<genus>.txt:
    adjusted_q_value < 0.01
    AMD class:    p_AMD > 0.6    & p_nonAMD < 0.1
    nonAMD class: p_nonAMD > 0.6 & p_AMD < 0.1
2. Read their presence/absence (0/1) from gene_presence_absence_<genus>.Rtab
3. Match genomes with group (Niche) and completeness
(checkm_completeness) from the Excel table
4. Run HAPPI on each cluster (AMD vs non-AMD, completeness as quality
variable); failed fits are retried with more starts (1, 5, 20)
5. Benjamini-Hochberg q-value within each genus

Inputs:  enrichment_results_<genus>.txt, gene_presence_absence_<genus>.Rtab,
         Acido_newniche_Matilde_Andrea1.xlsx
Outputs: results/happi_results_<genus>.tsv, results/happi_results_ALL.tsv
         (beta_AMD > 0 = more common in AMD, < 0 = more common in non-AMD)
Settings: epsilon = 0.1 (assumed false-positive rate; HAPPI default is 0),
         max_iterations = 1000, change_threshold = 0.001, seed = 13.
