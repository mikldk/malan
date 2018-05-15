## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(fig.width = 7)

## ---- message=FALSE------------------------------------------------------
library(malan)

## ------------------------------------------------------------------------
set.seed(1)

## ------------------------------------------------------------------------
population_sizes <- round(c(rep(100, 10), 100*1.02^(1:10)))
plot(population_sizes, type = "b")

## ------------------------------------------------------------------------
set.seed(1) # For reproducibility
sim_res_growth <- sample_geneology_varying_size(
  population_sizes = population_sizes, 
  
  # VRS = 0.2:
  enable_gamma_variance_extension = TRUE,
  gamma_parameter_shape = 5,
  gamma_parameter_scale = 1/5,
  
  # Live population: 
  # 3 generations
  generations_full = 3, 
  generations_return = 3,
  
  progress = FALSE)

## ------------------------------------------------------------------------
live_pop <- sim_res_growth$individuals_generations

## ------------------------------------------------------------------------
pedigrees <- build_pedigrees(sim_res_growth$population, progress = FALSE)
pedigrees
pedigrees_count(pedigrees)
pedigrees_table(pedigrees)
pedigree_size(pedigrees[[1]])

## ------------------------------------------------------------------------
g <- as_tbl_graph(pedigrees)
g

## ------------------------------------------------------------------------
if (requireNamespace("ggraph", quietly = TRUE)) {
  library(ggraph)
  p <- ggraph(g, layout = 'tree') +
    geom_edge_link() +
    geom_node_point(size = 8) +
    geom_node_text(aes(label = name), color = "white") +
    facet_nodes(~ ped_id) +
    theme_graph() 
  
  print(p)
}

## ------------------------------------------------------------------------
PED_ID <- 1

g_ped2 <- g %>% 
  activate(nodes) %>% 
  filter(ped_id == PED_ID)

if (requireNamespace("ggraph", quietly = TRUE)) {
  library(ggraph)
  p <- ggraph(g_ped2, layout = 'tree') +
      geom_edge_link() +
      geom_node_point(size = 8) +
      geom_node_text(aes(label = name), color = "white") +
      theme_graph() 
  print(p)
}

## ------------------------------------------------------------------------
ystr_markers

## ------------------------------------------------------------------------
ystr_kits
ystr_kits %>% count(Kit)

## ------------------------------------------------------------------------
partial_kit <- ystr_kits %>% 
  filter(Kit == "PowerPlex Y23") %>% 
  inner_join(ystr_markers, by = "Marker") %>% 
  filter(!(Marker %in% c("DYS437", "DYS448"))) %>% 
  rowwise() %>% # To work on each row
  mutate(IntegerAlleles = list(Alleles[Alleles == round(Alleles)]),
         MinIntAllele = min(IntegerAlleles),
         MaxIntAllele = max(IntegerAlleles)) %>% 
  ungroup() %>% 
  select(-Kit, -Alleles)
partial_kit

## ------------------------------------------------------------------------
mu <- partial_kit %>% pull(MutProb)
mu

## ------------------------------------------------------------------------
generate_random_haplotype <- function() {
  partial_kit %>% 
    rowwise() %>% 
    mutate(Allele = IntegerAlleles[sample.int(length(IntegerAlleles), 1)]) %>% 
    pull(Allele)
}

## ------------------------------------------------------------------------
generate_random_haplotype()
generate_random_haplotype()

## ------------------------------------------------------------------------
set.seed(1)
pedigrees_all_populate_haplotypes_custom_founders(
  pedigrees = pedigrees, 
  get_founder_haplotype = generate_random_haplotype,
  mutation_rates = mu,
  progress = FALSE)

## ------------------------------------------------------------------------
g_ped2 <- as_tbl_graph(pedigrees) %>% 
  activate(nodes) %>% 
  filter(ped_id == PED_ID) %>%
  group_by(name) %>% 
  mutate(haplotype_str = paste0(haplotype[[1]], collapse = ";"))
  #mutate(haplotype_str = map(haplotype, paste0, collapse = ";")[[1]])

if (requireNamespace("ggraph", quietly = TRUE)) {
  library(ggraph)
  p <- ggraph(g_ped2, layout = 'tree') +
      geom_edge_link() +
      geom_node_point(aes(color = haplotype_str), size = 8) +
      geom_node_text(aes(label = name), color = "white") +
      theme_graph() 
  print(p)
}

## ------------------------------------------------------------------------
set.seed(5)
Q_index <- sample.int(n = length(live_pop), size = 1)
Q <- live_pop[[Q_index]]
print_individual(Q)
Q_hap <- get_haplotype(Q)
Q_hap

## ------------------------------------------------------------------------
Q_ped <- get_pedigree_from_individual(Q)
ped_live_matches <- count_haplotype_occurrences_pedigree(
  pedigree = Q_ped, 
  haplotype = Q_hap, 
  generation_upper_bound_in_result = 2) # gen 0, 1, 2
pop_live_matches <- count_haplotype_occurrences_individuals(
  individuals = live_pop, 
  haplotype = Q_hap)

ped_live_matches
pop_live_matches

## ------------------------------------------------------------------------
path_details <- pedigree_haplotype_matches_in_pedigree_meiosis_L1_dists(
  suspect = Q, 
  generation_upper_bound_in_result = 2)

## ------------------------------------------------------------------------
nrow(path_details)
head(path_details)

