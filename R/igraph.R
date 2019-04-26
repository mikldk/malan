#' Convert igraph to population
#' 
#' @param x igraph, must be a forest of directed trees with unique positive integer names (as they will be pid's)
#' @param \dots Ignored
#' 
#' @return A population
#' 
#' @examples
#' g <- graph_from_literal( 2 +- 1 -+ 3, 4 -+ 5 )
#' plot(g)
#' pop <- from_igraph(g)
#' peds <- build_pedigrees(pop, progress = FALSE)
#' plot(peds)
#' infer_generations(peds)
#' get_generation(get_individual(pop, 1))
#' get_generation(get_individual(pop, 2))
#' get_generation(get_individual(pop, 3))
#' get_generation(get_individual(pop, 4))
#' get_generation(get_individual(pop, 5))
#' 
#' @importFrom igraph is_directed girth V get.edgelist
#' @export
from_igraph <- function(x, ...) {
  if (!is(x, "igraph")) stop("x must be an igraph")
  if (!igraph::is_directed(x)) stop("x must be a directed graph")
  if (igraph::girth(x, circle = FALSE)$girth != 0L) stop("x must be a tree (or a forest)")
  
  nms_chr <- names(igraph::V(x))
  nms <- strtoi(nms_chr)
  
  if (any(is.na(nms))) stop("x's vertices must have integer names (as they will pid's)")
  if (any(nms <= 0L)) stop("the vertex names must be positive integer")
  
  el_chr <- igraph::get.edgelist(x)
  el <- apply(el_chr, 2, strtoi)
  
  if (any(is.na(el))) stop("x's edges must be between positive integer named variables (as they will pid's)")
  if (any(el <= 0L)) stop("the edges must be between positive integers")
  
  pop <- from_igraph_rcpp(nms, el)
    
  return(pop)
}

