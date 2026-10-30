loc_map <- function(path = file.path("..", "Tables", "ferret_LOC_to_symbol_map.csv")) {
  m <- utils::read.csv(path, stringsAsFactors = FALSE)
  stats::setNames(m$new_symbol, m$loc_id)
}
rename_loc <- function(genes, map = loc_map()) {
  hit <- genes %in% names(map)
  new <- map[genes[hit]]
  clash <- new %in% genes[!hit] | duplicated(new)
  new[clash] <- genes[hit][clash]
  genes[hit] <- new
  genes
}

rename_loc_matrix <- function(mat, map = loc_map()) {
  rownames(mat) <- rename_loc(rownames(mat), map)
  mat
}
loc_report <- function(genes, panel, map = loc_map()) {
  after <- rename_loc(genes, map)
  data.frame(
    gene      = panel,
    before    = panel %in% genes,
    after     = panel %in% after,
    rescued   = !(panel %in% genes) & panel %in% after,
    stringsAsFactors = FALSE
  )
}
