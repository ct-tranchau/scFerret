root <- if (basename(getwd()) == "Analysis") normalizePath("..") else getwd()
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })
options(future.globals.maxSize = 64 * 1024^3, dplyr.summarise.inform = FALSE)
set.seed(1)

run <- Sys.getenv("GR0389_RUN", "v2_ncbi")
stopifnot(run %in% c("v1_ensembl", "v2_ncbi"))

dat <- file.path(root, "Data")
fig <- file.path(root, "Figures", run)
tab <- file.path(root, "Tables", run)
objdir <- file.path(root, "Working", paste0("objects_", run))
gfile <- file.path(objdir, "global.rds")
for (d in c(fig, tab, objdir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
ss <- read.csv(file.path(dat, "sample_sheet.csv"), stringsAsFactors = FALSE)
tp_levels <- intersect(c("PostV1", "PostV2", "PostChallenge"), unique(ss$timepoint))

lab_levels <- ss$label[order(factor(ss$group, c("Egg", "Cell")), ss$animal, ss$timepoint)]
ss$label <- factor(ss$label, lab_levels)
labs_of <- function(x) factor(ss$label[match(as.character(x), ss$sample)], lab_levels)

if (run == "v2_ncbi") {
  crdir <- function(s) file.path(dat, "cellranger", s)
  mtx_path <- function(s) file.path(crdir(s), "sample_filtered_feature_bc_matrix")
  vdj_dir <- function(s, ch) file.path(crdir(s), ch)
  source(file.path(root, "Analysis/rename_loc_genes.R"))
  gene_map <- loc_map(file.path(root, "Tables/ferret_LOC_to_symbol_map.csv"))
} else {
  crdir <- function(s) file.path(dat, s)
  mtx_path <- function(s) file.path(dat, "matrices", s)
  vdj_dir <- function(s, ch) file.path(dat, "vdj", s, ch)
  gene_map <- character()
}

vdj_path <- function(s, ch, f = "filtered_contig_annotations.csv") file.path(vdj_dir(s, ch), f)

mt_genes <- c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB")

ig_tr <- "^(IG[HKL][VDJC]|IGHM$|IGHA|IGHG|IGHD$|IGHE$|TR[ABDG][VDJC])"
ig_tr_loc <- local({
  f <- file.path(root, "Tables/ferret_ig_tr_genes.txt")
  if (run == "v2_ncbi" && file.exists(f)) readLines(f) else character()
})
drop_ig <- function(x) setdiff(grep(ig_tr, x, value = TRUE, invert = TRUE, ignore.case = TRUE), ig_tr_loc)

theme_set(theme_classic(base_size = 11) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold")))
rot <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
big <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 14), axis.text.y = element_text(size = 14), axis.title = element_text(size = 15), legend.text = element_text(size = 14))

ct_cols <- c("B cells" = "#1F77B4", "CD4 T cells" = "#2CA02C", "CD8 T cells" = "#98DF8A", "NK cells" = "#D62728", "Monocytes" = "#FF7F0E", "DCs" = "#E377C2",
             "Gamma delta T cells" = "#8C564B", "Neutrophils" = "#BCBD22", "Eosinophils" = "#17BECF", "Platelet-bound cells" = "#F7B6D2", "Doublets" = "#525252")
bd_cols <- c("Pre-B" = "#8C564B", "Immature B" = "#17BECF", "Naive B" = "#1F77B4", "Unswitched memory B" = "#7B2CBF", "Switched memory B" = "#2CA02C",
             "Exhausted B" = "#D62728", "Plasmablasts" = "#F781BF", "Plasma cells" = "#A6761D", "Damaged B" = "#C7C7C7", "Doublets" = "#525252")
t_cols <- c("CD4 T naive" = "#A1D99B", "CD4 T memory" = "#31A354", "Treg" = "#006D2C", "Tfh-like CD4 T" = "#E6550D", "CD8 T naive" = "#7FCDBB",
            "CD8 T effector memory" = "#225EA8", "Gamma delta T" = "#B15928", "NK cells" = "#E7298A", "Proliferating T/NK" = "#E7BA52", "Damaged T/NK" = "#BDBDBD", "Doublets" = "#525252")
iso_cols <- c(IgM = "#9ECAE1", IgD = "#C6DBEF", IgG = "#E6550D", IgA = "#31A354", IgE = "#756BB1")
state_cols <- c(bd_cols, t_cols, ct_cols[c("Monocytes", "DCs", "Neutrophils", "Eosinophils", "Platelet-bound cells")])
state_cols <- state_cols[!duplicated(names(state_cols))]

isg <- c("ISG15", "MX1", "MX2", "IFI44", "IFI44L", "IFIT2", "IFIT3", "OAS1", "RSAD2", "USP18")
global_dot <- list(B = c("MS4A1", "CD79A", "CD79B", "CD19", "PAX5", "IGHM"), Plasma = c("JCHAIN", "PRDM1", "MZB1", "TNFRSF17"),
                   T = c("CD3E", "CD3D", "CD2", "CD5", "IL7R"), CD4 = c("CD4", "CD40LG", "FOXP3"), CD8 = c("CD8A", "CD8B"), GD = "TRDC",
                   NK = c("NCR1", "GNLY", "KLRD1", "PRF1", "NKG7"), Mono = c("CD14", "CSF1R", "FCGR3A", "VCAN", "LYZ"),
                   Neut = c("S100A8", "S100A12", "CSF3R", "LTF"), Eos = c("IL5RA", "HRH4"), DC = c("CLEC10A", "CLEC9A", "PLD4"), Plt = c("GNG11", "NRGN", "CLU"))
b_dot <- list(Naive = c("IGHM", "IGHD", "FCER2", "SELL", "IL4R", "BACH2"), Transitional = c("SOX4", "VPREB3", "CD38"), Precursor = c("RAG1", "RAG2", "VPREB1"),
              IFN = c("ISG15", "MX1", "IFI44L"), Activated = c("CD69", "CD83", "NR4A1", "EGR1", "MYC", "NME1", "PSAT1"),
              Memory = c("CD27", "TNFRSF13B", "ITGB1", "CXCR3", "CD86", "IGHG1", "IGHA2"), Atypical = c("TBX21", "FCRL5", "ZEB2"), GC = c("AICDA", "BCL6"),
              Plasmablast = c("PRDM1", "XBP1", "IRF4", "MZB1", "JCHAIN", "TNFRSF17"), Cycling = c("TOP2A", "TYMS"), QC = c("CD3E", "CD14"))
t_dot <- list(T = c("CD3E", "CD4", "CD8A", "CD8B", "TRDC"), Naive = c("CCR7", "SELL", "LEF1", "TCF7"), Memory = c("IL7R", "S100A4", "ANXA1", "CD40LG", "ICOS", "ITGB1"),
              Treg = c("FOXP3", "IL2RA", "CTLA4"), Tfh = c("CXCR5", "PDCD1", "IL21", "BCL6"), Cytotoxic = c("GZMK", "GZMB", "GZMA", "PRF1", "NKG7", "CCL5"),
              NK = c("NCR1", "GNLY", "KLRD1"), Cycling = c("TOP2A", "TYMS"), QC = c("MS4A1", "CD14"))

gfind <- function(o, x) {
  r <- if (is.character(o)) o else rownames(o)
  unique(r[match(toupper(x), toupper(r), nomatch = 0)])
}

pct_expr <- function(o, group, genes) {
  x <- GetAssayData(o, layer = "counts")[genes, , drop = FALSE] > 0
  g <- factor(group)
  m <- sapply(levels(g), function(k) Matrix::rowMeans(x[, g == k, drop = FALSE]))
  matrix(m, nrow = length(genes), dimnames = list(genes, levels(g)))
}

pg <- function(pct, gene) {
  i <- match(toupper(gene), toupper(rownames(pct)))
  if (is.na(i)) setNames(rep(0, ncol(pct)), colnames(pct)) else pct[i, ]
}

pscore <- function(pct, genes) {
  g <- rownames(pct)[toupper(rownames(pct)) %in% toupper(genes)]
  if (!length(g)) return(pg(pct, "none"))
  colMeans(pct[g, , drop = FALSE])
}

rz <- function(v) (v - median(v)) / max(mad(v), .02)

recluster <- function(o, dims = 1:20, res = 1, nfeat = 2000) {
  suppressPackageStartupMessages(library(harmony))
  o <- CreateSeuratObject(GetAssayData(o, layer = "counts"), meta.data = o@meta.data)
  o[["RNA"]] <- split(o[["RNA"]], f = o$sample)
  o <- NormalizeData(o, verbose = FALSE)
  o <- FindVariableFeatures(o, nfeatures = nfeat + 200, verbose = FALSE)
  VariableFeatures(o) <- head(drop_ig(VariableFeatures(o)), nfeat)
  o <- ScaleData(o, verbose = FALSE)
  o <- RunPCA(o, npcs = 30, verbose = FALSE)
  o <- IntegrateLayers(o, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
  o <- JoinLayers(o)
  o <- FindNeighbors(o, reduction = "harmony", dims = dims, verbose = FALSE)
  o <- FindClusters(o, resolution = res, verbose = FALSE)
  RunUMAP(o, reduction = "harmony", dims = dims, verbose = FALSE)
}

global_rules <- function(o) {
  cl <- droplevels(factor(o$seurat_clusters))
  panels <- list(
    b    = c("CD79A", "CD79B", "MS4A1", "CD19", "PAX5", "BANK1"),
    pc   = c("JCHAIN", "MZB1", "TNFRSF17", "PRDM1"),
    t    = c("CD3E", "CD3D", "CD3G", "CD2", "CD5"),
    cd4  = c("CD4", "CD40LG"),
    cd8  = c("CD8A", "CD8B"),
    gd   = c("TRDC", "TRGC1", "TRGC2"),
    nk   = c("NCR1", "GNLY", "KLRD1", "PRF1", "NKG7", "EOMES"),
    mono = c("CD14", "CSF1R", "VCAN", "CD163", "MAFB", "SLC11A1"),
    neut = c("S100A8", "S100A12", "CSF3R", "LTF", "MMP8", "CAMP", "PGLYRP1", "CXCR2"),
    eos  = c("IL5RA", "HRH4", "CCR3", "EPX"),
    dc   = c("ZBTB46", "FLT3", "CADM1"),
    plt  = c("GNG11", "NRGN", "CLU", "PPBP", "PF4", "ITGA2B"))
  pct <- pct_expr(o, cl, gfind(o, unlist(panels)))
  s <- as.data.frame(lapply(panels, function(g) pscore(pct, g)))

  s$cell_type <- with(s, case_when(
    plt >= .4 ~ "Platelet-bound cells",
    pc >= .3 & pmax(t, mono, neut) < .5 ~ "B cells",
    pmin(b, t) >= .3 ~ "Doublets",
    pmin(b, pmax(mono, neut)) >= .3 ~ "Doublets",
    pmin(t, pmax(mono, neut)) >= .3 ~ "Doublets",
    pc >= .3 | b >= .4 ~ "B cells",
    t >= .5 & gd >= .4 ~ "Gamma delta T cells",
    t >= .5 & cd8 > cd4 ~ "CD8 T cells",
    t >= .5 ~ "CD4 T cells",
    nk >= .4 ~ "NK cells",
    eos >= .3 ~ "Eosinophils",
    neut >= .3 & neut > mono ~ "Neutrophils",
    mono >= .3 ~ "Monocytes",
    dc >= .3 ~ "DCs",
    t >= .3 & cd8 > cd4 ~ "CD8 T cells",
    t >= .3 ~ "CD4 T cells",
    b >= .25 ~ "B cells",
    TRUE ~ "Doublets"))
  data.frame(cluster = colnames(pct), n_cells = as.vector(table(cl)), cell_type = s$cell_type,
             genes = as.vector(tapply(o$nFeature_RNA, cl, median)), round(s[, names(s) != "cell_type"], 3))
}

b_rules <- function(o) {
  cl <- droplevels(factor(o$seurat_clusters))
  pct <- pct_expr(o, cl, gfind(o, c(unlist(b_dot), isg, "CD3D", "CD2", "CSF1R", "SHMT2", "WDR77")))
  switched <- as.vector(tapply(o$isotype %in% c("IgG", "IgA", "IgE"), cl, sum) / pmax(tapply(!is.na(o$isotype), cl, sum), 1))
  s <- data.frame(cd3 = pg(pct, "CD3E"), mye = (pg(pct, "CD14") + pg(pct, "CSF1R")) / 2, pc = pscore(pct, c("JCHAIN", "PRDM1", "TNFRSF17", "XBP1", "MZB1")),
                  top2a = pg(pct, "TOP2A"), rag = pg(pct, "RAG1"), tbx21 = pg(pct, "TBX21"), taci = pg(pct, "TNFRSF13B"), sox4 = pg(pct, "SOX4"), vpreb3 = pg(pct, "VPREB3"),
                  myc = pscore(pct, c("NME1", "PSAT1", "SHMT2", "WDR77")), act = pscore(pct, c("CD69", "CD83", "NR4A1", "EGR1", "MYC")),
                  switched = switched, genes = as.vector(tapply(o$nFeature_RNA, cl, median)))
  s$myc_z <- rz(s$myc)
  s$act_z <- rz(s$act)
  s$detail <- with(s, case_when(
    pc >= .3 ~ "Plasmablasts", cd3 >= .3 | mye >= .25 ~ "Doublets", rag >= .3 ~ "Pre-B",
    tbx21 >= .3 ~ "Exhausted B", taci >= .5 & switched >= .3 ~ "Switched memory B", taci >= .5 ~ "Unswitched memory B", sox4 >= .3 & vpreb3 >= .4 ~ "Immature B",
    act_z >= 2.5 ~ "Activated B", myc_z >= 2 ~ "Immature B",
    taci >= .2 ~ "Unswitched memory B", genes < .7 * median(o$nFeature_RNA) ~ "Damaged B", TRUE ~ "Naive B"))
  data.frame(cluster = colnames(pct), n_cells = as.vector(table(cl)), detail = s$detail, round(s[, names(s) != "detail"], 3))
}

t_rules <- function(o) {
  cl <- droplevels(factor(o$seurat_clusters))
  pct <- pct_expr(o, cl, gfind(o, c(unlist(t_dot), "CSF1R")))
  s <- data.frame(cd3 = pg(pct, "CD3E"), ms4a1 = pg(pct, "MS4A1"), mye = (pg(pct, "CD14") + pg(pct, "CSF1R")) / 2, ncr1 = pg(pct, "NCR1"), trdc = pg(pct, "TRDC"),
                  cd8 = as.numeric(pg(pct, "CD8B") >= .3 | (pg(pct, "CD8A") >= .5 & pg(pct, "CD4") < .5)), foxp3 = pg(pct, "FOXP3"), cxcr5 = pg(pct, "CXCR5"), pdcd1 = pg(pct, "PDCD1"),
                  icos = pg(pct, "ICOS"), itgb1 = pg(pct, "ITGB1"), top2a = pg(pct, "TOP2A"), naive = pscore(pct, c("CCR7", "SELL", "LEF1", "TCF7")),
                  mem = pscore(pct, c("ICOS", "CTLA4", "S100A4", "ANXA1", "PDCD1", "CD40LG")), cyto = pscore(pct, c("GZMK", "GZMB", "GZMA", "PRF1", "NKG7", "CCL5")))
  s$T_subtype <- with(s, case_when(
    ms4a1 >= .3 | mye >= .25 ~ "Doublets", top2a >= .15 ~ "Proliferating T/NK", cd3 < .3 & ncr1 >= .3 ~ "NK cells", cd3 < .7 ~ "Damaged T/NK", trdc >= .4 ~ "Gamma delta T",
    cd8 > 0 & (cyto >= .4 | rz(cyto) >= 1.5) ~ "CD8 T effector memory", cd8 > 0 ~ "CD8 T naive", foxp3 >= .15 ~ "Treg", cxcr5 >= .15 & (pdcd1 >= .15 | icos >= .25) ~ "Tfh-like CD4 T",
    icos >= .25 | rz(itgb1) >= 1 | rz(mem - naive) >= 1 ~ "CD4 T memory", TRUE ~ "CD4 T naive"))
  data.frame(cluster = colnames(pct), n_cells = as.vector(table(cl)), T_subtype = s$T_subtype, round(s[, names(s) != "T_subtype"], 3))
}

dotp <- function(o, l, group.by) {
  seen <- character()
  l2 <- list()
  for (n in names(l)) {
    g <- setdiff(gfind(o, l[[n]]), seen)
    if (length(g)) l2[[n]] <- g
    seen <- c(seen, g)
  }
  DotPlot(o, features = l2, group.by = group.by, cols = c("grey90", "#B2182B")) + rot + labs(x = NULL, y = NULL) + theme(strip.text.x = element_text(angle = 90, hjust = 0, size = 8), panel.spacing = unit(2, "pt"))
}

sub_umap <- function(o, cells, x, y) {
  s <- subset(o, cells = cells)
  e <- as.matrix(s@meta.data[, c(x, y)])
  colnames(e) <- c("subumap_1", "subumap_2")
  s[["umap"]] <- CreateDimReducObject(e, key = "subumap_", assay = DefaultAssay(s))
  s
}

count_frac <- function(md, var) {
  md |> mutate(g = droplevels(factor(.data[[var]]))) |> dplyr::count(sample, g, .drop = FALSE) |> group_by(sample) |> mutate(frac = n / sum(n)) |> ungroup() |> dplyr::rename(!!var := g) |> mutate(sample = as.character(sample)) |> left_join(ss, by = "sample") |> mutate(timepoint = factor(timepoint, tp_levels))
}

lab_bars <- function(cnt, var, cols, min_frac = .03) {
  txt_on <- function(hex) ifelse(colSums(col2rgb(hex) * c(.299, .587, .114)) > 150, "black", "white")
  if (!"label" %in% names(cnt)) cnt$label <- labs_of(cnt$sample)
  cnt <- cnt |> arrange(label, .data[[var]]) |> group_by(label) |> mutate(y_n = sum(n) - cumsum(n) + n / 2, y_p = 100 * (sum(frac) - cumsum(frac) + frac / 2)) |> ungroup()
  cnt <- cnt |> mutate(lab_n = ifelse(frac >= min_frac, n, ""), lab_p = ifelse(frac >= min_frac, paste0(round(100 * frac), "%"), ""), txt = txt_on(cols[as.character(.data[[var]])]))
  p_n <- ggplot(cnt, aes(label, n, fill = .data[[var]])) + geom_col(colour = "black", linewidth = .1) + geom_text(aes(y = y_n, label = lab_n, colour = txt), size = 3.4, fontface = "bold", show.legend = FALSE) + scale_colour_identity() + scale_fill_manual(values = cols) + rot + labs(x = NULL, y = "Number of cells", fill = NULL) + theme(legend.position = "none")
  p_p <- ggplot(cnt, aes(label, 100 * frac, fill = .data[[var]])) + geom_col(colour = "black", linewidth = .1) + geom_text(aes(y = y_p, label = lab_p, colour = txt), size = 3.4, fontface = "bold", show.legend = FALSE) + scale_colour_identity() + scale_fill_manual(values = cols) + rot + labs(x = NULL, y = "Percentage", fill = NULL)
  p_n | p_p
}

read_contigs <- function(ch) {
  truthy <- function(x) toupper(x) %in% c("TRUE", "T")
  one <- function(i) {
    f <- vdj_path(ss$sample[i], ch)
    if (!file.exists(f)) return(NULL)
    read.csv(f, colClasses = "character") |> filter(truthy(is_cell), truthy(high_confidence), truthy(productive)) |> mutate(umis = as.numeric(umis), cell = paste0(ss$sample[i], "_", barcode), animal = ss$animal[i])
  }
  bind_rows(lapply(seq_len(nrow(ss)), one))
}

best_chain <- function(a, chains, prefix) {
  a |> filter(chain %in% chains) |> group_by(cell) |> slice_max(umis, n = 1, with_ties = FALSE) |> ungroup() |> select(cell, animal, chain, contig_id, v_gene, j_gene, c_gene, cdr3, umis) |> rename_with(~ paste0(prefix, .x), -c(cell, animal))
}

clone_key <- function(animal, v, j, cdr3) {
  grp <- ifelse(is.na(cdr3), NA, paste(animal, v, j, nchar(cdr3), sep = "|"))
  out <- rep(NA_character_, length(cdr3))
  for (g in unique(na.omit(grp))) {
    i <- which(grp == g)
    u <- unique(cdr3[i])
    cl <- if (length(u) == 1) 1 else cutree(hclust(as.dist(outer(seq_along(u), seq_along(u), Vectorize(function(x, y) sum(strsplit(u[x], "")[[1]] != strsplit(u[y], "")[[1]]))) / nchar(u[1])), "single"), h = .15)
    out[i] <- paste(g, cl[match(cdr3[i], u)], sep = "|")
  }
  out
}
