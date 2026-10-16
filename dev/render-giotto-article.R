# Precompute the outputs shown in vignettes/articles/giotto.Rmd.
#
# CI does not run that article: it needs Giotto (>= 4.3.2), GiottoVisuals
# (>= 0.2.16, for spatDeconvPlot() with numeric cell type names), GiottoData
# and Giotto's environment. This script runs the article's code locally and
# saves what the article displays: printed tables as .txt and figures as
# .png, under vignettes/articles/figures/. Re-run it, and commit the outputs,
# when the Giotto side of the workflow changes.
#
#   Rscript dev/render-giotto-article.R
#
# Environment variables (optional), to render against unreleased code:
#   GIOTTO_PATH         a Giotto checkout to pkgload::load_all() instead of the
#                       installed Giotto
#   GIOTTOVISUALS_PATH  likewise for GiottoVisuals
suppressMessages(library(arrow))  # before Giotto: see GiottoDisk notes
load_or_lib <- function(pkg, path_var) {
    path <- Sys.getenv(path_var)
    if (nzchar(path)) {
        suppressMessages(pkgload::load_all(path, quiet = TRUE))
    } else {
        suppressMessages(library(pkg, character.only = TRUE))
    }
}
load_or_lib("GiottoVisuals", "GIOTTOVISUALS_PATH")
load_or_lib("Giotto", "GIOTTO_PATH")
library(spatialDWLS)

out <- "vignettes/articles/figures"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
save_txt <- function(x, name) {
    writeLines(utils::capture.output(print(x)), file.path(out, name))
}

options(giotto.use_conda = FALSE, giotto.verbose = FALSE,
        giotto.no_python_warn = TRUE)

## --- the article's code, chunk by chunk ---------------------------------

g <- GiottoData::loadGiottoMini("visium", python_path = NA)

expr <- getExpression(g, values = "normalized", output = "matrix")
meta <- pDataDT(g)
markers <- findMarkers_one_vs_all(g, method = "scran",
    cluster_column = "leiden_clus", min_feats = 20)
sign_gene <- unique(markers[, head(feats, 20), by = "cluster"]$V1)

sign_matrix <- makeSignMatrixDWLS(g,
    sign_gene = sign_gene,
    cell_type_vector = paste0("type_", meta$leiden_clus))
save_txt(dim(sign_matrix), "sign-matrix-dim.txt")
save_txt(round(sign_matrix[1:5, 1:5], 2), "sign-matrix-head.txt")

g <- runDWLSDeconv(g, sign_matrix = sign_matrix,
    cluster_column = "leiden_clus", n_cell = 10)

dwls <- getSpatialEnrichment(g, name = "DWLS", output = "data.table")
save_txt(head(dwls[, c("cell_ID", sort(setdiff(names(dwls), "cell_ID"))),
    with = FALSE], 4), "dwls-head.txt")

type_cols <- sort(setdiff(names(dwls), "cell_ID"))
palette <- stats::setNames(
    c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02",
      "#a6761d", "#666666")[seq_along(type_cols)], type_cols)

p <- spatDeconvPlot(g, deconv_name = "DWLS", pie_scale = 0.4,
    cell_color_code = palette,
    show_plot = FALSE, return_plot = TRUE, save_plot = FALSE)
ggplot2::ggsave(file.path(out, "giotto-deconv-pies.png"), p,
    width = 7, height = 5.5, dpi = 110)

dominant <- type_cols[max.col(as.matrix(dwls[, type_cols, with = FALSE]))]
g <- addCellMetadata(g, new_metadata = data.frame(cell_ID = dwls$cell_ID,
    dominant_type = dominant), by_column = TRUE, column_cell_ID = "cell_ID")
p2 <- spatPlot2D(g, cell_color = "dominant_type", point_size = 2.5,
    cell_color_code = palette,
    show_plot = FALSE, return_plot = TRUE, save_plot = FALSE)
ggplot2::ggsave(file.path(out, "giotto-dominant-type.png"), p2,
    width = 6.5, height = 5, dpi = 110)

message("wrote: ", paste(list.files(out), collapse = ", "))
