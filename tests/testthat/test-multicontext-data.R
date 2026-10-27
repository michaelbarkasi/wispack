# Port of the data-only checks from tutorial_multicontext.Rmd (checks 1-2). These just filter/count cell types in
# the shipped CSV -- no model fit -- so they're fast. The model-fit checks (3-13) live in
# test-slow-multicontext.R since they require a real wisp(fit_only = TRUE) fit on this data.

g <- .wisp_golden$multicontext

test_that("ACx laminar + celltype count data matches the vignette narrative", {
  path <- skip_unless_extdata("ACx_laminar_countdata_demo_with_celltype.csv.gz")
  countdata <- read.csv(path)

  unique_celltypes <- unique(countdata$celltype_MMC)
  expect_between(length(unique_celltypes), bounds = g$n_celltypes_raw)

  Glut_mask <- grepl("Glut", countdata$celltype_MMC)
  GABA_mask <- grepl("GABA", countdata$celltype_MMC)
  countdata$celltype_MMC[Glut_mask] <- "Glut"
  countdata$celltype_MMC[GABA_mask] <- "GABA"
  countdata <- countdata[Glut_mask | GABA_mask, ]

  unique_celltypes <- unique(countdata$celltype_MMC)
  expect_between(length(unique_celltypes), bounds = g$n_celltypes_filtered)
})
