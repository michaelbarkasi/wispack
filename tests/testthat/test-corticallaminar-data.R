# Port of the data-only checks from tutorial_corticallaminar.Rmd (checks 1-3). These just read the shipped CSV and
# count/index into it -- no model fit -- so they're fast. The model-fit checks (4-5) live in
# test-slow-corticallaminar.R since they require a real (default-settings) wisp() fit on this data.

g <- .wisp_golden$corticallaminar

test_that("S1 laminar count data matches the vignette narrative", {
  path <- skip_unless_extdata("S1_laminar_countdata_demo.csv.gz")
  countdata <- read.csv(path)

  num_genes <- length(unique(countdata$gene))
  expect_between(num_genes, bounds = g$n_genes)
  expect_equal(unique(countdata$gene)[2], g$second_gene_name)

  num_cells <- nrow(countdata) / num_genes
  expect_between(num_cells, bounds = g$n_cells)
})
