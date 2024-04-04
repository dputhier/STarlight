library(testthat)

test_that("STCompR...", {

  library(stcompr)
  set_verb_level(0)
  example_dataset("10819270/files/cmp_xen")
  test_data <- cmp_xen
  set.seed(123)
  m <- matrix(sample(1:50, 100, rep=TRUE), nc=4)

  expect_equal(round(sum(estimSf(m)), 2), 4.05)
  p <- heatmap_cmp(cmp_xen)
  expect_equal(sum(dim(p[[1]]$data)), 52)
  expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 37)

  p <- heatmap_cmp(cmp_xen, what = "c1")
  expect_equal(sum(dim(p[[1]]$data)), 52)
  expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 36.2)

  p <- heatmap_cmp(cmp_xen, what = "c2")
  expect_equal(sum(dim(p[[1]]$data)), 52)
  expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 73.2)

  p <- heatmap_cmp(cmp_xen, only_feat = c("Necab2", "Chat") )
  expect_equal(sum(dim(p[[1]]$data)), 7)
  expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 1.8)

  p <- heatmap_cmp(cmp_xen, only_feat = c("Necab2", "Chat"),  what = "c1")
  expect_equal(sum(dim(p[[1]]$data)), 7)
  expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 2)

  p <- heatmap_cmp(cmp_xen, only_feat = c("Necab2", "Chat"),  what = "c2")
  expect_equal(sum(dim(p[[1]]$data)), 7)
  expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 3.8)

  p <- heatmap_cmp(cmp_xen, del_feat = c("Necab2", "Chat") )
  expect_equal(sum(dim(p[[1]]$data)), 28)
  expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 16.8)

  p <- heatmap_cmp(cmp_xen, del_feat = c("Necab2", "Chat"),  what = "c1")
  expect_equal(sum(dim(p[[1]]$data)), 28)
  expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 17.6)

  p <- heatmap_cmp(cmp_xen, del_feat = c("Necab2", "Chat"),  what = "c2")
  expect_equal(sum(dim(p[[1]]$data)), 28)
  expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 34.4)

  expect_no_error(heatmap_cmp(cmp_xen, dist_method = "maximum"))
  expect_no_error(heatmap_cmp(cmp_xen, dist_method = "manhattan"))
  expect_no_error(heatmap_cmp(cmp_xen, dist_method = "canberra"))

  expect_no_error(heatmap_cmp(cmp_xen, hclust_method = "ward.D2"))
  expect_no_error(heatmap_cmp(cmp_xen, hclust_method = "single"))
  expect_no_error(heatmap_cmp(cmp_xen, hclust_method = "complete"))


  expect_equal(nb_feat(cmp_xen), length(feat_names(cmp_xen)))
  expect_equal(nb_feat(cmp_xen), 7)

  expect_equal(round(sum(-log10(stat_test(cmp_xen)$p_values))), 971)
  expect_equal(round(sum(abs(stat_test(cmp_xen)$log2_ratio))), 7)


})
