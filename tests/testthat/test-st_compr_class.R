library(testthat)

test_that("STCompR...", {

  library(STarlight)
  set_verb_level(0)
  example_dataset("11284296/files/cmp_xen")
  test_data <- cmp_xen
  set.seed(123)
  m <- matrix(sample(1:50, 100, rep=TRUE), nc=4)


  #p <- heatmap_cmp(cmp_xen)
  #expect_equal(sum(dim(p[[1]]$data)), 52)
  #expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 37)

  #p <- heatmap_cmp(cmp_xen, what = "c1")
  #expect_equal(sum(dim(p[[1]]$data)), 52)
  #expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 36.2)

  #p <- heatmap_cmp(cmp_xen, what = "c2")
  #expect_equal(sum(dim(p[[1]]$data)), 52)
  #expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 73.2)

  #p <- heatmap_cmp(cmp_xen, only_feat = c("Necab2", "Chat") )
  #expect_equal(sum(dim(p[[1]]$data)), 7)
  #expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 1.8)

  #p <- heatmap_cmp(cmp_xen, only_feat = c("Necab2", "Chat"),  what = "c1")
  #expect_equal(sum(dim(p[[1]]$data)), 7)
  #expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 2)

  #p <- heatmap_cmp(cmp_xen, only_feat = c("Necab2", "Chat"),  what = "c2")
  #expect_equal(sum(dim(p[[1]]$data)), 7)
  #expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 3.8)

  #p <- heatmap_cmp(cmp_xen, del_feat = c("Necab2", "Chat") )
  #expect_equal(sum(dim(p[[1]]$data)), 28)
  #expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 16.8)

  #p <- heatmap_cmp(cmp_xen, del_feat = c("Necab2", "Chat"),  what = "c1")
  #expect_equal(sum(dim(p[[1]]$data)), 28)
  #expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 17.6)

  #p <- heatmap_cmp(cmp_xen, del_feat = c("Necab2", "Chat"),  what = "c2")
  #expect_equal(sum(dim(p[[1]]$data)), 28)
  #expect_equal(round(sum(abs(p[[1]]$data$expression)), 1), 34.4)

  #expect_no_error(heatmap_cmp(cmp_xen, dist_method = "maximum"))
  #expect_no_error(heatmap_cmp(cmp_xen, dist_method = "manhattan"))
  #expect_no_error(heatmap_cmp(cmp_xen, dist_method = "canberra"))

  #expect_no_error(heatmap_cmp(cmp_xen, hclust_method = "ward.D2"))
  #expect_no_error(heatmap_cmp(cmp_xen, hclust_method = "single"))
  #expect_no_error(heatmap_cmp(cmp_xen, hclust_method = "complete"))


  expect_equal(nb_feat(cmp_xen), length(feat_names(cmp_xen)))
  expect_equal(nb_feat(cmp_xen), 7)


  expect_equal(round(sum(abs(stat_test(cmp_xen)$log2_ratio))), 7)

  example_dataset()
  xen <- Xenium_Mouse_Brain_Coronal_7g
  x_bins <-  bin_x(xen)[181:nbin_x(xen)]
  y_bins <-  bin_y(xen)[101:nbin_y(xen)]
  xen_r1.1 <- xen[x_bins, y_bins]
  x_bins <-  bin_x(xen)[61:101]
  y_bins <-  bin_y(xen)[101:nbin_y(xen)]
  xen_r2.1 <- xen[x_bins, y_bins]

  x_bins <-  bin_x(xen)[40:80]
  y_bins <-  bin_y(xen)[50:90]
  xen_r2.2 <- xen[x_bins, y_bins]

  expect_no_error(stcompr(list(xen_r1.1=xen_r1.1), list(xen_r2.1=xen_r2.1)))
  expect_equal(round(sum(-log10(stat_test(stcompr(list(xen_r1.1=xen_r1.1), list(xen_r2.1=xen_r2.1)))$padj))), 979)
  expect_no_error(stcompr(list(xen_r1.1=xen_r1.1), list(xen_r2.1=xen_r2.1, xen_r2.2=xen_r2.2), fit_type="mean"))
  expect_no_error(stcompr(list(xen_r1.1=xen_r1.1), list(xen_r2.1=xen_r2.1, xen_r2.2=xen_r2.2), fit_type="mean", p_adj_method = "holm"))

  expect_no_error(stcompr(xen_r1.1, xen_r2.1))
  expect_no_error(stat_test(cmp_xen))
  expect_equal(round(sum(stat_test(cmp_xen)$xen_r1), 0), 24418)
  expect_equal(nrow(stat_test(cmp_xen)), 7)
  expect_equal(nrow(stat_test(cmp_xen, melted_count = TRUE)), 7)

})
