library(testthat)
test_that("STGrid...", {

  library(stcompr)
  set_verb_level(0)


  example_dataset()
  test_data <- Xenium_Mouse_Brain_Coronal_7g
  x_bins <-  bin_x(test_data)[181:nbin_x(test_data)]
  y_bins <-  bin_y(test_data)[101:nbin_y(test_data)]
  xen_1 <- test_data[x_bins, y_bins]

  expect_equal(sum(dim(test_data)), 360)
  expect_equal(ncol(test_data), 142)
  expect_equal(nrow(test_data), 218)
  expect_equal(length(col_names(test_data)), 218)
  expect_equal(length(row_names(test_data)), 142)
  expect_equal(length(bin_x(test_data)), 218)
  expect_equal(length(bin_y(test_data)), 142)
  expect_equal(nbin_x(test_data), 218)
  expect_equal(nbin_y(test_data), 142)
  expect_equal(sum(round(coord(test_data)$x, 0)), 584994259)
  expect_equal(sum(round(coord(test_data)$y, 0)), 405922647)
  expect_equal(length(unique(coord(test_data)$feature)), 7)
  expect_equal(sum(bin_mat(test_data, del_bin = TRUE)), 152958)
  expect_equal(ncol(bin_mat(test_data, del_bin = TRUE)), 7)
  expect_equal(ncol(bin_mat(test_data, del_bin = FALSE)), 9)
  expect_equal(nb_items(test_data), 152958)
  expect_equal(nb_feat(test_data), 7)
  expect_equal(feat_names(test_data)[1], "Ano1")
  expect_equal(feat_names(test_data)[7], "Nwd2")
  expect_equal(bin_size(test_data), 25)
  t_2 <- compute_k_ripley(test_data[c("Ano1", "Ebf3"),])
  expect_equal(sum(round(ripley_k_function(t_2)$border),0), 248528293)
  expect_equal(nrow(coord(test_data)), 152958)

  ## coord

  ## [,]
  expect_equal(ncol(test_data) == ncol(test_data["Ano1",]), TRUE)
  expect_equal(nrow(test_data) == nrow(test_data["Ano1",]), TRUE)
  expect_equal(nbin_x(test_data) == nbin_x(test_data["Ano1",]), TRUE)
  expect_equal(nbin_y(test_data) == nbin_y(test_data["Ano1",]), TRUE)
  expect_equal(nrow(coord(test_data[bin_x(test_data)[1:2],])), 1)
  expect_equal(nrow(coord(test_data[bin_x(test_data)[1:3],])), 2)
  expect_equal(nrow(coord(test_data["Ano1",])) < 152958, TRUE)
  expect_equal(nrow(coord(test_data[feat_names(test_data),])) == 152958, TRUE)
  expect_equal(nrow(coord(test_data[,])) == 152958, TRUE)
  expect_equal(nb_items(test_data[bin_x(test_data), bin_y(test_data)]), nb_items(test_data))
  expect_equal(unique(coord(xen_1[bin_x(xen_1)[1], bin_y(xen_1)[1]])$feature), "Necab2")
  expect_equal(feat_names(xen_1[bin_x(xen_1)[1], bin_y(xen_1)[1]]), "Necab2")

  ## compute_k_ripley

  expect_equal(sum(round(compute_k_ripley(xen_1["Ano1",], verbose = FALSE)@ripley_k_function$border, 0)), 57780387)
  expect_equal(sum(round(compute_k_ripley(xen_1["Ano1",], verbose = FALSE)@ripley_k_function$iso, 0)), 51996873)

  ## load_spatial
  fp <- file.path(system.file("extdata", package = "stcompr"), "tyni_xenium.txt")
  expect_no_error(load_spatial(fp, method = "coordinates"))

})
