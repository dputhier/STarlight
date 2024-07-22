library(testthat)

test_that("STGrid...", {

  library(STarlight)
  set_verb_level(0)


  example_dataset()
  test_data <- Xenium_Mouse_Brain_Coronal_7g
  x_bins <-  bin_x(test_data)[181:nbin_x(test_data)]
  y_bins <-  bin_y(test_data)[101:nbin_y(test_data)]
  xen_1 <- test_data[x_bins, y_bins]

  expect_equal(sum(dim(test_data)), 358)
  expect_equal(ncol(test_data), 141)
  expect_equal(nrow(test_data), 217)
  expect_equal(length(col_names(test_data)), 217)
  expect_equal(length(row_names(test_data)), 141)
  expect_equal(length(bin_x(test_data)), 217)
  expect_equal(length(bin_y(test_data)), 141)
  expect_equal(nbin_x(test_data), 217)
  expect_equal(nbin_y(test_data), 141)
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
  expect_equal(sum(round(ripley_k_function(t_2)$border),0), 55633130)
  expect_equal(nrow(coord(test_data)), 152958)

  ## coord

  ## [,]
  expect_equal(ncol(test_data) == ncol(test_data["Ano1",]), TRUE)
  expect_equal(nrow(test_data) == nrow(test_data["Ano1",]), TRUE)
  expect_equal(nbin_x(test_data) == nbin_x(test_data["Ano1",]), TRUE)
  expect_equal(nbin_y(test_data) == nbin_y(test_data["Ano1",]), TRUE)
  expect_equal(nrow(coord(test_data[bin_x(test_data)[1:2],])), 2)
  expect_equal(nrow(coord(test_data[bin_x(test_data)[1:3],])),2)
  expect_equal(nrow(coord(test_data["Ano1",])) < 152958, TRUE)
  expect_equal(nrow(coord(test_data[feat_names(test_data),])) == 152958, TRUE)
  expect_equal(nrow(coord(test_data[,])) == 152958, TRUE)
  expect_equal(nb_items(test_data[bin_x(test_data), bin_y(test_data)]), nb_items(test_data))
  expect_equal(unique(coord(xen_1[bin_x(xen_1)[1], bin_y(xen_1)[1]])$feature), c("Necab2", "Nrp2"))
  expect_equal(feat_names(xen_1[bin_x(xen_1)[1], bin_y(xen_1)[1]]), c("Necab2", "Nrp2"))

  ## compute_k_ripley

  expect_equal(sum(round(compute_k_ripley(xen_1["Ano1",], verbose = FALSE)@ripley_k_function$border, 0)), 9655744)
  expect_equal(sum(round(compute_k_ripley(xen_1["Ano1",], verbose = FALSE)@ripley_k_function$iso, 0)), 9121974)

  ## load_spatial
  fp <- file.path(system.file("extdata", package = "STarlight"), "tyni_xenium.txt")
  expect_no_error(load_spatial(fp, method = "coordinates"))

  ## sum_over_contiguity_mat
  z <- structure(c(41L, 20L, 45L, 69L, 57L, 25L, 38L, 38L, 5L, 50L,
  5L, 92L, 87L, 67L, 18L, 69L, 43L, 18L, 86L, 80L, 91L, 82L, 85L,
  72L, 76L, 17L, 100L, 39L, 40L, 85L, 56L, 47L, 58L, 6L, 62L, 50L,
  7L, 14L, 4L, 46L, 4L, 97L, 2L, 100L, 74L, 10L, 76L, 57L, 96L,
  59L, 39L, 93L, 7L, 45L, 21L, 55L, 37L, 97L, 46L, 26L, 79L, 68L,
  44L, 88L, 40L, 17L, 43L, 55L, 29L, 56L, 7L, 21L, 5L, 4L, 94L,
  50L, 29L, 16L, 62L, 65L, 51L, 66L, 66L, 18L, 88L, 58L, 24L, 84L,
  56L, 13L, 88L, 36L, 86L, 47L, 55L, 56L, 91L, 46L, 73L, 61L), dim = c(10L,
                                                                      10L))
  u <- sum_over_contiguity_mat(z)

  expect_equal(dim(u)[1], 10)
  expect_equal(dim(u)[2], 10)

  expect_equal(u[1,1], 117)
  expect_equal(u[1,2], 326)
  expect_equal(u[1,3], 282)
  expect_equal(u[2,2], 456)
  expect_equal(u[2,9], 360)
  expect_equal(u[2,10], 357)
  expect_equal(u[6,1], 225)
  expect_equal(u[10,1], 171)
  expect_equal(u[10,10], 142)
  expect_equal(u[7,6], 410)

  u <- sum_over_contiguity_mat(z, include_center = TRUE)


  expect_equal(dim(u)[1], 10)
  expect_equal(dim(u)[2], 10)

  expect_equal(u[1,1], 117 + z[1,1])
  expect_equal(u[1,2], 326 + z[1,2])
  expect_equal(u[1,3], 282 + z[1,3])
  expect_equal(u[2,2], 456 + z[2,2])
  expect_equal(u[2,9], 360 + z[2,9])
  expect_equal(u[2,10], 357 + z[2,10])
  expect_equal(u[6,1], 225 + z[6,1])
  expect_equal(u[10,1], 171 + z[10,1])
  expect_equal(u[10,10], 142 + z[10,10])
  expect_equal(u[7,6], 410 + z[7,6])

  u <- sum_over_contiguity_mat(z, include_center = FALSE, ns=1, method = "rock")

  expect_equal(dim(u)[1], 10)
  expect_equal(dim(u)[2], 10)

  expect_equal(u[1,1], 25)
  expect_equal(u[1,2], 224)
  expect_equal(u[1,3], 143)
  expect_equal(u[2,2], 194)
  expect_equal(u[2,9], 174)
  expect_equal(u[2,10], 240)
  expect_equal(u[6,1], 164)
  expect_equal(u[10,1], 85)
  expect_equal(u[10,10], 86)
  expect_equal(u[7,6], 271)

  u <- sum_over_contiguity_mat(z, include_center = TRUE, ns=1, method = "rock")

  expect_equal(dim(u)[1], 10)
  expect_equal(dim(u)[2], 10)

  expect_equal(u[1,1], 25 + z[1,1])
  expect_equal(u[1,2], 224 + z[1,2])
  expect_equal(u[1,3], 143 + z[1,3])
  expect_equal(u[2,2], 194 + z[2,2])
  expect_equal(u[2,9], 174 + z[2,9])
  expect_equal(u[2,10], 240 + z[2,10])
  expect_equal(u[6,1], 164 + z[6,1])
  expect_equal(u[10,1], 85 + z[10,1])
  expect_equal(u[10,10], 86 + z[10,10])
  expect_equal(u[7,6], 271 + z[7,6])

  u <- sum_over_contiguity_mat(rbind(z,1), include_center = TRUE, ns=1, method = "rock")

  expect_equal(dim(u)[1], 11)
  expect_equal(dim(u)[2], 10)

  expect_equal(u[1,1], 25 + z[1,1])
  expect_equal(u[1,2], 224 + z[1,2])
  expect_equal(u[1,3], 143 + z[1,3])
  expect_equal(u[2,2], 194 + z[2,2])
  expect_equal(u[2,9], 174 + z[2,9])
  expect_equal(u[2,10], 240 + z[2,10])
  expect_equal(u[6,1], 164 + z[6,1])
  expect_equal(u[10,1], 86 + z[10,1])
  expect_equal(u[10,10], 87 + z[10,10])
  expect_equal(u[7,6], 271 + z[7,6])

  u <- sum_over_contiguity_mat(cbind(rbind(z,1), 1), include_center = TRUE, ns=1, method = "rock")

  expect_equal(dim(u)[1], 11)
  expect_equal(dim(u)[2], 11)

  expect_equal(u[1,1], 25 + z[1,1])
  expect_equal(u[1,2], 224 + z[1,2])
  expect_equal(u[1,3], 143 + z[1,3])
  expect_equal(u[2,2], 194 + z[2,2])
  expect_equal(u[2,9], 174 + z[2,9])
  expect_equal(u[2,10], 241 + z[2,10])
  expect_equal(u[6,1], 164 + z[6,1])
  expect_equal(u[10,1], 86 + z[10,1])
  expect_equal(u[10,10], 88 + z[10,10])
  expect_equal(u[7,6], 271 + z[7,6])

  u <- sum_over_contiguity_mat(z, ns=2)

  expect_equal(dim(u)[1], 10)
  expect_equal(dim(u)[2], 10)

  expect_equal(u[1,1], 507)
  expect_equal(u[1,2], 704)
  expect_equal(u[1,3], 721)
  expect_equal(u[2,2], 831)
  expect_equal(u[2,9], 708)
  expect_equal(u[2,10], 459)
  expect_equal(u[6,1], 721)
  expect_equal(u[10,1], 391)
  expect_equal(u[10,10], 415)
  expect_equal(u[7,6], 1104)

  u <- sum_over_contiguity_mat(z, ns=2, include_center = TRUE)

  expect_equal(dim(u)[1], 10)
  expect_equal(dim(u)[2], 10)

  expect_equal(u[1,1], 507 + z[1,1])
  expect_equal(u[1,2], 704 + z[1,2])
  expect_equal(u[1,3], 721 + z[1,3])
  expect_equal(u[2,2], 831 + z[2,2])
  expect_equal(u[2,9], 708 + z[2,9])
  expect_equal(u[2,10], 459 + z[2,10])
  expect_equal(u[6,1], 721 + z[6,1])
  expect_equal(u[10,1], 391 + z[10,1])
  expect_equal(u[10,10], 415 + z[10,10])
  expect_equal(u[7,6], 1104 + z[7,6])
})
