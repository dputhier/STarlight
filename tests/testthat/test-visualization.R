library(testthat)
test_that("spatial_image...", {


  library(STarlight)
  set_verb_level(0)
  example_dataset()
  test_data <- Xenium_Mouse_Brain_Coronal_7g

  # Test that spatial_image does not throw an error
  expect_no_error(spatial_image(test_data, feat=feat_names(test_data)[1]))

  # Test that spatial_image creates a ggplot object
  expect_true("ggplot" %in% class(spatial_image(test_data, feat=feat_names(test_data)[1:2])))

  # Test case: Check if features argument is NULL
  expect_error(spatial_image(test_data, features = NULL), "Please provide a feature name")

  # Test case: Check if specified features exist in the object
  expect_error(spatial_image(test_data, features = c("NonExistentFeature")), "The feature was not found")

  # Test case: Check if saturation argument is within range
  expect_error(spatial_image(test_data, feat=feat_names(test_data)[1], saturation = 2), "Saturation should be between 0 and 1")

  # Test case: Check if overlay_feature exists in the object
  expect_error(spatial_image(test_data, overlay_feature = "NonExistentFeature"), "The feature to overlay was not found in the object")

  # Test case: Check if grid_by argument is char
  expect_error(spatial_image(test_data, grid_by = "bla", feat=feat_names(test_data)[1]), "grid_by  should be an integer")

  # Test case: Check if color_grid argument is NULL
  expect_no_error(spatial_image(test_data, color_grid = NULL, feat=feat_names(test_data)[3]))

  # Test case: Check if color_grid argument is int
  expect_no_error(spatial_image(test_data, color_grid = 20, feat=feat_names(test_data)[3]))

  expect_no_error(spatial_image(test_data, feat=feat_names(test_data)[3], color_grid = "white"))

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 2, scale = F)$data$value), 2)), 7.52)

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 2, saturation = 0.8, scale = F)$data$value), 2)), 3)

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 10,  scale = F)$data$value), 2)), 2.26)

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 10, saturation = 0.8, scale = F)$data$value), 2)), 0.9)

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 10, saturation = 0.8, scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 2, saturation = 0.8, scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 2, saturation = 0.8, scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 10,  scale = T)$data$value), 2)), 1)

})






test_that("spatial_plot...", {

  library(STarlight)
  set_verb_level(0)
  example_dataset()
  test_data <- Xenium_Mouse_Brain_Coronal_7g
  test_colors <- c("red", "blue", "green")

  # Test case 1: Check if object is provided
  expect_error(spatial_plot())

  # Test case 2: Check if feat_list is provided
  expect_error(spatial_plot(test_data), "Please provide feature names")

  # Test case 3: Check if feat_list contains valid feature names
  expect_error(spatial_plot(test_data, feat_list = c("InvalidFeature")), "One or several features were not found")

  # Test case 4: Check if size argument works
  expect_no_error(spatial_plot(test_data, feat_list = feat_names(test_data), size = 0.05))

  # Test case 5: Check if coord_fixed argument works
  expect_no_error(spatial_plot(test_data, feat_list = feat_names(test_data), coord_fixed = FALSE))

  # Test case 6: Check if colors argument works
  expect_error(spatial_plot(test_data, feat_list = feat_names(test_data), colors = test_colors), "More colors are needed")

  # Test case 7: Check if x-axis label is correct
  expect_equal(spatial_plot(test_data, feat_list = feat_names(test_data))$labels$x, "x")

  # Test case 9: Check if y-axis label is correct
  expect_equal(spatial_plot(test_data, feat_list = feat_names(test_data))$labels$y, "y")

  # Test case 11: Check if axis text size is correct
  expect_equal(spatial_plot(test_data, feat_list = feat_names(test_data))$theme$axis.text$size, 6)

  # Test case 12: Check if scale_color_manual is added when colors argument is provided
  expect_error(spatial_plot(test_data, feat_list = feat_names(test_data), colors = test_colors), "More colors are needed")

  # Test case 13: Check if ggplot object is returned
  expect_true("ggplot" %in% class(spatial_plot(test_data[c("Chat", "Ano1"), ], feat_list = c("Chat", "Ano1"))))

})



test_that("cmp_bar_plot...", {

  library(STarlight)
  set_verb_level(0)
  example_dataset("11284296/files/cmp_xen")
  test_data <- cmp_xen
  feat <- feat_names(cmp_xen)

  expect_no_error(cmp_bar_plot(test_data, feat=feat))

  expect_error(cmp_bar_plot(test_data, feat=NULL), "Please provide some features")

  expect_equal(round(max(cmp_bar_plot(test_data, feat=feat, transform="log10")$data$Counts), 2), 4.47)

  test <- round(max(cmp_bar_plot(test_data, feat=feat, transform="log2")$data$Counts), 2)
  expect_true( test >= 14.80 & test < 14.90)

  expect_equal(round(max(cmp_bar_plot(test_data, feat=feat, transform="None")$data$Counts), 2), 29602)

  test <- round(max(cmp_bar_plot(test_data, feat=feat, transform="log")$data$Counts), 1)
  expect_true(test == 10.3)

})




test_that("cmp_boxplot...", {

  library(STarlight)
  set_verb_level(0)
  example_dataset("11284296/files/cmp_xen")
  test_data <- cmp_xen
  feat <- feat_names(cmp_xen)

  expect_no_error(cmp_boxplot(test_data))

  expect_no_error(cmp_boxplot(test_data, normalized = TRUE))

  expect_no_error(cmp_boxplot(test_data, normalized = FALSE))

  test <- max(cmp_boxplot(test_data, normalized = FALSE)$data$Counts)
  expect_true(test > 29100 & test < 29700)

  test <- round(max(cmp_boxplot(test_data, normalized = TRUE)$data$Counts), 2)
  expect_true(test > 10800 & test < 11050)

  expect_equal(round(max(cmp_bar_plot(test_data, feat=feat, transform="log10")$data$Counts), 1), 4.5)

  expect_equal(round(max(cmp_bar_plot(test_data, feat=feat, transform="log2")$data$Counts), 1), 14.9)

  test <- round(max(cmp_bar_plot(test_data, feat=feat, transform="None")$data$Counts), 2)

  expect_true(test > 29100 & test < 29700)

  expect_equal(round(max(cmp_bar_plot(test_data, feat=feat, transform="log")$data$Counts), 1), 10.3)

})



test_that("cmp_volcano...", {

  library(STarlight)
  set_verb_level(0)
  example_dataset("11284296/files/cmp_xen")
  test_data <- cmp_xen
  feat <- feat_names(cmp_xen)

  test <- round(sum(cmp_volcano(cmp_xen)$data$x), 1)
  expect_true( test > -1.7 & test < -1.2)

  expect_true(round(sum(-log10(cmp_volcano(cmp_xen)$data$y)), 0) >= 970)

  test <- round(sum(cmp_volcano(cmp_xen)$data$mean_counts), 2)
  expect_true(test > 20800 & test < 21250)

  expect_true(round(sum(-log10(cmp_volcano(cmp_xen, y_axis="p_values")$data$y)), 0) >= 970)


  expect_no_error(cmp_volcano(cmp_xen, color=c("red", "blue"),
                              text_y_lim = 200,
                              text_x_lim=2))
})



test_that("plot_rip_k", {

  library(STarlight)
  set_verb_level(0)
  example_dataset()
  test_data <- Xenium_Mouse_Brain_Coronal_7g

  test_data <-  compute_k_ripley(test_data, verbose = FALSE)

  expect_no_error(plot_rip_k(test_data))

  expect_equal(round(log10(sum(plot_rip_k(test_data)$data$border+1)), 1), 8.7 )

  expect_equal(round(log10(sum(plot_rip_k(test_data)$data$iso+1)), 1), 8.6)

  expect_no_error(plot_rip_k(test_data, max_feat_label = 1, color = "red"))
})




test_that("cmp_images...", {

  library(STarlight)
  set_verb_level(0)
  example_dataset()
  test_data <- Xenium_Mouse_Brain_Coronal_7g
  t_1 <- test_data[bin_x(test_data)[180:200], bin_y(test_data)[100:140]]
  t_2 <- test_data[bin_x(test_data)[101:141], bin_y(test_data)[61:101]]
  feat <- feat_names(test_data)

  expect_no_error(cmp_images(test_data, feat_list = feat[1:2]))

  # Test that cmp_images does not throw an error
  expect_no_error(cmp_images(test_data, feat_list=feat_names(test_data)[1]))

  # Test that cmp_images creates a ggplot object
  expect_true("ggplot" %in% class(cmp_images(test_data, feat_list=feat_names(test_data)[1:2])))

  # Test case: Check if features argument is NULL
  expect_error(cmp_images(test_data, feat_list = NULL), "Please provide a feature list")

  # Test case: Check if specified features exist in the object
  expect_error(cmp_images(test_data, features = c("NonExistentFeature")), "Please provide a feature list")

  # Test case: Check if saturation argument is within range
  expect_error(cmp_images(test_data, feat_list=feat_names(test_data)[1], saturation = 2), "Saturation should be between 0 and 1")

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 2, scale = F)$data$value), 2)), 7.52)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 2, saturation = 0.8, scale = F)$data$value), 2)), 3)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 10,  scale = F)$data$value), 2)), 2.26)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 10, saturation = 0.8, scale = F)$data$value), 2)), 0.9)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 10, saturation = 0.8, scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 2, saturation = 0.8, scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 2, saturation = 0.8, scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 10,  scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = NULL,  scale = F)$data$value), 2)), 183)

  expect_no_error(cmp_images(t_1, t_2, feat_list = feat_names(t_1)))
  expect_no_error(cmp_images(t_1, t_2,
                             feat_list = feat_names(t_1),
                             color_x_strip = rainbow(length(feat_names(t_1))),
                             color_y_strip = rainbow(2),
                             color_strip_text_x="black",
                             color_strip_text_y="white")
                             )

  expect_error(cmp_images(t_1, t_2,
                             feat_list = feat_names(t_1),
                             color_x_strip = rainbow(2),
                             color_y_strip = rainbow(2),
                             color_strip_text_x="black",
                             color_strip_text_y="white"))

  expect_error(cmp_images(t_1, t_2,
                             feat_list = feat_names(t_1),
                             color_x_strip = rainbow(length(feat_names(t_1))),
                             color_y_strip = rainbow(3),
                             color_strip_text_x="black",
                             color_strip_text_y="white"))

})



test_that("dist_st", {
  library(STarlight)
  set_verb_level(0)
  example_dataset()

  xen <- Xenium_Mouse_Brain_Coronal_7g
  x_bins <-  bin_x(xen)[180:nbin_x(xen)]
  y_bins <-  bin_y(xen)[100:nbin_y(xen)]
  xen_r1 <- xen[x_bins, y_bins]
  x_bins <-  bin_x(xen)[60:100]
  y_bins <-  bin_y(xen)[100:nbin_y(xen)]
  xen_r2 <- xen[x_bins, y_bins]
  x_bins <-  bin_x(xen)[20:60]
  y_bins <-  bin_y(xen)[20:60]
  xen_r3 <- xen[x_bins, y_bins]
  expect_no_error(dist_st(xen_r1, xen_r2, xen_r3,
                        fill_color=c("red", "black", "green"),
                        type="density",
                      transform="log10"))
  expect_no_error(dist_st(xen_r1, xen_r2, xen_r3,
                          fill_color=c("red", "black", "green"),
                          type="hist",
                          transform="log10"))
  expect_no_error(dist_st(xen_r1, xen_r2, xen_r3,
                          fill_color=c("red", "black", "green"),
                          type="boxplot",
                          transform="log10"))
  expect_no_error(dist_st(xen_r1, xen_r2, xen_r3,
                          fill_color=c("red", "black", "green"),
                          type="boxjitter",
                          transform="log10"))
})




