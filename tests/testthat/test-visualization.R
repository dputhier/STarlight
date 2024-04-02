library(stcompr)
set_verb_level(0)
library(testthat)

example_dataset()
test_data <- Xenium_Mouse_Brain_Coronal_7g

test_that("spatial_image...", {
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

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 2, saturation = 0.8, scale = F)$data$value), 2)), 2.81)

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 10,  scale = F)$data$value), 2)), 2.26)

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 10, saturation = 0.8, scale = F)$data$value), 2)), 0.85)

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 10, saturation = 0.8, scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 2, saturation = 0.8, scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 2, saturation = 0.8, scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(spatial_image(test_data, feat=feat_names(test_data), logb = 10,  scale = T)$data$value), 2)), 1)

})

context("spatial_plot function")


example_dataset()
test_data <- Xenium_Mouse_Brain_Coronal_7g
test_colors <- c("red", "blue", "green")

test_that("spatial_plot...", {

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

example_dataset("10819270/files/cmp_xen")
test_data <- cmp_xen
feat <- feat_names(cmp_xen)

test_that("cmp_bar_plot...", {

  expect_no_error(cmp_bar_plot(test_data, feat=feat))

  expect_error(cmp_bar_plot(test_data, feat=NULL), "Please provide some features")

  expect_equal(round(max(cmp_bar_plot(test_data, feat=feat, transform="log10")$data$Counts), 2), 4.47)

  expect_equal(round(max(cmp_bar_plot(test_data, feat=feat, transform="log2")$data$Counts), 2), 14.85)

  expect_equal(round(max(cmp_bar_plot(test_data, feat=feat, transform="None")$data$Counts), 2), 29602)

  expect_equal(round(max(cmp_bar_plot(test_data, feat=feat, transform="log")$data$Counts), 2), 10.3)

})


example_dataset("10819270/files/cmp_xen")
test_data <- cmp_xen
feat <- feat_names(cmp_xen)

test_that("cmp_boxplot...", {

  expect_no_error(cmp_boxplot(test_data))

  expect_no_error(cmp_boxplot(test_data, normalized = TRUE))

  expect_no_error(cmp_boxplot(test_data, normalized = FALSE))

  expect_equal(max(cmp_boxplot(test_data, normalized = FALSE)$data$Counts), 29602)

  expect_equal(round(max(cmp_boxplot(test_data, normalized = TRUE)$data$Counts), 2), 11046.53)

  expect_equal(round(max(cmp_bar_plot(test_data, feat=feat, transform="log10")$data$Counts), 2), 4.47)

  expect_equal(round(max(cmp_bar_plot(test_data, feat=feat, transform="log2")$data$Counts), 2), 14.85)

  expect_equal(round(max(cmp_bar_plot(test_data, feat=feat, transform="None")$data$Counts), 2), 29602)

  expect_equal(round(max(cmp_bar_plot(test_data, feat=feat, transform="log")$data$Counts), 2), 10.3)

})

example_dataset("10819270/files/cmp_xen")
test_data <- cmp_xen
feat <- feat_names(cmp_xen)

test_that("cmp_volcano...", {

  expect_equal(round(sum(cmp_volcano(cmp_xen)$data$x), 2), -1.56)

  expect_equal(round(sum(-log10(cmp_volcano(cmp_xen)$data$y)), 0), 971)

  expect_equal(round(sum(cmp_volcano(cmp_xen)$data$mean_counts), 2), 21247.58)

  expect_equal(round(sum(cmp_volcano(cmp_xen)$data$mean_counts), 2), 21247.58)

  expect_equal(round(sum(-log10(cmp_volcano(cmp_xen, y_axis="holm")$data$y)), 0), 968)

  expect_equal(round(sum(-log10(cmp_volcano(cmp_xen, y_axis="BH")$data$y)), 0), 969)

  expect_no_error(cmp_volcano(cmp_xen, y_axis="BH", x_lim = c(-1, 1)))

  expect_no_error(cmp_volcano(cmp_xen, color=c("red", "blue")))
})

example_dataset()
test_data <- Xenium_Mouse_Brain_Coronal_7g

test_that("plot_rip_k", {

  test_data <-  compute_k_ripley(test_data)

  expect_no_error(plot_rip_k(test_data))

  expect_equal(round(log10(sum(plot_rip_k(test_data)$data$border+1)), 2), 9.32)

  expect_equal(round(log10(sum(plot_rip_k(test_data)$data$iso+1)), 2), 9.31)

  expect_no_error(plot_rip_k(test_data, max_feat_label = 1, color = "red"))
})


example_dataset()
test_data <- Xenium_Mouse_Brain_Coronal_7g
t_1 <- test_data[bin_x(test_data)[180:200], bin_x(test_data)[100:140],]
t_2 <- test_data[bin_x(test_data)[101:141], bin_x(test_data)[61:101],]
feat <- feat_names(test_data)

test_that("cmp_images...", {

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

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 2, saturation = 0.8, scale = F)$data$value), 2)), 2.81)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 10,  scale = F)$data$value), 2)), 2.26)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 10, saturation = 0.8, scale = F)$data$value), 2)), 0.85)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 10, saturation = 0.8, scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 2, saturation = 0.8, scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 2, saturation = 0.8, scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = 10,  scale = T)$data$value), 2)), 1)

  expect_equal(max(round(max(cmp_images(test_data, feat_list=feat_names(test_data), logb = NULL,  scale = F)$data$value), 2)), 182)

  expect_no_error(cmp_images(t_1, t_2, feat_list = feat_names(t_1)))

})

