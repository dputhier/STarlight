# Changelog

* Added a report for a single sample (i.e. STGrid object) that can be called with reporting().
* plot_rip_k() can now be used with theorical values for a single feature.
* Changed argument return_tree to return_list in hc_tree().
* Added computation of L-function (normalized version of K-function).
* Implemented moran_index(), order_feat_by_moran_index(), plot_moran_i(), compute_moran_index().
* Moran's Index is now stored in the STGrid objects.
* in hc_tree() argument return_tree is now return_list
* In cmp_counts_st() change axis names from gene to feature

## v0.1.0

* multi_gene_contrast() has been renamed multi_feat_contrast()
* gene contrast() has been renamed feature_contrast().

## v0.0.8

* Added control_list() function.
* Added meta() method for STGrid.
* Added as_matrix() method for STGrid.
* Added connected_components(), cc_neighborhood(), create_hull(), find_contiguous().
* Added as_factor arg to spatial_image()
* Added axis_as_number arg to bin_mat()

## v0.0.7

* "Added function gene_contrast and multi_feat_contrast".

## v0.0.4

* Added new filtering method (cv) for heatmap_cmp(). This is the default now.


## v0.0.3

* Significant improvement of load_spatial() performances.

## v0.0.2

* load_spatial() now accept several constrains

## Wed May  1 15:46:15 CEST 2024

* Add meta_names() (to get access to meta_features)
* cmp_images() and plot_spatial now support meta_names.
* .DollarNames() is now exported and should allow completion using '$' on a STGrid object.
* Added compute_module_score() (see help files)
* Added dist_st() to draw count feature distributions of a set of STGrid object
* load_spatial has now a 'constrain' argument.
* "ward.D" is default argument for hc_tree().
