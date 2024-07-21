# Changelog

## v0.0.7

* "Added function gene_contrast and multi_gene_contrast".

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
