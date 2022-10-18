# GMOlabeler
Python/R code for cross-analyzing hyperspectral regression and (RGB) semantic segmentation results to calculate statistics for development of specific transgenic plant tissues

## Description
This repository contains R and Python scripts used in the final module of the [GMOdetector workflow](https://github.com/naglemi/GMOnotebook) for studying plant transformation with complementary RGB and fluorescent hyperspectral datasets. The scripts contained here are used to cross-reference semantic segmentation masks (with tissues segmented according to their type, e.g. callus and shoot) with weights from regression of fluorescent hyperspectral images over known components (e.g. chlorophyll and GFP) obtained from (CubeGLM)[https://github.com/naglemi/gmodetector_py]. Outputs produced include summary statistics of fluorescent protein signal in specific tissues, and p-values for the effects on these statistics by experimental treatment, genotype and interactions thereof.

These scripts are not intended for use outside of the [GMOdetector workflow](https://github.com/naglemi/GMOnotebook). This readme will provide a basic explanation of how these scripts are used in this workflow.

## Computation of summary statistics of fluorescent protein signal in specific tissues
The `main` function is used to integrate outputs from semantic segmentation and (CubeGLM)[https://github.com/naglemi/gmodetector_py]. Below is an example of how this function is used:

```
python main.py \
<complementary_image_paths.csv> \ # string
<my_aligned_grid.jpg> \ # string
<my_reporter_significance_threshold> \ # integer
<my_reporter_protein_ID> \ # string
<grid_type_integer> \ # integer, "12" or "20"; 12-explant or 20-explant grids supported
<output_format_string> \ #string, "hdf" or "csv"
<output_directory_path> \ # string
<my_grid_borders> # string
```
The parameters for this function are described in documentation for the [GMOdetector workflow](https://github.com/naglemi/GMOnotebook):
- [Preparing alignment for grid and finding grid borders](https://github.com/naglemi/GMOnotebook/tree/master/1_Decide_parameters/2_Align_and_crop_parameters)
- [Fluorescent protein settings](https://github.com/naglemi/GMOnotebook/blob/master/1_Decide_parameters/3_Other_parameters/3_Hyperspectral_settings.ipynb)
- Producing spreadsheet of corresponding RGB images, segmentation masks and hyperspectral weights ([Main workflow](https://github.com/naglemi/GMOnotebook/blob/master/2a_Deploy_workflow/GMOdetector_template_v0.62.ipynb), section "Prepare sample datasheet input")

## Plotting summary statistics and performing statistical tests for significance of effects
The `grid_item_plots.R` script is used to take outputs from the `main` function, produce plots and perform tests. Below is an example of basic use:

```
Rscript grid_item_plots.R \
-d <my_data_folder>" \ # string
-r <my_metadata_spreadsheet.xlsx> \ # string
--pixel_threshold <my_reporter_significance_threshold> \ # integer
-missing <if_some_explants_missing_or_contaminated> \ # logical
-MissingList <missing_explant_dataframe> \ # output from missing/contaminated explant script or user manual input
--grid_type <grid_type_integer> \ # integer, "12" or "20"; 12-explant or 20-explant grids supported
--samples-pre-labeling ${data}/samples_pre_labeling.csv \
--height <my_plot_height> \ # integer
--width <my_plot_width> \ # integer
--Reporter <my_reporter_protein_ID> \ # string
```
The additional parameters used here are also described in [GMOdetector workflow](https://github.com/naglemi/GMOnotebook) documentation.
- [Preparing metadata spreadsheet](https://github.com/naglemi/GMOnotebook/blob/master/1_Decide_parameters/1_Metadata_and_randomization/1-Generate_randomization_scheme.ipynb)
- [Missing/contaminated explant detection or manual user inputs](https://github.com/naglemi/GMOnotebook/blob/master/1_Decide_parameters/3_Other_parameters/2_Missing_or_contaminated_explants.ipynb)

## Acknowledgements
We thank the National Science Foundation Plant Genome Research Program for support (IOS #1546900, Analysis of genes affecting plant regeneration and transformation in poplar), and members of GREAT TREES Research Cooperative at OSU for its support of the Strauss laboratory.
