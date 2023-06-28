# FA Analyzer

This package contains a single Python/Jython script that can aid in the analysis of focal adhesion stainings.

The script can be opened and run directly in FIJI.

To run this script, you will have to first install the [MorphoLibJ](https://imagej.net/plugins/morpholibj) plugin. The easiest way to do this is to simply subscribe to the IJPB-Plugins site via the [FIJI update sites](https://imagej.net/update-sites/).

To run this script, open `fa_analyzer.py` in FIJI. This will bring up an editor window showing the script. Then simply click **Run**. Upon running the script, you will be prompted to select a directory. This directory should contain all of the focal adhesion stained image you want to analyze (for example, all images from a particular condition/experiment). The FA image should be single channel, 8- or 16-bit, any image size is fine. For each focal adhesion image, you will need to have an additional binary mask image with the segmented cell you would like to analyze. This should be an 8-bit-image with only values of 0 and 255. To create such an image, you can, for example, do the following:

1. Open your FA image or another channel that allows you to better see the cell boundaries
2. Using the Freehand Selection Tool or Polygon Selection Tool, draw a line around the cell boundary
3. *Edit > Selection > Add To Manager*
4. Create a new image: *File > New > Image* (*Type:* 8-bit, *Fill with:* Black, *Width/Height:* same as FA image), **OK**
5. Select your manually-drawn cell boundary from the ROI manager
6. Ensure that your foreground color is white using the Color Picker Tool
7. *Edit > Draw*
8. Save using the same name as your FA image, but with "_mask" appended to the end (for "cell2.tif", the mask image should be named "cell2_mask.tif")

As long as your mask image is an 8-bit binary image, you can create it any way you would like, for example using an automated segmentation in a previous analysis step. Just keep in mind that the naming convention must remain (always use the "_mask" suffix). If you would like to change this naming convention or the file type (in case you are using something other than Tiff images), you are free to do so in `fa_analyzer.py` in the `main` function at the bottom of the script. For more details on how the FA image and corresponding mask should be named and what the images should look like, refer to the included sample data.

The script will automatically create new directories for storing data. These new directories be named according to your FA images. Within these directories, there will be csv files saved with data for each individual focal adhesion and data related to focal adhesions at the level of the whole cell ("image data").

The focal adhesion data includes the following parameters:

- basename : the name of the image file
- id : the unique identifier for the focal adhesion
- fa_area_px : the area of the focal adhesion (in pixels)
- perim_px : the perimeter of the focal adhesion (in pixels)
- x0_px : the x-coordinate of the focal adhesion centroid (in pixels)
- y0_px : the x-coordinate of the focal adhesion centroid (in pixels)
- orientation_deg : the orientation angle of the focal adhesion (from -90 to 90 degrees)
- ar : the aspect ratio of the focal adhesion (length of principal axis / normal axis)
- sf : the shape factor of the focal adhesion (another measure of elongation, defined as sf = perimeter^2 / area)
- dist_px : the minimal distance from the focal adhesion centroid to the cell border
- angle_deg : the angle of the focal adhesion centroid with respect to the cell centroid (in degrees)
- rel_orientation_deg : the relative orientation of the focal adhesion orientation with respect to the angle to the centroid (0 degrees means orientation is parallel with the angle to the centroid, 90 degrees means the orientation is perpendicular with the angle to the centroid)

The cell level ("image") data includes the following additional parameters:

- mask_area_px : the area of the cell mask (in pixels)
- no_fas : the total number of focal adhesions within the cell mask area
- area_fraction : the fraction of total focal adhesion area to cell area
- frac_periph_fas : the fraction of "peripheral" focal adhesions (by default, the distance of "peripheral" focal adhesions to the cell border is less than 10% of the maximum distance from the border to the center of the cell; this 10% fraction is defined in the ```main``` function at the end of the script and can be changed, but it should be kept consistent for all images that are to be compared)

In addition to these additional cell level data, the mean and standard deviations (std) of most of the focal adhesion parameters are also output to the image data csv file.
For convenience, all of the focal adhesion and image data are collected together in csv files that are saved to the directory containing the images.

The script will also generate a 16-bit "label" image showing each focal adhesion that was found and analyzed.
Note that the label (intensity value) for each focal adhesion corresponds to the focal adhesion "id" in the csv fa_data file.

Please note that the area and perimeter data are always given in pixels. This avoids errors due to  inconsistencies in images already having pixel sizes defined in the image metadata. For your analysis/data presentation, it is recommended to convert such values to real physical values (microns, e.g.). Also please be aware that if you are using different imaging parameters (zoom/objectives, e.g.), this conversion will be different for different images. For this reason, it is recommended to be consistent with your imaging parameters for a given data set that will be analyzed together.

For more information on how to incorporate these scripts into your analysis pipeline for focal adhesion analysis, please see our accompanying article:

van Stalborch A.-M.D.,  Clark A.G., Sonnenberg A. and Margadant C. (2023) Imaging and quantitative analysis of integrin-dependent cell-matrix adhesions. *STAR Protocols*.

If you found these scripts to be useful, please feel free to cite our publication.