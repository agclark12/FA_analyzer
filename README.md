# FA Analyzer

This package contains a single Python/Jython script that can aid in the analysis of focal adhesion stainings.

The script can be opened and run directly in FIJI.

To run this script, you will have to first install the [MorphoLibJ](https://imagej.net/plugins/morpholibj) plugin. The easiest way to do this is to simply subscribe to the IJPB-Plugins site via the [FIJI update sites](https://imagej.net/update-sites/).

To run this script, open ```fa_analyzer.py``` in FIJI. This will bring up an editor window showing the script. Then simply click **Run**. Upon running the script, you will be prompted to select a directory. This directory should contain all of the focal adhesion stained image you want to analyze (for example, all images from a particular condition/experiment). The FA image should be single channel, 8- or 16-bit, any image size is fine. For each focal adhesion image, you will need to have an additional binary mask image with the segmented cell you would like to analyze. This should be an 8-bit-image with only values of 0 and 255. To create such an image, you can, for example, do the following:

1. Open your FA image or another channel that allows you to better see the cell boundaries
2. Using the Freehand Selection Tool or Polygon Selection Tool, draw a line around the cell boundary
3. Edit -> Selection -> Add To Manager
4. Create a new image: File -> New -> Image (Type: 8-bit, Fill with: Black, Width/Height: same as FA image), **OK**
5. Select your manually-drawn cell boundary from the ROI manager
6. Ensure that your foreground color is white using the Color Picker Tool
7. Edit -> Draw
8. Save using the same name as your FA image, but with "_mask" appended to the end (for "cell2.tif", the mask image should be named "cell2_mask.tif")

As long as your mask image is an 8-bit binary image, you can create it any way you would like, for example using an automated segmentation in a previous analysis step. Just keep in mind that the naming convention must remain (anyways use the "_mask" suffix). If you would like to change this naming convention or the file type (in case you are using something other than Tiff images), you are free to do so in ```fa_analyzer.py``` in the ```main``` function at the bottom of the script. For more details on how the FA image and corresponding mask should be named and what the images should look like, refer to the included sample data.

For more information on how to incorporate these scripts into your analysis pipline for focal adhesion analysis, please see our accompanying article:

XXX.

If you found these scripts to be useful, please feel free to cite our publication.