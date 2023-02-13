"""

This is a jython script to analyze focal adhesion data from a directory of focal adhesion images with corresponding masks.

This script can be run directly in FIJI and requires installation of morpholibJ:

https://imagej.net/plugins/morpholibj

The easiest is to simply subscribe to the IJPB-Plugins site via the FIJI update sites.

When running ths script, you will first be prompted to select a directory containing images.

The only important modifiable parameter here is the fraction counted as 'peripheral' vs. 'central' adhesions.

This paramter is 'periph_frac', which should have a value 0-1 can be found in the main() function at the bottom.

"""

from __future__ import print_function, division

__author__ = "Andrew Clark"
__copyright__ = "Copyright 2022, FA Analyzer"

from exceptions import Exception as PythonException
from re import split as re_split
from os import listdir, sep, mkdir
from os.path import join, splitext, isdir
from collections import OrderedDict
from math import sqrt, atan2, pi, radians, degrees

from ij import IJ, ImagePlus
from ij.plugin import ImageCalculator
from ij.plugin.filter import MaximumFinder, ImageMath
from ij.measure import ResultsTable
from ij.gui import WaitForUserDialog, Overlay, TextRoi, MessageDialog
#from java.awt import Font, Color

from inra.ijpb.binary import BinaryImages
from inra.ijpb.label import LabelImages

def raise_error(message):

	window = IJ.showMessage("Error",message)
	raise PythonException(message)

def mean(ell):
	"""Returns the mean value from a list of numbers
	
	Parameters
    ----------
    ell : list of numbers
        the list over which you would like to calculate the mean
	Returns
    ----------
    mean : float
        the mean of the list	
	"""

	return float(sum(ell)) / float(len(ell))
	
def std(ell):
	"""Returns the standard deviation from a list of numbers
	
	Parameters
    ----------
    ell : list of numbers
        the list over which you would like to calculate the standard deviation
	Returns
    ----------
    std_dev : float
        the standard deviation of the list	
	"""
	
	ell_mean = mean(ell)
	dev = [None] * len(ell)
	for i in range(len(ell)):
		dev[i] = abs(float(ell[i]) - ell_mean)**2
	
	return sqrt(mean(dev))

def write_csv(array,path):
    """Writes a csv file from a 2D array
    
    Parameters
    ----------
    array : 2D list of lists
        the data you would like to save
    path : string
        path where the .csv file should be saved
    """

    ofile = open(path,'w')
    for i in range(len(array)):
        line = ','.join([str(_) for _ in array[i]]) + '\r\n'
        ofile.write(line)
    ofile.close()

def od_to_matrix(od):
	"""Converts an OrderedDictionary object to a 2D matrix (list of lists) with the keys as headers
	
	Parameters
    ----------
    od : OrderedDictionary
        the ordered dictionary you want to convert
	Returns
    ----------
    matrix : list (2D)
        a 2D list of lists with the keys as headers and data for each key in columns
	"""
	
	matrix = [[_] for _ in od.keys()]
	for i in range(len(matrix)):
		data_to_append = od[matrix[i][0]]
		
		if type(data_to_append) is list:
			matrix[i].extend(data_to_append)
		else:
			matrix[i].append(data_to_append)
		
	matrix = list(zip(*matrix))
	return matrix

def natural_sort(ell):
    """Returns a naturally sorted list (taking numbers into account)
    
    Parameters
    ----------
    ell : list
        the list to be sorted
    Returns
    -------
    ell_sorted : list
        the sorted list
    """

    convert = lambda text: int(text) if text.isdigit() else text.lower()
    key = lambda key: [convert(c) for c in re_split('([0-9]+)', key)]
    ell_sorted = sorted(ell, key = key)
    return ell_sorted

def get_unique_list(seq):
    """Returns a list sorted to have only unique values (in the same order as original list)
    
    Parameters
    ----------
    ell : list
        the list to be sorted
    Returns
    -------
    ell_sorted : list
        the sorted list
    """
    seen = set()
    seen_add = seen.add
    ell_sorted = [x for x in seq if not (x in seen or seen_add(x))]
    return ell_sorted

def concat_ods(od1,od2):
	"""Concatenates two ordered dictionaries according to the keys (preserves all keys from both dictionaries)
	Parameters
    ----------
    od1 : OrderedDictionary
        the OrderedDict to be appended to
    od2 : OrderedDictionary
        the OrderedDict to append
	Returns
    ----------
    od_full : OrderedDictionary
        the concatenated ordered dictionaries
	"""
	
	# gets keys and initializes a new OrderedDict
	keys = get_unique_list(od1.keys() + od2.keys())
	od_full = OrderedDict()
	
	# goes through all the keys
	for key in keys:
		if key in od1.keys():
			if type(od1[key]) is list:
				od_full[key] = od1[key]
			else:
				od_full[key] = [od1[key]]
		else:
			od_full[key] = []
		if key in od2.keys():
			if type(od2[key]) is list:
				od_full[key].extend(od2[key])
			else:
				od_full[key].append(od2[key])
	return od_full

def get_angular_diff(angle1, angle2):
	"""Gets the angular difference between two angles or orientations (in degrees)
	
    Parameters
    ----------
    angle1 : float
    	first angle (in degrees)
    angle2 : float
    	second angle (in degrees)
    Returns
    ----------
    diff : float
        angular difference (in degrees)
	"""

 	diff = abs(angle1 - angle2)
 	if diff > 180:
 		raise PythonException("Angular difference cannot exceed 180 degrees")
 	if diff >90:
 		diff = 180 - diff
	return diff

def angle_to_orientation(angle):
	"""Converts an angle (in degrees) to an orientation between -90 and 90 degrees
	
	Parameters
    ----------
    angle : float
    	an angle from -180 to 180
    Returns
    ----------
    new_angle : float
        orientation between -90 and 90 degrees 
	
	"""
	if not all((angle <= 180, angle >= -180)):
		raise_exception("Input angle must be between -180 and 180 degrees")

	if angle > 90:
		new_angle = angle - 180
	elif angle <= -90:
		new_angle = angle + 180
	else:
		new_angle = angle
	return new_angle

def analyze_image(file_path,mask_suffix,file_ending,periph_frac):
	"""Analyzes a focal adhesion and corresponding mask image. A label image and the FA and image data
	are saved automatically to automatically-generated directories with the same name as the image.
	
	Parameters
    ----------
    file_path : string
    	path to the focal adhesion image
    mask_suffix : string
    	suffix of the mask file with respect to the basename of the image file
    file_ending : string
    	the file ending (usually an image extension like ".tif")
    periph_frac : float
    	the distance fraction to be considered a peripheral adhesion (fraction of max distance of the mask from the border)
    Returns
    ----------
    fa_data : OrderedDict
        ordered dictionary containing the individual focal adhesion data
    image_data : OrderedDict
    	ordered dictionary containing the image-level data and means/sds for focal adhesion data
	
	"""

	print("Analyzing Image: " + file_path)

	# opens image and corresponding mask image
	img = IJ.openImage(file_path)
	mask_path = file_path.replace(file_ending, mask_suffix + file_ending)
	mask = IJ.openImage(mask_path)
	
	# checks some requirements that the FA image and mask image match
	if not mask.getProcessor().isBinary():
		raise_error("Mask image must be binary image (8-bit-image with only 0 and 255)")
	if not all((img.getNSlices() == 1, mask.getNSlices() == 1)):
		raise_error("FA image and mask must have only one slice")
	if not all((img.getWidth() == mask.getWidth(), img.getHeight() == mask.getHeight())):
		raise_error("FA image and mask must have same dimensions")
	
	# removes scale from image and mask
	IJ.run(img, "Properties...", "channels=1 slices=1 frames=1 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
	IJ.run(mask, "Properties...", "channels=1 slices=1 frames=1 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
	
	# creates distance map
	dist_map = BinaryImages.distanceMap(mask)
	
	# converts mask to a a [0,1] binary
	IJ.run(mask, "Divide...", "value=255")
	IJ.run(mask, "16-bit", "")
	
	# gets the centroid of the cell mask
	IJ.run(mask, "Analyze Regions", "pixel_count centroid");
	rt_mask = ResultsTable.getActiveTable()
	area_mask = float(rt_mask.getValue('PixelCount',0))
	x0_mask = float(rt_mask.getValue('Centroid.X',0))
	y0_mask = float(rt_mask.getValue('Centroid.Y',0))
	IJ.selectWindow(rt_mask.title)
	IJ.run("Close")
	
	# make sure image requirements are met
	##TODO: check sizes are same, only single slices, mask should be 8-bit binary, 0-255
			
	# makes a new directory for saving the data
	basename = splitext(file_path.split(sep)[-1])[0]
	save_dir = file_path.replace(file_ending,"")
	if not isdir(save_dir):
		mkdir(save_dir)
		
	# thresholds and segments the focal adhesions
	IJ.run(img, "Gaussian Blur...", "sigma=1")
	IJ.run(img, "Auto Threshold", "method=Otsu white")
	IJ.run(img, "Watershed", "")
	bin = ImageCalculator.run(img, mask, "Multiply create");
	bin_filt = BinaryImages.sizeOpening(bin, 10)
	labels = BinaryImages.componentsLabeling(bin_filt, 4, 16)
	
	# saves image of labeled focal adhesions
	IJ.run(labels, "glasbey", "")
	#labels.show()
	save_path = join(save_dir,basename+"_fa_labels.tif")
	IJ.saveAs(labels, "Tiff", save_path);
	
	# gets morphometric data for labels
	# (this could probably be done better by using the morpholibj API directly instead of the plugin, but it's much simpler with the plugin)
	IJ.run(labels, "Analyze Regions", "pixel_count perimeter equivalent_ellipse ellipse_elong.");
	rt = ResultsTable.getActiveTable()
	no_rows = rt.size()
		
	# makes a data dictionary for keeping track of the morphology data
	fa_data = OrderedDict()
	fa_data["basename"] = [basename] * no_rows #basename of image
	fa_data["id"] = [None] * no_rows  #focal adhesion id
 	fa_data["fa_area_px"] = [None] * no_rows #area in pixels
	fa_data["perim_px"] = [None] * no_rows #perimeter in pixels
 	fa_data["x0_px"] = [None] * no_rows #centroid x in pixels #####make sure these centroid values are not affected by scale (otherwise remove scale)
	fa_data["y0_px"] = [None] * no_rows #centriod y in pixels #####make sure these centroid values are not affected by scale (otherwise remove scale)
	fa_data["orientation_deg"] = [None] * no_rows #absolute orientation angle
	fa_data["ar"] = [None] * no_rows #aspect ratio
	fa_data["sf"] = [None] * no_rows #shape factor
	fa_data["dist_px"] = [None] * no_rows #distance from focal adhesion to boundary
	fa_data["angle_deg"] = [None] * no_rows #angle of the FA centroid with respect to the mask center of mass
	fa_data["rel_orientation_deg"] = [None] * no_rows #relative orientation of the FA with respect to the mask center of mass
	
#	#for testing
#	IJ.run(labels, "Specify...", "width=10 height=10 x="+str(int(round(x0_mask)))+" y="+str(int(round(y0_mask)))+" oval centered")
#	IJ.setForegroundColor(255, 0, 0)
#	IJ.run(labels, "Fill", "slice")

	# makes lists for keeping track of peripheral and central focal adhesions
	periph_count = 0
	central_count = 0
	
	# extracts the FA data and calculates some additional parameters
	if no_rows > 0:
		for i in range(no_rows):
		
			fa_data["id"][i] = i + 1
			fa_data["fa_area_px"][i] = float(rt.getValue('PixelCount',i))
			fa_data["perim_px"][i] = float(rt.getValue('Perimeter',i))
			fa_data["x0_px"][i] = float(rt.getValue('Ellipse.Center.X',i))
			fa_data["y0_px"][i] = float(rt.getValue('Ellipse.Center.Y',i))
			fa_data["orientation_deg"][i] = -float(rt.getValue('Ellipse.Orientation',i)) #negative because origin is top left
			fa_data["ar"][i] = float(rt.getValue('Ellipse.Elong',i))
			fa_data["sf"][i] = float(rt.getValue('Perimeter',i)) / sqrt(float(rt.getValue('PixelCount',i)))
			
			#gets distance to nearest contour pixel using distance transform
			x0_round = int(round(float(rt.getValue('Ellipse.Center.X',i))))
			y0_round = int(round(float(rt.getValue('Ellipse.Center.Y',i))))
			fa_data["dist_px"][i] = dist_map.getPixel(x0_round,y0_round)[0]
			periph_thresh = dist_map.getProcessor().getMax() * periph_frac
			if fa_data["dist_px"][i] < periph_thresh:
				periph_count += 1
			else:
				central_count += 1
		
			#calculates the relative orientation with respect a vector to the mask centroid
			fa_data["angle_deg"][i] = degrees(-atan2(fa_data["y0_px"][i] - y0_mask, fa_data["x0_px"][i] - x0_mask)) #negative because origin is top left
			fa_pos_orient = angle_to_orientation(fa_data["angle_deg"][i])
			fa_data["rel_orientation_deg"][i] = get_angular_diff(fa_pos_orient,fa_data["orientation_deg"][i])

	#closes the results window
	IJ.selectWindow(rt.title)
	IJ.run("Close")
	
	#starts an ordered dict for collecting the image data
	image_data = OrderedDict()
	image_data["basename"] = [] #basename of image
	image_data["mask_area_px"] = [] #area of the mask in pixels
	image_data["no_fas"] = [] #total number of focal adhesions
	image_data["area_fraction"] = [] #fraction of area covered by focal adhesions
	image_data["frac_periph_fas"] = [] #fraction of peripheral fas compared to total fas
	keys = ["fa_area_px","perim_px","orientation_deg","ar","sf","dist_px","rel_orientation_deg"]
	for key in keys:
		image_data[key + "_mean"] = []
		image_data[key + "_std"] = []
	
	# calculates some parameters on the level of the image and records
	if no_rows > 0:
	
		image_data["basename"] = basename
		image_data["mask_area_px"] = area_mask 
		image_data["no_fas"] = len(fa_data["id"])
		image_data["area_fraction"] = sum(fa_data["fa_area_px"]) / area_mask
		image_data["frac_periph_fas"] = periph_count / (central_count + periph_count)

		# adds means/SDs of the fa data to the image data
		for key in keys:
			image_data[key + "_mean"] = mean(fa_data[key])
			image_data[key + "_std"] = std(fa_data[key])
			
	#saves the focal adhesion and image data
	print("Saving Focal Adhesion and Image Data")
	fa_data_to_save = write_csv(od_to_matrix(fa_data),join(save_dir,basename + "_fa_data.csv"))
	image_data_to_save = write_csv(od_to_matrix(image_data),join(save_dir,basename + "_image_data.csv"))
	
	# close remaining images
	img.close()
	bin.close()
	labels.close()
			
	return fa_data, image_data

def analyze_image_dir(parent_dir,mask_suffix="_mask",file_ending=".tif",periph_frac=0.1):
	"""Goes through a directory of images for analysis of focal adhesions.
	
	Parameters
    ----------
    parent_dir : string
    	path to the directory containing images
    mask_suffix : string
    	suffix of the mask file with respect to the basename of the image file
    file_ending : string
    	the file ending (usually an image extension like ".tif")
    periph_frac : float
    	the distance fraction to be considered a peripheral adhesion (fraction of max distance of the mask from the border)
    Returns
    ----------
    fa_data : OrderedDict
        ordered dictionary containing the individual focal adhesion data, compiled for all images in the directory
    image_data : OrderedDict
    	ordered dictionary containing the image-level data and means/sds for focal adhesion data, compiled for all images in the directory
	
	"""


	print("Analyzing Directory: " + parent_dir)

	# initializes dictionaries for saving all data together
	fa_data = OrderedDict()
	image_data = OrderedDict()
	
	# gets a list of image files in the directory and sorts
	file_list = [_ for _ in listdir(parent_dir) if _.endswith(file_ending) and not mask_suffix in _]
	file_list = natural_sort(file_list)

	# goes through the parent directory and finds the data directories and data files to work on
	for file_name in file_list:
		
		# gets the file path
		file_path = join(parent_dir,file_name)
			
		# analyzes the image
		fa_data_tmp, image_data_tmp = analyze_image(file_path,mask_suffix,file_ending,periph_frac)
		
		# appends the data to keep track of everything together
		fa_data = concat_ods(fa_data,fa_data_tmp)
		image_data = concat_ods(image_data,image_data_tmp)
												
	return fa_data, image_data

def main():
	"""Main function for the focal adhesion analyzer. 
	   Any adjustable parameters can be found here.
	   Upon running the script, you will be prompted to open a directory containing images.
	
	"""

	#ensures black background for binary images
	IJ.run("Options...", "iterations=1 count=1 black")
	IJ.run("Colors...", "foreground=white background=black selection=yellow")

	# specifies the parent directory and gets the FA data
	parent_dir = IJ.getDirectory("Please select directory containing FA and mask images")
	periph_frac = 0.1 #the fraction of FAs considered 'peripheral' with respect to the max distance from the edge to the center of the cell
	mask_suffix = "_mask" #the suffix at the end of the mask images (note: the rest of the file name should be the same between the fa image and mask
	file_ending = ".tif" #the file extension for your images (usually '.tif')
	fa_data, image_data = analyze_image_dir(parent_dir,mask_suffix="_mask",file_ending=".tif",periph_frac=periph_frac)
	
	#writes csv files to output data
	print("Saving Compiled Data")
	write_csv(od_to_matrix(fa_data), join(parent_dir,"fa_data.csv"))
	write_csv(od_to_matrix(image_data), join(parent_dir,"image_data.csv"))
	
	print("All Done!")
	
if __name__ in ["__builtin__","__main__"]:
    main()