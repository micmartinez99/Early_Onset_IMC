#!/usr/bin/python

from pathlib import Path
import napari
import pandas as pd
import tifffile

panel_file = "/Users/michaelmartinez/Desktop/Combined_EOLO/Data/combined_steinbock/panel.csv"
img_file = "/Users/michaelmartinez/Desktop/Combined_EOLO/Data/combined_steinbock/img/Slide2-5_CRCTMA3_CRCTMA4_001.tiff"
mask_file = "/Users/michaelmartinez/Desktop/Combined_EOLO/Data/combined_steinbock/masks_deepcell/Slide2-5_CRCTMA3_CRCTMA4_001.tiff"

channel_names = None
if Path(img_file).exists():
	panel = pd.read_csv(panel_file)
	channel_names = panel.loc[panel["keep"] == 1, "name"]
	print(len(channel_names), "channels")
	
img = None
if Path(img_file).exists():
	img = tifffile.imread(img_file)
	
mask = None
if Path(mask_file).exists():
	mask = tifffile.imread(mask_file)
	
viewer = napari.Viewer()

viewer.axes.visible = True
viewer.dims.axis_labels = ("y", "x")
viewer.scale_bar.visible = True
viewer.scale_bar.unit = "um"

if img is not None:
	img_layers = viewer.add_image(
	data = img,
	channel_axis = 0, 
	colormap = "gray",
	name = channel_names,
	blending = "additive",
	visible = False)
		
if mask is not None:
	mask_layer = viewer.add_labels(
	data = mask,
	name = "Cells",
	blending = "translucent",
	visible = False)