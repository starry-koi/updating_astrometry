#!/usr/bin/env python
import numpy as np
import reproject as rp
import copy
import sys
import warnings
import os
import subprocess as sp
from astropy.io import fits
from astropy import wcs
from astropy.utils.exceptions import AstropyWarning
from scipy import signal
from os import listdir
from os.path import isfile,join

#Author: Lilikoi Latimer
#April 2018
#Updated September 2022

'''
=========================================================================================


### Quirks and possible issues ###

- Only works on .fits files

- Uses ONLY ONE REFERENCE IMAGE with multiple broken images. The path to the folder for
  the reference file is purely for convenience.

- Have to MANUALLY set the tolerance between methods (multiplying arrays versus 
  correlating them) and the approximate maximum distance the astrometry can be off by.
  Though keeping these at the default has mostly worked out fine in my experience.

- When calculating maxima the code averages if there are multiple maxima - this works
  well if the maxima are close together, but very poorly if they are from two separate,
  distant sources.
  
- This program assumes that the broken image is pretty close (astrometrically speaking)
  to the reference image - the astrometry is mostly correct, just slightly off. While
  this is a reasonable assumption, if the broken image is too far off, the multiplication
  method won't work, and every image will print both multiplication and correlated
  results - and you won't be able to check the correlation method against the 
  multiplication method, so you might get some error there.
  
- Correlation can also mess up if there's a source brighter than the source we care
  about - as the correlation can sometimes pick that as the max, leading to wildly 
  incorrect answers. That's why we check the multiplication method and the correlation
  method against eachother. While they both have flaws, the flaws usually result from
  different scenarios.




### Broad Overview ###

This program reads in a single reference image and a folder containing the 'broken' 
images that need astrometry corrections. From there we cycle through the broken images,
correcting astrometry as we go. There are two methods employed here - correlation between
the images, and a more simple multiplication of arrays. 

First we get the reference image in the 'pixel coordinates' of the broken image, giving
us two arrays of the same size containing our relevant data. There are two ways to do
this projection - exact and interpolated. The default is exact, although it takes a bit
longer (how much exactly depends on your computer).

For the multiplication method: we take our two arrays and simply multiply them together.
We take the resulting array and find the location of the maximum, which tells us roughly
where the two arrays are most similar. We then 'draw a box' around that point of side 
lengths 2*rough_max_separation, and search for the maximum in the broken data and the
reference data. Then we simply take the difference and append the header accordingly
(Note that we also average if more than one maximum is found).

For the correlation method: we auto-correlate the reference image to get a baseline, then
cross-correlate the reference and broken image. We take the difference of the coordinates
of the maximum of the auto-correlated array and the cross-correlated array, and append
the header accordingly.

From here we check whether the results from the multiplication and correlation methods
agree to within a user-specified tolerance, and if they do we pick one of the methods
to go with (default is correlation). If they disagree, then we save the output of both
methods.

The final files are written out with the filename appended with _fixed_mult or 
_fixed_corr as appropriate.

Also, throughout the rest of the program:
ref -> reference image
off or broken -> broken image


=========================================================================================
'''
#################################################
################### VARIABLES ###################
#################################################


#full path of reference folder. Note that it should only contain ONE reference file -
#defining the folder is just for convenience.
ref_path = os.getcwd() + "/reference/"

#full path of the folder containing the broken images
broken_path = os.getcwd() + "/broken/"

tol 		= 3             #in pixels - tolerance for difference in corr. and mult. methods

sep 		= 2             #in arcsec - rough maximum of initial astometrical differences

default 	= 'mult'      	#define what to take when both options agree within tolerance
				#options of "corr" or "mult" for correlation and multiplication
				#method, respectively
				
projection 	= "exact"  	#define what projection to use
				#options of "exact" or "interp", for exact and interpolated
				
save_opt 	= 'sep'  	#where to put the corrected files
				#options of 'sep' and 'mix'. 'sep' puts the corrected files
				#in a new folder (labelled /results) and 'mix' puts the
				#corrected files in with the to-be-fixed files (in the
				#/broken folder) (with different filenames, obviously).


#Note: defaults are
#tol = 3
#sep = 2
#default = 'corr'
#projection = 'exact'
#save_opt = 'mix'

#we'll print to the screen what all of these are - both here and at the end (ish, we'll
#only print the paths here).
print('Tolerance of',tol,'pixels')
print("Rough separation of",sep,"arcseconds")
print("Default method:",default)
print("Projection method:",projection)

#========================================================================================
#========================================================================================

#################################################
############ DECLARING FOR LATER USE ############
#################################################


#suppresses the annoying FITSFixedWarning about the RADECSYS keyword being deprecated,
#which as far as I can tell doesn't actually affect anything.
warnings.simplefilter('ignore',AstropyWarning)

#getting the reference file
ref_file_temp = [f for f in listdir(ref_path) if isfile(join(ref_path,f))]
ref_file = ref_path + ref_file_temp[0]

#getting the list of broken filenames from the folder
broken_files = [f for f in listdir(broken_path) if isfile(join(broken_path,f))]

#opening the reference file before the for loop starts
ref_full = fits.open(ref_file)

#declaring some lists we'll use to print the files where the methods agreed and where
#they disagreed. Also get the separation distance into degrees.
disagree_files = []
disagree_delta_list = []
agree_files = []
agree_correction_list = []
arcsec = sep*2/3600

#Making the results folder and storing a few different filepaths for convenience
current_dir = os.getcwd()
results_dir = os.getcwd() + "/results/"
if save_opt == 'sep':
	sp.run(["mkdir", "-p", results_dir])

#========================================================================================
#========================================================================================

#################################################
################ FIXING THE FILES ###############
#################################################


for n in range(0,len(broken_files)):
	
	
	#================================================================================
	#================================================================================
	###############################
	##### DECLARING VARIABLES #####
	###############################
	'''
	The first part of this is just getting the strings for the filenames and for
	printing information to the screen all together in one place for convenience.
	This has to go inside the for loop because it changes for each new file.
	
	The remaining part is dealing with the data from the broken image - importing
	it, and projecting the reference image's data into the frame (pixel coords) of
	the broken image.
	'''
	
	#getting the filenames and filepaths all together, as well as the names of the
	#fixed files that we'll be writing out to.
	broken_file = broken_path + broken_files[n]
	if save_opt == 'sep':
		fixed_file_mult = str.replace(broken_files[n],".fits","_fixed_mult.fits")
		fixed_file_corr = str.replace(broken_files[n],".fits","_fixed_corr.fits")
	elif save_opt == 'mix':
		fixed_file_mult = str.replace(broken_file,".fits","_fixed_mult.fits")
		fixed_file_corr = str.replace(broken_file,".fits","_fixed_corr.fits")
	else:
		print("Please specify save option (save_opt variable)!")
		sys.exit()
	
	#defining some strings for convenience, just used for printing the information to
	#the screen.
	current_num = '(' + str(n+1) + '/' + str(len(broken_files)) + ')'
	placeholder = '     '
	print(current_num,"Fixing",broken_files[n],"...")
	
	#get the data and header from the broken file
	data_off,header_off = fits.getdata(broken_file,header = True)
	
	#simple if statement to pick which reprojection we will use.
	if projection == "exact":
		reference,footprint = rp.reproject_exact(ref_full,header_off)
	elif projection == "interp":
		reference,footprint = rp.reproject_interp(ref_full,header_off)
	else:
		print(placeholder,"Default projection method not specified - check code")
		sys.exit()
	
	#fixing the projected reference image - if the broken image doesn't lie
	#completely within the reference image, there will be some cutoff, and the
	#reproject function replaces those values of the reference image that don't
	#actually exist with NaNs. This just sets those to zero.
	#Also, new problem that popped up like a year after writing this: apparently,
	#the broken image can also have NaN values (seems that along the edge of the
	#image we can get a few NaNs), so we'll go ahead and set those to zero as well.
	ref_fixed = reference
	ref_fixed[np.isnan(ref_fixed)] = 0
	data_off[np.isnan(data_off)] = 0
	
	#================================================================================
	#================================================================================
	
	###############################
	#### MULTIPLICATION METHOD ####
	###############################
	'''
	Since the two images are supposedly reasonably close together (as we'll assume
	that the astrometry is mostly correct) we can just multiply them and we should
	get a somewhat reasonable approximation of where they overlap. Note that if some
	noise happens to be unfortunately placed, this might give some strange
	corrections. Usually the noise between the reference image and the broken image
	is different enough that this doesn't happen.
	
	We also account for if there are multiple maxima (like if the original pixel size
	of the reference image is much bigger than the broken image's pixel size, which
	results in multiple pixels having the same value). We simply average over the 
	coordinates, rounding to the nearest integer. This only really helps if the 
	maxima are close together, though.
	
	We then take this location of the maximum as the center of our box that we search
	within to find the location of the maxima of the reference image and the maxima
	of the broken image.
	
	We take the difference between the two maxima and call that our correction.
	'''
	
	#multiplying the broken and reference data together and finding the location of
	#the maximum.
	multiplied = ref_fixed*data_off
	mult_temp = np.where(multiplied == np.amax(multiplied))
	mult_max_coords = np.array([mult_temp[0].astype(int),mult_temp[1].astype(int)])
	
	#this gives us the means to transform between the broken image's pixel coords
	#and the WCS (RA, dec).
	w_off = wcs.WCS(header_off,broken_file)
	
	#take the coordinates of the maximum from the mult.array, transform them into
	#WCS coords (averaging to account for multiple maxima).
	center_init = w_off.all_pix2world(mult_max_coords[0],mult_max_coords[1],0)
	ra_avg = np.mean(center_init[0])
	dec_avg = np.mean(center_init[1])
	center = np.array([ra_avg,dec_avg])
	
	#defining the corners of our box (in WCS coords) that we'll search.
	wcs_corner1 = center + np.array([-arcsec,-arcsec])
	wcs_corner2 = center + np.array([-arcsec,arcsec])
	wcs_corner3 = center + np.array([arcsec,-arcsec])
	wcs_corner4 = center + np.array([arcsec,arcsec])
	
	#transforming the box into the pixel coords of the broken image
	corner1 = w_off.all_world2pix(wcs_corner1[0],wcs_corner1[1],0)
	corner2 = w_off.all_world2pix(wcs_corner2[0],wcs_corner2[1],0)
	corner3 = w_off.all_world2pix(wcs_corner3[0],wcs_corner3[1],0)
	corner4 = w_off.all_world2pix(wcs_corner4[0],wcs_corner4[1],0)
	
	#gathering the row,column values of our box
	corner_row = np.array([corner1[0],corner2[0],corner3[0],corner4[0]])
	corner_col = np.array([corner1[1],corner2[1],corner3[1],corner4[1]])
	
	#setting any negative values(indicating a corner of the box lies outside the
	#array) to zero.
	corner_row[corner_row<0] = 0
	corner_col[corner_col<0] = 0
	
	#finding the minimum and maximum row and column values, which will make up the 
	#corners of our box. The +1 is so the beginnings and endings are nice and 
	#symmetric (as x[a:b] goes from x[a] to x[b-1]). We also round everything to
	#the nearest integer as arrays only take ints as arguments.
	row_start = int(round(np.amin(corner_row)))
	col_start = int(round(np.amin(corner_col)))
	row_end = int(round(np.amax(corner_row))) + 1
	col_end = int(round(np.amax(corner_col))) + 1
	
	#essentially just making an array that has our data from the broken image inside
	#our little box and zeros everywhere else.
	interim_broken = data_off[row_start:row_end, col_start:col_end]
	search_broken = np.zeros(np.shape(data_off))
	search_broken[row_start:row_end,col_start:col_end] = interim_broken
	
	#same thing but for the reference image, now.
	interim_ref = ref_fixed[row_start:row_end, col_start:col_end]
	search_ref = np.zeros(np.shape(ref_fixed))
	search_ref[row_start:row_end,col_start:col_end] = interim_ref
	
	#finding the location of the brightest pixel in our box in the broken image, 
	#and averaging over multiple maxima if we have them.
	broken_max = np.amax(search_broken)
	broken_max_coords = np.where(search_broken == broken_max)
	broken_row_avg = int(round(np.mean(broken_max_coords[0])))
	broken_col_avg = int(round(np.mean(broken_max_coords[1])))
	broken_avg = np.array([broken_row_avg,broken_col_avg])
	
	#same thing but for the reference image.
	ref_max = np.amax(search_ref)
	ref_max_coords = np.where(search_ref == ref_max)
	ref_row_avg = int(round(np.mean(ref_max_coords[0])))
	ref_col_avg = int(round(np.mean(ref_max_coords[1])))
	ref_avg = np.array([ref_row_avg,ref_col_avg])
	
	#evaluating the difference between the maxima - essentially the difference
	#between our sources, as with correct astrometry they would be in the same
	#place.
	delta_row_mult = ref_avg[0] - broken_avg[0]
	delta_col_mult = ref_avg[1] - broken_avg[1]
	delta_vec_mult = [delta_row_mult,delta_col_mult]
	
	#making a copy (have to use the copy.deepcopy or else the original changes too)
	#of the header and modifying it with our corrections. crpix1 moves left/right (- for right)
	#and crpix2 moves up/down (+ for down)
	header_off_mult = copy.deepcopy(header_off)
	header_off_mult["CRPIX1"] = header_off["CRPIX1"] - delta_col_mult
	header_off_mult["CRPIX2"] = header_off["CRPIX2"] - delta_row_mult
	
	#printing the corrections to the screen.
	print(placeholder,"Mult. correction of",delta_vec_mult,"pixels.")
	
	#================================================================================
	#================================================================================
	
	###############################
	###### CORRELATION METHOD #####
	###############################
	'''
	Much, much simpler. We take the auto-correlation of the reference image, and find
	the maximum - this is our baseline. Then we cross-correlate the reference and 
	broken image, and find the maximum of that. We account for possible multiple
	maxima by taking the average.
	
	The difference between the maxima is simply the correction, which we apply to the
	header.
	
	Note that if there's another source inside the array that is brighter than the 
	source we care about, this can give wildly incorrect answers as the cross-
	correlation might pick out that point as the maximum. That's why we check this
	with the multiplication method, which shouldn't be as susceptible to that flaw.
	'''
	
	#define the auto-correlated and the cross-correlated arrays.
	auto_corr = signal.correlate(ref_fixed,ref_fixed)
	corr = signal.correlate(ref_fixed,data_off)
	
	#find the location of the brightest pixel for the auto-correlated array
	#(accounting for possible multiple maxima).
	auto_temp = np.where(auto_corr == np.amax(auto_corr))
	auto_corr_max_coords = np.array([int(auto_temp[0]),int(auto_temp[1])])
	auto_row_avg = int(round(np.mean(auto_corr_max_coords[0])))
	auto_col_avg = int(round(np.mean(auto_corr_max_coords[1])))
	auto_corr_coords_avg = np.array([auto_row_avg,auto_col_avg])
	
	#same thing but for the cross-correlated array
	corr_temp = np.where(corr == np.amax(corr))
	corr_max_coords = np.array([int(corr_temp[0]),int(corr_temp[1])])
	corr_row_avg = int(round(np.mean(corr_max_coords[0])))
	corr_col_avg = int(round(np.mean(corr_max_coords[1])))
	corr_coords_avg = np.array([corr_row_avg,corr_col_avg])
	
	#find the difference between the auto-correlated and cross-correlated arrays
	difference = auto_corr_coords_avg - corr_coords_avg
	
	#making a vector for aesthetic similarity to multiplication method
	delta_row_corr = -difference[0]
	delta_col_corr = -difference[1]
	delta_vec_corr = [delta_row_corr,delta_col_corr]
	
	#making a copy (same as in the multiplication method) and modifying the header
	header_off_corr = copy.deepcopy(header_off)
	header_off_corr["CRPIX1"] = header_off["CRPIX1"] - delta_col_corr
	header_off_corr["CRPIX2"] = header_off["CRPIX2"] - delta_row_corr
	
	#printing the corrections to the screen
	print(placeholder,"Corr. correction of",delta_vec_corr,"pixels.")
	
	#================================================================================
	#================================================================================
	
	###############################
	###### COMPARING METHODS ######
	###############################
	
	#finding the difference between the results of the multiplication method and the
	#correlation method, and printing to the screen. 
	delta = np.abs(np.array(delta_vec_mult) - np.array(delta_vec_corr))
	print(placeholder,"Difference is:",delta,"pixels.")
	
	#if statement detailing what happens if the difference exceeds the tolerance,
	#in which case we save and write out the results of both methods, and if the 
	#differences are within tolerance, in which case we take the results of the 
	#default method.
	if np.any(delta > tol) == True:
		
		print(placeholder,placeholder,"Difference outside tolerance bounds!")
		print(placeholder,placeholder,"Saving both mult. and corr. fixed files.")
		
		#Change directories to results folder, save images, then change back
		if save_opt == 'sep':
			os.chdir(results_dir)
		fits.writeto(fixed_file_mult,data_off,header_off_mult)
		fits.writeto(fixed_file_corr,data_off,header_off_corr)
		os.chdir(current_dir)
		
		disagree_files.append(broken_files[n])
		disagree_delta_list.append(delta)
		
	else:
		agree_files.append(broken_files[n])
		
		#Change directories to results folder, save images, then change back
		if save_opt == 'sep':
			os.chdir(results_dir)
		if default == "mult":
			fits.writeto(fixed_file_mult,data_off,header_off_mult)
			agree_correction_list.append(delta_vec_mult)
		elif default == "corr":
			fits.writeto(fixed_file_corr,data_off,header_off_corr)
			agree_correction_list.append(delta_vec_corr)
		else:
			sys.exit("Default preference not specified - check code")
		os.chdir(current_dir)

#========================================================================================
#========================================================================================

#################################################
############# DISPLAYING THE RESULTS ############
#################################################

#printing the currently selected options to the screen (again):
print('-------------------------')
print('Tolerance of',tol,'pixels')
print("Rough separation of",sep,"arcseconds")
print("Default method:",default)
print("Projection method:",projection)

#printing to the screen how many files had both methods in agreement (to within the 
#specified tolerance).
agreements = '(' + str(len(agree_files)) + '/' + str(len(broken_files)) + ')'
print("Agreements:",agreements)
for m in range(0,len(agree_files)):
	print(placeholder,agree_files[m],"- change of",agree_correction_list[m])

#printing to the screen how many file had the methods in disagreement (outside of the
#specified tolerance).
disagreements = '(' + str(len(disagree_files)) + '/' + str(len(broken_files)) + ')'
print("Disagreements:",disagreements)
for m in range(0,len(disagree_files)):
	print(placeholder,disagree_files[m],"- diff. of",disagree_delta_list[m])


print("Done.")

#========================================================================================
#========================================================================================
#########################################################################################
########################################## END ##########################################
#########################################################################################
