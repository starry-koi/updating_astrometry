Updating Astrometry
===

This code functions as a first pass to correcting astrometry on fits images. It takes a reference image assumed to have correct astrometry (stored in the /reference folder) and a list of 'broken' images that have astrometry in need of correcting (stored in the /broken folder). The reference image is then used to correct the astrometry of the broken images. Two methods are used to provide some measure of redundancy. It 

The corrected images should still be viewed and compared manually with the reference image (e.g. using an RGN frame in DS9). While this code can make the astrometry-correcting process significantly less tedious, the results are not always perfect, especially if there is a large difference in pixel sizes between reference and broken images. Some manual fine-tuning might still be necessary.

For more depth than the ReadMe, see `outline.pdf`.


Setup
---

You will need to delete the `placeholder.fits` files in the /reference and /broken folders and replace them with the appropriate images. The /reference folder should contain the image with correct astrometry (in the past I've used images from a large survey such as SDSS), and the /broken folder should contain all the images you want to correct the astrometry of.


Brief Rundown of Code
---

A more in-depth walkthrough is available in the `astrom_final.py` code as well as the `outline.pdf` file (the latter probably being easier to read). The overview is as follows:

* We project the reference image into the WCS frame of each broken image (effectively getting them both in the same pixel/array coordinates). 
* The first method used to find the astrometric correction is 'mult' - this involves multiplying the two arrays (reference and broken) together. This tends to reduce the effect of noise.
* The second method is 'corr' - correlating the two image arrays. When this works, it tends to be more accurate than the multiplication method, but it's more susceptible to being thrown off by noise.
* We then compare the astrometric corrections found by the two methods. 
  * If they're the same to within tolerance (the `tol` variable), then we use the results of method specified by the `default` variable and save the image with updated astrometry with the name format `originalfilename_fixed.fits`.
  * If they differ by more than the tolerance, we save two copies of the image with the astrometry updated via each method, with the name formats as `originalfilename_fixed_mult.fits` and `originalfilename_fixed_corr.fits`. 
* We save the resulting images either in the /broken folder or in a separate /results folder, depending on what the variable `save_opt` is set to. We also print out to the screen the filenames and the final astrometry corrections (in pixel/array coordinates) and whether or not the two methods agreed on the correction (see `outline.pdf` for an example of the terminal output).


Variables the User Will Have to Care About
---

* `save_opt` 
  * Where to save the resulting images. Options are 'sep' and 'mix', for storing them in a separate /results folder or putting them in the /broken folder with the starting images.

There are a few other variables (see `outline.pdf` or the code file `astrom_final.py` for details on them) but these, at least in my experience, don't need much adjustment from their defaults.


A Sample Run
---

Let's say you have 3 HST images in differing bands that you'd like to correct the astrometry of and combine into a composite RGB image. You first need to determine what image you will be matching to; if you want to correct the overall astrometry, then you likely download an image from a larger survey, such as SDSS, to serve as a reference image. If all you care about is that the final RGB composite look good and the three HST images match each other, then perhaps you pick one of the HST images as a reference.

Whichever you choose, you copy the reference image into the /reference folder and the rest of the images into the /broken folder. You then check that the `save_opt` variable is what you want it to be and run the code. Afterwards, you check the corrected images in DS9 and make small adjustments to the astrometry, if needed.


Rerunning the Code
---

Before rerunning the code, the user should make sure that the /broken folder only contains the original images to be corrected, as if `save_opt` was set to 'mix' the corrected images will have been saved to the /broken folder as well. If not, the code will attempt to correct the already-corrected images. 

The code will also throw an error if it detects a file with a name like the corrected images from the first run (e.g. `originalfilename_fixed.fits`), so if `save_opt` was set to 'sep', with the corrected images from the first run stored in a separate results folder, the user should rename this folder to avoid collisions.


Quirks
---

* The pixel corrections shown in the terminal output are in the frame of the 'broken' image, not the reference image.
* This code functions as a first pass - some manual fine-tuning will likely still be needed. Hopefully this code makes the entire process less tedious, though.





