# STEPUP Image Analysis 

<img src="https://github.com/helenarichie/helenarichie/blob/master/images/STEPUP_logo.png" alt="logo" width="200"/>


STEPUP Image Analysis (SIA) is a pipeline written for easy photometric analysis to extact light curves from image data using differential aperture photometry. SIA runs in three main (independent) steps: instrument signature removal (**ISR**), generating plate solutions (**ASTROM**), and differential aperture photometry (**PHOT**). Please feel free to contact the author, [Helena Richie](https://helenarichie.github.io/helenarichie/), with any questions!

This respository contains the SIA pipeline and a user manual written explicity for use by University of Pittsburgh students that are either members of the STEPUP team or students in ASTRON 1263. The guide can also be applied for general use, omitting location-specific steps. A manual for general use is coming soon.

## Installation
To run SIA_II, user must have Python 3 with standard libraries and the following Python packages: 
- Matplotlib
- NumPy
- AstroPy
- PhotUtils

These packages are all available for installation using `pip`. Additionally, the user must install and compile the [WCSTools software package](http://tdc-www.harvard.edu/wcstools/) and add the executables to their device's `PATH`.

## Outside Resources
SIA_II uses [Astrometry.net](http://astrometry.net/) and [WCSTools software package](http://tdc-www.harvard.edu/wcstools/) to plate solve inputted datasets.

## Sample Data Products
SIA outputs plots that depict the behavior of comarison stars, the target's light curve behavior, and a summary of image centroiding. In addition to this, SIA outputs data files with absolute and net count values of the target star for each image. Examples of plots can be seen below and in the SIA user manual.

Target light curve:  
<img src="https://github.com/helenarichie/STEPUP_image_analysis_II/blob/master/user_manual/wasplightcurve.png" alt="lc" width="450"/>

Centroid corrections:  
<img src="https://github.com/helenarichie/STEPUP_image_analysis_II/blob/master/user_manual/centroid.png" alt="centroid" width="450"/>

Unscaled comparison star magnitudes:  
<img src="https://github.com/helenarichie/STEPUP_image_analysis_II/blob/master/user_manual/comp1r.png" alt="complc" width="450"/>
