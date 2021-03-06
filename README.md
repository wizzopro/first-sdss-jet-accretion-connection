# first-sdss-jet-accretion-connection
"A deeper look at the jet-accretion connection in AGN (Active Galactic
Nuclei)."

This project page will feature code I developed for attaining my Ph.D. in
black hole astrophysics. The paper with the above title that will
summarise the results obtained with this code is in preparation and
will be submitted.

Practically speaking this uses Python 2.7/astropy to combine archival data from the VLA (Very Large Array) FIRST (Faint Images of the Radio Sky at Twenty Centimeters) survey, with optical spectroscopic data from the SDSS (Sloan Digital Sky Survey).

As the final paper has not been published yet, currently only the first main module is online (first_imtab). Please see the description below. 


|FILE | Description |
| --- | --- |
|README.md           | # This file|
|first_imtab.py       |# Module exploring from which original radio image to get radio data that matches SDSS data|
|imfuncs.py           |# Helper functions for cutouts of radio data|
|requirements.txt     |# List of required packages from pip freeze|
|settings.py          |# Constants and settings that should be accessible in multiple modules|



README for first_imtab module:

To start, download all available data from the VLA (Very Large
Array) FIRST (Faint Images of the Radio Sky at Twenty-Centimeters)
survey (http://sundog.stsci.edu/). These images are available through anonymous ftp from
ftp://archive.stsci.edu/pub/vla_first/data and comprise a total
size of about 340 GB. Set the folder where the data is downloaded in
the folders.firstdata_folder variable. This module consequently selects the last version of every available
FIRST image, gathers and constructs a table with important necessary
of these images (most importantly the ra/dec covered by all useable first images as
well as the RMS as estimated by the survey itself). Also does a new estimate for RMS from multiple areas in a FIRST image
(needed later for statistics). This module can be run independently to create the table, but is also called from main module (first.py).
