# first-sdss-jet-accretion-connection
"A deeper look at the jet-accretion connection in AGN (Active Galactic
Nuclei)."

This project page will feature code I developed for attaining my Ph.D. in
black hole astrophysics. The paper with the above title that will
summarise the results obtained with this code is in preparation and
will be submitted end of summer 2017.

Practically speaking this uses Python/astropy to combine archival data from the VLA (Very Large Array) FIRST (Faint Images of the Radio Sky at Twenty Centimeters) survey, with optical spectroscopic data from the SDSS (Sloan Digital Sky Survey).

As the final paper has not been published yet, currently only the first main module is online (first_imtab). Please see the description in this file. 


|FILE | Description |
| --- | --- |
|README.md           | # This file|
|first_imtab.py       |# Module exploring from which original radio image to get radio data that matches SDSS data|
|imfuncs.py           |# Helper functions for cutouts of radio data|
|requirements.txt     |# List of required packages from pip freeze|
|settings.py          |# Constants and settings that should be accessible in multiple modules|
