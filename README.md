# DiPhot

https://github.com/viyh/diphot

Scripts to simplify differential photometry data reduction using PyRAF.

# Overview

These scripts can be used to reduce data for differential photometry. The general order of usage is as follows:

* process.py - Fully automatic processing of your FITS images to a light curve. This includes all of the steps from the other scripts below.
* reduce.py - Organizes and creates master zero, dark, and flat images, and removes these from object images.
* curveofgrowth.py - Create a curve of growth which helps you choose the currect aperture to use for photometry and gives an idea of the SNR of the object images.
* photometry.py - Find stars in images and collect magnitude data and created a data dump file.
* lightcurve.py - Parses the dump file to align stars by ID and dump to CSV. Create a differential light curve.

# Quick Start

DiPhot can be run in a completely automatic mode where it will process raw FITS files from start to light curve.

* Setup your diphot.yml file (copy the diphot_sample.yml and start with that). Make sure the path to your raw FITS files and the output directory are correct. You will need to specify the `target_x` and `target_y` of your target star in the first object FITS image in order to generate a light curve. Add any other settings that you want.

* Run diphot:
`python process.py -c diphot.yml`

# Configuration

To configure DiPhot, a YAML file needs to be created. The included `diphot_sample.yml` can be copied and edited as needed. It's a good idea to create a separate configuration file for each of your FITS data sets, that way you know exactly what settings you used for processing to reference later.

DiPhot can be invoked with the `-c` argument (or `--config_file`) to specify which configuration file to use. By default, it will look for a file named `diphot.yml`.

The default settings are in `diphot_defaults.yml` **which should never be edited**. This file is for reference only. A description of each DiPhot setting is in this file. The PyRAF settings show are the ones that DiPhot uses, however any additional PyRAF settings can be set in a similar fashion. Consult the [PyRAF documentation](http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?apphot.hlp) for more details about PyRAF specific settings.

# Scripts

## Usage

All of the scripts have the same command line arguments.

For example:

`python process.py [-h] [--debug] [--ignore_id IGNORE_IDS] [--comp] [--config_file CONFIG_FILE]`

    optional arguments:
      -h, --help            show this help message and exit
      --debug               turn on debug messages
      --ignore_id IGNORE_ID, -i IGNORE_ID
                            star to ignore (use multiple times to ignore stars)
      --comp                show individual star graphs, "comparison mode"
      --config_file CONFIG_FILE, -c CONFIG_FILE
                            configuration file (default: diphot.yml)

## process.py
This is for automatic processing of the FITS files. This will do the following:
* reduce the raw FITS files
* create a curve of growth and determine FWHM and optimal aperture for photometry
* find target and comparison stars in each FITS image
* sort stars based on small changes in X, Y coordinates between images
* create differential photometry light curve

Each of these steps can be run using the individual scripts which are described below in more detail.

## reduce.py
Use this script to create a master bias, dark, and flat, and apply them to a set of FITS cubes.

The script will expect all of your FITS cubes to be in a source directory, and will create output directory containing the masters and processed FITS cubes.

This will organize the FITS files into the following directory structure in the <OUTPUT DIRECTORY>:
    /<OUTPUT DIRECTORY>/
        zero/
            zero (bias) FITS files
        dark/
            dark FITS files
        flat/
            flat FITS files
        object/
            object FITS files
        master/
            master FITS files that are created for zero, dark, and flat images
        tmp/
            temporary directory used by these scripts

The object directory contains the final object science images with the master zero, dark, and flat subtracted.

## curveofgrowth.py
This calculates the FWHM and aperture with the maximum signal-to-noise ratio. An output curve of growth graph can optionally be displayed and/or saved if specified. By default, it chooses the FITS image halfway through the set to calculate these values from. A contour plot of the PSF can also be created for the target star.

## photometry.py
This runs the star finding algoritm, calculates the instrumental magnitude for each found star, and creates the txdump file containing these values. Since the star IDs that PyRAF assigns vary for a given star from image to image, this code attempts to identify the stars from image to image across the entire data set and assign them consistent IDs.

## lightcurve.py
This uses the txdump and creates the differential light curve. Additionally, it can create a TSV file which can be uploaded to the [TRESCA database](http://var2.astro.cz/EN/tresca/index.php).
