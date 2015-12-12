# pyraf-scripts
Scripts to simplify data reduction with PyRAF

## pre-process.py
Use this script to create a master bias, dark, and flat, and apply them to a set of FITS cubes.

The script will expect all of your FITS cubes to be in a source directory, and will create output directory containing the masters and processed FITS cubes.

To run the code:

`pre-process.py [-h] [--src_dir SRC_DIR] [--dst_dir DST_DIR]`

    optional arguments:
      -h, --help         show this help message and exit
      --src_dir SRC_DIR  source directory of raw FITS files
      --dst_dir DST_DIR  destination directory for processed and master FITS files

## txdump-parse.py
Use this script to parse the output of PyRAF's txdump tool and create a consistent set of star data for photometry use. This script attempts to find "the same" star from image to image and output the data set for each one. Daofind is not good at providing consistently numbered stars, so this attempts to mitigate that issue.

The txdump output needs to have the following columns in this order:
    IMAGE,ID,XCENTER,YCENTER,OTIME,MAG,MERR

These settings at the top of the code can adjust the parser:

* px_threshold (int) - The maximum pixel change (in the x or y direction) between images for a star to be considered "the same".
* mag_threshold (int) - The maximum magnitude change between images for the star to be considered "the same"
* skip_px_threshold (int) - The maximum pixel change (in the x or y direction) between images for a star to be considered "different".
* skip_mag_threshold (int) - The maximum magnitude change between images for a star to be considered "different".
* assume (True or False) - This can force the script to assume stars between the px_threshold and the skip_px_threshold, or the mag_threshold and the skip_mag_threshold are always the same or never the same stars. If this is not set (or set to None), all stars between these values will be confirmed via the user.

The script outputs a "data.csv" file with the cleaned data.

To run the code:

`txdump-parse.py <filename>`

    arguments:
        <filename>  The file from txdump to be parsed
