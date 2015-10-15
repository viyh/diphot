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
