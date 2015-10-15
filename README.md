# pyraf-scripts
Scripts to simplify data reduction with PyRAF

## pre-process.py
Use this script to create a master bias, dark, and flat, and apply them to a set of FITS cubes.

The script will expect all of your FITS cubes to be in a "data" directory, and will create an "output" directory containing the masters and processed FITS cubes. It does not modify the original FITS files.

To run the code:

`python pre-process.py`
