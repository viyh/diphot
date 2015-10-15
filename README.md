# pyraf-scripts
Scripts to simplify data reduction with PyRAF

## pre-process.py
Use this script to create a master bias, dark, and flat, and apply them to a set of FITS cubes.

The script will expect all of your FITS cubes to be in a "data" directory, and will create an "output" directory containing the masters and processed FITS cubes. It does not modify the original FITS files.

Pre-requisites include a working PyRAF/IRAF installation complete with a login.cl file, and a Cal Poly instrument "dat" file. This can be setup by adding the following lines to these two files in your IRAF "uparm" directory:

* ccdsetint.par:
`instrument,s,a,'cp',,,'Instrument ID (type ? for a list)'`

* imdccdred.par:
`instrument,s,h,'cp.dat',,,'CCD instrument file'`

Create a "cp.dat" file which contains the following:

    subset          FILTER

    darktime        EXPTIME

    'Dark Frame'    dark
    'Bias Frame'    zero
    'Light Frame'   object
    'Flat Field'    flat

Once this is setup, you can run the code:

`python pre-process.py`
