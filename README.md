# pyraf-scripts
Scripts to simplify data reduction with PyRAF

## pre-process.py
Use this script to create a master bias, dark, and flat, and apply them to a set of FITS cubes.

The script will expect all of your FITS cubes to be in a source directory, and will create output directory containing the masters and processed FITS cubes.

To run the code:

`python pre-process.py <--src_dir=data> <--dst_dir=output>`

* --src_dir
Default: data
This is the source directory with your raw FITS files. These will not be modified.

* --dst_dir
Default: output
The destination directory which will contain the master bias, dark, flat, and processed FITS files once the script has been run.
