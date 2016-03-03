To Do:

* Find out if greater precision for merr is possible
* Add "skip" functionality for reduction (skip-flat, skip-zero, etc.)
* Display average delta star comps
* Rework data point/star objects

Done:

v0.2
* Fix trailing slash failures on inputted directory names
* Implement configuration file for settings
* Fix duplicate log messages from running multiple operations
* Add x,y arguments
* Implement target guess based of user provided coords of star in first image
* Make compatible with dirac4
    - dirac4 needs upgrade: numpy, matplotlib
* Add PSF radial plots
* Add tresca dump functionality
* Fix argument parser
* Document config settings in comments

v0.1
* Implement "use_new_imt = no" workaround for IRAF bug
* Fix broken dir creation when type dirs already exist
