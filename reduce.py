#!/usr/bin/python

#
# DiPhot: Image Reduction - 2016-01-20
# https://github.com/viyh/diphot
#
# Create master zero, dark, and flat images, and apply them to a set of FITS cubes
# to prep for science use.
#

import diphot

if __name__ == "__main__":
    p = diphot.Reduce()
    p.auto_reduce()

