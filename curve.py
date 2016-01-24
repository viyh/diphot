#!/usr/bin/python

#
# DiPhot: Curve of Growth - 2016-01-20
# https://github.com/viyh/diphot
#
# Create a curve of growth that analyzes reduced FITS files.
#

import diphot

if __name__ == "__main__":
    p = diphot.CurveOfGrowth()
    p.process()
