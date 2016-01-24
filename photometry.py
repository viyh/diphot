#!/usr/bin/python

#
# DiPhot: Photometry - 2016-01-20
# https://github.com/viyh/diphot
#
# Create coordinate, relative magnitude, and data dump files from FILE images
#

import diphot

if __name__ == "__main__":
    l = diphot.Photometry()
    l.process()
