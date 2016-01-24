#!/usr/bin/python

#
# DiPhot: Full Image Processing - 2016-01-23
# https://github.com/viyh/diphot
#
# Reduce images, generate curve of growth, create CSV dump
#

import diphot

if __name__ == "__main__":
    r = diphot.Reduce()
    r.auto_reduce()

    g= diphot.CurveOfGrowth()
    g.create_curve()

    l = diphot.Photometry()
    l.fwhm = g.fwhm
    l.aperture = g.max_snr_aperture
    l.run_photometry()

    t = diphot.TxdumpParse()
    t.parse()
