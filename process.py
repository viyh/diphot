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
    r.process()

    g= diphot.CurveOfGrowth()
    g.process()

    p = diphot.Photometry()
    p.fwhm = g.fwhm
    p.aperture = g.max_snr_aperture
    p.process()

    t = diphot.TxdumpParse()
    t.process()

    for star in t.data:
        print "\n" + str(star)
        print star.data[0]

    l = diphot.LightCurve(target_id=t.found_target, data=t.data)
    l.process()
