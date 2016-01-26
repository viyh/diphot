#!/usr/bin/python

#
# DiPhot: Light Curve - 2016-01-25
# https://github.com/viyh/diphot
#
# Generate a light curve
#

import diphot

if __name__ == '__main__':
    t = diphot.TxdumpParse()
    t.process()

    for star in t.data:
        print "\n" + str(star)
        print star.data[0]

    l = diphot.LightCurve(target_id=7, data=t.data)
    l.process()
