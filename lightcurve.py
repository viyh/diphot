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

    l = diphot.LightCurve(target_id=t.found_target, data=t.data)
    l.process()
