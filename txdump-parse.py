#!/usr/bin/python

#
# DiPhot: Txdump Parse - 2016-01-20
# https://github.com/viyh/diphot
#
# Parse a dump file and attempt to align stars by ID
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
