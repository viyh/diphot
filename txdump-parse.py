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
    t.parse()
