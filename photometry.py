#!/usr/bin/python

#
# DiPhot: Photometry - 2016-01-20
# https://github.com/viyh/diphot
#
# Create coordinate, relative magnitude, and data dump files from FILE images
#

import os, sys, shutil, argparse, glob
from pyraf import iraf
from collections import defaultdict
import matplotlib.pyplot as plt
import tempfile

import diphot

logger = diphot.logger_init('photometry')

class Photometry:
    def __init__(self, src_dir='data', fwhm='8', aperture='15', target_id=None):
        """
        @param src_dir: source directory of reduced FITS cubes
        @type src_dir: str
        """
        self.src_dir = src_dir
        diphot.initialize_instrument()
        self.fwhm = fwhm
        self.aperture = aperture
        self.target_id = target_id

    def run_photometry(self):
        logger.info('Creating light curve...')
        diphot.set_datapars(params=[('fwhmpsf', self.fwhm)])
        duphot.set_centerpars(params=[('cbox', self.fwhm * 2.0)])
        diphot.set_photpars(params=[('apertures', self.aperture)])
        diphot.set_fitskypars()
        diphot.set_findpars()
        diphot.cleanup_tmp(self.src_dir)
        self.create_data_files()

    def run_daofind(self, filename):
        iraf.noao.digiphot.apphot(_doprint=0)
        iraf.noao.digiphot.apphot.daofind(
            image = filename,
            output = self.src_dir + '/tmp/default',
            starmap = '',
            skymap = '',
            datapars = '',
            findpars = '',
            boundary = 'nearest',
            constant = 0,
            interactive = 'no',
            icommands = '',
            gcommands = '',
            verify = 'no'
        )

    def get_txdump(self, filemask, fields):
        iraf.noao.digiphot.ptools(_doprint=0)
        return iraf.noao.digiphot.ptools.txdump(
            textfiles = filemask,
            fields = fields,
            expr = 'yes',
            headers = 'no',
            parameters = 'yes',
            Stdout=1
        )

    def write_file_from_array(self, filename, contents):
        with open(filename,'w') as f:
            f.write("\n".join(contents))
            f.write("\n")

    def create_filelist(self, files, page):
        page_files = files[page*240:(page+1)*240]
        self.write_file_from_array('images.list', page_files)

    def do_photometry(self, filemask):
        files = sorted(glob.glob(filemask))
        for i in range(0, len(files)/240 + 1):
            self.create_filelist(files, i)
            self.run_phot('@images.list')

    def run_phot(self, filemask):
        iraf.noao.digiphot.apphot.phot(
            image = filemask,
            coords = self.src_dir + '/tmp/default',
            output = self.src_dir + '/tmp/default',
            skyfile = '',
            plotfile = '',
            datapars = '',
            centerpars = '',
            fitskypars = '',
            photpars = '',
            interactive = 'no',
            radplots = 'no',
            icommands = '',
            gcommands = '',
            verify = 'no'
        )

    def create_data_files(self):
        """Create coord file (daofind) and mag file (phot) for science images."""
        logger.info('Creating coords, mag for science images.')

        files = diphot.get_files_of_type(self.src_dir + '/object', 'object')
        self.run_daofind(self.src_dir + '/object/*.fits')
        self.do_photometry(self.src_dir + '/object/*.fits')
        mag_files = self.src_dir + '/tmp/*.mag.1'
        mag_dump = self.get_txdump(mag_files, 'IMAGE,ID,XCENTER,YCENTER,OTIME,MAG,MERR')
        self.write_file_from_array(self.src_dir + '/txdump.txt', mag_dump)

def parse_args():
    parser = argparse.ArgumentParser(description='Create growth curve')
    parser.add_argument('--src_dir', type=str, default='data',
        help='source directory of reduced FITS files')
    parser.add_argument('--fwhm', type=float, default=8.0,
        help='FWHM')
    parser.add_argument('--aperture', type=float, default=15.0,
        help='aperture size')
    parser.add_argument('--id', type=int,
        help='star number of the primary target from the coordinate file')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    os.environ.get('iraf','/usr/local/iraf')
    os.environ.get('IRAFARCH','linux64')
    l = Photometry(src_dir=args.src_dir, fwhm=args.fwhm, aperture=args.aperture, target_id=args.id)
    l.run_photometry()
