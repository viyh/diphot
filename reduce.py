#!/usr/bin/python

#
# DiPhot: Curve of Growth - 2016-01-20
# https://github.com/viyh/diphot
#
# Create master zero, dark, and flat images, and apply them to a set of FITS cubes
# to prep for science use.
#

import os, sys, shutil, argparse
from pyraf import iraf
import diphot

logger = diphot.logger_init('reduce')

class Reduce:
    def __init__(self, src_dir='data', dst_dir='output', dry_run=False):
        """
        @param src_dir: source directory of original FITS cubes
        @type src_dir: str
        @param dst_dir: destination directory where master and calibrated FITS cubes are written
        @type dst_dir: str
        """
        self.src_dir = src_dir
        self.dst_dir = dst_dir
        self.master_dir = dst_dir + '/master'
        self.tmp_dir = dst_dir + '/tmp'
        self.filetypes = ['zero', 'dark', 'flat', 'object']
        self.dry_run = dry_run
        self.initialize_dirs()
        diphot.initialize_instrument()

    def auto_reduce(self):
        logger.info('Creating IRAF FITS images...')
        self.convert_to_fits()
        self.organize_files()
        self.create_zero()
        self.apply_zero()
        self.create_dark()
        self.apply_dark()
        self.create_flat()
        self.apply_flat()

    def initialize_dirs(self):
        if not os.path.exists(self.src_dir):
            logger.info("Source directory '{}' does not exist!".format(self.src_dir))
            sys.exit(1)
        diphot.mkdir(self.dst_dir)
        diphot.mkdir(self.master_dir)
        diphot.mkdir(self.tmp_dir)
        for filetype in self.filetypes:
            diphot.mkdir(self.dst_dir + '/' + filetype)

    def convert_to_fits(self):
        iraf.dataio(_doprint=0)
        iraf.dataio.rfits(
            fits_file = self.src_dir + '/*',
            file_list = '',
            iraf_file = self.tmp_dir + '/',
            make_image = 'yes',
            long_header = 'no',
            short_header = 'yes',
            datatype = '',
            blank = 0.0,
            scale = 'yes',
            oldirafname = 'no',
            offset = 0,
            mode = 'ql',
            Stderr='dev$null'
        )

    def organize_files(self):
        for filetype in self.filetypes:
            logger.info("Moving {} files to {}".format(filetype, self.dst_dir + '/' + filetype))
            files = diphot.get_files_of_type(self.dst_dir + '/tmp', filetype)
            diphot.move_files(files, self.dst_dir + '/' + filetype)

    def create_zero(self):
        """Creates a master zero image from a set of zero FITS cubes."""
        logger.info('Creating master zero...')
        iraf.noao.imred(_doprint=0)
        iraf.noao.imred.ccdred(_doprint=0)
        iraf.noao.imred.ccdred.ccdproc.unlearn()
        iraf.noao.imred.ccdred.zerocombine.unlearn()
        iraf.noao.imred.ccdred.zerocombine(
            input = self.dst_dir + '/zero/*.fits',
            output = self.master_dir + '/masterzero',
            combine = 'average',
            reject = 'minmax',
            ccdtype = 'zero',
            process = 'no',
            delete = 'no',
            clobber = 'no',
            scale = 'none',
            rdnoise = 8.8,
            gain = 1.3
        )

    def apply_zero(self):
        """Apply master zero image to a set of FITS cubes."""
        logger.info('Applying master zero to data files.')
        iraf.noao.imred(_doprint=0)
        iraf.noao.imred.ccdred(_doprint=0)
        iraf.noao.imred.ccdred.ccdproc.unlearn()
        iraf.noao.imred.ccdred.ccdproc.output = ''
        iraf.noao.imred.ccdred.ccdproc.ccdtype = ''
        iraf.noao.imred.ccdred.ccdproc.max_cac = 0
        iraf.noao.imred.ccdred.ccdproc.noproc = 'no'
        iraf.noao.imred.ccdred.ccdproc.fixpix = 'no'
        iraf.noao.imred.ccdred.ccdproc.oversca = 'no'
        iraf.noao.imred.ccdred.ccdproc.trim = 'no'
        iraf.noao.imred.ccdred.ccdproc.zerocor = 'yes'
        iraf.noao.imred.ccdred.ccdproc.darkcor = 'no'
        iraf.noao.imred.ccdred.ccdproc.flatcor = 'no'
        iraf.noao.imred.ccdred.ccdproc.illumco = 'no'
        iraf.noao.imred.ccdred.ccdproc.fringec = 'no'
        iraf.noao.imred.ccdred.ccdproc.readcor = 'no'
        iraf.noao.imred.ccdred.ccdproc.scancor = 'no'
        iraf.noao.imred.ccdred.ccdproc.readaxis = 'line'
        iraf.noao.imred.ccdred.ccdproc.fixfile = ''
        iraf.noao.imred.ccdred.ccdproc.biassec = ''
        iraf.noao.imred.ccdred.ccdproc.trimsec = ''
        iraf.noao.imred.ccdred.ccdproc.zero = self.master_dir + '/masterzero'
        iraf.noao.imred.ccdred.ccdproc.dark = ''
        iraf.noao.imred.ccdred.ccdproc.flat = ''

        filemasks = [
            self.dst_dir + '/dark/*.fits',
            self.dst_dir + '/flat/*.fits',
            self.dst_dir + '/object/*.fits'
        ]
        for filemask in filemasks:
            iraf.noao.imred.ccdred.ccdproc(images=filemask)

    def create_dark(self):
        """Creates a master dark image from a set of dark FITS cubes."""
        logger.info('Creating master dark...')
        iraf.noao.imred(_doprint=0)
        iraf.noao.imred.ccdred(_doprint=0)
        iraf.noao.imred.ccdred.ccdproc.unlearn()
        iraf.noao.imred.ccdred.darkcombine.unlearn()
        iraf.noao.imred.ccdred.ccdproc.output = ''
        iraf.noao.imred.ccdred.ccdproc.ccdtype = ''
        iraf.noao.imred.ccdred.ccdproc.max_cac = 0
        iraf.noao.imred.ccdred.ccdproc.noproc = 'no'
        iraf.noao.imred.ccdred.ccdproc.fixpix = 'no'
        iraf.noao.imred.ccdred.ccdproc.oversca = 'no'
        iraf.noao.imred.ccdred.ccdproc.trim = 'no'
        iraf.noao.imred.ccdred.ccdproc.zerocor = 'yes'
        iraf.noao.imred.ccdred.ccdproc.darkcor = 'no'
        iraf.noao.imred.ccdred.ccdproc.flatcor = 'no'
        iraf.noao.imred.ccdred.ccdproc.illumco = 'no'
        iraf.noao.imred.ccdred.ccdproc.fringec = 'no'
        iraf.noao.imred.ccdred.ccdproc.readcor = 'no'
        iraf.noao.imred.ccdred.ccdproc.scancor = 'no'
        iraf.noao.imred.ccdred.ccdproc.readaxis = 'line'
        iraf.noao.imred.ccdred.ccdproc.fixfile = ''
        iraf.noao.imred.ccdred.ccdproc.biassec = ''
        iraf.noao.imred.ccdred.ccdproc.trimsec = ''
        iraf.noao.imred.ccdred.ccdproc.zero = self.master_dir + '/masterzero'
        iraf.noao.imred.ccdred.ccdproc.dark = ''
        iraf.noao.imred.ccdred.ccdproc.flat = ''
        iraf.noao.imred.ccdred.darkcombine(
            input = self.dst_dir + '/dark/*.fits',
            output = self.master_dir + '/masterdark',
            combine = 'median',
            reject = 'minmax',
            ccdtype = 'dark',
            process = 'yes',
            delete = 'no',
            clobber = 'no',
            scale = 'exposure',
            rdnoise = 8.8,
            gain = 1.3
        )

    def apply_dark(self):
        """Apply master dark image to a set of FITS cubes."""
        logger.info('Applying master dark to data files.')
        iraf.noao.imred(_doprint=0)
        iraf.noao.imred.ccdred(_doprint=0)
        iraf.noao.imred.ccdred.ccdproc.unlearn()
        iraf.noao.imred.ccdred.ccdproc.output = ''
        iraf.noao.imred.ccdred.ccdproc.ccdtype = ''
        iraf.noao.imred.ccdred.ccdproc.max_cac = 0
        iraf.noao.imred.ccdred.ccdproc.noproc = 'no'
        iraf.noao.imred.ccdred.ccdproc.fixpix = 'no'
        iraf.noao.imred.ccdred.ccdproc.oversca = 'no'
        iraf.noao.imred.ccdred.ccdproc.trim = 'no'
        iraf.noao.imred.ccdred.ccdproc.zerocor = 'no'
        iraf.noao.imred.ccdred.ccdproc.darkcor = 'yes'
        iraf.noao.imred.ccdred.ccdproc.flatcor = 'no'
        iraf.noao.imred.ccdred.ccdproc.illumco = 'no'
        iraf.noao.imred.ccdred.ccdproc.fringec = 'no'
        iraf.noao.imred.ccdred.ccdproc.readcor = 'no'
        iraf.noao.imred.ccdred.ccdproc.scancor = 'no'
        iraf.noao.imred.ccdred.ccdproc.readaxis = 'line'
        iraf.noao.imred.ccdred.ccdproc.fixfile = ''
        iraf.noao.imred.ccdred.ccdproc.biassec = ''
        iraf.noao.imred.ccdred.ccdproc.trimsec = ''
        iraf.noao.imred.ccdred.ccdproc.zero = ''
        iraf.noao.imred.ccdred.ccdproc.dark = self.master_dir + '/masterdark'
        iraf.noao.imred.ccdred.ccdproc.flat = ''

        filemasks = [
            self.dst_dir + '/flat/*.fits',
            self.dst_dir + '/object/*.fits'
        ]
        for filemask in filemasks:
            iraf.noao.imred.ccdred.ccdproc(images=filemask)

    def create_flat(self):
        """Creates a master flat image from a set of flat FITS cubes."""
        logger.info('Creating master flat...')
        iraf.noao.imred(_doprint=0)
        iraf.noao.imred.ccdred(_doprint=0)
        iraf.noao.imred.ccdred.ccdproc.unlearn()
        iraf.noao.imred.ccdred.flatcombine.unlearn()
        iraf.noao.imred.ccdred.ccdproc.output = ''
        iraf.noao.imred.ccdred.ccdproc.ccdtype = ''
        iraf.noao.imred.ccdred.ccdproc.max_cac = 0
        iraf.noao.imred.ccdred.ccdproc.noproc = 'no'
        iraf.noao.imred.ccdred.ccdproc.fixpix = 'no'
        iraf.noao.imred.ccdred.ccdproc.oversca = 'no'
        iraf.noao.imred.ccdred.ccdproc.trim = 'no'
        iraf.noao.imred.ccdred.ccdproc.zerocor = 'no'
        iraf.noao.imred.ccdred.ccdproc.darkcor = 'no'
        iraf.noao.imred.ccdred.ccdproc.flatcor = 'no'
        iraf.noao.imred.ccdred.ccdproc.illumco = 'no'
        iraf.noao.imred.ccdred.ccdproc.fringec = 'no'
        iraf.noao.imred.ccdred.ccdproc.readcor = 'no'
        iraf.noao.imred.ccdred.ccdproc.scancor = 'no'
        iraf.noao.imred.ccdred.ccdproc.readaxis = 'line'
        iraf.noao.imred.ccdred.ccdproc.fixfile = ''
        iraf.noao.imred.ccdred.ccdproc.biassec = ''
        iraf.noao.imred.ccdred.ccdproc.trimsec = ''
        iraf.noao.imred.ccdred.ccdproc.zero = ''
        iraf.noao.imred.ccdred.ccdproc.dark = ''
        iraf.noao.imred.ccdred.ccdproc.flat = ''
        iraf.noao.imred.ccdred.flatcombine(
            input = self.dst_dir + '/flat/*.fits',
            output = self.master_dir + '/masterflat',
            combine = 'median',
            reject = 'avsigclip',
            ccdtype = 'flat',
            subsets = 'yes',
            process = 'yes',
            delete = 'no',
            clobber = 'no',
            statsec = '[700:1400:400:1100]',
            scale = 'median',
            rdnoise = 8.8,
            gain = 1.3
        )

    def apply_flat(self):
        """Apply master flat image to a set of FITS cubes."""
        logger.info('Applying master flat to data files.')
        iraf.noao.imred(_doprint=0)
        iraf.noao.imred.ccdred(_doprint=0)
        iraf.noao.imred.ccdred.ccdproc.unlearn()
        iraf.noao.imred.ccdred.ccdproc.output = ''
        iraf.noao.imred.ccdred.ccdproc.ccdtype = ''
        iraf.noao.imred.ccdred.ccdproc.max_cac = 0
        iraf.noao.imred.ccdred.ccdproc.noproc = 'no'
        iraf.noao.imred.ccdred.ccdproc.fixpix = 'no'
        iraf.noao.imred.ccdred.ccdproc.oversca = 'no'
        iraf.noao.imred.ccdred.ccdproc.trim = 'no'
        iraf.noao.imred.ccdred.ccdproc.zerocor = 'no'
        iraf.noao.imred.ccdred.ccdproc.darkcor = 'no'
        iraf.noao.imred.ccdred.ccdproc.flatcor = 'yes'
        iraf.noao.imred.ccdred.ccdproc.illumco = 'no'
        iraf.noao.imred.ccdred.ccdproc.fringec = 'no'
        iraf.noao.imred.ccdred.ccdproc.readcor = 'no'
        iraf.noao.imred.ccdred.ccdproc.scancor = 'no'
        iraf.noao.imred.ccdred.ccdproc.readaxis = 'line'
        iraf.noao.imred.ccdred.ccdproc.fixfile = ''
        iraf.noao.imred.ccdred.ccdproc.biassec = ''
        iraf.noao.imred.ccdred.ccdproc.trimsec = ''
        iraf.noao.imred.ccdred.ccdproc.zero = ''
        iraf.noao.imred.ccdred.ccdproc.dark = ''
        iraf.noao.imred.ccdred.ccdproc.flat = self.master_dir + '/masterflat*'

        filemasks = [
            self.dst_dir + '/object/*.fits'
        ]
        for filemask in filemasks:
            iraf.noao.imred.ccdred.ccdproc(images=filemask)

def parse_args():
    parser = argparse.ArgumentParser(description='Pre-process FITS cubes')
    parser.add_argument('--src_dir', type=str, default='data',
        help='source directory of raw FITS files')
    parser.add_argument('--dst_dir', type=str, default='output',
        help='destination directory for processed and master FITS files')
    parser.add_argument('--dry_run', action='store_true',
        help='don\'t actually do anything (for testing)')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    os.environ.get('iraf','/usr/local/iraf')
    os.environ.get('IRAFARCH','linux64')
    p = Reduce(src_dir=args.src_dir, dst_dir=args.dst_dir, dry_run=args.dry_run)
    p.auto_reduce()
