#!/usr/bin/python

#
# PyRAF FITS pre-processing script - 2015-10-15
# https://github.com/viyh/pyraf-scripts
#
# PyRAF script to create master bias, dark, and flat, and apply them to a set of FITS cubes
# to prep FITS files for science use.
#

import os, sys, shutil
from pyraf import iraf

class PreProcess:
    def __init__(self, src_dir='data', dst_dir='output'):
        """
        @param src_dir: source directory of original FITS cubes
        @type src_dir: str
        @param dst_dir: destination directory where master and calibrated FITS cubes are written
        @type dst_dir: str
        """
        self.src_dir = src_dir
        self.dst_dir = dst_dir
        self.master_dir = dst_dir + '/master'
        self.initialize_dirs()
        self.initialize_instrument()
        iraf.noao.imred(_doprint=0)
        iraf.noao.imred.ccdred(_doprint=0)
        iraf.dataio(_doprint=0)

    def mkdir(self, dirname):
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        else:
            answer = raw_input("Directory '{}' already exists. Delete all files? (y/N) ".format(dirname)).lower() == 'y'
            if answer:
                shutil.rmtree(dirname, ignore_errors=True)

    def initialize_dirs(self):
        if not os.path.exists(self.src_dir):
            print("Source directory '{}' does not exist!".format(self.src_dir))
            sys.exit(1)
        self.mkdir(self.dst_dir)
        self.mkdir(self.master_dir)

    def initialize_instrument(self):
        inst_file_name = 'cp.dat'
        if not os.path.exists(inst_file_name):
            with open(inst_file_name, 'w+') as inst_file:
                inst_file.write('subset          FILTER\n\n')
                inst_file.write('darktime        EXPTIME\n\n')
                inst_file.write('\'Dark Frame\'    dark\n')
                inst_file.write('\'Bias Frame\'    zero\n')
                inst_file.write('\'Light Frame\'   object\n')
                inst_file.write('\'Flat Field\'    flat\n')
        iraf.noao.imred.ccdred.setinst(
            instrument='cp',
            review='no',
            mode='h',
            dir='',
            site=''
        )

    def convert_to_fits(self):
        iraf.dataio.rfits(
            fits_file = self.src_dir + '/*',
            file_list = '',
            iraf_file = self.dst_dir + '/t',
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

    def auto_prep(self):
        print('Creating IRAF FITS images...')
        self.convert_to_fits()
        self.create_zero()
        self.apply_zero()
        self.create_dark()
        self.apply_dark()
        self.create_flat()
        self.apply_flat()

    def create_zero(self):
        """Creates a master bias image from a set of bias FITS cubes."""
        print('Creating master bias...')
        iraf.noao.imred.ccdred.zerocombine.unlearn()
        iraf.noao.imred.ccdred.zerocombine(
            input = self.dst_dir + '/t*.fits',
            output = self.master_dir + '/masterbias',
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
        """Apply master bias image to a set of FITS cubes."""
        print('Applying master bias to data files.')
        iraf.noao.imred.ccdred.ccdproc(
            images = self.dst_dir + '/t*.fits',
            output = '',
            ccdtype = '',
            max_cac = 0,
            noproc = 'no',
            fixpix = 'no',
            oversca = 'no',
            trim = 'no',
            zerocor = 'yes',
            darkcor = 'no',
            flatcor = 'no',
            illumco = 'no',
            fringec = 'no',
            readcor = 'no',
            scancor = 'no',
            readaxis = 'line',
            fixfile = '',
            biassec = '',
            trimsec = '',
            zero = self.master_dir + '/masterbias',
            dark = '',
            flat = ''
        )

    def create_dark(self):
        """Creates a master dark image from a set of dark FITS cubes."""
        print('Creating master dark...')
        iraf.noao.imred.ccdred.darkcombine.unlearn()
        iraf.noao.imred.ccdred.darkcombine(
            input = self.dst_dir + '/t*.fits',
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
        print('Applying master dark to data files.')
        iraf.noao.imred.ccdred.ccdproc(
            images = self.dst_dir + '/t*.fits',
            output = '',
            ccdtype = '',
            max_cac = 0,
            noproc = 'no',
            fixpix = 'no',
            oversca = 'no',
            trim = 'no',
            zerocor = 'no',
            darkcor = 'yes',
            flatcor = 'no',
            illumco = 'no',
            fringec = 'no',
            readcor = 'no',
            scancor = 'no',
            readaxis = 'line',
            fixfile = '',
            biassec = '',
            trimsec = '',
            zero = '',
            dark = self.master_dir + '/masterdark',
            flat = ''
        )

    def create_flat(self):
        """Creates a master flat image from a set of flat FITS cubes."""
        print('Creating master flat...')
        iraf.noao.imred.ccdred.flatcombine.unlearn()
        iraf.noao.imred.ccdred.flatcombine(
            input = self.dst_dir + '/t*.fits',
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
        print('Applying master flat to data files.')
        iraf.noao.imred.ccdred.ccdproc(
            images = self.dst_dir + '/t*.fits',
            output = '',
            ccdtype = '',
            max_cac = 0,
            noproc = 'no',
            fixpix = 'no',
            oversca = 'no',
            trim = 'no',
            zerocor = 'no',
            darkcor = 'no',
            flatcor = 'yes',
            illumco = 'no',
            fringec = 'no',
            readcor = 'no',
            scancor = 'no',
            readaxis = 'line',
            fixfile = '',
            biassec = '',
            trimsec = '',
            zero = '',
            dark = '',
            flat = self.master_dir + '/masterflat*',
        )

if __name__ == "__main__":
    p = PreProcess()
    p.auto_prep()
