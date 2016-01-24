#!/usr/bin/python

#
# DiPhot: Library - 2016-01-20
# https://github.com/viyh/diphot
#
# Library of common functions
#

import logging, tempfile, csv
import os, shutil, glob, argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from pyraf import iraf
from collections import defaultdict, OrderedDict

class DiPhot():
    def __init__(self, name):
        os.environ.get('iraf','/usr/local/iraf')
        os.environ.get('IRAFARCH','linux64')
        self.name = name
        self.args = self.init_parse_args()
        self.debug = self.args.debug
        self.raw_dir = self.args.raw_dir
        self.output_dir = self.args.output_dir
        self.dry_run = self.args.dry_run
        self.logger = self.logger_init(self.name)
        self.pyraf = Pyraf(self.logger, self.debug)

    def init_parse_args(self):
        self.parser = argparse.ArgumentParser(description='DiPhot - Differential Photometry')
        self.parser.add_argument('--raw_dir', type=str, default='input',
            help='source directory of raw FITS files')
        self.parser.add_argument('--output_dir', type=str, default='output',
            help='destination directory for processed and master FITS files')
        self.parser.add_argument('--dry_run', action='store_true',
            help='don\'t actually do anything (for testing)')
        self.parser.add_argument('--debug', action='store_true',
            help='turn on debug messages')
        self.parse_args()
        return self.parser.parse_known_args()[0]

    def parse_args(self):
        raise("Argument parser not implemented!")

    def logger_init(self, log_name):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console = logging.StreamHandler()
        console.setFormatter(formatter)
        handler = logging.FileHandler(log_name + '.log')
        handler.setFormatter(formatter)
        logger.addHandler(console)
        logger.addHandler(handler)
        return logger

    def cleanup_tmp(self, target_dir):
        files = glob.glob(target_dir + '/tmp/*')
        for f in files:
            os.remove(f)

    def write_file_from_array(self, filename, contents):
        with open(filename,'w') as f:
            f.write("\n".join(contents))
            f.write("\n")

    def move_files(self, files, target_dir):
        for f in files:
            shutil.move(f, target_dir)

    def mkdir(self, dirname):
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        else:
            answer = raw_input("Directory '{}' already exists. Delete all files? (y/N) ".format(dirname)).lower() == 'y'
            if answer:
                shutil.rmtree(dirname, ignore_errors=True)

class Pyraf():
    def __init__(self, logger, debug=False):
        self.debug = debug
        self.logger = logger
        self.initialize_instrument()

    def initialize_instrument(self):
        inst_file_name = 'cp.dat'
        iraf.noao.imred(_doprint=0)
        iraf.noao.imred.ccdred(_doprint=0)
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

    def get_files_of_type(self, src_dir, filetype):
        iraf.noao(_doprint=0)
        return iraf.noao.imred.ccdred.ccdlist(
            images = src_dir + '/*.fits',
            ccdtype = filetype,
            names = 'yes',
            Stdout=1
        )

    def set_datapars(self, params=[]):
        """Set datapars parameters for photometry."""
        iraf.noao.digiphot(_doprint=0)
        iraf.noao.digiphot.apphot(_doprint=0)
        iraf.noao.digiphot.apphot.datapars.setParam('scale', '1.0')
        iraf.noao.digiphot.apphot.datapars.setParam('fwhmpsf', '8')
        iraf.noao.digiphot.apphot.datapars.setParam('emission', 'yes')
        iraf.noao.digiphot.apphot.datapars.setParam('sigma', '60')
        iraf.noao.digiphot.apphot.datapars.setParam('datamin', '0')
        iraf.noao.digiphot.apphot.datapars.setParam('datamax', '60000')
        iraf.noao.digiphot.apphot.datapars.setParam('noise', 'poisson')
        iraf.noao.digiphot.apphot.datapars.setParam('ccdread', '')
        iraf.noao.digiphot.apphot.datapars.setParam('gain', '')
        iraf.noao.digiphot.apphot.datapars.setParam('readnoise', '8.8')
        iraf.noao.digiphot.apphot.datapars.setParam('epadu', '1.3')
        iraf.noao.digiphot.apphot.datapars.setParam('exposur', 'EXPTIME')
        iraf.noao.digiphot.apphot.datapars.setParam('airmass', 'AIRMASS')
        iraf.noao.digiphot.apphot.datapars.setParam('filter', 'FILTER')
        iraf.noao.digiphot.apphot.datapars.setParam('obstime', 'DATE-OBS')
        iraf.noao.digiphot.apphot.datapars.setParam('itime', 'INDEF')
        iraf.noao.digiphot.apphot.datapars.setParam('xairmass', 'INDEF')
        iraf.noao.digiphot.apphot.datapars.setParam('ifilter', 'INDEF')
        iraf.noao.digiphot.apphot.datapars.setParam('otime', 'INDEF')
        for param in params:
            iraf.noao.digiphot.apphot.datapars.setParam(param[0], param[1])

    def set_centerpars(self, params=[]):
        """Set centerpars parameters for photometry."""
        iraf.noao.digiphot(_doprint=0)
        iraf.noao.digiphot.apphot(_doprint=0)
        iraf.noao.digiphot.apphot.centerpars.setParam('calgorithm', 'centroid')
        iraf.noao.digiphot.apphot.centerpars.setParam('cbox', '16')
        iraf.noao.digiphot.apphot.centerpars.setParam('cthreshold', '0')
        iraf.noao.digiphot.apphot.centerpars.setParam('minsnratio', '1')
        iraf.noao.digiphot.apphot.centerpars.setParam('cmaxiter', '10')
        iraf.noao.digiphot.apphot.centerpars.setParam('maxshift', '3')
        iraf.noao.digiphot.apphot.centerpars.setParam('clean', 'no')
        iraf.noao.digiphot.apphot.centerpars.setParam('mkcente', 'yes')
        for param in params:
            iraf.noao.digiphot.apphot.centerpars.setParam(param[0], param[1])

    def set_photpars(self, params=[]):
        """Set photpars parameters for photometry."""
        iraf.noao.digiphot(_doprint=0)
        iraf.noao.digiphot.apphot(_doprint=0)
        iraf.noao.digiphot.apphot.photpars.setParam('weighting', 'constant')
        iraf.noao.digiphot.apphot.photpars.setParam('apertures', '15')
        iraf.noao.digiphot.apphot.photpars.setParam('zmag', '25')
        iraf.noao.digiphot.apphot.photpars.setParam('mkapert', 'yes')
        for param in params:
            iraf.noao.digiphot.apphot.photpars.setParam(param[0], param[1])

    def set_fitskypars(self, params=[]):
        """Set fitskypars parameters for photometry."""
        iraf.noao.digiphot(_doprint=0)
        iraf.noao.digiphot.apphot(_doprint=0)
        iraf.noao.digiphot.apphot.fitskypars.setParam('salgorithm', 'centroid')
        iraf.noao.digiphot.apphot.fitskypars.setParam('annulus', '25')
        iraf.noao.digiphot.apphot.fitskypars.setParam('dannulus', '5')
        iraf.noao.digiphot.apphot.fitskypars.setParam('skyvalue', '0')
        iraf.noao.digiphot.apphot.fitskypars.setParam('smaxiter', '10')
        iraf.noao.digiphot.apphot.fitskypars.setParam('sloclip', '0')
        iraf.noao.digiphot.apphot.fitskypars.setParam('shiclip', '0')
        iraf.noao.digiphot.apphot.fitskypars.setParam('snreject', '50')
        iraf.noao.digiphot.apphot.fitskypars.setParam('sloreject', '3')
        iraf.noao.digiphot.apphot.fitskypars.setParam('shireject', '3')
        iraf.noao.digiphot.apphot.fitskypars.setParam('khist', '3')
        iraf.noao.digiphot.apphot.fitskypars.setParam('binsize', '0.1')
        iraf.noao.digiphot.apphot.fitskypars.setParam('smooth', 'no')
        iraf.noao.digiphot.apphot.fitskypars.setParam('rgrow', '0')
        iraf.noao.digiphot.apphot.fitskypars.setParam('mksky', 'yes')
        for param in params:
            iraf.noao.digiphot.apphot.fitskypars.setParam(param[0], param[1])

    def set_findpars(self, params=[]):
        """Set findpars parameters for photometry."""
        iraf.noao.digiphot(_doprint=0)
        iraf.noao.digiphot.apphot(_doprint=0)
        iraf.noao.digiphot.apphot.findpars.setParam('threshold', '10')
        iraf.noao.digiphot.apphot.findpars.setParam('nsigma', '1.5')
        iraf.noao.digiphot.apphot.findpars.setParam('ratio', '1')
        iraf.noao.digiphot.apphot.findpars.setParam('theta', '0')
        iraf.noao.digiphot.apphot.findpars.setParam('sharplo', '-4')
        iraf.noao.digiphot.apphot.findpars.setParam('sharphi', '4')
        iraf.noao.digiphot.apphot.findpars.setParam('roundlo', '-4')
        iraf.noao.digiphot.apphot.findpars.setParam('roundhi', '4')
        iraf.noao.digiphot.apphot.findpars.setParam('mkdetections', 'yes')
        for param in params:
            iraf.noao.digiphot.apphot.findpars.setParam(param[0], param[1])

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

    def run_rfits(self, raw_dir, output_dir):
        iraf.dataio(_doprint=0)
        iraf.dataio.rfits(
            fits_file = raw_dir + '/*',
            file_list = '',
            iraf_file = output_dir + '/tmp/',
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

    def run_daofind(self, output_dir, filename):
        iraf.noao.digiphot.apphot(_doprint=0)
        iraf.noao.digiphot.apphot.daofind(
            image = filename,
            output = output_dir + '/tmp/default',
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

    def run_phot(self, output_dir, filemask):
        iraf.noao.digiphot.apphot.phot(
            image = filemask,
            coords = output_dir + '/tmp/default',
            output = output_dir + '/tmp/default',
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

    def run_psfmeasure(self, filename, coords):
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        for coord in coords:
            c = coords[coord]
            if self.debug:
                self.logger.debug('\t{}: {}, {}'.format(coord, c['x'][0], c['y'][0]))
            tmp_file.write('{}, {}, {}\n'.format(c['x'][0], c['y'][0], coord))
        tmp_file.close()
        iraf.noao.obsutil(_doprint=0)
        psf_out = iraf.noao.obsutil.psfmeasure(
            images = filename,
            iterations = 1,
            logfile = '',
            imagecur = tmp_file.name,
            display = 'no',
            Stdout = 1,
            Stderr = '/dev/null',
            StdoutG = '/dev/null'
        )
        os.unlink(tmp_file.name)
        fwhm = float(psf_out[-1].split('width at half maximum (FWHM) of ')[-1])
        self.logger.info('Found average FWHM of {} stars in image {}: {}'.format(len(coords.keys()), filename, fwhm))
        return fwhm

    def create_zero(self, output_dir):
        """Creates a master zero image from a set of zero FITS cubes."""
        self.logger.info('Creating master zero...')
        iraf.noao.imred(_doprint=0)
        iraf.noao.imred.ccdred(_doprint=0)
        iraf.noao.imred.ccdred.ccdproc.unlearn()
        iraf.noao.imred.ccdred.zerocombine.unlearn()
        iraf.noao.imred.ccdred.zerocombine(
            input = output_dir + '/zero/*.fits',
            output = output_dir + '/master/masterzero',
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

    def apply_zero(self, output_dir):
        """Apply master zero image to a set of FITS cubes."""
        self.logger.info('Applying master zero to data files.')
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
        iraf.noao.imred.ccdred.ccdproc.zero = output_dir + '/master/masterzero'
        iraf.noao.imred.ccdred.ccdproc.dark = ''
        iraf.noao.imred.ccdred.ccdproc.flat = ''

        filemasks = [
            output_dir + '/dark/*.fits',
            output_dir + '/flat/*.fits',
            output_dir + '/object/*.fits'
        ]
        for filemask in filemasks:
            iraf.noao.imred.ccdred.ccdproc(images=filemask)

    def create_dark(self, output_dir):
        """Creates a master dark image from a set of dark FITS cubes."""
        self.logger.info('Creating master dark...')
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
        iraf.noao.imred.ccdred.ccdproc.zero = output_dir + '/master/masterzero'
        iraf.noao.imred.ccdred.ccdproc.dark = ''
        iraf.noao.imred.ccdred.ccdproc.flat = ''
        iraf.noao.imred.ccdred.darkcombine(
            input = output_dir + '/dark/*.fits',
            output = output_dir + '/master/masterdark',
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

    def apply_dark(self, output_dir):
        """Apply master dark image to a set of FITS cubes."""
        self.logger.info('Applying master dark to data files.')
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
        iraf.noao.imred.ccdred.ccdproc.dark = output_dir + '/master/masterdark'
        iraf.noao.imred.ccdred.ccdproc.flat = ''

        filemasks = [
            output_dir + '/flat/*.fits',
            output_dir + '/object/*.fits'
        ]
        for filemask in filemasks:
            iraf.noao.imred.ccdred.ccdproc(images=filemask)

    def create_flat(self, output_dir):
        """Creates a master flat image from a set of flat FITS cubes."""
        self.logger.info('Creating master flat...')
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
            input = output_dir + '/flat/*.fits',
            output = output_dir + '/master/masterflat',
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

    def apply_flat(self, output_dir):
        """Apply master flat image to a set of FITS cubes."""
        self.logger.info('Applying master flat to data files.')
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
        iraf.noao.imred.ccdred.ccdproc.flat = output_dir + '/master/masterflat*'

        filemasks = [
            output_dir + '/object/*.fits'
        ]
        for filemask in filemasks:
            iraf.noao.imred.ccdred.ccdproc(images=filemask)

class Reduce(DiPhot):
    def __init__(self, raw_dir='data', output_dir='output', dry_run=False):
        """
        @param raw_dir: source directory of original FITS cubes
        @type raw_dir: str
        @param output_dir: destination directory where master and calibrated FITS cubes are written
        @type output_dir: str
        """
        DiPhot.__init__(self, 'reduce')
        self.filetypes = ['zero', 'dark', 'flat', 'object']
        self.initialize_dirs()

    def parse_args(self):
        pass

    def auto_reduce(self):
        self.logger.info('Creating IRAF FITS images...')
        self.pyraf.run_rfits(self.raw_dir, self.output_dir)
        self.organize_files()
        self.pyraf.create_zero(self.output_dir)
        self.pyraf.apply_zero(self.output_dir)
        self.pyraf.create_dark(self.output_dir)
        self.pyraf.apply_dark(self.output_dir)
        self.pyraf.create_flat(self.output_dir)
        self.pyraf.apply_flat(self.output_dir)

    def initialize_dirs(self):
        if not os.path.exists(self.raw_dir):
            self.logger.info("Source directory '{}' does not exist!".format(self.raw_dir))
            sys.exit(1)
        self.mkdir(self.output_dir)
        self.mkdir(self.output_dir + '/master')
        self.mkdir(self.output_dir + '/tmp')
        for filetype in self.filetypes:
            self.mkdir(self.output_dir + '/' + filetype)

    def organize_files(self):
        for filetype in self.filetypes:
            self.logger.info("Moving {} files to {}".format(filetype, self.output_dir + '/' + filetype))
            files = self.pyraf.get_files_of_type(self.output_dir + '/tmp', filetype)
            self.move_files(files, self.output_dir + '/' + filetype)

class CurveOfGrowth(DiPhot):
    def __init__(self):
        DiPhot.__init__(self, 'curveofgrowth')
        self.display = self.args.display
        self.graph_file = self.args.graph_file

    def create_curve(self):
        self.logger.info('Creating curve of growth...')
        self.cleanup_tmp(self.output_dir)
        self.pyraf.set_datapars()
        self.pyraf.set_centerpars()
        self.pyraf.set_photpars()
        self.pyraf.set_fitskypars()
        self.pyraf.set_findpars()
        self.create_data_files()
        self.cleanup_tmp(self.output_dir)

    def parse_txdump(self, dump, fields):
        fields.pop(0)
        parsed = defaultdict(list)
        for line in dump:
            arr = line.split()
            data_id = arr.pop(0)
            line_parsed = {}
            if len(arr) % len(fields) != 0: continue
            data_points = len(arr) / len(fields)
            for i in range(0, (len(fields))):
                line_parsed[fields[i]] = arr[i*data_points:(i+1)*data_points]
            parsed[data_id] = line_parsed
        return parsed

    def get_data_points(self, mags):
        max_snr_x, max_snr = None, 0
        for i in sorted(mags, key=int):
            v = mags[str(i)]
            if 'INDEF' in v['merr'] or 'INDEF' in v['flux']: continue
            x = map(lambda x: float(x), v['aperture'])
            y1 = map(lambda y1: 1.0/float(y1), v['merr'])
            y2 = map(lambda y2: float(y2), v['flux'])
            max_snr = 1.0 / float(max(v['merr']))
            max_snr_x = v['merr'].index(min(v['merr'])) + 1
            self.logger.info('Max SNR [{:.2f}], aperture {}px'.format(max_snr, max_snr_x))
            self.max_snr = max_snr
            self.max_snr_aperture = max_snr_x
            return max_snr_x, x, y1, y2
        self.logger.info('Could not find stars with complete data set!'.format())
        sys.exit()

    def create_plot(self, x, y1, y2, max_x):
        fig, ax1 = plt.subplots()

        max_x = float(max_x)
        x_np = np.array(x)
        y1_np = np.array(y1)
        y2_np = np.array(y2)

        x_smooth = np.linspace(x_np.min(), x_np.max(), 200)
        y1_smooth = interp1d(x_np, y1_np, kind='cubic')(x_smooth)
        y2_smooth = interp1d(x_np, y2_np, kind='cubic')(x_smooth)

        ax1.plot(x_smooth, y1_smooth, 'b-', lw=2)
        ax1.set_xlabel('Aperture (px)')
        ax1.set_ylabel('SNR', color='b')
        ax1.set_ylim([0, max(y1) + max(y1)/20])

        ax1.axvline(max_x, color='k', linestyle='--')
        ax1.text(max_x + 1, 2, 'Aperture with max SNR: ' + str(max_x) + 'px')

        ax2 = ax1.twinx()
        ax2.plot(x_smooth, y2_smooth, 'r-', lw=2)
        ax2.set_ylabel('Flux', color='r')
        ax2.set_ylim([0, max(y2) + max(y2)/20])

        plt.subplots_adjust(bottom=.13, left=.13, right=.85, top=.95)
        if self.display:
            plt.show()
        if self.graph_file:
            self.logger.info('Creating curve of growth graph image: {}'.format(self.graph_file))
            plt.savefig(self.graph_file)

    def create_data_files(self):
        """Create coord file (daofind) and mag file (phot) for first science image."""
        self.logger.info('Creating coords, mag for first science image.')

        files = self.pyraf.get_files_of_type(self.output_dir + '/object', 'object')
        file_num = int(len(files) / 2)
        self.pyraf.run_daofind(self.output_dir, files[file_num])

        coord_files = self.output_dir + '/tmp/' + os.path.basename(files[file_num]) + '.coo.1'
        coord_dump = self.pyraf.get_txdump(coord_files, 'ID,XCENTER,YCENTER')
        coords = self.parse_txdump(coord_dump, ['id', 'x', 'y'])

        self.fwhm = self.pyraf.run_psfmeasure(files[file_num], coords)
        iraf.noao.digiphot.apphot.datapars.setParam('fwhmpsf', self.fwhm)
        iraf.noao.digiphot.apphot.centerpars.setParam('cbox', self.fwhm*2.0)

        iraf.noao.digiphot.apphot.photpars.setParam('apertures', '1:50:1')
        self.pyraf.run_phot(self.output_dir, files[file_num])

        mag_files = self.output_dir + '/tmp/' + os.path.basename(files[file_num]) + '.mag.*'
        mag_dump = self.pyraf.get_txdump(mag_files, 'ID,RAPERT,MAG,MERR,FLUX')
        mags = self.parse_txdump(mag_dump, ['id', 'aperture', 'mag', 'merr', 'flux'])

        max_snr_aperture, x, y1, y2 = self.get_data_points(mags)
        self.create_plot(x, y1, y2, max_snr_aperture)

    def parse_args(self):
        self.parser.add_argument('--display', action='store_true',
            help='display graph')
        self.parser.add_argument('--graph_file', type=str, default=None,
            help='output graph file name')

class Photometry(DiPhot):
    def __init__(self):
        DiPhot.__init__(self, 'photometry')
        self.fwhm = self.args.fwhm
        self.aperture = self.args.aperture

    def run_photometry(self):
        self.logger.info('Creating light curve...')
        self.pyraf.set_datapars(params=[('fwhmpsf', self.fwhm)])
        self.pyraf.set_centerpars(params=[('cbox', self.fwhm * 2.0)])
        self.pyraf.set_photpars(params=[('apertures', self.aperture)])
        self.pyraf.set_fitskypars()
        self.pyraf.set_findpars()
        self.cleanup_tmp(self.output_dir)
        self.create_data_files()

    def create_filelist(self, files, page):
        page_files = files[page*240:(page+1)*240]
        self.write_file_from_array(self.output_dir + '/images.list', page_files)

    def do_photometry(self, filemask):
        files = sorted(glob.glob(filemask))
        for i in range(0, len(files)/240 + 1):
            self.create_filelist(files, i)
            self.pyraf.run_phot(self.output_dir, '@' + self.output_dir + '/images.list')

    def create_data_files(self):
        """Create coord file (daofind) and mag file (phot) for science images."""
        self.logger.info('Creating coords, mag for science images.')
        self.pyraf.run_daofind(self.output_dir, self.output_dir + '/object/*.fits')
        self.do_photometry(self.output_dir + '/object/*.fits')
        mag_files = self.output_dir + '/tmp/*.mag.1'
        mag_dump = self.pyraf.get_txdump(mag_files, 'IMAGE,ID,XCENTER,YCENTER,OTIME,MAG,MERR')
        self.write_file_from_array(self.output_dir + '/txdump.txt', mag_dump)

    def parse_args(self):
        self.parser.add_argument('--src_dir', type=str, default='data',
            help='source directory of reduced FITS files')
        self.parser.add_argument('--fwhm', type=float, default=8.0,
            help='FWHM')
        self.parser.add_argument('--aperture', type=float, default=15.0,
            help='aperture size')

class TxdumpParse(DiPhot):
    def __init__(self):
        DiPhot.__init__(self, 'curveofgrowth')
        self.px_threshold = 75.0
        self.mag_threshold = 0.8
        self.skip_px_threshold = 90.0
        self.skip_mag_threshold = 1.0
        self.assume = False

    def parse_args(self):
        pass

    def parse(self):
        dump_file = self.output_dir + '/txdump.txt'
        dump = self.read_dump(dump_file)
        dump = self.sort_dump(dump)

        print
        images = set()
        for star in dump:
            images |= set(dump[star].keys())

        final_dump = defaultdict(OrderedDict)

        for star, values in dump.iteritems():
            star_images = set(values.keys())
            if set(images).difference(star_images):
                print "Not a complete data set: {} [missing {} datapoints out of {}]".format(star, len(set(images).difference(star_images)), len(images))
            else:
                print "Complete data set for star [{}]!".format(star)

            if float( len(set(images).difference(star_images)) ) / float ( len(images) ) * 100.0 < 3:
                print "Saving star [{}]!".format(star)
                final_dump[star] = values

        self.write_csv(images, final_dump)

    def read_dump(self, dump_filename):
        dump = defaultdict(list)
        with open(dump_filename) as dumpfile:
            for line in dumpfile.readlines():
                image, id, x, y, time, mag, merr = line.split()
                dump[image].append({'time': time, 'x': float(x), 'y': float(y), 'mag': mag, 'merr': merr})
        return OrderedDict(sorted(dump.items()))

    def px_test(self, c1, c2):
        return abs(c1 - c2) < self.px_threshold

    def skip_px_test(self, c1, c2):
        return abs(c1 - c2) > self.skip_px_threshold

    def mag_test(self, star1, star2):
        if not star1.has_key('mag') or not star2.has_key('mag'):
            return True
        return star1['mag'].strip() == 'INDEF' or \
            star2['mag'].strip() == 'INDEF' or \
            abs(float(star1['mag']) - float(star2['mag'])) < self.mag_threshold

    def skip_mag_test(self, star1, star2):
        if not star1.has_key('mag') or not star2.has_key('mag'):
            return False
        return star1['mag'].strip() != 'INDEF' and \
            star2['mag'].strip() != 'INDEF' and \
            abs(float(star1['mag']) - float(star2['mag'])) > self.skip_mag_threshold

    def find_similar_star(self, sorted_dump, star, image):
        if not sorted_dump:
            return False
        for star_id, star_data in sorted_dump.iteritems():
            last_image, last_data = self.last(star_data)
            if self.px_test(last_data['x'], star['x']) and \
            self.px_test(last_data['y'], star['y']) and \
            self.mag_test(last_data, star):
                return star_id

        # no match, lets try manual
        print '\nimage: {}, time: {}, x: {}, y: {}, mag: {}, merr: {}'.format(image, star['time'], star['x'], star['y'], star['mag'], star['merr'])
        for star_id, star_data in sorted_dump.iteritems():
            last_image, last_data = self.last(star_data)
            if self.skip_px_test(last_data['x'], star['x']) or \
            self.skip_px_test(last_data['y'], star['y']) or \
            self.skip_mag_test(last_data, star):
                continue
            print '\timage: {}, time: {}, x: {}, y: {}, mag: {}, merr: {}'.format(last_image, last_data['time'], last_data['x'], last_data['y'], last_data['mag'], last_data['merr']),
            print "\t\tIs this the above star? ",
            if self.assume == False:
                continue
            if self.assume == True or raw_input().lower() == 'y':
                return star_id
            print '\n'
        return False

    def get_new_star_id(self, sorted_dump):
        if sorted_dump:
            return max(sorted_dump.keys()) + 1
        else:
            return 1

    def last(self, ordered_dict):
        key = next(reversed(ordered_dict))
        return (key, ordered_dict[key])

    def sort_dump(self, dump):
        sorted_dump = defaultdict(OrderedDict)
        for image, data in dump.iteritems():
            for star in data:
                star_id = self.find_similar_star(sorted_dump, star, image)
                if not star_id:
                    star_id = self.get_new_star_id(sorted_dump)
                    print 'New star found - id: {}, image: {}, time: {}, x: {}, y: {}, mag: {}, merr: {}'.format(star_id, image, star['time'], star['x'], star['y'], star['mag'], star['merr'])
                sorted_dump[star_id][image] = star
        return sorted_dump

    def write_csv(self, images, sorted_dump):
        with open('data.csv', 'w') as csvfile:
            csvhandle = csv.writer(csvfile, delimiter=',')
            csvhandle.writerow(['id', 'image', 'time', 'x', 'y', 'mag', 'merr'])
            for star, datapoints in sorted_dump.iteritems():
                for image in sorted(images):
                    if datapoints.has_key(image):
                        datapoint = datapoints[image]
                    else:
                        datapoint = {'time': 0, 'x': 0, 'y': 0, 'mag': 'INDEF', 'merr': 'INDEF'}
                    row = [star, image, datapoint['time'], '{:.3f}'.format(datapoint['x']), '{:.3f}'.format(datapoint['y']), '{}'.format(datapoint['mag']), '{}'.format(datapoint['merr'])]
                    csvhandle.writerow(row)
