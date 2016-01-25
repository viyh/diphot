#!/usr/bin/python

#
# DiPhot: Library - 2016-01-20
# https://github.com/viyh/diphot
#
# Library of common functions
#

import logging, tempfile, csv, sys, math
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
        self.filetypes = ['zero', 'dark', 'flat', 'object']
        self.initialize_dirs()
        self.logger = self.logger_init(self.name)
        self.cleanup_tmp(self.output_dir)
        self.pyraf = Pyraf(self.logger, self.debug)
        self.pyraf.initialize_instrument(self.output_dir)

    def arguments(self):
        raise("Arguments function not implemented!")

    def process(self):
        raise("Process function not implemented!")

    def initialize_dirs(self):
        self.mkdir(self.output_dir, force=True)
        self.mkdir(self.output_dir + '/tmp', force=True)
        self.mkdir(self.output_dir + '/logs', force=True)

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
        self.arguments()
        return self.parser.parse_known_args()[0]

    def logger_init(self, log_name):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console = logging.StreamHandler()
        console.setFormatter(formatter)
        handler = logging.FileHandler(self.output_dir + '/logs/' + log_name + '.log')
        handler.setFormatter(formatter)
        logger.addHandler(console)
        logger.addHandler(handler)
        return logger

    def cleanup_tmp(self, target_dir):
        self.logger.info("Cleaning up tmp directory.".format(target_dir + '/tmp'))
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

    def mkdir(self, dirname, force=False):
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        elif not force:
            answer = raw_input("Directory '{}' already exists. Delete all files? (y/N) ".format(dirname)).lower() == 'y'
            if answer:
                shutil.rmtree(dirname, ignore_errors=True)

class Pyraf():
    def __init__(self, logger, debug=False):
        self.debug = debug
        self.logger = logger

    def initialize_instrument(self, output_dir):
        self.logger.info("Initializing instrument.")
        inst_file_name = output_dir + '/tmp/cp.dat'
        if self.debug:
            iraf.set(debug=1)
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
            instrument = 'cp',
            review = 'no',
            mode = 'h',
            dir = output_dir + '/tmp/',
            site = ''
        )

    def get_files_of_type(self, src_dir, filetype):
        iraf.noao(_doprint=0)
        return iraf.noao.imred.ccdred.ccdlist(
            images = src_dir + '/*.fits',
            ccdtype = filetype,
            names = 'yes',
            Stdout = 1
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
            Stdout = 1
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
            Stderr = 'dev$null'
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
            verify = 'no',
            cache = 0
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
            verify = 'no',
            cache = 0
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
    def __init__(self):
        """
        @param raw_dir: source directory of original FITS cubes
        @type raw_dir: str
        @param output_dir: destination directory where master and calibrated FITS cubes are written
        @type output_dir: str
        """
        DiPhot.__init__(self, 'reduce')
        self.filetypes = ['zero', 'dark', 'flat', 'object']
        self.initialize_type_dirs()

    def arguments(self):
        pass

    def process(self):
        self.logger.info('Creating IRAF FITS images...')
        self.pyraf.run_rfits(self.raw_dir, self.output_dir)
        self.organize_files()
        self.pyraf.create_zero(self.output_dir)
        self.pyraf.apply_zero(self.output_dir)
        self.pyraf.create_dark(self.output_dir)
        self.pyraf.apply_dark(self.output_dir)
        self.pyraf.create_flat(self.output_dir)
        self.pyraf.apply_flat(self.output_dir)

    def initialize_type_dirs(self):
        if not os.path.exists(self.raw_dir):
            self.logger.info("Source directory '{}' does not exist!".format(self.raw_dir))
            sys.exit(1)
        self.mkdir(self.output_dir + '/master')
        for filetype in self.filetypes:
            self.mkdir(self.output_dir + '/' + filetype)

    def organize_files(self):
        for filetype in self.filetypes:
            type_dir = self.output_dir + '/' + filetype + '/'
            self.logger.info("Moving {} files to {}".format(filetype, type_dir))
            files = self.pyraf.get_files_of_type(self.output_dir + '/tmp', filetype)
            self.move_files(files, type_dir)

class CurveOfGrowth(DiPhot):
    def __init__(self):
        DiPhot.__init__(self, 'curveofgrowth')
        self.display = self.args.display
        self.graph_file = self.args.graph_file

    def process(self):
        self.logger.info('Creating curve of growth...')
        self.pyraf.set_datapars()
        self.pyraf.set_centerpars()
        self.pyraf.set_photpars()
        self.pyraf.set_fitskypars()
        self.pyraf.set_findpars()
        self.create_data_files()

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
            self.x = map(lambda x: float(x), v['aperture'])
            self.y1 = map(lambda y1: 1.0/float(y1), v['merr'])
            self.y2 = map(lambda y2: float(y2), v['flux'])
            self.max_snr = 1.0 / float(max(v['merr']))
            self.max_snr_aperture = float(v['merr'].index(min(v['merr'])) + 1)
            self.logger.info('Max SNR [{:.2f}], aperture {}px'.format(self.max_snr, self.max_snr_aperture))
            return True
        self.logger.info('Could not find stars with complete data set!'.format())
        sys.exit()

    def create_plot(self):
        x_np = np.array(self.x)
        y1_np = np.array(self.y1)
        y2_np = np.array(self.y2)
        x_smooth = np.linspace(x_np.min(), x_np.max(), 200)
        y1_smooth = interp1d(x_np, y1_np, kind='cubic')(x_smooth)
        y2_smooth = interp1d(x_np, y2_np, kind='cubic')(x_smooth)

        fig, ax1 = plt.subplots()
        ax1.plot(x_smooth, y1_smooth, 'b-', lw=2)
        ax1.set_xlabel('Aperture (px)')
        ax1.set_ylabel('SNR', color='b')
        ax1.set_ylim([0, max(self.y1) + max(self.y1)/20])
        ax1.axvline(self.max_snr_aperture, color='k', linestyle='--')
        ax1.text(self.max_snr_aperture + 1, 2, 'Aperture with max SNR: ' + str(self.max_snr_aperture) + 'px')
        ax2 = ax1.twinx()
        ax2.plot(x_smooth, y2_smooth, 'r-', lw=2)
        ax2.set_ylabel('Flux', color='r')
        ax2.set_ylim([0, max(self.y2) + max(self.y2)/20])
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
        iraf.noao.digiphot.apphot.centerpars.setParam('cbox', self.fwhm * 2.0)
        iraf.noao.digiphot.apphot.photpars.setParam('apertures', '1:50:1')
        self.pyraf.run_phot(self.output_dir, files[file_num])

        mag_files = self.output_dir + '/tmp/' + os.path.basename(files[file_num]) + '.mag.*'
        mag_dump = self.pyraf.get_txdump(mag_files, 'ID,RAPERT,MAG,MERR,FLUX')
        mags = self.parse_txdump(mag_dump, ['id', 'aperture', 'mag', 'merr', 'flux'])

        self.get_data_points(mags)
        self.create_plot()

    def arguments(self):
        self.parser.add_argument('--display', action='store_true',
            help='display graph')
        self.parser.add_argument('--graph_file', type=str, default=None,
            help='output graph file name')

class Photometry(DiPhot):
    def __init__(self):
        DiPhot.__init__(self, 'photometry')
        self.fwhm = self.args.fwhm
        self.aperture = self.args.aperture

    def arguments(self):
        self.parser.add_argument('--fwhm', type=float, default=8.0,
            help='FWHM')
        self.parser.add_argument('--aperture', type=float, default=15.0,
            help='aperture size')

    def process(self):
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

class TxdumpParse(DiPhot):
    def __init__(self):
        DiPhot.__init__(self, 'curveofgrowth')
        self.px_threshold = 75.0
        self.mag_threshold = 0.8
        self.skip_px_threshold = 90.0
        self.skip_mag_threshold = 1.0
        self.assume = False
        self.missing_tolerance_percent = 17
        self.data = []

    def arguments(self):
        pass

    def process(self):
        dump_file = self.output_dir + '/txdump.txt'
        points = self.read_dump(dump_file)
        dump = self.sort_dump(points)
        images = self.get_full_image_list(dump)
        final_dump = self.clean_dump(images, dump)
        self.organize_star_data(images, final_dump)
        self.write_csv()

    def read_dump(self, dump_filename):
        dump = defaultdict(list)
        with open(dump_filename) as dumpfile:
            for line in dumpfile.readlines():
                image, id, x, y, time, mag, merr = line.split()
                dump[(image, time)].append({'x': float(x), 'y': float(y), 'mag': mag, 'merr': merr})
        return OrderedDict(sorted(dump.items()))

    def sort_dump(self, dump):
        sorted_dump = defaultdict(OrderedDict)
        for (image, time), data in dump.iteritems():
            for star in data:
                star_id = self.find_similar_star(sorted_dump, star, image, time)
                if not star_id:
                    star_id = self.get_new_star_id(sorted_dump)
                    print 'New star found [{}]'.format(star_id)
                    self.show_star(image, time, star)
                sorted_dump[star_id][(image, time)] = star
        return sorted_dump

    def find_similar_star(self, sorted_dump, star, image, time):
        if not sorted_dump:
            return False
        for star_id, star_data in sorted_dump.iteritems():
            last_image, last_data = self.last(star_data)
            if self.match_last(last_data, star):
                return star_id
        return self.manual_match(sorted_dump, star, image, time)

    def get_new_star_id(self, sorted_dump):
        if sorted_dump:
            return max(sorted_dump.keys()) + 1
        else:
            return 1

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

    def match_last(self, last_data, star):
        return self.px_test(last_data['x'], star['x']) and \
            self.px_test(last_data['y'], star['y']) and \
            self.mag_test(last_data, star)

    def match_last_skip(self, last_data, star):
        return self.skip_px_test(last_data['x'], star['x']) or \
            self.skip_px_test(last_data['y'], star['y']) or \
            self.skip_mag_test(last_data, star)

    def show_star(self, image, time, star):
        print('Image:\t{}'.format(image))
        print('Time:\t{}'.format(time))
        print('X:\t{}'.format(star['x']))
        print('Y:\t{}'.format(star['y']))
        print('Mag:\t{}'.format(star['mag']))
        print('MErr:\t{}'.format(star['merr']))

    def manual_match(self, sorted_dump, star, image, time):
        print "\n\n" + "=" * 50
        self.show_star(image, time, star)
        print "=" * 50 + "\n"
        for star_id, star_data in sorted_dump.iteritems():
            last_image, last_data = self.last(star_data)
            if self.match_last_skip(last_data, star): continue
            self.show_star(last_image, time, last_data)
            print "\nIs this the above star (y/N)?",
            if self.assume == False:
                print '\n'
                continue
            if self.assume == True or raw_input().lower() == 'y':
                print '\n'
                return star_id
            print '\n'
        return False

    def last(self, ordered_dict):
        key = next(reversed(ordered_dict))
        return (key, ordered_dict[key])

    def get_row(self, star, point):
        return [
            star.star_id,
            point.image,
            point.time,
            '{:.3f}'.format(point.x),
            '{:.3f}'.format(point.y),
            '{}'.format(point.mag),
            '{}'.format(point.merr)
        ]

    def normalize_star_data(self, images, star, datapoints):
        for image, time in sorted(images):
            if datapoints.has_key((image, time)):
                datapoint = self.point(image=image, time=time, **datapoints[(image, time)])
            else:
                datapoint = self.point(image=image, time=time, x=0, y=0, mag='INDEF', merr='INDEF')
            star.data.append(datapoint)

    def organize_star_data(self, images, sorted_dump):
        for star, datapoints in sorted_dump.iteritems():
            s = self.star(star_id=star, data=[])
            data = self.normalize_star_data(images, s, datapoints)
            self.data.append(s)

    def write_csv(self):
        with open(self.output_dir + '/data.csv', 'w') as csvfile:
            csvhandle = csv.writer(csvfile, delimiter=',')
            csvhandle.writerow(['id', 'image', 'time', 'x', 'y', 'mag', 'merr'])
            for star in self.data:
                for point in star.data:
                    csvhandle.writerow(self.get_row(star, point))

    def get_full_image_list(self, dump):
        images = set()
        for star in dump:
            images |= set(dump[star].keys())
        return images

    def clean_dump(self, images, dump):
        final_dump = defaultdict(OrderedDict)
        for star, values in dump.iteritems():
            star_images = set(values.keys())
            diff = set(images).difference(star_images)
            if diff:
                print "[Star {}] - Incomplete data set: [missing {} datapoints out of {}]".format(star, len(diff), len(images))
            else:
                print "[Star {}] - Complete data set!".format(star)
            if float( len(diff) ) / float ( len(images) ) * 100.0 < self.missing_tolerance_percent:
                print "Saving star [{}]!".format(star)
                final_dump[star] = values
        return final_dump

    class star():
        def __init__(self, star_id=None, data=[]):
            self.star_id = star_id
            self.data = data

        def __str__(self):
            return 'Star [{}]'.format(self.star_id)

    class point():
        def __init__(self, time, image, x , y, mag, merr):
            self.time = time
            self.image = image
            self.x = x
            self.y = y
            self.mag = mag
            self.merr = merr

        def __str__(self):
            return 'Image: {}, Time: {}, X: {:.3f}, Y: {:.3f}, Mag: {}, MErr: {}'.format(
                self.image, self.time, self.x, self.y, self.mag, self.merr)

class LightCurve(DiPhot):
    def __init__(self, target_id=None, data=[]):
        DiPhot.__init__(self, 'lightcurve')
        self.target_id = target_id
        self.data = data
        self.x = []
        self.y1 = []
        self.y2 = []

    def arguments(self):
        self.parser.add_argument('--target_id', type=int, help='target star ID')

    def process(self):
        self.logger.info('Creating light curve...')
        self.separate_stars()
        self.calculate_differential()
        self.create_plot()

    def separate_stars(self):
        stars = defaultdict(list)
        for star in self.data:
            for point in star.data:
                stars[point.time].append((star.star_id, point.mag, point.merr))
        for time in sorted(stars.keys()):
            if 'INDEF' in [p[1] for p in stars[time]] or 'INDEF' in [p[2] for p in stars[time]]:
                stars.pop(time, None)
        self.target_data = self.get_target_data(stars)
        self.comp_data = self.get_comp_data(stars)

    def get_comp_data(self, stars):
        times = sorted(stars.keys())
        comp = OrderedDict()
        for time in times:
            for i, (star_id, mag, merr) in enumerate(stars[time]):
                if star_id == self.target_id:
                    del(stars[time][i])
        for time in times:
            mags = [float(p[1]) for p in stars[time]]
            avg_mag = sum(mags) / len(mags)
            sq_merrs = [math.pow(float(p[2]),2) for p in stars[time]]
            avg_merr = math.sqrt(sum(sq_merrs))
            comp[time] = (avg_mag, avg_merr)
        return comp

    def get_target_data(self, stars):
        times = sorted(stars.keys())
        target = OrderedDict()
        for time in times:
            for i, (star_id, mag, merr) in enumerate(stars[time]):
                if star_id == self.target_id:
                    target[time] = (float(mag), float(merr))
        return target

    def calculate_differential(self):
        for time in self.comp_data.keys():
            from datetime import datetime
            date = datetime.strptime(time, '%H:%M:%S.%f')
            print date, self.comp_data[time][0], self.target_data[time][0], self.comp_data[time][1], self.target_data[time][1]
            self.x.append(date)
            self.y1.append(self.comp_data[time][0] - self.target_data[time][0])
            self.y2.append(self.comp_data[time][1] - self.target_data[time][1])

    def create_plot(self):
        fig, ax1 = plt.subplots()
        ax1.plot(self.x, self.y1, 'b-')
        ax1.set_xlabel('time')
        ax1.set_ylabel('diff mag', color='b')
        ax1.set_ylim([min(self.y1) - max(self.y1)/20, max(self.y1) + max(self.y1)/20])
        plt.show()

    # def create_data_files(self):
    #     """Create coord file (daofind) and mag file (phot) for science images."""
    #     logger.info('Creating light curve')

    #     avg_mags, target_mags = {}, {}
    #     for frame in mags:
    #         avg, n = 0, 0
    #         from datetime import datetime
    #         # date = datetime.strptime(mags[frame][0]['date'] + ' ' + mags[frame][0]['otime'] + '000', '%Y-%m-%d %H:%M:%S.%f')
    #         date = datetime.strptime(mags[frame][0]['otime'] + '000', '%H:%M:%S.%f')
    #         for mag in mags[frame]:
    #             if mag['id'] == str(self.target_id):
    #                 try:
    #                     target_mags[date] = float(mag['mag'])
    #                     continue
    #                 except:
    #                     continue
    #             try:
    #                 avg += float(mag['mag'])
    #                 n += 1
    #             except ValueError:
    #                 pass
    #         try:
    #             avg /= n
    #             avg_mags[date] = avg
    #         except ZeroDivisionError:
    #             pass

    #     x, y_raw, y= [], [], []
    #     for date in sorted(set(avg_mags.keys() + target_mags.keys())):
    #         if not target_mags.has_key(date) or not avg_mags.has_key(date):
    #             continue
    #         x.append(date)
    #         y_raw.append(target_mags[date] - avg_mags[date])

    #     for mag in y_raw:
    #         y.append(mag + min(y_raw))

    #     self.create_plot(x, y)

    class star():
        def __init__(self, star_id=None, data=[]):
            self.star_id = star_id
            self.data = data

        def __str__(self):
            return 'Star [{}]'.format(self.star_id)

    class point():
        def __init__(self, time, image,x , y, mag, merr):
            self.time = time
            self.image = image
            self.x = x
            self.y = y
            self.mag = mag
            self.merr = merr
