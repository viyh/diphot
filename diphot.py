#!/usr/bin/python

#
# DiPhot: Library - 2016-01-20
# https://github.com/viyh/diphot
#
# Library of common functions
#

__version__ = '0.2'

import logging, tempfile, csv, sys, math
import os, shutil, glob, argparse
import yaml
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import dates
import pyfits
import astropy.time
from scipy.interpolate import interp1d, griddata
from pyraf import iraf
from collections import defaultdict, OrderedDict
from datetime import datetime

class _Singleton(type):
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(_Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

class Singleton(_Singleton('SingletonMeta', (object,), {})): pass

class Logger(Singleton):
    def __init__(self, diphot):
        self.logger = self.logger_init(diphot.output_dir, 'diphot')

    def logger_init(self, output_dir, log_name):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console = logging.StreamHandler()
        console.setFormatter(formatter)
        handler = logging.FileHandler(output_dir + '/logs/' + log_name + '.log')
        handler.setFormatter(formatter)
        logger.addHandler(console)
        logger.addHandler(handler)
        return logger

class DiPhot():
    def __init__(self, name):
        os.environ.get('iraf','/usr/local/iraf')
        os.environ.get('IRAFARCH','linux64')
        self.name = name
        self.config = self.init_parse_config()
        if hasattr(self.config['diphot'], 'debug'):
            self.debug = self.config['diphot']['debug']
        self.init_parse_args()
        self.set_attributes()
        self.raw_dir = os.path.expanduser(self.raw_dir.rstrip('/'))
        self.output_dir = os.path.expanduser(self.output_dir.rstrip('/'))
        self.initialize_dirs()
        self.logger = Logger(self).logger
        self.cleanup_tmp(self.output_dir)
        self.pyraf = PyRAF(self.config, self.logger, self.debug)
        self.pyraf.initialize_instrument(self.output_dir)
        self.pyraf.initialize_parameters()

    def process(self):
        raise("Process function not implemented!")

    def init_parse_config(self):
        config = yaml.safe_load(open('diphot_defaults.yml'))
        custom_config = {}
        if os.path.exists('diphot.yml'):
            custom_config = yaml.safe_load(open('diphot.yml'))
        config = self.config_merge(custom_config, config)
        return config

    def config_merge(self, custom, default):
        if isinstance(custom, dict) and isinstance(default, dict):
            for k,v in default.iteritems():
                if k not in custom:
                    custom[k] = v
                else:
                    custom[k] = self.config_merge(custom[k], v)
        return custom

    def set_attributes(self):
        for k, v in self.config['diphot'].iteritems():
            if k == 'debug': continue
            if not hasattr(self, k):
                setattr(self, k, v)
            if self.debug:
                print('Arg - {}: {}'.format(k, getattr(self, k)))

    def initialize_dirs(self):
        self.mkdir(self.output_dir, force=True)
        self.mkdir(self.output_dir + '/tmp', force=True)
        self.mkdir(self.output_dir + '/logs', force=True)

    def init_parse_args(self):
        self.parser = argparse.ArgumentParser(description='DiPhot - Differential Photometry')
        self.parser.add_argument('--debug', action='store_true', help='turn on debug messages')
        self.parser.add_argument('--ignore_id', '-i', dest='ignore_ids', type=int, action='append', help='star to ignore')
        self.parser.add_argument('--comp', action='store_true', help='show individual star graphs')
        self.parser.parse_known_args(namespace=self)

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
                os.makedirs(dirname)

class PyRAF(Singleton):
    def __init__(self, config, logger, debug=False):
        self.config = config
        self.logger = logger
        self.debug = debug

    def initialize_instrument(self, output_dir):
        self.logger.info("Initializing instrument.")
        inst_file_name = output_dir + '/tmp/cp.dat'
        if self.debug:
            iraf.set(debug=1)
        iraf.set(use_new_imt='no')
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
        iraf.noao.digiphot(_doprint=0)
        iraf.noao.digiphot.apphot(_doprint=0)
        iraf.dataio(_doprint=0)
        iraf.noao.digiphot.ptools(_doprint=0)

    def get_files_of_type(self, src_dir, filetype):
        return iraf.noao.imred.ccdred.ccdlist(
            images = src_dir + '/*.fits',
            ccdtype = filetype,
            names = 'yes',
            Stdout = 1
        )

    def initialize_parameters(self):
        sections = [
            'iraf.noao.digiphot.apphot.datapars',
            'iraf.noao.digiphot.apphot.centerpars',
            'iraf.noao.digiphot.apphot.photpars',
            'iraf.noao.digiphot.apphot.fitskypars',
            'iraf.noao.digiphot.apphot.findpars'
        ]
        for section_name in sections:
            if section_name not in self.config: continue
            self.set_params(section_name, params=self.config[section_name].items())

    def set_params(self, ns, params=[]):
        """Set PyRAF parameters"""
        instance = eval(ns)
        for param in params:
            instance.setParam(param[0], param[1])

    def get_txdump(self, filemask, fields):
        return iraf.noao.digiphot.ptools.txdump(
            textfiles = filemask,
            fields = fields,
            expr = 'yes',
            headers = 'no',
            parameters = 'yes',
            Stdout = 1
        )

    def run_rfits(self, raw_dir, output_dir):
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

    def run_daofind(self, output_dir, filename, coords='/tmp/default'):
        daofind = iraf.noao.digiphot.apphot.daofind
        daofind.image = filename
        daofind.output = output_dir + coords
        daofind.starmap = ''
        daofind.skymap = ''
        daofind.datapars = ''
        daofind.findpars = ''
        daofind.boundary = 'nearest'
        daofind.constant = 0
        daofind.interactive = 'no'
        daofind.icommands = ''
        daofind.gcommands = ''
        daofind.verify = 'no'
        daofind.cache = 0
        daofind(filename)

    def run_phot(self, output_dir, filemask, coords='/tmp/default', mags='/tmp/default'):
        phot = iraf.noao.digiphot.apphot.phot
        phot.image = filemask
        phot.coords = output_dir + coords
        phot.output = output_dir + mags
        phot.skyfile = ''
        phot.plotfile = ''
        phot.datapars = ''
        phot.centerpars = ''
        phot.fitskypars = ''
        phot.photpars = ''
        phot.interactive = 'no'
        phot.radplots = 'no'
        phot.icommands = ''
        phot.gcommands = ''
        phot.verify = 'no'
        phot.cache = 0
        phot(filemask, Stdout = 1, Stderr = '/dev/null', StdoutG = '/dev/null')

    def run_psfmeasure(self, filename, coords):
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        for coord in coords:
            c = coords[coord]
            if self.debug:
                self.logger.debug('\t{}: {}, {}'.format(coord, c['x'][0], c['y'][0]))
            tmp_file.write('{}, {}, {}\n'.format(c['x'][0], c['y'][0], coord))
        tmp_file.close()
        iraf.noao.obsutil(_doprint=0)
        psfmeasure = iraf.noao.obsutil.psfmeasure
        psfmeasure.images = filename
        psfmeasure.iterations = 1
        psfmeasure.logfile = ''
        psfmeasure.imagecur = tmp_file.name
        psfmeasure.display = 'no'
        psf_out = psfmeasure(filename, Stdout = 1, Stderr = '/dev/null', StdoutG = '/dev/null')
        os.unlink(tmp_file.name)
        fwhm = float(psf_out[-1].split('width at half maximum (FWHM) of ')[-1])
        self.logger.info('Found average FWHM of {} stars in image {}: {}'.format(len(coords.keys()), filename, fwhm))
        return fwhm

    def create_zero(self, output_dir):
        """Creates a master zero image from a set of zero FITS cubes."""
        self.logger.info('Creating master zero...')
        zerocombine = iraf.noao.imred.ccdred.zerocombine
        zerocombine.input = output_dir + '/zero/*.fits'
        zerocombine.output = output_dir + '/master/masterzero'
        zerocombine.combine = 'average'
        zerocombine.reject = 'minmax'
        zerocombine.ccdtype = 'zero'
        zerocombine.process = 'no'
        zerocombine.delete = 'no'
        zerocombine.clobber = 'no'
        zerocombine.scale = 'none'
        zerocombine.rdnoise = iraf.noao.digiphot.apphot.datapars.readnoise
        zerocombine.gain = iraf.noao.digiphot.apphot.datapars.epadu
        zerocombine._runCode()

    def apply_zero(self, output_dir):
        """Apply master zero image to a set of FITS cubes."""
        self.logger.info('Applying master zero to data files.')
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
        darkcombine = iraf.noao.imred.ccdred.darkcombine
        darkcombine.input = output_dir + '/dark/*.fits'
        darkcombine.output = output_dir + '/master/masterdark'
        darkcombine.combine = 'median'
        darkcombine.reject = 'minmax'
        darkcombine.ccdtype = 'dark'
        darkcombine.process = 'yes'
        darkcombine.delete = 'no'
        darkcombine.clobber = 'no'
        darkcombine.scale = 'exposure'
        darkcombine.rdnoise = iraf.noao.digiphot.apphot.datapars.readnoise
        darkcombine.gain = iraf.noao.digiphot.apphot.datapars.epadu
        darkcombine._runCode()

    def apply_dark(self, output_dir):
        """Apply master dark image to a set of FITS cubes."""
        self.logger.info('Applying master dark to data files.')
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
        flatcombine = iraf.noao.imred.ccdred.flatcombine
        flatcombine.input = output_dir + '/flat/*.fits'
        flatcombine.output = output_dir + '/master/masterflat'
        flatcombine.combine = 'median'
        flatcombine.reject = 'avsigclip'
        flatcombine.ccdtype = 'flat'
        flatcombine.subsets = 'yes'
        flatcombine.process = 'yes'
        flatcombine.delete = 'no'
        flatcombine.clobber = 'no'
        flatcombine.statsec = '[700:1400:400:1100]'
        flatcombine.scale = 'median'
        flatcombine.rdnoise = iraf.noao.digiphot.apphot.datapars.readnoise
        flatcombine.gain = iraf.noao.digiphot.apphot.datapars.epadu
        flatcombine._runCode()

    def apply_flat(self, output_dir):
        """Apply master flat image to a set of FITS cubes."""
        self.logger.info('Applying master flat to data files.')
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

    def process(self):
        self.logger.info('Creating curve of growth...')
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
            self.max_snr = 1.0 / float(min(v['merr']))
            self.max_snr_aperture = float(v['merr'].index(min(v['merr'])) + 1)
            self.logger.info('Max SNR [{:.2f}], aperture {}px'.format(self.max_snr, self.max_snr_aperture))
            return True
        self.logger.info('Could not find stars with complete data set!'.format())
        sys.exit()

    def generate_psf_data(self, filename, radius):
        hdulist = pyfits.open(filename)
        scidata = hdulist[0].data
        x_data, y_data, z_data = [], [], []
        for y in range( int(self.target_y) - radius, int(self.target_y) + radius ):
            for x in range( int(self.target_x) - radius, int(self.target_x) + radius ):
                x_data.append(x)
                y_data.append(y)
                z_data.append(scidata[y][x])
        self.logger.debug('Generating grid data...')
        X, Y = np.meshgrid(x_data, y_data)
        Z = griddata((x_data, y_data), z_data, (X, Y), method='linear')
        return X, Y, Z

    def generate_psf_graphs(self, X, Y, Z, radius):
        plt.rcParams.update({'font.size': 10})
        fig, ax = plt.subplots(1, 3, figsize=(9, 3))
        cm = plt.get_cmap("jet")
        stride = radius * 0.25
        ax[0].contourf(Y, Z, X, zdir='x', cstride=stride, rstride=stride, offset=self.target_x - radius, cmap=cm, hold='on')
        ax[0].set_xlabel('y')
        ax[0].set_ylabel('counts')
        ax[1].contourf(X, Z, Y, zdir='y', cstride=stride, rstride=stride, offset=self.target_y + radius, cmap=cm, hold='on')
        ax[1].set_xlabel('x')
        ax[1].set_ylabel('counts')
        ax[2].contourf(X, Y, Z, zdir='z', cstride=stride, rstride=stride, offset=0, cmap=cm, hold='on')
        ax[2].set_xlabel('x')
        ax[2].set_ylabel('y')

    def create_psf_plot(self, filename, radius=50):
        if not self.display_psf or self.psf_graph_file: return
        self.logger.info('Creating contour PSF plots (this can take a minute)...')
        X, Y, Z = self.generate_psf_data(filename, radius)
        self.logger.debug('Generating graphs...')
        self.generate_psf_graphs(X, Y, Z, radius)
        if self.display_psf:
            plt.show()
        if self.psf_graph_file:
            self.logger.info('Creating curve of growth graph image: {}'.format(self.psf_graph_file))
            plt.savefig(self.psf_graph_file)
        plt.show()

    def generate_cog_data(self):
        x_np = np.array(self.x)
        y1_np = np.array(self.y1)
        y2_np = np.array(self.y2)
        x_smooth = np.linspace(x_np.min(), x_np.max(), 200)
        y1_smooth = interp1d(x_np, y1_np, kind='cubic')(x_smooth)
        y2_smooth = interp1d(x_np, y2_np, kind='cubic')(x_smooth)
        return x_smooth, y1_smooth, y2_smooth

    def generate_cog_graph(self, x_smooth, y1_smooth, y2_smooth):
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

    def create_cog_plot(self):
        if not self.display_cog or self.cog_graph_file: return
        self.logger.info('Creating curve of growth graph...')
        x_smooth, y1_smooth, y2_smooth = self.generate_cog_data()
        self.logger.debug('Generating graphs...')
        self.generate_cog_graph(x_smooth, y1_smooth, y2_smooth)
        if self.display_cog:
            plt.show()
        if self.cog_graph_file:
            self.logger.info('Creating curve of growth graph image: {}'.format(self.cog_graph_file))
            plt.savefig(self.cog_graph_file)

    def create_data_files(self):
        """Create coord file (daofind) and mag file (phot) for first science image."""
        self.logger.info('Creating coords, mag for first science image.')

        files = self.pyraf.get_files_of_type(self.output_dir + '/object', 'object')
        file_num = int(len(files) / 2)
        cog_image = files[file_num]
        if self.cog_image:
            cog_image = self.cog_image

        self.pyraf.run_daofind(self.output_dir, cog_image)

        coord_files = self.output_dir + '/tmp/' + os.path.basename(cog_image) + '.coo.1'
        coord_dump = self.pyraf.get_txdump(coord_files, 'ID,XCENTER,YCENTER')
        coords = self.parse_txdump(coord_dump, ['id', 'x', 'y'])

        self.fwhm = self.pyraf.run_psfmeasure(cog_image, coords)
        self.pyraf.set_params(
            ns='iraf.noao.digiphot.apphot.datapars',
            params=[('fwhmpsf', self.fwhm)]
        )
        self.pyraf.set_params(
            ns='iraf.noao.digiphot.apphot.centerpars',
            params=[('cbox', self.fwhm * 2.0)]
        )
        self.pyraf.set_params(
            ns='iraf.noao.digiphot.apphot.photpars',
            params=[('apertures', '1:50:1')]
        )
        self.pyraf.run_phot(self.output_dir, cog_image)

        mag_files = self.output_dir + '/tmp/' + os.path.basename(cog_image) + '.mag.*'
        mag_dump = self.pyraf.get_txdump(mag_files, 'ID,RAPERT,MAG,MERR,FLUX')
        mags = self.parse_txdump(mag_dump, ['id', 'aperture', 'mag', 'merr', 'flux'])

        self.get_data_points(mags)
        self.create_cog_plot()
        self.create_psf_plot(files[0])

class Photometry(DiPhot):
    def __init__(self):
        DiPhot.__init__(self, 'photometry')

    def process(self):
        self.logger.info('Creating light curve...')
        self.pyraf.set_params(
            ns='iraf.noao.digiphot.apphot.datapars',
            params=[('fwhmpsf', self.fwhm)]
        )
        self.pyraf.set_params(
            ns='iraf.noao.digiphot.apphot.centerpars',
            params=[('cbox', self.fwhm * 2.0)]
        )
        self.pyraf.set_params(
            ns='iraf.noao.digiphot.apphot.photpars',
            params=[('apertures', self.aperture)]
        )
        self.initialize_phot_dirs()
        self.create_data_files()

    def initialize_phot_dirs(self):
        self.mkdir(self.output_dir + '/mag')
        self.mkdir(self.output_dir + '/coord')

    def create_filelist(self, files, page):
        page_files = files[page*240:(page+1)*240]
        self.write_file_from_array(self.output_dir + '/images.list', page_files)

    def do_photometry(self, filemask):
        files = sorted(glob.glob(filemask))
        for i in range(0, len(files)/240 + 1):
            self.create_filelist(files, i)
            self.pyraf.run_phot(self.output_dir, '@' + self.output_dir + '/images.list', coords='/coord/default', mags='/mag/default')

    def create_data_files(self):
        """Create coord file (daofind) and mag file (phot) for science images."""
        self.logger.info('Creating coords, mag for science images.')
        self.pyraf.run_daofind(self.output_dir, self.output_dir + '/object/*.fits', coords='/coord/default')
        self.do_photometry(self.output_dir + '/object/*.fits')
        mag_files = self.output_dir + '/mag/*.mag.1'
        mag_dump = self.pyraf.get_txdump(mag_files, 'IMAGE,ID,XCENTER,YCENTER,OTIME,MAG,MERR')
        self.write_file_from_array(self.output_dir + '/txdump.txt', mag_dump)

class TxdumpParse(DiPhot):
    def __init__(self):
        DiPhot.__init__(self, 'txdumpparse')
        self.data = []

    def process(self):
        dump_file = self.output_dir + '/txdump.txt'
        self.create_dump(dump_file)
        points = self.read_dump(dump_file)
        dump = self.sort_dump(points)
        images = self.get_full_image_list(dump)
        final_dump = self.clean_dump(images, dump)
        self.organize_star_data(images, final_dump)
        self.write_csv()
        self.found_target = self.find_target()
        # self.display_image()

    def find_target(self):
        if not self.target_x or not self.target_y: return False
        for star in self.data:
            x = star.data[0].x
            y = star.data[0].y
            if abs(x - self.target_x) < 10 and abs(y - self.target_y) < 10:
                self.logger.info('Found target ID [{}]: ({}, {})'.format(star.star_id, x, y))
                return star.star_id
        return False

    def create_dump(self, dump_filename):
        mag_files = self.output_dir + '/mag/*.mag.1'
        mag_dump = self.pyraf.get_txdump(mag_files, 'IMAGE,ID,XCENTER,YCENTER,OTIME,MAG,MERR')
        self.write_file_from_array(self.output_dir + '/txdump.txt', mag_dump)

    def read_dump(self, dump_filename):
        dump = defaultdict(list)
        with open(dump_filename) as dumpfile:
            for line in dumpfile.readlines():
                print line
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
            if last_image == (image, time): continue
            if self.match_last_skip(last_data, star): continue
            if self.assume == False: continue
            if self.assume == True: return star_id
            self.show_star(last_image[0], last_image[1], last_data)
            print "\nIs this the above star (y/N)?",
            if raw_input().lower() == 'y':
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
                if datapoint.mag == 'INDEF':
                    datapoint.mag = np.nan
                if datapoint.merr == 'INDEF':
                    datapoint.merr = np.nan
            else:
                datapoint = self.point(image=image, time=time, x=0, y=0, mag=np.nan, merr=np.nan)
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

    def display_image(self):
        fitsfile = self.output_dir + '/object/' + self.data[0].data[0].image
        hdulist = pyfits.open(fitsfile)
        tbdata = hdulist[0].data
        plt.imshow(tbdata, cmap='gray')
        plt.colorbar()

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
    def __init__(self, target_id=None, data=[], ignore_ids=[]):
        DiPhot.__init__(self, 'lightcurve')
        self.raw_data = data
        self.points = {}
        self.filtered_points = {}
        self.binned_points = {}
        if not hasattr(self, 'target_id') or not self.target_id:
            self.target_id = target_id
        if not hasattr(self, 'ignore_ids') or not self.ignore_ids:
            self.ignore_ids = ignore_ids

    def process(self):
        self.logger.info('Creating light curve...')
        self.remove_ignored()
        self.separate_stars()
        if (hasattr(self, 'comp') and (self.comp) or not self.target_id):
            self.create_comp_plots()
            sys.exit(0)
        self.calculate_differential()
        self.remove_outliers()
        self.bin_data()
        self.create_diff_plot()
        self.tresca_dump()

    def quad(self, arr):
        sqrs = [math.pow(float(a),2) for a in arr]
        return math.sqrt(sum(sqrs))

    def avg(self, arr):
        return np.average(arr)

    def med(self, arr):
        return np.median(arr)

    def numpy_zip(self, x, y):
        return np.array(zip(x, y), dtype=[('x',float),('y',float)])

    def remove_ignored(self):
        self.raw_data = [s for s in self.raw_data if s.star_id not in self.ignore_ids]

    def separate_stars(self):
        stars = defaultdict(list)
        for star in self.raw_data:
            for point in star.data:
                stars[point.time].append((star.star_id, point.mag, point.merr))
        for time in sorted(stars.keys()):
            if np.nan in [p[1] for p in stars[time]] or np.nan in [p[2] for p in stars[time]]:
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
            avg_mag = self.avg(mags)
            avg_merr = self.quad([p[2] for p in stars[time]])
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
            date = datetime.strptime(time, '%H:%M:%S.%f')
            diff_mag = abs(self.target_data[time][0] - self.comp_data[time][0])
            diff_merr = self.quad([self.comp_data[time][1], self.target_data[time][1]])
            self.points[dates.date2num(date)] = {'mag': diff_mag, 'merr': diff_merr}

    def remove_outliers(self):
        filtered = {}
        mags = [p['mag'] for p in self.points.values()]
        median = np.median(mags)
        for date, value in self.points.iteritems():
            if abs(value['mag'] - median) < self.lightcurve_sigma * np.std(mags):
                filtered[date] = value
        self.points = filtered

    def bin_data(self):
        binned = {}
        for i in xrange(0, len(self.points), self.lightcurve_bin):
            chunk = self.points.keys()[i:i + self.lightcurve_bin]
            chunk_date = self.avg(chunk)
            mag = self.med([self.points[date]['mag'] for date in chunk])
            merr = self.quad([self.points[date]['merr'] for date in chunk])
            binned[chunk_date] = {'mag': mag, 'merr': merr}
        self.points = binned

    def get_comp_average(self):
        avg_comp = []
        for time in self.comp_data.keys():
            # avg_mag = self.avg([float(d[time][0]) for d in self.comp_data])
            # avg_merr = self.quad([float(d[time][1]) for d in self.comp_data])
            avg_comp.append(self.datapoint(0, time, self.comp_data[time][0], self.comp_data[time][1]))
        self.avg_comp = avg_comp

    def create_comp_ax(self, axs, i, dim_x, dim_y, data, desc):
        if dim_x == 1:
            ax = axs
        elif dim_y == 1:
            ax = axs[int(i%dim_x)]
        else:
            ax = axs[int(i/dim_x), int(i%dim_x)]
        x = [datetime.strptime(p.time, '%H:%M:%S.%f') for p in data]
        y1 = [float(p.mag) for p in data]
        y2 = [float(p.merr) for p in data]
        self.create_plot(ax, x, y1)
        ax.set_title(desc)

    def create_comp_plots(self):
        num_stars = len(self.raw_data)
        if num_stars == 0:
            self.logger.info('No stars found! Try adjusting the tolerence percentage and max/px thresholds.')
            sys.exit(0)
        num_stars += 1
        dim_x = int(math.ceil(math.sqrt(num_stars)))
        dim_y = int(math.ceil(math.sqrt(num_stars)))
        fig, axs = plt.subplots(dim_y, dim_x)
        for i, star in enumerate(self.raw_data):
            self.create_comp_ax(axs, i, dim_x, dim_y, star.data, "Star ID " + str(star.star_id))

        self.get_comp_average()
        self.create_comp_ax(axs, i+1, dim_x, dim_y, self.avg_comp, "Average")
        fig.autofmt_xdate()

        for i in range(num_stars, dim_x * dim_y):
            fig.delaxes(axs[int(i/dim_x), int(i%dim_x)])
        plt.show()

    def create_diff_plot(self):
        fig, (ax1, ax2) = plt.subplots(1, 2)
        x = self.points.keys()
        y1 = [v['mag'] for v in self.points.values()]
        y2 = [v['merr'] for v in self.points.values()]
        self.create_plot(ax1, x, y1, y2)
        fig.autofmt_xdate()
        ax2.hexbin(x, y1, gridsize=100, bins=20)
        time_fmt = dates.DateFormatter('%H:%M:%S')
        ax2.xaxis.set_major_formatter(time_fmt)
        ax2.set_xlim([min(x), max(x)])
        ax2.set_ylim([min(y1) - 0.05, max(y1) + 0.05])
        plt.show()

    def get_type_filelist(self, filetype):
        filemask = self.output_dir + '/' + filetype + '/*.fits'
        return sorted(glob.glob(filemask))

    def get_julian_dates(self):
        files = self.get_type_filelist('object')
        dts = {}
        for f in files:
            hdulist = pyfits.open(f)
            header = hdulist[0].header
            dt = datetime.strptime(header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
            jd = astropy.time.Time(dt).jd
            dts[dt.strftime("%H:%M:%S")] = jd
        return dts

    def tresca_dump(self):
        dts = self.get_julian_dates()
        x = [dts[dates.num2date(k).strftime("%H:%M:%S")] for k in self.points.keys()]
        y1 = [v['mag'] for v in self.points.values()]
        y2 = [v['merr'] for v in self.points.values()]
        rows = zip(x, y1, y2)
        f = open(self.output_dir + '/tresca.csv', 'wb')
        writer = csv.writer(f, delimiter = '\t')
        writer.writerow(['JD', 'MAG', 'ERROR'])
        for row in rows:
            writer.writerow(row)

    def create_plot(self, ax, x, y1, y2=None):
        ax.plot_date(x, y1, 'b.')
        if y2:
            ax.errorbar(x, y1, yerr=y2, linestyle='None', color='gray')
        time_fmt = dates.DateFormatter('%H:%M:%S')
        ax.xaxis.set_major_formatter(time_fmt)
        ax.set_xlabel('time')
        ax.set_ylabel('mag', color='b')
        ax.set_xlim([min(x), max(x)])

    class datapoint():
        def __init__(self, star_id, time, mag, merr):
            self.star_id = star_id
            self.time = time
            self.mag = mag
            self.merr = merr

        def __str__(self):
            return 'Time: {}, Mag: {}, MErr: {}'.format(
                self.time, self.mag, self.merr)
