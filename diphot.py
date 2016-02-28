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
import matplotlib.pyplot as plt
from matplotlib import dates
from scipy.interpolate import interp1d
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
        self.args = self.init_parse_args()
        if hasattr(self.config['diphot'], 'debug'):
            self.debug = self.config['diphot']['debug']
        if hasattr(self.args, 'debug'):
            self.debug = self.args.debug
        self.set_attributes()
        self.raw_dir = os.path.expanduser(self.raw_dir.rstrip('/'))
        self.output_dir = os.path.expanduser(self.output_dir.rstrip('/'))
        self.initialize_dirs()
        self.logger = Logger(self).logger
        self.cleanup_tmp(self.output_dir)
        self.pyraf = PyRAF(self.config, self.logger, self.debug)
        self.pyraf.initialize_instrument(self.output_dir)

    def arguments(self):
        raise("Arguments function not implemented!")

    def process(self):
        raise("Process function not implemented!")

    def init_parse_config(self):
        config = yaml.safe_load(open('diphot_defaults.yml'))
        custom_config = {}
        if os.path.exists('diphot.yml'):
            custom_config = yaml.safe_load(open('diphot.yml'))
        for k in config.keys():
            if k not in custom_config: continue
            config[k].update(custom_config[k])
        return config

    def set_attributes(self):
        for k, v in self.config['diphot'].iteritems():
            if k == 'debug': continue
            setattr(self, k, v)
            if hasattr(self.args, k) and getattr(self.args, k):
                setattr(self, k, getattr(self.args, k))
            if self.debug:
                print('Arg - {}: {}'.format(k, getattr(self, k)))

    def initialize_dirs(self):
        self.mkdir(self.output_dir, force=True)
        self.mkdir(self.output_dir + '/tmp', force=True)
        self.mkdir(self.output_dir + '/logs', force=True)

    def init_parse_args(self):
        self.parser = argparse.ArgumentParser(description='DiPhot - Differential Photometry')
        self.parser.add_argument('--raw_dir', type=str,
            help='source directory of raw FITS files')
        self.parser.add_argument('--output_dir', type=str,
            help='destination directory for processed and master FITS files')
        self.parser.add_argument('--target_x', '-x', type=float, dest='target_x',
            help='approximate X pixel of target in first object image')
        self.parser.add_argument('--target_y', '-y', type=float, dest='target_y',
            help='approximate Y pixel of target in first object image')
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

    def initialize_parameters(self, section_name, instance):
        if section_name not in self.config: return
        for k, v in self.config[section_name].iteritems():
            instance.setParam(k, v)

    def set_datapars(self, params=[]):
        """Set datapars parameters for photometry."""
        datapars = iraf.noao.digiphot.apphot.datapars
        self.initialize_parameters('iraf.noao.digiphot.apphot.datapars', datapars)
        for param in params:
            datapars.setParam(param[0], param[1])

    def set_centerpars(self, params=[]):
        """Set centerpars parameters for photometry."""
        centerpars = iraf.noao.digiphot.apphot.centerpars
        self.initialize_parameters('iraf.noao.digiphot.apphot.centerpars', centerpars)
        for param in params:
            centerpars.setParam(param[0], param[1])

    def set_photpars(self, params=[]):
        """Set photpars parameters for photometry."""
        photpars = iraf.noao.digiphot.apphot.photpars
        self.initialize_parameters('iraf.noao.digiphot.apphot.photpars', photpars)
        for param in params:
            photpars.setParam(param[0], param[1])

    def set_fitskypars(self, params=[]):
        """Set fitskypars parameters for photometry."""
        fitskypars = iraf.noao.digiphot.apphot.fitskypars
        self.initialize_parameters('iraf.noao.digiphot.apphot.fitskypars', fitskypars)
        for param in params:
            fitskypars.setParam(param[0], param[1])

    def set_findpars(self, params=[]):
        """Set findpars parameters for photometry."""
        findpars = iraf.noao.digiphot.apphot.findpars
        self.initialize_parameters('iraf.noao.digiphot.apphot.findpars', findpars)
        for param in params:
            findpars.setParam(param[0], param[1])

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
        self.initialize_parameters('iraf.noao.imred.ccdred.zerocombine', zerocombine)
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
        self.initialize_parameters('iraf.noao.imred.ccdred.darkcombine', darkcombine)
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
        self.initialize_parameters('iraf.noao.imred.ccdred.flatcombine', flatcombine)
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
        self.cog_image = self.args.cog_image

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
            self.max_snr = 1.0 / float(min(v['merr']))
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
        cog_image = files[file_num]
        if self.cog_image:
            cog_image = self.cog_image

        self.pyraf.run_daofind(self.output_dir, cog_image)

        coord_files = self.output_dir + '/tmp/' + os.path.basename(cog_image) + '.coo.1'
        coord_dump = self.pyraf.get_txdump(coord_files, 'ID,XCENTER,YCENTER')
        coords = self.parse_txdump(coord_dump, ['id', 'x', 'y'])

        self.fwhm = self.pyraf.run_psfmeasure(cog_image, coords)
        iraf.noao.digiphot.apphot.datapars.setParam('fwhmpsf', self.fwhm)
        iraf.noao.digiphot.apphot.centerpars.setParam('cbox', self.fwhm * 2.0)
        iraf.noao.digiphot.apphot.photpars.setParam('apertures', '1:50:1')
        self.pyraf.run_phot(self.output_dir, cog_image)

        mag_files = self.output_dir + '/tmp/' + os.path.basename(cog_image) + '.mag.*'
        mag_dump = self.pyraf.get_txdump(mag_files, 'ID,RAPERT,MAG,MERR,FLUX')
        mags = self.parse_txdump(mag_dump, ['id', 'aperture', 'mag', 'merr', 'flux'])

        self.get_data_points(mags)
        self.create_plot()

    def arguments(self):
        self.parser.add_argument('--display', action='store_true',
            help='display graph')
        self.parser.add_argument('--graph_file', type=str, default=None,
            help='output graph file name')
        self.parser.add_argument('--cog_image', type=str, default=None,
            help='image to use for curve of growth measurement')

class Photometry(DiPhot):
    def __init__(self):
        DiPhot.__init__(self, 'photometry')

    def arguments(self):
        self.parser.add_argument('--fwhm', type=float,
            help='FWHM')
        self.parser.add_argument('--aperture', type=float,
            help='aperture size')

    def process(self):
        self.logger.info('Creating light curve...')
        self.pyraf.set_datapars(params=[('fwhmpsf', self.fwhm)])
        self.pyraf.set_centerpars(params=[('cbox', self.fwhm * 2.0)])
        self.pyraf.set_photpars(params=[('apertures', self.aperture)])
        self.pyraf.set_fitskypars()
        self.pyraf.set_findpars()
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

    def arguments(self):
        self.parser.add_argument('--px_threshold', '-p', type=float,
            help='Pixel shift in X or Y direction to consider this the same star')
        self.parser.add_argument('--mag_threshold', '-m', type=float,
            help='Magnitude delta to consider this the same star')
        self.parser.add_argument('--skip_px_threshold', '-sp', type=float,
            help='Pixel shift in X or Y direction to not consider this the same star')
        self.parser.add_argument('--skip_mag_threshold', '-sm', type=float,
            help='Magnitude delta to not consider this the same star')
        self.parser.add_argument('--tolerance', dest='missing_tolerance_percent', type=float,
            help='Percentage of missing data points allowed to include in light curve')
        self.parser.add_argument('--assume-true', action='store_true', dest='assume',
            help='Assume stars more than mag/px thesholds but less than skip threshold are the same star')
        self.parser.add_argument('--assume-false', action='store_false', dest='assume',
            help='Assume stars more than mag/px thesholds but less than skip threshold are _not_ the same star')
        pass

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
        import pyfits
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
        self.bin_size = self.args.bin
        self.sigma = self.args.sigma
        self.target_id = self.args.target_id
        if not self.target_id:
            self.target_id = target_id
        self.ignore_ids = self.args.ignore_id
        if not self.ignore_ids:
            self.ignore_ids = ignore_ids

    def arguments(self):
        self.parser.add_argument('-t', '--target_id', type=int, help='target star ID')
        self.parser.add_argument('--bin', type=int, default=1, help='data binning size')
        self.parser.add_argument('--sigma', type=float, default=3.5, help='stddev to remove outlier data points')
        self.parser.add_argument('-i', '--ignore_id', type=int, action='append', help='star to ignore')
        self.parser.add_argument('--comp', action='store_true', help='show individual star graphs')

    def process(self):
        self.logger.info('Creating light curve...')
        self.remove_ignored()
        self.separate_stars()
        # self.get_full_average()
        if self.args.comp:
            self.create_comp_plots()
            sys.exit(0)
        self.calculate_differential()
        self.remove_outliers()
        self.bin_data()
        self.create_diff_plot()

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
            if abs(value['mag'] - median) < self.sigma * np.std(mags):
                filtered[date] = value
        self.points = filtered

    def bin_data(self):
        binned = {}
        for i in xrange(0, len(self.points), self.bin_size):
            chunk = self.points.keys()[i:i + self.bin_size]
            chunk_date = self.avg(chunk)
            mag = self.med([self.points[date]['mag'] for date in chunk])
            merr = self.quad([self.points[date]['merr'] for date in chunk])
            binned[chunk_date] = {'mag': mag, 'merr': merr}
        self.points = binned

    def get_full_average(self):
        avg_comp  = defaultdict(list)
        for time in self.comp_data.keys():
            date = datetime.strptime(time, '%H:%M:%S.%f')
            avg_mag = self.avg([float(d[0]) for d in self.comp_data])
            avg_merr = self.quad([float(d[1]) for d in self.comp_data])
            avg_comp[date] = {'mag': avg_mag, 'merr': avg_merr}
        self.avg_comp = avg_comp

    def create_comp_plots(self):
        self.get_full_average()
        avgs = [self.avg_comp[a]['mag'] for a in self.avg_comp]
        num_stars = len(self.raw_data)
        dim_x = int(math.ceil(math.sqrt(num_stars)))
        dim_y = int(math.ceil(math.sqrt(num_stars)))
        fig, axs = plt.subplots(dim_y, dim_x)
        for i, star in enumerate(self.raw_data):
            ax = axs[int(i/dim_x), i%dim_x]
            x = [datetime.strptime(p.time, '%H:%M:%S.%f') for p in star.data]
            y1 = [float(p.mag) for p in star.data]
            # y1 = [float(p.mag) - float(q) for p, q in zip(star.data, avgs)]
            y2 = [float(p.merr) for p in star.data]
            self.create_plot(ax, x, y1)
            ax.set_title("Star ID " + str(star.star_id))
            # ax.set_ylim([min(y1) - 0.05, max(y1) + 0.05])
        fig.autofmt_xdate()
        for i in range(num_stars, dim_x * dim_y):
            fig.delaxes(axs[int(i/dim_x), i%dim_x])
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

    def create_plot(self, ax, x, y1, y2=None):
        ax.plot_date(x, y1, 'b.')
        if y2:
            ax.errorbar(x, y1, yerr=y2, linestyle='None', color='gray')
        time_fmt = dates.DateFormatter('%H:%M:%S')
        ax.xaxis.set_major_formatter(time_fmt)
        ax.set_xlabel('time')
        ax.set_ylabel('diff mag', color='b')
        ax.set_xlim([min(x), max(x)])
