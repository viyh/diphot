#!/usr/bin/python

#
# DiPhot: Library - 2016-01-20
# https://github.com/viyh/diphot
#
# Library of common functions
#

import logging
import os, shutil, glob
from pyraf import iraf

def logger_init(log_name):
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

def cleanup_tmp(src_dir):
    files = glob.glob(src_dir + '/tmp/*')
    for f in files:
        os.remove(f)

def write_file_from_array(filename, contents):
    with open(filename,'w') as f:
        f.write("\n".join(contents))
        f.write("\n")

def move_files(files, dst_dir):
    for f in files:
        shutil.move(f, dst_dir)

def mkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    else:
        answer = raw_input("Directory '{}' already exists. Delete all files? (y/N) ".format(dirname)).lower() == 'y'
        if answer:
            shutil.rmtree(dirname, ignore_errors=True)

def initialize_instrument():
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

def get_files_of_type(src_dir, filetype):
    iraf.noao(_doprint=0)
    return iraf.noao.imred.ccdred.ccdlist(
        images = src_dir + '/*.fits',
        ccdtype = filetype,
        names = 'yes',
        Stdout=1
    )

def set_datapars(params=[]):
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

def set_centerpars(params=[]):
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

def set_photpars(params=[]):
    """Set photpars parameters for photometry."""
    iraf.noao.digiphot(_doprint=0)
    iraf.noao.digiphot.apphot(_doprint=0)
    iraf.noao.digiphot.apphot.photpars.setParam('weighting', 'constant')
    iraf.noao.digiphot.apphot.photpars.setParam('apertures', '15')
    iraf.noao.digiphot.apphot.photpars.setParam('zmag', '25')
    iraf.noao.digiphot.apphot.photpars.setParam('mkapert', 'yes')
    for param in params:
        iraf.noao.digiphot.apphot.photpars.setParam(param[0], param[1])

def set_fitskypars(params=[]):
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

def set_findpars(params=[]):
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

def get_txdump(filemask, fields):
    iraf.noao.digiphot.ptools(_doprint=0)
    return iraf.noao.digiphot.ptools.txdump(
        textfiles = filemask,
        fields = fields,
        expr = 'yes',
        headers = 'no',
        parameters = 'yes',
        Stdout=1
    )

def run_daofind(src_dir, filename):
    iraf.noao.digiphot.apphot(_doprint=0)
    iraf.noao.digiphot.apphot.daofind(
        image = filename,
        output = src_dir + '/tmp/default',
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

def run_phot(src_dir, filemask):
    iraf.noao.digiphot.apphot.phot(
        image = filemask,
        coords = src_dir + '/tmp/default',
        output = src_dir + '/tmp/default',
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
