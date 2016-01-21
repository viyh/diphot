#!/usr/bin/python

#
# DiPhot: Curve of Growth - 2016-01-20
# https://github.com/viyh/diphot
#
# Create a curve of growth that analyzes reduced FITS files.
#

import os, sys, shutil, argparse, glob
from pyraf import iraf
from collections import defaultdict
import matplotlib.pyplot as plt
import tempfile
import numpy as np
from scipy.interpolate import spline, interp1d

import diphot

logger = diphot.logger_init('curve')

class GrowthCurve:
    def __init__(self, src_dir='data', debug=False, display=False, graph_file=None):
        """
        @param src_dir: source directory of pre-processed FITS cubes
        @type src_dir: str
        """
        self.src_dir = src_dir
        self.debug = debug
        self.display = display
        self.graph_file = graph_file
        diphot.initialize_instrument()

    def create_curve(self):
        logger.info('Creating curve of growth...')
        diphot.cleanup_tmp(self.src_dir)
        diphot.set_datapars()
        diphot.set_centerpars()
        diphot.set_photpars()
        diphot.set_fitskypars()
        diphot.set_findpars()
        self.create_data_files()
        diphot.cleanup_tmp(self.src_dir)

    def parse_txdump(self, dump, fields):
        parsed = defaultdict(list)
        for line in dump:
            arr = line.split()
            if len(arr) != len(fields): continue
            parsed[arr[0]].append(dict(zip(fields, arr)))
        return parsed

    def get_data_points(self, mags):
        i = 1
        while i < len(mags.keys()):
            if any([ m for m in mags[str(i)] if m['merr'] == 'INDEF' or m['flux'] == 'INDEF']):
                i += 1
                continue
            x, y1, y2 = [], [], []
            max_snr_x, max_snr = None, 0
            sorted_mags = sorted(mags[str(i)], key=lambda m: float(m['aperture']))
            for m in sorted_mags:
                x.append(float(m['aperture']))
                y1.append(1.0 / float(m['merr']))
                y2.append(float(m['flux']))
                if 1.0 / float(m['merr']) > max_snr:
                    max_snr = 1.0 / float(m['merr'])
                    max_snr_x = m['aperture']
            logger.info('Max SNR [{:.2f}], aperture {}px'.format(max_snr, max_snr_x))
            return max_snr_x, x, y1, y2
        return False

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
            logger.info('Creating curve of growth graph image: {}'.format(self.graph_file))
            plt.savefig(self.graph_file)

    def run_daofind(self, filename):
        iraf.noao.digiphot(_doprint=0)
        iraf.noao.digiphot.apphot(_doprint=0)
        iraf.noao.digiphot.apphot.daofind(
            image = filename,
            output = self.src_dir + '/tmp/',
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
        iraf.noao.digiphot(_doprint=0)
        iraf.noao.digiphot.ptools(_doprint=0)
        return iraf.noao.digiphot.ptools.txdump(
            textfiles = filemask,
            fields = fields,
            expr = 'yes',
            headers = 'no',
            parameters = 'yes',
            Stdout=1
        )

    def run_phot(self, filemask):
        iraf.noao.digiphot(_doprint=0)
        iraf.noao.digiphot.apphot.phot(
            image = filemask,
            coords = self.src_dir + '/tmp/',
            output = self.src_dir + '/tmp/',
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
            c = coords[coord][0]
            if self.debug:
                logger.debug('\t{}: {}, {}'.format(c['id'], c['x'], c['y']))
            tmp_file.write('{}, {}, {}\n'.format(c['x'], c['y'], c['id']))
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
        logger.info('Found average FWHM of {} stars in image {}: {}'.format(len(coords.keys()), filename, fwhm))
        return fwhm

    def create_data_files(self):
        """Create coord file (daofind) and mag file (phot) for first science image."""
        logger.info('Creating coords, mag for first science image.')

        files = diphot.get_files_of_type(self.src_dir + '/object', 'object')
        file_num = int(len(files) / 2)
        self.run_daofind(files[file_num])

        coord_files = self.src_dir + '/tmp/' + os.path.basename(files[file_num]) + '.coo.1'
        coord_dump = self.get_txdump(coord_files, 'ID,XCENTER,YCENTER')
        coords = self.parse_txdump(coord_dump, ['id', 'x', 'y'])

        fwhm = self.run_psfmeasure(files[file_num], coords)
        iraf.noao.digiphot.apphot.datapars.setParam('fwhmpsf', fwhm)
        iraf.noao.digiphot.apphot.centerpars.setParam('cbox', fwhm*2.0)

        for i in range(1, 50):
            iraf.noao.digiphot.apphot.photpars.setParam('apertures', i)
            self.run_phot(files[file_num])

        mag_files = self.src_dir + '/tmp/' + os.path.basename(files[file_num]) + '.mag.*'
        mag_dump = self.get_txdump(mag_files, 'ID,RAPERT,MAG,MERR,FLUX')
        mags = self.parse_txdump(mag_dump, ['id', 'aperture', 'mag', 'merr', 'flux'])

        max_snr_aperture, x, y1, y2 = self.get_data_points(mags)
        self.create_plot(x, y1, y2, max_snr_aperture)

def parse_args():
    parser = argparse.ArgumentParser(description='Create growth curve')
    parser.add_argument('--src_dir', type=str, default='data',
        help='source directory of FITS files')
    parser.add_argument('--debug', action='store_true',
        help='turn on debug messages')
    parser.add_argument('--display', action='store_true',
        help='display graph')
    parser.add_argument('--graph_file', type=str, default=None,
        help='output graph file name')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    os.environ.get('iraf','/usr/local/iraf')
    os.environ.get('IRAFARCH','linux64')
    p = GrowthCurve(src_dir=args.src_dir, debug=args.debug, display=args.display, graph_file=args.graph_file)
    p.create_curve()
