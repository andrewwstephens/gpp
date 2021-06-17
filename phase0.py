#!/usr/bin/env python3
import argparse
from   astropy.io import ascii
import logging
import numpy
import sys

# Data at https://docs.google.com/spreadsheets/d/1-Rz-uqEU3AM3LiDMLYS2rbh5PrSDp3cdyjF7flV_Ztg/

__version__ = '2019-Feb-12'  # Bryan Miller, original version
__version__ = '2021-Jun-16'  # astephens, separate imaging and spectroscopy spreadsheets

def main(args):
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S', level=getattr(logging, args.loglevel.upper()))
    logger = logging.getLogger()

    if args.mode == 'imaging':
        configfile = 'phase0.imaging.csv'
    else:
        configfile = 'phase0.spectroscopy.csv'
    instmodes = ascii.read(configfile, format='csv')
    nrows = len(instmodes)
    logger.debug('Read %d instrument modes from %s', nrows, configfile)
    logger.debug('instmodes:\n%s', instmodes)

    if args.mode == 'imaging':
        instmodes['filters'] = [s.split(',') for s in instmodes['filters']]  # convert string to list
        filtermatch = numpy.array([len(set(args.filter) & set(f)) > 0 for f in instmodes['filters']])
        logger.debug('filtermatch: %s', filtermatch)
        aomatch = (instmodes['AO'] == 'yes') if (args.iq < 0.2) else numpy.full(nrows, True)  # Allow AO even if not needed
        logger.debug('aomatch: %s', aomatch)
        specklematch = (instmodes['capabilities'] == 'speckle').filled(False) if ('speckle' in args.capabilities) else (instmodes['capabilities'] != 'speckle').filled(True)
        logger.debug('specklematch: %s', specklematch)
        idx = numpy.where((instmodes['FoV'] >= args.fov) & filtermatch & aomatch & specklematch)[0]

    elif args.mode == 'spectroscopy':
        instmodes['Focal Plane'] = [s.split(',') for s in instmodes['Focal Plane']]  # convert string to list
        fpmatch = [args.fpu in f for f in instmodes['Focal Plane']]
        logger.debug('fpmatch: %s', fpmatch)
        nsmatch = (instmodes['capabilities'] == 'Nod&Shuffle').filled(False) if ('nodshuffle' in args.capabilities) else (instmodes['capabilities'] != 'Nod&Shuffle').filled(True)
        logger.debug('nsmatch: %s', nsmatch)
        aomatch = (instmodes['AO'] == 'yes') if (args.iq < 0.2) else numpy.full(nrows, True)
        logger.debug('aomatch: %s', aomatch)
        coronamatch = (instmodes['capabilities'] == 'coronagraph').filled(False) if ('coronagraph' in args.capabilities) else (instmodes['capabilities'] != 'coronagraph').filled(True)
        logger.debug('coronamatch: %s', coronamatch)
        idx = numpy.where(
            (instmodes['wave min'] <= args.wave) &
            (instmodes['wave max'] >= args.wave) &
            (instmodes['resolution'] >= args.res) &
            (instmodes['wave range'] >= args.range) &
            (instmodes['slit length'] >= args.fov) &
            fpmatch & nsmatch & aomatch & coronamatch)[0]

    logger.debug('idx: %s', idx)
    nmatch = len(idx)
    logger.debug('Found %d matches!\n' % nmatch)

    if args.mode == 'imaging':
        print('Field of View:', args.fov)
        print('Filters:', args.filter)
    elif args.mode == 'spectroscopy':
        print('Focal Plane: %s' % args.fpu)
        print('Central Wavelength: %0.3f um' % args.wave)
        print('Spectral Resolution: %d' % args.res)
        print('Slit Length: %0.1f arcsec' % args.fov)
    print('Special Capabilities: %s' % args.capabilities)

    print('\nMatching configurations:')
    if args.mode == 'imaging':  # Just print the matching configurations
        for i in idx:
            print(instmodes['instrument'][i])
    elif args.mode == 'spectroscopy':  # prioritize the configurations since there may be many
        delta_wav = numpy.abs(instmodes['wave optimal'] - args.wave)     # Difference in wavelength
        delta_slit_width = numpy.abs(instmodes['slit width'] - args.iq)  # Difference in slit width
        delta_res = numpy.abs(instmodes['resolution'] - args.res)        # Difference in resolution
        score = numpy.zeros(nrows)
        for i in idx:

            if args.iq > 0.2 and instmodes['AO'][i] == 'no':  # give a bump to non-AO modes (but don't discount them)
                score[i] += 0.5

            score[i] += args.wave/(args.wave + delta_wav[i])  # Wavelength match

            # If wavelength > 0.65mu, then prefer settings with a filter to avoid 2nd order contamination
            if args.wave > 0.65 and instmodes['filter'][i] != 'none':
                score[i] += 0.5

            score[i] += args.res/(args.res + delta_res[i])  # Resolution match

            if 'ifu' in instmodes['Focal Plane'][i]:  # Slit width match to the seeing (the IFU always matches)
                score[i] += 1.0
            else:
                score[i] += args.iq / (args.iq + delta_slit_width[i])

        print('Instrument  FPU                 Disperser  Filter Resolution  Score')
        for i in numpy.flip(numpy.argsort(score)):
            if score[i] > 0:
                print('%-11s %-19s %-10s %-8s %8d %6.3f' %
                      (instmodes['instrument'][i], instmodes['fpu'][i], instmodes['disperser'][i],
                       instmodes['filter'][i], instmodes['resolution'][i], score[i]))

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description='Phase 0 experimentation.\n\n' +
        'examples:\n' +
        '  phase0.py  imaging  --filter r \n' +
        '  phase0.py  imaging  --filter K  --iq 0.1\n' +
        '  phase0.py  imaging  --filter g  --capabilities speckle\n' +
        '  phase0.py  spectroscopy  --wave 0.61  --res 1000  --iq 0.8\n' +
        '  phase0.py  spectroscopy  --wave 0.55  --range 0.3 --capabilities nodshuffle\n' +
        '  phase0.py  spectroscopy  --fpu multislit --wave 0.48  --res 2000\n' +
        '  phase0.py  spectroscopy  --wave 1.65  --range 2.0\n' +
        '  phase0.py  spectroscopy  --fpu ifu --wave 2.2\n'
        '  phase0.py  spectroscopy  --fpu ifu --wave 2.2  --capabilities coronagraph',
        epilog='Version: ' + __version__)

    parser.add_argument('mode', default=None, type=str,
                        choices=['imaging', 'spectroscopy'])

    parser.add_argument('--filter', default=None, type=str, action='append',
                        help='Imaging filter')

    parser.add_argument('--fpu', default='singleslit', type=str, action='store',
                        choices=['singleslit', 'multislit', 'ifu'],
                        help='Spectroscopy focal plane unit [singleslit]')

    parser.add_argument('--fov', default=1.0, type=float, action='store',
                        help='Field of view or slit length (arcsec)')

    parser.add_argument('--res', action='store', type=int, default=1,
                        help='Spectral resolution')

    parser.add_argument('--wave', action='store', type=float, default=0.5,
                        help='Spectroscopy central wavelength (microns)')

    parser.add_argument('--range', action='store', type=float, default=0.0,
                        help='Spectroscopy wavelength range (microns)')

    parser.add_argument('--iq', action='store', type=float, default=1.0,
                        help='Image quality (arcseconds)')

    parser.add_argument('--capabilities', default=[None], type=str, action='append',
                        choices=['nodshuffle', 'speckle', 'coronagraph'])

    parser.add_argument('--loglevel', type=str, default='info',
                        choices=['debug', 'info', 'warning', 'error'])

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    if args.mode == 'imaging' and args.filter is None:
        args.filter = ['r']  # parser.error('Please specify an imaging filter')

    main(args)

# --------------------------------------------------------------------------------------------------
