#!/usr/bin/env python

# Experiments with a Phase 0 tool
# Bryan Miller

from __future__ import print_function
# from astropy.table import Table
from astropy.io import ascii
import numpy as np

instmodes = ascii.read("phase0_Instrument_Matrix.tsv",format='tab')

# For ToOs instmodes could be limited by what is on the telescope and installed in the instruments

# print(instmodes)
# nrows = instmodes.__sizeof__()
nrows = len(instmodes['mode'])
print(nrows)

#################
# Imaging

mode = 'Imaging'
wlen = 2.15
dwmin = 0.0
dwmax = 'Any'
iqmax = 0.1
fov = 60.
coron = 'No'
# Minimum exposure time (to select a fast read mode?)
minexp = 1.0

if dwmax == 'Any':
    l_dwmax = 10.0
else:
    l_dwmax = dwmax

if dwmin == 'Any':
    l_dwmin = 0.0
else:
    l_dwmin = dwmin

if fov == 'Any':
    l_fov = 1.0
else:
    l_fov = fov

if iqmax == 'Any':
    l_iqmax = 10.0
else:
    l_iqmax = iqmax

# Difference in central wavelength
dcw = np.abs(instmodes['wavelength'] - wlen)
# print(dw)

# Select all optioms that meet the requirements
ii = np.where(np.logical_and(instmodes['mode'] == mode.lower(),
              np.logical_and(instmodes['band_width'] >= l_dwmin,
              np.logical_and(instmodes['band_width'] <= l_dwmax,
              np.logical_and(instmodes['coronagraph'] == coron.lower(),
              np.logical_and(instmodes['minexp'] <= minexp,
              np.logical_and(instmodes['iq_min'] <= l_iqmax, np.logical_and(instmodes['fov'] >= l_fov,
              np.logical_and(instmodes['grcwlen_min'] <= wlen, instmodes['grcwlen_max'] >= wlen)))))))))[0]

nmatch = len(ii)
score = np.zeros(nrows)
irec = []

print()
print(mode)
if len(ii) == 0:
    print('No options found.')
else:
    # print()
    # print('Options:')
    # Try to make a recommendation if multiple options, try to choose the least restrictive
    # Could also use target visibility to decide if options are from different sites

    for i in ii:
        # Configuration match
        score[i] += 1

        # Prefer non-AO if not required for IQ
        if 'no' in instmodes['ao'][i].lower():
            score[i] += 1

        # Central wavelength match
        # if dcw[i] == dcwmin:
        #     score[i] += 1
        dcwscore = wlen/(wlen + dcw[i])
        score[i] += dcwscore

        # print(instmodes['instrument'][i], instmodes['filter'][i], score[i])

    isort = np.argsort(score)
    # Reverse order, high to low
    # ir = np.array([isort[-(j+1)] for j in range(nmatch)])
    ir = np.flip(isort,axis=0)

    print()
    print('Recommendations:')
    for i in range(nmatch):
        irec.append(ir[i])
        print('{:12s} {:10s} {:5.2f}'.format(instmodes['instrument'][ir[i]], instmodes['filter'][ir[i]], score[ir[i]]))

#################
# Spectroscopy

mode = 'Spec'
wlen = 0.8
dwmin = 0.0
dwmax = 'Any'
rmin = 300.
iqmax = 1.0  # proxy for slit width
fov = 'Any'
dims = 1
coron = 'No'
mos = 'No'
# Normal/high quality sky subtraction (to select N&S if target brightness
# not available?
skysub = 'High' # normal, high
minexp = 1.0

if dwmin == 'Any':
    l_dwmin = 0.0
else:
    l_dwmin = dwmin

if dwmax == 'Any':
    l_dwmax = 10.0
else:
    l_dwmax = dwmax

if fov == 'Any':
    l_fov = 1.0
else:
    l_fov = fov

if iqmax == 'Any':
    l_iqmax = 10.0
else:
    l_iqmax = iqmax

# High sky subtraction (N&S) only for the optical, for now
if wlen >= 1.0:
    l_skysub = 'normal'
else:
    l_skysub = skysub.lower()

# Difference in central wavelength
dcw = np.abs(instmodes['wavelength'] - wlen)
# print(dw)

# Difference in slit width
ds = np.abs(instmodes['slit_width'] - l_iqmax)
# print(dw)

# Difference in resolution
dr = np.abs(instmodes['resolution'] - rmin)

# Criteria matching
ii = np.where(np.logical_and(instmodes['mode'] == mode.lower(),
              np.logical_and(instmodes['resolution'] >= rmin,
              np.logical_and(instmodes['spatial_dims'] >= dims,
              np.logical_and(instmodes['band_width'] >= l_dwmin,
              np.logical_and(instmodes['band_width'] <= l_dwmax,
              np.logical_and(instmodes['coronagraph'] == coron.lower(),
              np.logical_and(instmodes['minexp'] <= minexp,
              np.logical_and(instmodes['mos'] == mos.lower(),
              np.logical_and(instmodes['skysub'] == l_skysub.lower(),
              np.logical_and(instmodes['iq_min'] <= l_iqmax, np.logical_and(instmodes['fov'] >= l_fov,
              np.logical_and(instmodes['grcwlen_min'] <= wlen, instmodes['grcwlen_max'] >= wlen)))))))))))))[0]

nmatch = len(ii)
score = np.zeros(nrows)
# print(nrows, nmatch, len(ii))
# print(ii)
irec = []

print()
print(mode)
if nmatch == 0:
    print('No options found.')
else:
    # print()
    # print('Options:')
    for i in ii:
        # Configuration match
        score[i] += 1.

        # Prefer non-AO if not required for IQ
        if 'no' in instmodes['ao'][i].lower():
            score[i] += 1.

        # Wavelength match
        dcwscore = wlen/(wlen + dcw[i])
        score[i] += dcwscore

        # If wavelength > 0.65mu, then prefer settings with a filter (avoid 2nd order)
        if mode.lower() == 'spec' and wlen > 0.65 and instmodes['filter'][i] != 'none':
            score[i] += 0.5

        # Resolution match
        score[i] += rmin/(rmin + dr[i])

        # Slit width match, IFU and MOS always matches
        if int(instmodes['spatial_dims'][i]) == 2 or 'yes' in instmodes['mos'][i].lower():
            score[i] += 1.
        else:
            score[i] += l_iqmax/(l_iqmax + ds[i])

        # Could include a term for sensitivity (S/N) if info available

        # print(instmodes['instrument'][i], instmodes['resolution'][i], instmodes['disperser'][i],
        #       instmodes['filter'][i], instmodes['fpu'][i], dcwscore, score[i])

    isort = np.argsort(score)
    # Reverse, high to low
    # ir = np.array([isort[-(j+1)] for j in range(nmatch)])
    ir = np.flip(isort,axis=0)

    print()
    print('Recommendations:')
    for i in range(nmatch):
        irec.append(ir[i])
        print('{:12s} {:7.0f} {:6s} {:10s} {:18s} {:3s} {:5.2f}'.format(instmodes['instrument'][ir[i]], instmodes['resolution'][ir[i]], instmodes['disperser'][ir[i]],
              instmodes['filter'][ir[i]], instmodes['fpu'][ir[i]], instmodes['ao'][ir[i]], score[ir[i]]))


