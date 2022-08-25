#!/bin/env python3
#coding:utf-8
#
# raw image from fits data

import astropy.io.fits
import numpy as np
import sys

argc = len(sys.argv)
if argc < 5:
    sys.stderr.write('Usage: ' + sys.argv[0] + ' img.fits HDU Nout out.raw\n')
    sys.exit()


fitsfile = sys.argv[1]
HDU = int(sys.argv[2])
Nout = int(sys.argv[3])
OUTFILE = sys.argv[4]

## read image ##
img = astropy.io.fits.open(fitsfile)[HDU].data

## trimming ##
Ny, Nx = img.shape
xi = int( (Nx-1)//2 - (Nout-1)//2 )
yi = int( (Ny-1)//2 - (Nout-1)//2 )
img_trim = img[yi:yi+Nout,xi:xi+Nout]

## output ##
img_trim.astype(np.float32).tofile(OUTFILE)

