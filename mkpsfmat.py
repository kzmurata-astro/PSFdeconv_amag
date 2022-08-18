#!/bin/env python
#coding:utf-8
#

import numpy as np
from scipy import sparse
from scipy.interpolate import interp2d
import sys
from scipy import io

def main():
    import argparse
    parser = argparse.ArgumentParser(
    description='Produce PSF observation matrix.'
    )
    parser.add_argument("--psfimg", help="PSF image in raw-file in float32", type=str)
    parser.add_argument("--img_size", help="image size", type=int, default=64)
    parser.add_argument("--psfimg_size", help="PSF image size", type=int, default=40)
    parser.add_argument("--fl", help="fraction of light included in the output PSF", default=0.9999, type=float)
    parser.add_argument("--prefix", help="prefix of output matrix file", default='psf', type=str)
    args = parser.parse_args()

    ## parameters ##
    # img size, assuming the same x & y size #
    img_size = args.img_size
    origimg_size = img_size # assuming original size is the same as the observed one.
    psfimg_size = args.psfimg_size

    # fraction light in output PSF #
    f_light = args.fl

    #output#
    prefix = args.prefix
    outfile = prefix + '.mat'

    ## Reading PSF image in raw ##
    if args.psfimg != None:
        psfimg = np.fromfile(args.psfimg, np.float32).reshape(psfimg_size, psfimg_size)
    else:
        sys.stderr.write('Specify PSF file. Exit.\n')
        sys.exit()
    
    ## coordinates of orig, obs, psf matrices ##
    xi, yi = xymat(origimg_size)
    xo, yo = xymat(img_size) #observed image
    xp, yp = xymat(psfimg_size)

    # function for psf values #
    psffunc = interp2d(xp[0], yp[:,0], psfimg, kind='linear', bounds_error=False, fill_value=0)

    # max_radius up to which PSF matrix is calculated #
    dist = flight_radius(psfimg, f_light)
    
    ## prepare for observed sparse matrix ##
    A = sparse.lil_matrix((img_size*img_size, origimg_size*origimg_size),dtype=np.float32)
    
    # loop for original index; orig light is divided into obs img #
    sys.stderr.write('Converting {:} -> {:}\n'.format(args.psfimg, outfile))
    for i in range(origimg_size*origimg_size):
        dx = xo[0] - xi.flatten()[i]
        dy = yo[:,0] - yi.flatten()[i]
        # ignore too far pixels #
        wx = np.where( abs(dx) < dist)[0]
        wy = np.where( abs(dy) < dist)[0]
        index_out = (wy.reshape(wy.size,1) * img_size + wx.reshape(1,wx.size)).flatten()
        if (wx.size > 0) & (wy.size > 0):
            A[index_out, i] = psffunc(dx[wx], dy[wy]).flatten()


    # max singular value of A #
    sys.stderr.write('Calculating ATA for Lipschitz Constant beta\n')
    AA_norm = maxSingularValue(A.dot(A.T))
            
    ## output ##
    io.savemat(outfile,{"A":A, "AA_norm":AA_norm})
    
    # for confirmation output raw images #
    #A.toarray().astype(np.float32).tofile(prefix + '_check_{0}x{1}x{2}.raw'.format(origimg_size, origimg_size, img_size*img_size))


def xymat(Nsize, scale=1.0):
    x = (np.arange(Nsize) - (Nsize-1)//2) * scale
    y = (np.arange(Nsize) - (Nsize-1)//2) * scale
    x, y = np.meshgrid(x,y)
    return x, y


# calculate radius where fraction of light > f_light 
def flight_radius(psfimg, f_light, pixsize=1.0):
    flux_tot = np.sum(psfimg)
    Npsf = psfimg.shape[0]
    xp, yp = xymat(Npsf, scale=1.0) # calculate in pix
    rp = np.hypot(xp, yp)
    r_ref = np.arange(0, 0.5*Npsf, 0.5)
    flux = np.zeros_like(r_ref)
    for i in range(r_ref.size):
        w = np.where(rp < r_ref[i])
        flux[i] = np.sum(psfimg[w])
    dist = np.interp(f_light, flux/flux_tot, r_ref, left=np.min(r_ref), right=np.max(r_ref))
    dist *= pixsize # pixsize to in_pix
    return dist


# maximum value of singular value decompositon #
def maxSingularValue(L):
    from scipy.sparse.linalg import svds
    u, s, vt = svds(L, k=1)
    return s[0]


if __name__ == '__main__':
    main()


