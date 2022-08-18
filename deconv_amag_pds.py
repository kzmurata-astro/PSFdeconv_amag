#!/bin/env python
#coding:utf-8
#
## PSF Deconvolution with Primal-Dual splitting ##

import sys
import numpy as np 
from scipy import io
from scipy import sparse


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Perform PSF deconvolution via Primal-Dual splitting'
    )
    parser.add_argument("--psfmat", type=str, default='psf.mat', help="A: Observation matrix")
    parser.add_argument("--obsimg", type=str, default='obs.raw', help="Y: Observed image file (float32)")
    parser.add_argument("--varimg", type=str, default='var.raw', help="sig^2: Observed variance file (float32)")
    parser.add_argument("--outimg", type=str, default='out.raw', help="X: Output image file")
    parser.add_argument("--xnimg", type=str, default=None, help="If specified, output Xn (images at each iteration).")
    parser.add_argument("--lam", type=float, default=0.1, help="lambda: Regularization parameter")
    parser.add_argument("--bsig", type=float, default=1.0, help="b in amag to be bsig x sigma; default is b = 1.0 sigma")
    parser.add_argument("--Nite", type=int, default=1000, help="Number of Iterations")
    parser.add_argument("--costfile", type=str, default=None, help="output cost file")
    parser.add_argument("--costplot", default=False, action="store_true", help="plot cost func?")
    args = parser.parse_args()

    
    ## PSF matrix ##
    psfmat = args.psfmat
    A = io.loadmat(psfmat)["A"]
    # Calculating ATA #
    try:
        AA_norm = io.loadmat(psfmat)["AA_norm"][0,0]
    except:
        sys.stderr.write('Calculating ATA for Lipschitz Constant beta, i.e. for data fidelity term.\n')
        AA_norm = maxSingularValue(A.dot(A.T))

    ## observed image and variance ##
    y = np.fromfile(args.obsimg, np.float32)
    yvar = np.fromfile(args.varimg, np.float32)

    ## observed & output image size ##
    ysize, xsize = A.shape # 
    Ny = np.sqrt(ysize).astype(int) # obs image size
    Nx = np.sqrt(xsize).astype(int) # out image size
    
    ## regularization matrix ##
    L = gradmat(Nx)
    LL_norm = np.square(maxSingularValue(L)) #||L||
    
    ## hyper parameters ##
    lam = args.lam
    scale, b, tau, sig = calcHyperParam(AA_norm, LL_norm, yvar, lam)
    b = b * args.bsig
    
    ###########################
    ## Primal-dual splitting ##
    ###########################
    ximg, cost_chi2v, cost_reg = deconvAmagPds(y, yvar, A, L, b, lam, tau, sig, scale=scale, Nite=args.Nite)
    
    # output images #
    ximg[-1].tofile(args.outimg)
    if args.xnimg != None:
        ximg.astype(np.float32).tofile(args.xnimg)
    
    ## Save cost function ##
    if args.costfile !=None: 
        costs = np.stack((cost_chi2v, cost_reg))
        np.savetxt(args.costfile, costs.T, fmt='%e', header='chi2v, reg, lam={:}'.format(lam))

    ## plot ##
    if args.costplot:  
        import matplotlib.pyplot as plt
        plt.plot(0.5*cost_chi2v*y.size, label='0.5 * chi2')
        plt.plot(lam*cost_reg, label='reg term')
        plt.plot(0.5*cost_chi2v*y.size + lam*cost_reg, 'k-',label='total cost')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Iteration')
        plt.ylabel('Cost function')
        plt.legend()
        plt.show()
    

def calcHyperParam(AA_norm, LL_norm, yvar, lam):
    ## Scaling pixel value for fast & stable calculation ##
    # -> sig to be around lam
    beta_tmp = AA_norm / np.min(yvar) # Lipshitz constant of the chi2 term
    scale =  np.sqrt( 0.5 * beta_tmp / LL_norm) / np.sqrt(lam)
    ## Lipshitz constant of the chi2 term ##
    beta = beta_tmp / np.square(scale)
    #beta /= 3.0
    ## Hyper parameters ##
    b = np.median(np.sqrt(yvar)) * scale # for amag calculation
    tau = 1 / beta
    #sig = (1/tau - beta/2) / LL_norm
    sig = 0.5 * beta / LL_norm
    ## output ##
    return scale, b, tau, sig


def deconvAmagPds(y, yvar, A, L, b, lam, tau, sig, scale=1.0, Nite=100):
    a = 1.085736204 #-2.5 / np.log(10)
    magthresh = img2amag(0,b) #non-negativity in amag domain
    cost_chi2v = []
    cost_reg = []
    ximg = []

    ## Pre-process (scaling) ##
    y = y * scale
    yvar = yvar * np.square(scale)
    
    ## initiualize ##
    x0 = np.copy(y) # no sub-sampling is assumed
    xflux_n = x0
    xmag_n = img2amag(xflux_n, b)
    v_n = np.zeros_like(L.dot(xmag_n))
    
    ## loop ##
    for i in range(Nite):
        # gradient #
        chi2v, grad = calc_chi2v(x=xflux_n, y=y, yvar=yvar, A=A)
        cost_chi2v.append(chi2v)
        cost_reg.append( reg_amag(x=xflux_n, b=b, regmat=L) )
        cost = 0.5 * cost_chi2v[i] * y.size + lam * cost_reg[i]
        sys.stderr.write('\r{0:d} {1:.3f} {2:.3f} {3:.3f}'.format(i, cost_chi2v[i], lam*cost_reg[i], cost))

        ## Primal update ##
        xflux_n_tmp = xflux_n - tau * grad 
        xmag_n_tmp = xmag_n - tau * L.T.dot(v_n)

        # non-negativity #
        xflux_n_tmp[xflux_n_tmp < 0] = 0
        xmag_n_tmp[xmag_n_tmp < magthresh] = magthresh

        # fidelity in flux & acmag
        magwei = np.square(a / np.hypot(xflux_n, 2*b))
        xflux_np1 = (xflux_n_tmp + amag2img(xmag_n_tmp, b)*magwei) / (1 + magwei)
        xmag_np1 = img2amag(xflux_np1, b)

        ## Dual update ##
        alpha = 2*lam / (2*lam + sig)
        #alpha = 2 / (2 + sig)
        v_np1 = alpha  * (v_n + sig*L.dot(2*xmag_np1 - xmag_n) )

        ## preparation for next loop ##
        xflux_n = xflux_np1
        xmag_n = xmag_np1
        v_n = v_np1
        ximg.append(xflux_n)

    ## output #
    cost_chi2v = np.array(cost_chi2v)
    cost_reg = np.array(cost_reg)
    ximg = np.array(ximg) / scale # recover the scale
    return ximg, cost_chi2v, cost_reg


def gradmat(Nin):
    Bxy = sparse.lil_matrix((2*Nin*Nin, Nin*Nin),dtype=np.float32)
    #x,y#
    x = np.arange(Nin)
    y = np.arange(Nin)
    x2 = (x.repeat(repeats=Nin).reshape(Nin,Nin).T).flatten()
    y2 = (y.repeat(repeats=Nin).reshape(Nin,Nin)).flatten()
    ##grad x##
    w = np.where((x2>0) & (y2>0))
    index_x0 = (x2+y2*Nin)[w]
    index_0 = index_x0
    index_xm1 = index_0 - 1
    #substitution#
    Bxy[index_x0,index_0] = 1
    Bxy[index_x0,index_xm1] = -1
    #grad y#
    index_y0 = ((x2+y2*Nin) + Nin*Nin)[w]
    index_ym1 = index_0 - Nin
    #substitution#
    Bxy[index_y0,index_0] = 1
    Bxy[index_y0,index_ym1] = -1
    #Convert to csr#
    return Bxy.tocsr()

def maxSingularValue(L):
    from scipy.sparse.linalg import svds
    #maximum value of singular value decompositon #
    u, s, vt = svds(L, k=1)
    return s[0]

# convert flux image to asinh mag image #
def img2amag(x, b):
    # a -> -a
    a = 1.085736204 #-2.5 / np.log(10)
    amag = a * (np.arcsinh(x/(2*b)) + np.log(b))
    return amag

# convert asinh mag image to flux image#
def amag2img(amag, b):
    # a -> -a
    a = 1.085736204 #-2.5 / np.log(10)
    x = np.sinh( amag/a - np.log(b) ) * 2*b
    return x


## data fidelity term ##
def calc_chi2v(x, y, yvar, A):
    dof = y.size - 1
    # residuals #
    dy = A.dot(x) - y
    #chi2#
    chi2 = np.sum( np.square(dy) / yvar)
    chi2_v = chi2 / dof
    #gradient; a factor of 2 disappears because multiplied by 1/2#
    #grad = A.T.dot(dy / yvar) / dof
    grad = A.T.dot(dy / yvar)
    #output#
    return chi2_v, grad

## regularization term ##
def reg_amag(x, b, regmat):
    amag = img2amag(x, b)
    Lx = regmat.dot(amag)
    regvar = np.square(np.linalg.norm(Lx,ord=2)) # ||Lx||^2
    return regvar


if __name__ == '__main__':
    main()

