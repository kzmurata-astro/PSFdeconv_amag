# Deblurring galaxy images with Tikhonov Regularization on magnitude domain

<p align="center">
<img src="https://github.com/kzmurata-astro/git-tutorial/blob/main/img/spi2_before_after.png" width=40% height=40%>
<br> (c) Murata & Takeuchi in prep. </br>
</p>

There are two steps required to run the code:

## 1. producing PSF matrix
```py
python3 mkpsfmat.py --psfimg psf.raw --prefix psf --img_size 64 --psfimg_size 40 --fl 0.9999
```
- psgimg: PSF image (raw image with flaot32) 
- prefix: Output matrix file is [prefix].mat
- img_size: Obs / Deconv image size
- psfimg_size: PSF image size
- fl: Fraction of light included in the output PSF


## 2. Perform PSF deconvlution

```py
python3 deconv_amag_pds.py --obsimg observed.raw --varimg variance.raw --psfmat psf.mat --outimg out_deconv.raw --lam 1.0 --bsig 1.0 --Nite 1000 [--xnimg xn.raw --costfile cost.txt --costplot]
```
Options

- obsimg: Observed image (raw image with flaot32)
- varimg: Variance image (raw image with flaot32)
- psfmat: PSF matrix
- outimg: Output PSF-deconvolved image (raw image with flaot32)
- lam: Regulatization parameter lambda
- bsig: scale parameter, where b in asinh magnitude is bsig x sqrt(<variance>)
- Nite: Number of iterations
- xnimg: If specified, images at each iteration are output.
- costfile: if specified, the data-fidelity and regulatization terms normalized by d.o.f. are output.
- costplot: if set, plot the cost sequence. 


## Required Environment
- `numpy`
- `scipy`
- `matplotlib` if you want to show the sequence of cost function



