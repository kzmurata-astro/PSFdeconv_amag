To get data sample used in our paper,

- Go to HSC-SSP data release site.
https://hsc-release.mtk.nao.ac.jp/

- Upload `gal.lst`, at Data Access / Image cutout
https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr2/

- Upload `psf.lst`, at Data Access / PSF Picker
https://hsc-release.mtk.nao.ac.jp/psf/pdr2/


- Obtain raw & PSF images from the fits data with the following command. 

### Spiral-2, DUD
```py
python3 fits2raw.py 2-cutout-HSC-I-9813-pdr2_dud.fits 1 64 spi2_dud_img_64x64.raw
python3 fits2raw.py 2-cutout-HSC-I-9813-pdr2_dud.fits 3 64 spi2_dud_var_64x64.raw
python3 fits2raw.py 2-psf-calexp-pdr2_dud-HSC-I-9813-4,6-150.17282-2.60907.fits 0 40 spi2_dud_psf_40x40.raw
```

### Spiral-2, Wide
```py
python3 fits2raw.py 4-cutout-HSC-I-9813-pdr2_wide.fits 1 64 spi2_wide_img_64x64.raw
python3 fits2raw.py 4-cutout-HSC-I-9813-pdr2_wide.fits 3 64 spi2_wide_var_64x64.raw
python3 fits2raw.py 4-psf-calexp-pdr2_wide-HSC-I-9813-4,6-150.17282-2.60907.fits 0 40 spi2_wide_psf_40x40.raw
```

### Elliptical-2, DUD
```py
python3 fits2raw.py 3-cutout-HSC-I-9813-pdr2_dud.fits 1 64 ell2_dud_img_64x64.raw
python3 fits2raw.py 3-cutout-HSC-I-9813-pdr2_dud.fits 3 64 ell2_dud_var_64x64.raw
python3 fits2raw.py 3-psf-calexp-pdr2_dud-HSC-I-9813-2,6-150.68784-2.61864.fits 0 40 ell2_dud_psf_40x40.raw
```

### Elliptical-2, Wide
```py
python3 fits2raw.py 5-cutout-HSC-I-9813-pdr2_wide.fits 1 64 ell2_wide_img_64x64.raw
python3 fits2raw.py 5-cutout-HSC-I-9813-pdr2_wide.fits 3 64 ell2_wide_var_64x64.raw
python3 fits2raw.py 5-psf-calexp-pdr2_wide-HSC-I-9813-2,6-150.68784-2.61864.fits 0 40 ell2_wide_psf_40x40.raw
```

The sample data can be analyzed as follows (Spiral-2 as an example).

```py
python3 mkpsfmat.py --psfimg spi2_dud_psf_40x40.raw --prefix spi2_dud_psf --img_size 64 --psfimg_size 40 --fl 0.9999
```

```py
python3 deconv_amag_pds.py --obsimg spi2_dud_64x64.raw --varimg spi2_dud_var_64x64.raw --psfmat spi2_dud_psf.mat --outimg spi2_dud_deconv.raw --lam 2.0 --bsig 1.0 --Nite 1000 --xnimg spi2_dud_xn.raw --costfile cost.txt
```


