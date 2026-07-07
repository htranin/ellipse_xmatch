# ellipse-xmatch

Crossmatch a list of sky coordinates against a galaxy catalog with
elliptical extents (R1 = semi-major axis, R2 = semi-minor axis, PA =
position angle), returning every point that falls inside a galaxy's
ellipse.

Efficient on large catalogs: galaxies are binned by R1 (geomspace bins),
and `astropy.coordinates.SkyCoord.search_around_sky` is run per bin using
that bin's max R1 as the search radius. This is guaranteed not to miss
any true match (the ellipse is always inscribed in a circle of radius R1),
and an exact ellipse test then removes circle-but-not-ellipse false
positives.

## Install

Local editable install (recommended while developing):

```bash
cd ellipse_xmatch
pip install -e .
```

Or directly from a folder/tarball without cloning:

```bash
pip install /path/to/ellipse_xmatch
```

Or straight from a git repo, once you push this to GitHub:

```bash
pip install git+https://github.com/yourname/ellipse-xmatch.git
```

## Use as a library

```python
from ellipse_xmatch import crossmatch_fits, crossmatch_ellipses

# quickest path: FITS in, FITS out
result = crossmatch_fits(
    "galaxies.fits", "points.fits", "matches.fits",
    gal_ra_col="gal_ra", gal_dec_col="gal_dec",
    gal_r1_col="R1", gal_r2_col="R2", gal_pa_col="PA",
    pt_ra_col="ra", pt_dec_col="dec",
    nbins=20, verbose=True,
)

# or work directly with SkyCoord arrays you already have in memory
gal_idx, pt_idx = crossmatch_ellipses(gal_coord, pt_coord, R1, R2, PA, nbins=20)
```

## Use from the command line

```bash
ellipse-xmatch galaxies.fits points.fits matches.fits \
    --gal-ra gal_ra --gal-dec gal_dec --r1 R1 --r2 R2 --pa PA \
    --pt-ra ra --pt-dec dec --nbins 20 -v
```

## PA convention

Position angle is assumed to be in degrees, East of North (the standard
astronomical convention). If your PA is defined the opposite way, negate
it before passing it in (or negate it in your input file).
