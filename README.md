# ellipse-xmatch

Crossmatch a list of sky coordinates against a galaxy catalog with
elliptical extents (R1 = semi-major axis, R2 = semi-minor axis, PA =
position angle, degrees East of North), returning every point that falls
inside a galaxy's ellipse.

Efficient on large catalogs: galaxies are binned by R1 (geomspace bins),
and `astropy.coordinates.SkyCoord.search_around_sky` is run per bin using
that bin's max R1 as the search radius. This can't miss a true match (the
ellipse is always inscribed in a circle of radius R1), and an exact
ellipse test then removes circle-but-not-ellipse false positives.

## Install

```bash
cd ellipse_xmatch
pip install -e .
```

## Galaxy catalog

Need a `gal_fits` input? Download the REGALADE catalog (v2, July 6th 2026,
reduced number of columns, 4GB):
[Google Drive link](https://drive.google.com/file/d/1dvnjt-0Y4zN0hcy67QIlwGXMx-tHNgoP/view?usp=sharing)

(Full documentation available in [this repo](https://github.com/htranin/regalade))

## Use as a library

```python
from ellipse_xmatch import crossmatch_fits, crossmatch_ellipses

# quickest path: FITS in, FITS out
result = crossmatch_fits(
    "galaxies.fits", "points.fits", "matches.fits",
    gal_ra_col="gal_ra", gal_dec_col="gal_dec",
    gal_r1_col="R1", gal_r2_col="R2", gal_pa_col="PA",
    pt_ra_col="ra", pt_dec_col="dec",
    dlr_factor=1.5, nbins=20, verbose=True,
)

# or work directly with SkyCoord arrays you already have in memory
gal_idx, pt_idx, separation = crossmatch_ellipses(gal_coord, pt_coord, R1, R2, PA,
                                                   dlr_factor=1.5, nbins=20)
```

## Use from the command line

```bash
ellipse-xmatch galaxies.fits points.fits matches.fits \
    --gal-ra gal_ra --gal-dec gal_dec --r1 R1 --r2 R2 --pa PA \
    --pt-ra ra --pt-dec dec --dlr-factor 1.5 --nbins 20 -v
```

## Notes

- **DLR factor**: `dlr_factor` (`--dlr-factor` on the CLI) scales both R1
  and R2 before matching, growing (>1) or shrinking (<1) every galaxy's
  ellipse uniformly. Default is 1.0; ~1.5 is recommended to avoid missing
  transients hosted with large offsets.
- **`separation` column**: for each match, the point's angular distance
  from the galaxy center divided by the ellipse's radius in that
  direction — 0 at the center, 1 on the boundary. 
- **PA convention**: degrees East of North. If your PA is defined the
  opposite way, negate it before passing it in.
- **File formats**: `gal_fits`, `pts_fits`, and `output_fits` accept any
  format `astropy.table.Table.read`/`.write` supports (FITS, CSV, ECSV,
  HDF5, ...), inferred from the file extension — not just FITS.
