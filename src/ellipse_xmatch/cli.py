"""
Command-line interface for ellipse_xmatch.

Usage:
    ellipse-xmatch galaxies.fits points.fits output.fits \\
        --gal-ra gal_ra --gal-dec gal_dec --r1 R1 --r2 R2 --pa PA \\
        --pt-ra ra --pt-dec dec --dlr-factor 1.5 --nbins 20 -v
"""

import argparse
from .core import crossmatch_fits


def main():
    p = argparse.ArgumentParser(
        prog="ellipse-xmatch",
        description="Crossmatch a list of sky coordinates against a galaxy "
                     "catalog with elliptical extents (R1, R2, PA)."
    )
    p.add_argument("gal_fits", help="FITS table of galaxies (ra, dec, R1, R2, PA)")
    p.add_argument("pts_fits", help="FITS table of query points (ra, dec)")
    p.add_argument("output_fits", help="Output FITS table of associations")

    p.add_argument("--gal-ra", default="gal_ra")
    p.add_argument("--gal-dec", default="gal_dec")
    p.add_argument("--r1", default="R1", help="Semi-major axis column (arcsec)")
    p.add_argument("--r2", default="R2", help="Semi-minor axis column (arcsec)")
    p.add_argument("--pa", default="PA", help="Position angle column (deg, East of North)")
    p.add_argument("--pt-ra", default="ra")
    p.add_argument("--pt-dec", default="dec")

    p.add_argument("--dlr-factor", type=float, default=1.0,
                    help="Directional light radius (DLR) scale factor applied to "
                         "both R1 and R2 before matching, growing (>1) or shrinking "
                         "(<1) every galaxy's ellipse (default: 1.0). A factor of "
                         "~1.5 is recommended to avoid missing transients hosted "
                         "with large offsets.")
    p.add_argument("--nbins", type=int, default=20,
                    help="Number of geomspace R1 bins (default: 20)")
    p.add_argument("-v", "--verbose", action="store_true")

    args = p.parse_args()

    crossmatch_fits(
        args.gal_fits, args.pts_fits, args.output_fits,
        gal_ra_col=args.gal_ra, gal_dec_col=args.gal_dec,
        gal_r1_col=args.r1, gal_r2_col=args.r2, gal_pa_col=args.pa,
        pt_ra_col=args.pt_ra, pt_dec_col=args.pt_dec,
        dlr_factor=args.dlr_factor, nbins=args.nbins, verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
