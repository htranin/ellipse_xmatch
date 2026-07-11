"""
Core crossmatch logic: points vs galaxy ellipses, binned by R1 with
search_around_sky for speed + completeness on large catalogs.
"""

import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
import time

__all__ = ["in_ellipse", "ellipse_separation", "crossmatch_ellipses", "crossmatch_fits"]


def ellipse_separation(gal_c, pt_c, r1, r2, pa_deg):
    """
    Normalized elliptical separation (vectorized).

    Computes the angular distance from gal_c to pt_c divided by the
    galaxy ellipse's radius in that same direction, i.e. how far along
    the ray from the galaxy center to the point the ellipse boundary
    falls. Values <= 1 are inside the ellipse, > 1 are outside.

    Parameters
    ----------
    gal_c, pt_c : SkyCoord
        Equal-length arrays of paired candidate (galaxy, point) coordinates.
    r1, r2 : array-like
        Semi-major / semi-minor axes, in arcsec.
    pa_deg : array-like
        Position angle in degrees, East of North convention.
        If your PA is defined the other way round, negate it before calling.

    Returns
    -------
    numpy.ndarray of float
        Normalized separation: 0 at the galaxy center, 1 on the ellipse
        boundary, >1 outside.
    """
    dlon, dlat = gal_c.spherical_offsets_to(pt_c)  # dlon ~ East, dlat ~ North
    dE = dlon.to(u.arcsec).value
    dN = dlat.to(u.arcsec).value

    pa = np.deg2rad(np.asarray(pa_deg, dtype=float))
    x_major = dE * np.sin(pa) + dN * np.cos(pa)
    y_minor = dE * np.cos(pa) - dN * np.sin(pa)

    return np.sqrt((x_major / r1) ** 2 + (y_minor / r2) ** 2)


def in_ellipse(gal_c, pt_c, r1, r2, pa_deg):
    """
    Exact ellipse membership test (vectorized).

    Parameters
    ----------
    gal_c, pt_c : SkyCoord
        Equal-length arrays of paired candidate (galaxy, point) coordinates.
    r1, r2 : array-like
        Semi-major / semi-minor axes, in arcsec.
    pa_deg : array-like
        Position angle in degrees, East of North convention.
        If your PA is defined the other way round, negate it before calling.

    Returns
    -------
    numpy.ndarray of bool
        True where pt_c falls inside the ellipse centered on gal_c.
    """
    return ellipse_separation(gal_c, pt_c, r1, r2, pa_deg) <= 1.0


def crossmatch_ellipses(gal_coord, pt_coord, R1, R2, PA, dlr_factor=1.0,
                         nbins=20, verbose=False, t0=time.time()):
    """
    Find every (point, galaxy) pair where the point falls inside the
    galaxy's ellipse, using geomspace R1 bins + search_around_sky.

    Parameters
    ----------
    gal_coord : SkyCoord
        Galaxy centers.
    pt_coord : SkyCoord
        Query points.
    R1, R2 : array-like
        Semi-major / semi-minor axes, in arcsec. Same length as gal_coord.
    PA : array-like
        Position angle, degrees, East of North. Same length as gal_coord.
    dlr_factor : float, optional
        Directional light radius (DLR) scale factor applied to both R1 and
        R2 before matching, growing or shrinking every galaxy's ellipse
        (default 1.0, i.e. unchanged). A factor of ~1.5 is recommended to
        avoid missing transients hosted with large offsets.
    nbins : int, optional
        Number of geomspace bins in R1 (default 20).
    verbose : bool, optional
        Print per-bin progress.

    Returns
    -------
    gal_idx, pt_idx : numpy.ndarray of int
        Indices into gal_coord / pt_coord for every valid association.
    separation : numpy.ndarray of float
        Normalized elliptical separation for each association (0 at the
        galaxy center, 1 on the ellipse boundary); see `ellipse_separation`.
    """
    R1 = np.asarray(R1, dtype=float) * dlr_factor
    R2 = np.asarray(R2, dtype=float) * dlr_factor
    PA = np.asarray(PA, dtype=float)

    r1_pos = R1[R1 > 0]
    if r1_pos.size == 0:
        return (np.array([], dtype=int), np.array([], dtype=int),
                np.array([], dtype=float))

    bin_edges = np.geomspace(r1_pos.min(), R1.max(), nbins + 1)

    match_gal_idx = []
    match_pt_idx = []
    match_sep = []

    for i in range(nbins):
        lo, hi = bin_edges[i], bin_edges[i + 1]
        if i < nbins - 1:
            sel = (R1 >= lo) & (R1 < hi)
        else:
            sel = (R1 >= lo) & (R1 <= hi)

        if not np.any(sel):
            continue

        sub_idx = np.where(sel)[0]
        sub_coord = gal_coord[sub_idx]
        seplimit = hi * u.arcsec

        # NOTE: self.search_around_sky(other, seplimit) returns
        # (idx_into_other, idx_into_self, sep2d, dist3d) -- the order is
        # the reverse of what you might expect from the call syntax.
        idx_p, idx_g, sep2d, _ = sub_coord.search_around_sky(pt_coord, seplimit)

        if verbose:
            print(f"bin {i:2d}  R1 in [{lo:8.2f}, {hi:8.2f}] arcsec  "
                  f"n_gal={sel.sum():7d}  candidate pairs={len(idx_g)}")
            print(f"Elapsed time: {time.time()-t0:.0f}s")
            
        if len(idx_g) == 0:
            continue

        gal_idx = sub_idx[idx_g]
        sep = ellipse_separation(gal_coord[gal_idx], pt_coord[idx_p],
                                  R1[gal_idx], R2[gal_idx], PA[gal_idx])
        keep = sep <= 1.0

        match_gal_idx.append(gal_idx[keep])
        match_pt_idx.append(idx_p[keep])
        match_sep.append(sep[keep])

    if match_gal_idx:
        return (np.concatenate(match_gal_idx), np.concatenate(match_pt_idx),
                np.concatenate(match_sep))
    return np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=float)


def crossmatch_fits(gal_fits, pts_fits, output_fits,
                     gal_ra_col="gal_ra", gal_dec_col="gal_dec",
                     gal_r1_col="R1", gal_r2_col="R2", gal_pa_col="PA",
                     pt_ra_col="ra", pt_dec_col="dec",
                     dlr_factor=1.0, nbins=20, verbose=False):
    """
    Convenience wrapper: read two FITS tables, crossmatch, write result.

    Parameters
    ----------
    dlr_factor : float, optional
        Directional light radius (DLR) scale factor applied to both R1 and
        R2 before matching (default 1.0). A factor of ~1.5 is recommended
        to avoid missing transients hosted with large offsets.

    Returns the result Table (also written to output_fits).
    """
    t0 = time.time()
    gal = Table.read(gal_fits)
    pts = Table.read(pts_fits)
    # reduce search area if possible
    if pts[pt_ra_col].min()>60:
        gal = gal[gal[gal_ra_col]>pts[pt_ra_col].min()-10]
    if pts[pt_ra_col].max()<300:
        gal = gal[gal[gal_ra_col]<pts[pt_ra_col].max()+10]
    if pts[pt_dec_col].min()>-50:
        gal = gal[gal[gal_dec_col]>pts[pt_dec_col].min()-10]
    if pts[pt_dec_col].max()<50:
        gal = gal[gal[gal_dec_col]<pts[pt_dec_col].max()+10]
    
    gal_coord = SkyCoord(ra=np.asarray(gal[gal_ra_col]) * u.deg,
                          dec=np.asarray(gal[gal_dec_col]) * u.deg)
    pt_coord = SkyCoord(ra=np.asarray(pts[pt_ra_col]) * u.deg,
                         dec=np.asarray(pts[pt_dec_col]) * u.deg)

    R1 = np.asarray(gal[gal_r1_col], dtype=float)
    R2 = np.asarray(gal[gal_r2_col], dtype=float)
    PA = np.asarray(gal[gal_pa_col], dtype=float)
    
    if verbose:
        print(f"{len(gal)} galaxies, {len(pts)} points")
        print(f"Elapsed time: {time.time()-t0:.0f}s")

    gal_idx, pt_idx, separation = crossmatch_ellipses(
        gal_coord, pt_coord, R1, R2, PA,
        dlr_factor=dlr_factor, nbins=nbins, verbose=verbose, t0=t0)

    result = Table()
    result["gal_idx"] = gal_idx
    result["point_idx"] = pt_idx
    for col in gal.colnames:
        if col == gal_ra_col:
            out_col = "gal_ra"
        elif col == gal_dec_col:
            out_col = "gal_dec"
        else:
            out_col = col
        result[out_col] = gal[col][gal_idx]
    # reflect the dlr_factor-scaled radii actually used for matching
    result[gal_r1_col] = R1[gal_idx] * dlr_factor
    result[gal_r2_col] = R2[gal_idx] * dlr_factor
    for col in pts.colnames:
        out_col = f"pts_{col}" if col in result.colnames else col
        result[out_col] = pts[col][pt_idx]
    result["separation"] = separation
    result.meta["DLRFACT"] = dlr_factor

    result.write(output_fits, overwrite=True)
    if verbose:
        print(f"Total associations: {len(result)} -> written to {output_fits}")
        print(f"Elapsed time: {time.time()-t0:.0f}s")
        
    return result
