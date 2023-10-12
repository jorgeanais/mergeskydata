from astropy.table import Table, hstack, vstack
from astropy.table.column import MaskedColumn, Column
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import numpy.typing as npt

def _table_to_array(t: Table, columns: list[str]) -> npt.NDArray:
    """Convert columns of an astropy table into a numpy array"""
    
    output_array = np.zeros(shape=(len(t), len(columns)))
    for i, col in enumerate(columns):
        output_array[:, i] = t[col].data.data
    
    return output_array

def _table_to_marray(t: Table, columns: list[str]) -> tuple[npt.NDArray, npt.NDArray]:
    """Convert columns of an astropy table into a numpy masked array"""
    
    output_array = np.zeros(shape=(len(t), len(columns)))
    output_mask = np.ma.empty(shape=(len(t), len(columns)))
    for i, col in enumerate(columns):
        output_array[:, i] = t[col].data.data
        output_mask[:, i] = t[col].data.mask
    
    return output_array, output_mask



def match_tables(
    table1: Table,
    table2: Table,
    tol_arcsec: float = 0.34,
) -> Table:
    """Match two tables based on the values of two columns.

    Parameters
    ----------
    table1 : Table
        First table.
    table2 : Table
        Second table.
    tol_arcsec : float, optional
        Tolerance in arcseconds, by default 0.34.

    Returns
    -------
    Table
        Table with the matched rows.

    """

    
    # Check if the tables have the same columns
    if set(table1.colnames) != set(table2.colnames):
        raise ValueError("The tables have different columns")
    
    magnitude_colnames = [col for col in table1.colnames if "mag" in col]
    error_colnames = [col for col in table1.colnames if "err" in col]

    c1 = SkyCoord(
        ra=table1["ra"],
        dec=table1["dec"],
        unit=(u.deg, u.deg)
    )

    c2 = SkyCoord(
        ra=table2["ra"],
        dec=table2["dec"],
        unit=(u.deg, u.deg)
    )

    idx, d2d, _ = c1.match_to_catalog_sky(c2)
    matched_indx = d2d < tol_arcsec * u.arcsec

    # Rows from c2 that match
    matched_indx_c2 = np.in1d(np.arange(len(c2)), idx[matched_indx])

    # Add unmatched rows from both catalgos to the output table
    unmatched_sources = vstack([table1[~matched_indx], table2[~matched_indx_c2]])

    # Combine matched rows
    aux_matched_sources = hstack(
        [
            table1[matched_indx],
            table2[idx[matched_indx]]
        ],
        table_names=["a", "b"]
    )

    # Missing values treatment. Magnitudes are computed using 
    # columns = ["mag_1_a", "mag_1_b"]
    # x = _table_to_array(aux_matched_sources, ["mag_1_a", "mag_1_b"])
    # w = 1./np.power(_table_to_array(aux_matched_sources, ["er_1_a", "er_1_a"]),2)
    # np.average(x, weights=w, axis=1)


    # Using masked arrays
    columns = ["mag_1_a", "mag_1_b"]

    x, mask = _table_to_marray(aux_matched_sources, ["mag_1_a", "mag_1_b"])
    x_masked = np.ma.masked_array(x, mask=mask)
    w = 1./np.power(_table_to_array(aux_matched_sources, ["er_1_a", "er_1_b"]),2)

    result = np.ma.average(x_masked, weights=w, axis=1)

    # Error calculation
    xer, xer_mask = _table_to_marray(aux_matched_sources, ["er_1_a", "er_1_b"])
    result_error = np.sqrt(1. / np.sum(np.ma.masked_array(np.power(1. / xer, 2), mask=xer_mask), axis=1))
    

    aux_matched_sources["result"] = result
    aux_matched_sources["result_error"] = result_error
    aux_matched_sources.write("partialresults.fits", overwrite=True)