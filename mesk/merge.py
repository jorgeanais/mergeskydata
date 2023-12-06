from astropy.table import Table, hstack, vstack
from astropy.table.column import MaskedColumn, Column
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import numpy.typing as npt


MATCHED_COLUMN_SUFFIX = ["a", "b"]


def _table_to_array(t: Table, columns: list[str]) -> npt.NDArray:
    """Convert columns of an astropy table into a numpy array"""

    output_array = np.zeros(shape=(len(t), len(columns)))
    for i, col in enumerate(columns):
        output_array[:, i] = t[col].data.data

    return output_array


def _table_to_marray(t: Table, columns: list[str]) -> tuple[npt.NDArray, npt.NDArray]:
    """Convert columns of an astropy table into a numpy array and a mask"""

    output_array = np.zeros(shape=(len(t), len(columns)))
    output_mask = np.ma.empty(shape=(len(t), len(columns)))
    for i, col in enumerate(columns):
        output_array[:, i] = t[col].data.data
        output_mask[:, i] = t[col].data.mask

    return output_array, output_mask


def compute_weighted_magnitude_and_error(
    table: Table, magnitude_cols: list[str, str], error_cols: list[str, str]
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Compute the weighted average of two magnitudes from a astropy table.

    | mag-A | error-mag-A | mag-B | error-mag-B |
    |-------|-------------|-------|-------------|
    |  ...  |     ...     |  ...  |     ...     |

    Parameters
    ----------
    table : Table
        Table with the magnitudes.
    magnitude_colnames : list[str, str]
        Names of the columns with the magnitudes.
    error_colnames : list[str, str]
        Names of the columns with the errors.

    Returns
    -------
    tuple[Column, Column]
        Weighted average of the magnitudes and the errors.

    """

    magnitudes, m_mask = _table_to_marray(table, magnitude_cols)
    magerrors, e_mask = _table_to_marray(table, error_cols)

    # First, compute the weighted average of the magnitudes
    masked_magnitudes = np.ma.masked_array(magnitudes, mask=m_mask)
    weights = 1.0 / np.power(magerrors, 2.0)
    weighted_magnitude = np.ma.average(masked_magnitudes, weights=weights, axis=1)

    # Second, compute the error of the weighted average
    weighted_magnitude_error = np.sqrt(
        1.0
        / np.sum(np.ma.masked_array(np.power(1.0 / magerrors, 2), mask=e_mask), axis=1)
    )

    return weighted_magnitude, weighted_magnitude_error


def match_tables(
    table1: Table,
    table2: Table,
    tol_arcsec: float = 0.34,
) -> tuple[Table, Table]:
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
    tuple[Table, Table]
        Tuple containing tables with unmatched sources and matched sources.

    """


    # Check if the tables have the same columns
    if set(table1.colnames) != set(table2.colnames):
        raise ValueError("The tables have different columns")

    c1 = SkyCoord(ra=table1["ra"], dec=table1["dec"], unit=(u.deg, u.deg))

    c2 = SkyCoord(ra=table2["ra"], dec=table2["dec"], unit=(u.deg, u.deg))

    idx, d2d, _ = c1.match_to_catalog_sky(c2)
    matched_indx = d2d < tol_arcsec * u.arcsec

    # Rows from c2 that match
    matched_indx_c2 = np.in1d(np.arange(len(c2)), idx[matched_indx])

    # Add unmatched rows from both catalgos to the output table
    unmatched_sources = vstack([table1[~matched_indx], table2[~matched_indx_c2]])

    # Combine matched rows
    matched_sources = hstack(
        [table1[matched_indx], table2[idx[matched_indx]]],
        table_names=MATCHED_COLUMN_SUFFIX,
    )

    return unmatched_sources, matched_sources


def avg_positions(
    table: Table,
    lon: str = "ra",
    lat: str = "dec",
) -> Table:
    """Average the positions (lon, lat) from matched sources."""

    def _get_average_value(
        table: Table,
        col: str,
    ):
        """Get the average value of two columns from a table (small angle aproximation)"""
        x_deg = _table_to_array(
            table,
            [f"{col}_{MATCHED_COLUMN_SUFFIX[0]}", f"{col}_{MATCHED_COLUMN_SUFFIX[1]}"],
        )
        x_rad = np.deg2rad(x_deg)
        x_rad_mean = np.average(x_rad, axis=1)
        return np.rad2deg(x_rad_mean)

    table = table.copy()
    table[lon] = _get_average_value(table, lon)
    table[lat] = _get_average_value(table, lat)

    table.remove_columns(
        [
            f"{lon}_{MATCHED_COLUMN_SUFFIX[0]}",
            f"{lon}_{MATCHED_COLUMN_SUFFIX[1]}",
            f"{lat}_{MATCHED_COLUMN_SUFFIX[0]}",
            f"{lat}_{MATCHED_COLUMN_SUFFIX[1]}",
        ]
    )

    return table


def avg_magnitudes(
    table: Table,
    band_names: list[str],
    error_names: list[str],
) -> Table:
    """
    Combine the magnitudes from matched sources for different bands using weighted average.

    Parameters
    ----------

    table : Table
        Table with the matched sources.
    band_names : list[str]
        Names of the bands.
    error_names : list[str]
        Names of the errors.

    Returns
    -------
    Table
        Table with the combined magnitudes.
    """

    for band_name, error_name in zip(band_names, error_names):
        w_mag, w_error = compute_weighted_magnitude_and_error(
            table,
            magnitude_cols=[
                f"{band_name}_{MATCHED_COLUMN_SUFFIX[0]}",
                f"{band_name}_{MATCHED_COLUMN_SUFFIX[1]}",
            ],
            error_cols=[
                f"{error_name}_{MATCHED_COLUMN_SUFFIX[0]}",
                f"{error_name}_{MATCHED_COLUMN_SUFFIX[1]}",
            ],
        )

        # Remove the columns with the magnitudes and errors
        table.remove_columns(
            [
                f"{band_name}_{MATCHED_COLUMN_SUFFIX[0]}",
                f"{band_name}_{MATCHED_COLUMN_SUFFIX[1]}",
            ]
        )
        table.remove_columns(
            [
                f"{error_name}_{MATCHED_COLUMN_SUFFIX[0]}",
                f"{error_name}_{MATCHED_COLUMN_SUFFIX[1]}",
            ]
        )

        # Add the new columns with the weighted magnitudes and errors
        table[f"{band_name}"] = w_mag
        table[f"{error_name}"] = w_error

    return table


def merge_tables(
    table1: Table,
    table2: Table,
) -> Table:

    # Get column names
    magnitude_colnames = [col for col in table1.colnames if "mag" in col]
    error_colnames = [col for col in table1.colnames if "er" in col]

    # Match and combine tables
    unmatched_sources_t, matched_sources_t = match_tables(table1, table2)

    matched_sources_t = avg_magnitudes(
        table=matched_sources_t,
        band_names=magnitude_colnames,
        error_names=error_colnames,
    )

    matched_sources_t = avg_positions(matched_sources_t)

    return vstack([unmatched_sources_t, matched_sources_t])
