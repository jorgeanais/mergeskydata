import click
from pathlib import Path

from astropy.table import Table

from mesk.merge import merge_tables

"""
Program that merges two fits tables containing photometrical data into one table.
Example:
    python main.py \
    -f1 /run/media/jorge/BLUE/DATA/phd/vvvx/jhk0512.fits \
    -f2 /run/media/jorge/BLUE/DATA/phd/vvvx/jhk0511.fits \
    -o restuls.fits
"""


@click.command()
@click.option("-f1", help="File of the input fits table 1", required=True)
@click.option("-f2", help="File of the input fits table 2", required=True)
@click.option("-o", help="Output file for the merged table", required=True)
def main(f1, f2, o):
    """
    Program that merges two fits tables containing photometrical data into one table.
    """

    file1 = Path(f1)
    file2 = Path(f2)
    output_file = Path(o)

    table1 = Table.read(file1)
    table2 = Table.read(file2)

    merged_table = merge_tables(table1, table2)
    merged_table.write(output_file, format="fits")


if __name__ == "__main__":
    main()
