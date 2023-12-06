import click
from pathlib import Path

from astropy.table import Table


"""
Simple script to transform a .cals file into a .fits file.
"""

@click.command()
@click.option("-f", help="File of the input cals table", required=True)
@click.option("-o", help="Output file for the merged table", required=False)
def main(f, o=None):
    """
    Program that merges two fits tables containing photometrical data into one table.
    """
    input_file = Path(f)

    if not o:
        o = input_file.with_suffix(".fits")
    output_file = Path(o)

    print(f"Reading table {input_file}...")
    table = Table.read(input_file, format="ascii")

    print(f"Writing table {output_file}...")
    table.write(output_file, format="fits", overwrite=True)


if __name__ == "__main__":
    main()