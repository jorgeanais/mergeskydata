import click
from pathlib import Path

from astropy.table import Table, hstack, vstack

"""
Program that merges two fits tables into one.
Example:
    python main.py \
    -f1 file1.fits \
    -f2 file2.fits \
    -o results.fits
"""

LMIN = 0.0
LMAX = 10.8
BMIN = -10.3
BMAX = -3.0

@click.command()
@click.option("-f1", help="File of the input fits table 1", required=True)
@click.option("-f2", help="File of the input fits table 2", required=True)
@click.option("-o", help="Output file for the merged table", required=True)
def main(f1, f2, o):
    file1 = Path(f1)
    file2 = Path(f2)
    output_file = Path(o)
    
    table1 = Table.read(file1)
    table2 = Table.read(file2)
    
    # remove objects outside roi
    mask_lmin = table2["l"] >= LMIN
    mask_lmax = table2["l"] < LMAX
    mask_bmin = table2["b"] >= BMIN
    mask_bmax = table2["b"] < BMAX
    mask = mask_lmin * mask_lmax * mask_bmin * mask_bmax
    table2 = table2[mask]
    
    
    print(f"Processing tables {file1} and {file2}...")
    merged_table = vstack([table1, table2])
    merged_table.write(output_file, format="fits", overwrite=True)
    

if __name__ == "__main__":
    main()



