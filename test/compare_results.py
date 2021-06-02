#!/usr/bin/env python3

import sys
import os
import argparse
import math

QUANTUM_ENERGY_IDX=4
INTENSITY_INDICES = {
    "einstein": 2,
    "nu": 3,
}

def find_start_end_block(lines, blockname):
    """Identify start and end of blocks in .chk files"""
    for i, line in enumerate(lines):
        if "Start " + blockname in line:
            start_idx = i
        elif "End " + blockname in line:
            end_idx = i
    return start_idx, end_idx

def find_log_block(lines, blockname):
    """Identify start and end of block in log file"""
    start_idx = end_idx = None
    for i, line in enumerate(lines):
        if blockname in line:
            start_idx = i
            break

    if start_idx is None:
        raise Exception(f"{blockname} not found")

     # find first "done" after we've found blockname
    for i, line in enumerate(lines[start_idx:]):
        if "done" in line:
            end_idx = i + start_idx
            break

    if end_idx is None:
        raise Exception(f"Could not find end of {blockname}")

    return start_idx, end_idx

def read_chk_file(fname):
    """Read checkpoint file as a list of lines"""
    with open(fname, 'r') as fp:
        lines = fp.readlines()
    return lines

def strip_newlines(lines):
    """Remove all lines which are just newlines"""
    return [line for line in lines if line != '\n']

def strip_comments(lines):
    """Remove all lines with comments"""
    return [line for line in lines if '<-' not in line]

def extract_column(lines, column_no):
    """Extract column of numbers from list of str lines"""
    return [float(line.split()[column_no]) for line in lines]

def read_energy_column(fname, column_no):
    """Extract energies from a column in file fname"""
    lines = read_chk_file(fname)
    lines = strip_newlines(strip_comments(lines))
    # remove last line (which is not part of the actual data)
    lines = lines[:-1]
    return extract_column(lines, column_no)

def read_quantum_energies(fname):
    """Extract quantum energies from fname"""
    lines = read_chk_file(fname)
    lines = strip_newlines(strip_comments(lines))
    start, end = find_start_end_block(lines, "Quantum")
    # take out first 4 lines and last line of block
    lines = lines[start+4:end]
    return extract_column(lines, QUANTUM_ENERGY_IDX)

def read_intensity_column(fname, column_name):
    """Extract quantum energies from fname"""
    assert column_name in INTENSITY_INDICES.keys(), f"Intensity column name must be one of {INTENSITY_INDICES.keys()}"

    lines = read_chk_file(fname)
    lines = strip_newlines(lines)
    start, end = find_log_block(lines, "Linestrength")
    # take out first and last line of block
    lines = lines[start+1:end]

    return extract_column(lines, INTENSITY_INDICES[column_name])

def compare_columns(col1, col2, abs_precision=0.0, rel_precision=1e-10):
    """Compare two columns of numbers to a given absolute or relative precision"""
    difference_exists = False
    for i, (e1, e2) in enumerate(zip(col1, col2)):
        if not math.isclose(e1, e2, abs_tol=abs_precision, rel_tol=rel_precision):
            difference_exists = True
            print(f"{e1} and {e2} differ by {abs(e1-e2)} at index {i}")

    assert difference_exists == False

def compare_energy_files(fname1, fname2, column_no, precision=1e-10):
    """Compare two energy files"""
    energies1 = read_energy_column(fname1, column_no)
    energies2 = read_energy_column(fname2, column_no)
    compare_columns(energies1, energies2, abs_precision=precision)

def compare_quantum_files(fname1, fname2, precision=1e-10):
    """Compare two files in quantum form"""
    energies1 = read_quantum_energies(fname1)
    energies2 = read_quantum_energies(fname2)
    # Note, this uses the more accurate rel_precision
    compare_columns(energies1, energies2, rel_precision=precision)

def compare_intensity_files(fname1, fname2, precision=1e-10):
    """Compare two files in quantum form"""
    for col_name in INTENSITY_INDICES.keys():
        col1 = read_intensity_column(fname1, col_name)
        col2 = read_intensity_column(fname2, col_name)
        compare_columns(col1, col2, abs_precision=precision)

def main():
    parser = argparse.ArgumentParser(description='Compare output files from TROVE')
    parser.add_argument('filenames', nargs='+',
                        help='names of files to compare in each folder')
    parser.add_argument('--folder1', required=True,
                        help='first folder to compare')
    parser.add_argument('--folder2', required=True,
                        help='second folder to compare')
    parser.add_argument('--kind', choices=['quantum', 'column', 'intensity'], required=True, help='type of file')
    parser.add_argument('--column', default=-1, type=int, help='column to compare when \'column\' is supplied to --kind (index starts at 0)')
    parser.add_argument('--precision', default=1e-10, type=float, help='relative precision of which two values must differ by to be considered nonequal')

    args = parser.parse_args()

    filelist=args.filenames

    folder1=args.folder1
    folder2=args.folder2

    if args.kind == 'column' and args.column == -1:
        raise argparse.ArgumentError('column number must be supplied when using --kind=column')

    if args.kind == 'column':
        for fname in filelist:
            compare_energy_files(os.path.join(folder1, fname), os.path.join(folder2, fname), args.column, precision=args.precision)
    elif args.kind == 'quantum':
        for fname in filelist:
            try:
                compare_quantum_files(folder1 + "/" + fname, folder2 + "/" + fname, precision=args.precision)
            except IndexError:
                print(folder1 + "/" + fname, "does not match regular quantum file format")
    elif args.kind == 'intensity':
        for fname in filelist:
            compare_intensity_files(folder1 + "/" + fname, folder2 + "/" + fname, precision=args.precision)

    exit(0)

if __name__ == '__main__':
    main()
