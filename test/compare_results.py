#!/usr/bin/env python3

import sys
import os
import argparse
from pytest import approx

QUANTUM_ENERGY_IDX=4

def find_start_end_block(lines, blockname):
    """Identify start and end of blocks"""
    for i, line in enumerate(lines):
        if "Start " + blockname in line:
            start_idx = i
        elif "End " + blockname in line:
            end_idx = i
    return start_idx, end_idx

def extract_quantum_block(lines):
    """Extract only the quantum block from a chk file"""
    idxs = find_start_end_block(lines, "Quantum")
    temp = lines[idxs[0]+4:idxs[-1]]
    return [line.split() for line in temp]

def extract_quantum_energies(block):
    """Extract quantum energies from entire block"""
    return [float(line[QUANTUM_ENERGY_IDX]) for line in block]

def read_chk_file(fname):
    """Read checkpoint file as a list of lines"""
    with open(fname, 'r') as fp:
        # Strip newlines and lines with comments
        lines = [line for line in fp.readlines() if line != '\n' and '<-' not in line]
    return lines

def read_energy_column(fname, column_no):
    """Extract energies from a column in file fname"""
    lines = read_chk_file(fname)
    lines = lines[:-1] # remove last line (which is not part of the actual data)
    return [float(line.split()[column_no]) for line in lines]

def read_quantum_block(fname):
    """Extract quantum energies from fname"""
    lines = read_chk_file(fname)
    return extract_quantum_block(lines)

def compare_columns(fname1, fname2, column_no, precision=1e-10):
    """Compare two energy files"""
    energies1 = read_energy_column(fname1, column_no)
    energies2 = read_energy_column(fname2, column_no)

    for e1, e2 in zip(energies1, energies2):
        assert e1 == approx(e2, abs=precision), \
            f"{e1} and {e2} differ by {abs(e1-e2)}"

def compare_quantum_files(fname1, fname2, precision=1e-10):
    """Compare two files in quantum form"""
    energy_block1 = read_quantum_block(fname1)
    energy_block2 = read_quantum_block(fname2)

    energies1 = extract_quantum_energies(energy_block1)
    energies2 = extract_quantum_energies(energy_block2)

    assert energies1 == approx(energies2, rel=precision)

def main():
    parser = argparse.ArgumentParser(description='Compare output files from TROVE')
    parser.add_argument('filenames', nargs='+',
                        help='names of files to compare in each folder')
    parser.add_argument('--folder1', required=True,
                        help='first folder to compare')
    parser.add_argument('--folder2', required=True,
                        help='second folder to compare')
    parser.add_argument('--kind', choices=['quantum', 'column'], required=True, help='type of file')
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
            compare_columns(os.path.join(folder1, fname), os.path.join(folder2, fname), args.column, precision=args.precision)
    elif args.kind == 'quantum':
        for fname in filelist:
            try:
                compare_quantum_files(folder1 + "/" + fname, folder2 + "/" + fname, precision=args.precision)
            except IndexError:
                print(folder1 + "/" + fname, "does not match regular quantum file format")

    exit(0)

if __name__ == '__main__':
    main()
