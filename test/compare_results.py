#!/usr/bin/env python3

import sys
import argparse
from pytest import approx

QUANTUM_ENERGY_IDX=4
ENERGY_DIFF_THRESHOLD=1e-10 # How different can two energies be?

def find_start_end_block(lines, blockname):
    for i, line in enumerate(lines):
        if "Start " + blockname in line:
            start_idx = i
        elif "End " + blockname in line:
            end_idx = i
    return start_idx, end_idx

def extract_quantum_block(lines, idxs):
    temp = lines[idxs[0]+4:idxs[-1]]
    return [line.split() for line in temp]

def extract_quantum_energies(block):
    return [float(line[QUANTUM_ENERGY_IDX]) for line in block]

def read_quantum_file(fname):
    with open(fname, 'r') as fp:
        lines = fp.readlines()
        if lines:
            lines  = list(filter( lambda x : x != '\n', lines )) # remove empty lines
    return lines

def read_quantum_block(fname):
    lines = read_quantum_file(fname)
    quantum_idxs = find_start_end_block(lines, "Quantum")
    return extract_quantum_block(lines, quantum_idxs)

def compare_quantum_files(fname1, fname2):
    energy_block1 = read_quantum_block(fname1)
    energy_block2 = read_quantum_block(fname2)

    energies1 = extract_quantum_energies(energy_block1)
    energies2 = extract_quantum_energies(energy_block2)

    assert energies1 == approx(energies2, rel=ENERGY_DIFF_THRESHOLD)

def main():
    parser = argparse.ArgumentParser(description='Compare output files from TROVE')
    parser.add_argument('filenames', nargs='+',
                        help='names of files to compare in each folder')
    parser.add_argument('--folder1', require=True,
                        help='first folder to compare')
    parser.add_argument('--folder2', require=True,
                        help='second folder to compare')

    args = parser.parse_args()

    quantum_filelist=args.filenames

    folder1=args.folder1
    folder2=args.folder2

    for fname in quantum_filelist:
        try:
            compare_quantum_files(folder1 + "/" + fname, folder2 + "/" + fname)
        except IndexError:
            print(folder1 + "/" + fname, "does not match regular quantum file format")
