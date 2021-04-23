#!/usr/bin/env python3

import sys

QUANTUM_ENERGY_IDX=4

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
    return [line[QUANTUM_ENERGY_IDX] for line in block]

def read_quantum_file(fname):
    with open(fname, 'r') as fp:
        lines = fp.readlines()
        if lines:
            lines  = list(filter( lambda x : x != '\n', lines )) # remove empty lines
    return lines

def read_quantum_energies(fname):
    lines = read_quantum_file(fname)
    fprint_idxs = find_start_end_block(lines, "Fingerprints")
    quantum_idxs = find_start_end_block(lines, "Quantum")
    return extract_quantum_energies(
        extract_quantum_block(lines, quantum_idxs)
    )


print(read_quantum_energies(sys.argv[1]))

