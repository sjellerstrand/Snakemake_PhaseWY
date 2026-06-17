#! /usr/bin/env python

# Version 2025-08-27
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

from sys import *
import os, time, argparse, re
import msprime, tskit, pyslim
import warnings
import numpy as np
from textwrap import wrap

parser = argparse.ArgumentParser(description='Simulates mutations on a tree sequences file')
parser.add_argument('-i', '--input', dest='i', help="input tree sequence file", required=True)
parser.add_argument('-s', '--seed', dest='s', type=int, help="seed", required=True)
parser.add_argument('-r', '--recrate', dest='r', type=float, help="recombination rate", required=True)
parser.add_argument('-m', '--mutrate', dest='m', type=float, help="mutation rate", required=True)
parser.add_argument('-N', '--popsize', dest='N', type=int, help="ancestral population size", required=True)
parser.add_argument('-g', '--gen', dest='g', type=int, required=True, help="Generation number to include in FASTA header")
parser.add_argument('-o', '--output', dest='o', help="output file [required]", required=True)
parser.add_argument('-x', '--fasta-template', dest='x', required=True, help="Input ancestral FASTA file")
parser.add_argument('-z', '--fasta-out', dest='z', required=True, help="Output FASTA file")
args = parser.parse_args()

input = args.i
seed = args.s
recrate = args.r
mutrate = args.m
popsize = args.N
gen = args.g
output = args.o
fasta_template = args.x
fasta_out = args.z

# Load the .trees file
ts = tskit.load(input)    # no simplify!

# Suppress warning about time-unit mismatch
warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

# Recapitate!
recap = pyslim.recapitate(ts, ancestral_Ne=popsize, recombination_rate=recrate, random_seed=seed)

# Remove old mutations, such as sex-determining locus (otherwise this might fail)
tables = recap.dump_tables()
tables.sites.clear()
tables.mutations.clear()
recap = tables.tree_sequence()

# Overlay mutations
mutated = msprime.sim_mutations(recap, rate=mutrate, random_seed=seed)
mutated.dump(output)

# Calculate nucleotide diversity
windows = np.linspace(0, mutated.sequence_length, num=3)
pi = mutated.diversity(windows=windows, mode="site")
pi_output_path = output.replace("overlaid.trees", "pi.txt")

with open(pi_output_path, "w") as f:
    f.write("start\tend\tpi\n")
    for i in range(len(pi)):
        f.write(f"{windows[i]}\t{windows[i+1]}\t{pi[i]}\n")

# Load ancestral FASTA and modify sequence
with open(fasta_template) as f:
    lines = f.readlines()
    header = lines[0].strip()
    if header.startswith(">"):
        header = header[1:]
    new_header = f">" + header + f"_gen_{gen}"
    sequence = "".join(line.strip() for line in lines[1:])
    sequence = list(sequence)

# Apply mutations
for site in mutated.sites():
    if site.mutations:
        pos = int(site.position)
        if 0 <= pos < len(sequence):
            ancestral = site.ancestral_state
            sequence[pos] = ancestral

# Write modified FASTA
with open(fasta_out, "w") as f:
    f.write(f"{new_header}\n")
    for line in wrap("".join(sequence), 60):
        f.write(f"{line}\n")
