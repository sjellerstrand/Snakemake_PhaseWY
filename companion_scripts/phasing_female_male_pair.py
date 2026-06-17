#! /usr/bin/env python

# Version 2025-06-05
# Author: Simon Jacobsen Ellerstrand
# Github: sjellerstrand

import os, time, argparse, re
import gzip

parser = argparse.ArgumentParser(description='phases sex-linked genotypes from one female-male pair in a vcf. First individual should be homogametic, second individual should be heterogametic.')
parser.add_argument('-i', '--input', dest='i', help="input file in freebayes vcf format [required]", required=True)
parser.add_argument('-o', '--output', dest='o', help="output file [required]", required=True)
args = parser.parse_args()

# Check if the input vcf file is gzipped
if args.i.endswith(".gz"):
    inputF = gzip.open(args.i, 'rt')
else:
    inputF = open(args.i, 'r')

outputF = open(args.o, 'w')

for Line in inputF:

    # Check if the header section
    if re.match('^#', Line) is None:

        # Get the columns of that line
        columns = Line.strip("\n").split("\t")

        # Add the info to the site
        result = columns[0:9]

        # Get the genotypes
        hetgam_gt0 = columns[9]
        homgam_gt0 = columns[10]

        # Check each genotype if it is diploid
        hetgam_gt = hetgam_gt0.split(":")
        homgam_gt = homgam_gt0.split(":")

        # Check genotype field
        if len(str(hetgam_gt[0])) != 1 and len(str(homgam_gt[0])) != 1:

            # Split genotypes
            if "/" in hetgam_gt[0]:
                hetgam_alleles = hetgam_gt[0].split("/")
            else:
                hetgam_alleles = hetgam_gt[0].split("|")

            if "/" in homgam_gt[0]:
                homgam_alleles = homgam_gt[0].split("/")
            else:
                homgam_alleles = homgam_gt[0].split("|")

            # Retrieve W and Z alleles
            homgamZ1 = homgam_alleles[0]
            homgamZ2 = homgam_alleles[1]

            if hetgam_alleles[0] == hetgam_alleles[1]:
                hetgamW = hetgam_alleles[0]
                hetgamZ = hetgam_alleles[0]

            elif hetgam_alleles[0] == homgamZ1 or hetgam_alleles[0] == homgamZ2 or hetgam_alleles[1] == homgamZ1 or hetgam_alleles[1] == homgamZ2:

                if hetgam_alleles[0] == homgamZ1 or hetgam_alleles[0] == homgamZ2:

                    if hetgam_alleles[1] == homgamZ1 or hetgam_alleles[1] == homgamZ2:
                        hetgamW = "."
                        hetgamZ = "."
                    else:
                        hetgamW = hetgam_alleles[1]
                        hetgamZ = hetgam_alleles[0]

                elif hetgam_alleles[1] == homgamZ1 or hetgam_alleles[1] == homgamZ2:
                    hetgamW = hetgam_alleles[0]
                    hetgamZ = hetgam_alleles[1]

            else:
                hetgamW = "."
                hetgamZ = "."

            # Extract FORMAT and genotype fields
            format_fields = columns[8].split(":")
            hetgam_fields = hetgam_gt0.split(":")

            # Initialize AD values
            ad_values = []
            if "AD" in format_fields:
                ad_index = format_fields.index("AD")
                try:
                    ad_values = hetgam_fields[ad_index].split(",")
                except:
                    ad_values = []

            # Default AD strings
            ad_stringW = ad_stringZ = "."

            # Try to parse allele indexes
            try:
                aW = int(hetgamW)
                aZ = int(hetgamZ)
            except:
                aW = aZ = None

            # Proceed only if AD values are valid
            if ad_values and aW is not None and aZ is not None:
                ad_len = len(ad_values)
                if aW == aZ and aW < ad_len:
                    # Homozygous: split count between both haploids
                    try:
                        count = float(ad_values[aW])
                        half = str(int(round(count / 2)))

                        ad_copyW = ["0"] * ad_len
                        ad_copyZ = ["0"] * ad_len
                        ad_copyW[aW] = half
                        ad_copyZ[aZ] = half

                        ad_stringW = ",".join(ad_copyW)
                        ad_stringZ = ",".join(ad_copyZ)
                    except:
                        pass
                elif aW < ad_len and aZ < ad_len:
                    # Heterozygous: assign only each allele's value
                    ad_copyW = ["0"] * ad_len
                    ad_copyZ = ["0"] * ad_len
                    ad_copyW[aW] = ad_values[aW]
                    ad_copyZ[aZ] = ad_values[aZ]

                    ad_stringW = ",".join(ad_copyW)
                    ad_stringZ = ",".join(ad_copyZ)

            # Compose haploid FORMAT strings
            hetgamW_str = f"{hetgamW}:.:.:{ad_stringW}"
            hetgamZ_str = f"{hetgamZ}:.:.:{ad_stringZ}"

            # Append phased hetgam genotypes and homgam
            result.extend([hetgamW_str, hetgamZ_str, homgam_gt0])

        else:
            # If haploid genotype
            hetgamW = "."
            hetgamZ = "."
            homgamZ1 = "."
            homgamZ2 = "."
            result.extend([".:.:.:.", ".:.:.:.", homgam_gt0])

        # Write line
        outputF.write('\t'.join(result) + "\n")

    elif re.match('^#CHROM', Line) is None:
        # Other header lines, write as-is
        outputF.write(Line)

    else:
        # Handle #CHROM line by inserting new sample names
        columns = Line.strip("\n").split("\t")
        result = columns[0:9]
        hetgam = columns[9]
        homgam = columns[10]
        hetgamW = hetgam + "_W"
        hetgamZ = hetgam + "_Z"
        result.extend([hetgamW, hetgamZ, homgam])
        outputF.write('\t'.join(result) + "\n")

inputF.close()
outputF.close()
