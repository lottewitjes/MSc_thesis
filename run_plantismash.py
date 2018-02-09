#!/usr/bin/evn/python

"""A Python script that runs plantiSMASH. You should run this script in the directory of plantiSMASH.

python witjes_run_plantismash.py  <gff3_plant_genome> <fasta_plant_genome> <expression_file> <output_directory>

Keyword arguments:
- gff3_plant_genome --> plant genome annotation in GFF3 format
- fasta_plant_genome --> plant genomic DNA (or just one chromosome) in FASTA format
- expression_file --> .soft file with expression data
- output_directory --> path to output directory, choose your directory name here
Returns:
- filled output_directory containing plantiSMASH results
"""

from __future__ import division
import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__version__ = "1.0"

def run_plantismash(gff3_plant_genome, fasta_plant_genome, expression_file, output_directory):
    """A function that runs plantiSMASH on the input files.

    Keyword arguments:
    - gff3_plant_genome --> plant genome annotation in GFF3 format
    - fasta_plant_genome --> plant genomic DNA (or just one chromosome) in FASTA format
    - expression_file --> .soft file with expression data
    - output_directory --> path to output directory, choose your directory name here
    Returns:
    - filled output_directory containing plantiSMASH results
    """
    cmd = "python run_antismash.py --taxon plants --gff3 {} --coexpress --coexpress-csv_file {} --clusterblast --knownclusterblast --min-domain-number 1 --cdh-cutoff 0.6 --outputfolder {} {}".format(gff3_plant_genome, expression_file, output_directory, fasta_plant_genome)
    if os.path.exists(output_directory):
        pass
    else:
        try:
            results = subprocess.check_call(cmd, shell=True)
            return results
        except subprocess.CalledProcessError as err:
            print err.output
            sys.exit()

if __name__ == "__main__":
    gff3_plant_genome = sys.argv[1]
    fasta_plant_genome = sys.argv[2]
    expression_file = sys.argv[3]
    output_directory = sys.argv[4]

    run_plantismash(gff3_plant_genome, fasta_plant_genome,expression_file, output_directory)





