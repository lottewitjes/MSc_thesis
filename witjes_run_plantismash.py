#!/usr/bin/evn/python

"""A Python script that runs plantiSMASH. You should run this script in the directory of plantiSMASH.

python witjes_run_plantismash.py  <gff3_plant_genome> <fasta_plant_genome> <output_directory >

Keyword arguments:
- gff3_plant_genome --> plant genome annotation in GFF3 format
- fasta_plant_genome --> plant genomic DNA (or just one chromosome) in FASTA format
- output_directory --> path to output directory, choose your directory name here
Returns:
- filled output_directory containing plantiSMASH results
"""

__author__ == "Lotte Witjes"
__version__ == "1.0"
__email__ == "lotte.witjes@wur.nl"

from __future__ import division
from sys import argv
import subprocess
import os.path

def run_plantismash(gff3_plant_genome, fasta_plant_genome, output_directory):
    """A function that runs plantiSMASH on the input files.

    Keyword arguments:
    - gff3_plant_genome --> plant genome annotation in GFF3 format
    - fasta_plant_genome --> plant genomic DNA (or just one chromosome) in FASTA format
    - output_directory --> path to output directory, choose your directory name here
    Returns:
    - filled output_directory containing plantiSMASH results
    """
    cmd = "python run_plantismash.py --taxon plants --gff3 {} --clusterblast --knownclusterblast --outputfolder {} {}".format(gff3_plant_genome, output_directory, fasta_plant_genome)
    if os.path.exists(output_directory):
        pass
    else:
        try:
            results = subprocess.check_call(cmd, shell=True)
            return results
        except subprocess.CalledProcessError as err:
            print err.strerror
            sys.exit()

if __name__ == "__main__":
gff3_plant_genome = argv[1]
fasta_plant_genome = argv[2]
output_directory = argv[3]

run_plantismash(gff3_plant_genome, fasta_plant_genome, output_directory)





