#! /usr/bin/evn python

"""A Python script that parses marker chromosomal position to an xQTL file.

python marker_position_parser.py <xQTL_file> <marker_position_file>

Keyword arguments:
    xQTL_file - a tab-separated file containing xQTL, chr, peak_mb, inf_mb, sup_mb, and lod-score.
    marker_position_file - a tab-separated file containing markerID, chr, pos, gene1, and gene2.
Returns:
    parsed_xQTLs.txt - a tab-separated file as xQTL_file but with parsed correct genomic locations.
"""

from __future__ import division
from sys import argv
import subprocess
import os.path
from ~/MSc_thesis/overlap_xQTL_BGC import xQTL_parser

__author__ = "Lotte Witjes"
__email__= "lottewitjes@outlook.com"
__version__ = "1.0"
__date__ = "15 Jan 2018"

def marker_position_parser(marker_position_file):
    """A function that parses marker's chromosome and position into a dictionary.

    Keyword arguments:
        marker_position_file - a tab-separated file containing markerID, chr, pos, gene1, and gene2.
    Returns:
        thedic - a dictionary containing marker's ID as key, and chr and pos as values.
    """
    thedic = {}
    with open(marker_position_file, "r") as thefile:
        next(thefile) #skip the header
        for line in thefile:
            elements = line.split("\t")
            key, chr, position = elements[0], elements[1], elements[2]
            thedic[key] = [chr, position]
    return thedic

if __name__ == "__main__":
    #Get files from command line
    mQTL_file = argv[1]
    marker_position_file = argv[2]

    #Parse files
    mQTL_list = xQTL_parser(mQTL_file)
    print mQTL_list
    marker_position_dic = marker_position_parser(marker_position_file)
    #print marker_position_dic

    #Substitute xQTL_file positions with correct positions

    #Write correct parsed output
