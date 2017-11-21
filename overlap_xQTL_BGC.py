#!/usr/bin/evn/python
"""A Python script that overlaps xQTL with BGCs found with plantiSMASH.

python overlap_xQTL_BGC.py <BGC_file> <eQTL_file> <mQTL_file>

Keyword arguments:
    BGC_file --> A file containing the BGCs
    eQTL_file --> A file containing the eQTLs
    mQTL_file --> A file containing the mQTLs

Returns:
    A file with xQTL per BGC
"""

from __future__ import division
from sys import argv
import subprocess
import os.path
import re

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "21 Nov 2017"
__version__ = "1.0"

def xQTL_parser(xQTL_file):
    thefile = open(xQTL_file, "r")
    thelist = []
    next(thefile) #skip the header
    for line in thefile:
        elements = line.split()
        thelist.append(elements)
    return thelist

def BGC_parser(BGC_file):
    thefile = open(BGC_file, "r")
    thedic = {}
    next(thefile) #skip the header
    for line in thefile:
        elements = line.split("\t")
        chr, cluster_id = elements[0].split("_c")
        type = elements[1]
        from_bp, to_bp = elements[3].split(";")
        genes = elements[4].split(";")
        genes = [re.sub(r"-.*", "", gene) for gene in genes]
        thedic[cluster_id] = [type, chr, from_bp, to_bp, genes]
    return thedic

if __name__ == "__main__":
    #Get files from command line
    BGC_file = argv[1]
    eQTL_file = argv[2]
    mQTL_file = argv[3]

    #Parse the files
    BGC_dic = BGC_parser(BGC_file)
    eQTL_list = xQTL_parser(eQTL_file)
    mQTL_list = xQTL_parser(mQTL_file)
