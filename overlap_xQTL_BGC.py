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

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "21 Nov 2017"
__version__ = "1.0"

def xQTL_parser():
    thefile = open(xQTL_file, "r")

def BGC_parser():
    thefile = open(BGC_file, "r")

if __name__ == "__main__":
BGC_file = argv[1]
eQTL_file = argv[2]
mQTL_file = argv[3]

BGC_list = BGC_parser(BGC_file)
eQTL_list = xQTL_parser(eQTL_file)
mQTL_list = xQTL_parser(mQTL_file) 
