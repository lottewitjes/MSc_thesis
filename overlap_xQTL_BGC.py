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

def xQTL_parser(xQTL_file):
    thefile = open(xQTL_file, "r")
    thelist = []
    for line in thefile:
        elements = line.split()
        thelist.append(elements)
    return thelist

def BGC_parser(BGC_file):
    thefile = open(BGC_file, "r")
    thedic = {}
    for line in thefile:
        elements = line.split()
        thedic[elements[0]] = [elements[1], elements[2], elements[3], elements[4], elements[5], elements[6], elements[7], elements[8], elements[9]]
    return thedic

if __name__ == "__main__":
    BGC_file = argv[1]
    eQTL_file = argv[2]
    mQTL_file = argv[3]

    #GC_list = BGC_parser(BGC_file)
    #print BGC_list
    eQTL_list = xQTL_parser(eQTL_file)
    #print eQTL_list
    mQTL_list = xQTL_parser(mQTL_file) 
    print mQTL_list
