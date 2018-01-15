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
from overlap_xQTL_BGC import xQTL_parser

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

def substitute_xQTL_position(xQTL_list, marker_position_dic):
    """A function that substitutes xQTL positions in the xQTL_list with the correct positions from the marker_position_dic.

    Keyword arguments:
        xQTL_list - a list of lists containing xQTL, chr, peak_mb, inf_mb, sup_mb, and lod-score per xQTL.
        marker_position_dic - a dictionary containing marker's ID as key, and chr and pos as value.
    Returns:
        parsed_xQTL_list - same as xQTL_list but with correct positions.
    """
    parsed_xQTL_list = xQTL_list[:]
    for xQTL in parsed_xQTL_list:
        if str(xQTL[1]) in marker_position_dic:
            chr = int(marker_position_dic[str(xQTL[1])][0])
            peak_mb = float(int(marker_position_dic[str(xQTL[1])][1])/1000000)
            inf_mb = float((int(marker_position_dic[str(xQTL[1])][1])-1)/1000000)
            sup_mb = float((int(marker_position_dic[str(xQTL[1])][1])+1)/1000000)
            xQTL[1], xQTL[2], xQTL[3], xQTL[4] = chr, peak_mb, inf_mb, sup_mb
        else:
            print "This is an unknown marker."
    return parsed_xQTL_list

def write_file(output_name, parsed_xQTL_list):
    """A function that write the output of substitute_xQTL_position() to a .tsv.

    Keyword arguments:
        output_name - name of output file. 
        parsed_xQTL_list - same as xQTL_list but with correct positions.
    Returns:
        output_name.txt - a tab-separated file containing the values of parsed_xQTL_list.
    """
    with open(output_name, "w") as thefile:
        thefile.write("metabolite\tchr\tpeak_mb\tinf_mb\tsup_mb\tlod\n")
        for xQTL in parsed_xQTL_list:
            xQTL = [str(element) for element in xQTL]
            line = "\t".join(xQTL)
            thefile.write(line + "\n")

if __name__ == "__main__":
    #Get files from command line
    mQTL_file = argv[1]
    marker_position_file = argv[2]
    output_name = "parsed_xQTLs.txt"

    #Parse files
    mQTL_list = xQTL_parser(mQTL_file)
    marker_position_dic = marker_position_parser(marker_position_file)

    #Substitute xQTL_file positions with correct positions
    parsed_mQTL_list = substitute_xQTL_position(mQTL_list, marker_position_dic)

    #Write correct parsed output
    write_file(output_name, parsed_mQTL_list)
