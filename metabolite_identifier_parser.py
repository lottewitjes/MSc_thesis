#! /usr/bin/evn python
"""A Python script that changes metabolite identifiers in a mQTL file based on another file into molecule name
   or previousID_molweight_rettime.

python metabolite_identifier_parser.py <mQTL_file> <metaboliteID_file>

Keyword arguments:
    mQTL_file - a .tsv containing metaboliteID, chr, peak_mb, inf_mb, sup_mb and lod_score
    metaboliteID_file - a .tsv containing metaboliteID, molecule name, molecular weight (Da) and retention time (min)
Returns:
    parsed_mQTLs.txt - same as mQTL_file but with parsed metaboliteIDs
"""

from __future__ import division
from sys import argv
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__version__ = "1.0"
__date__ = "23 Nov 2017"

def parse_mQTL_file(mQTL_file):
    """A function to parse a file which IDs have to be changed into a list of list containing the elements per row.

    Keyword arguments:
        mQTL_file - as described above
    Returns:
        thelist - as described above
    """
    with open(mQTL_file, "r") as thefile:
        thelist = []
        next(thefile) #skip the header
        for line in thefile:
            elements = line.split()
            thelist.append(elements)
        return thelist

def parse_metaboliteID_file(metaboliteID_file):
    """A function to parse metaboliteID_file into a dictionary with metaboliteID as key and the rest as values.

    Keyword arguments:
        metaboliteID_file - as described above
    Returns:
        thedic - as described above
    """
    with open(metaboliteID_file, "r") as thefile:
        thedic = {}
        next(thefile) #skip the header
        for line in thefile:
            elements = line.split()
            print elements
            thedic[elements[0]] = elements[1:]
        return thedic

def write_file(output_name, id_dic, mQTL_list):
    """A function to write a new file with values from mQTL_list but with altered IDs according to id_dic.

    Keyword arguments:
        output_name - the name of the result file
        id_dic - a dictionary with one ID as the key and the matching characteristics as values
        mQTL_list - as described above
    Returns:
        parsed_mQTL.txt - a file, like mQTL_file, with changed IDs
    """
    with open(output_name, "w") as thefile:
        for mQTL in mQTL_list:
            if id_dic[mQTL[0]][0] != "NaN":
                id = id_dic[mQTL[0]][0]
                line = "{}  {}  {}  {}  {}  {}".format(id, mQTL[1], mQTL[2], mQTL[3], mQTL[4], mQTL[5])
            elif id_dic[mQTL[0]][0] == "NaN":
                id = "{}_{}_{}".format(mQTL[0], id_dic[mQTL[0]][1], id_dic[mQTL[0]][2])
                line = "{}  {}  {}  {}  {}  {}".format(id, mQTL[1], mQTL[2], mQTL[3], mQTL[4], mQTL[5])
            thefile.write(line + "\n")

if __name__ == "__main__":
    #Get the files from the command line
    mQTL_file = argv[1]
    metaboliteID_file = argv[2]
    output_name = "parsed_mQTLs.txt"

    #Parse the files
    mQTL_list = parse_mQTL_file(mQTL_file)
    id_dic = parse_metaboliteID_file(metaboliteID_file)

    #Change IDs and write to output_name
    write_file(output_name, id_dic, mQTL_list)
