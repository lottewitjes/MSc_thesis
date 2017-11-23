#! /usr/bin/evn python

"""A Python scripts that changes gene identifiers in an eQTL file based on another file.

python identifier_parser.py <eQTL_file> <geneID_file>

Keyword arguments:
- eQTL_file
- geneID_file
Returns:
- parsed_eQTLs.txt
"""

from __future__ import division
from sys import argv
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__version__ = "1.0"
__date__ = "17 Nov 2017"


def parser_geneID_file(geneID_file):
    """A function to parse a file with IDs into a dictionary with one ID as the key and the other as the value.

    Keyword arguments:
        geneID_file - a file containing 2 columns with matching IDs
    Returns:
        dic - a dictionary with one ID as the key and the matching ID as value
    """
    with open(geneID_file, "r") as thefile:
        dic = {}
        for line in thefile:
            elements = line.split()
            if len(elements) != 1:
                dic[elements[1]] = elements[0]
        return dic

def parser_eQTL_file(eQTL_file):
    """A function to parse a file which IDs have to be changed into a list of list containing the elements per  row.

    Keyword arguments:
        eQTL_file - a file which IDs have to be changed
    Returns:
        thelist - a list of lists containing the values in the eQTL_file as a list per row.
    """
    with open(eQTL_file, "r") as thefile:
        thelist = []
        for line in thefile:
            elements = line.split()
            thelist.append(elements)
        return thelist

def write_file(output_name, id_dic, eQTL_list):
    """A function to write a new file with values from list_to_be_changed but with altered IDs according to id1_id2_dic.

    Keyword arguments:
        eQTL_list- a list of lists containing the values in the file_to_be_changed as a list per row.
        id_dic - a dictionary with one ID as the key and the matching ID as value.
        output_name - the name of the result file.
    Returns:
        parsed_eQTLs.txt - a file, like eQTL_file, with changed IDs.
    """
    with open(output_name, "w") as thefile:
        for eQTL in eQTL_list:
            if eQTL[0] in id_dic:
                line = "{}  {}  {}  {}  {}  {}".format(id_dic[eQTL[0]], eQTL[1], eQTL[2], eQTL[3], eQTL[4], eQTL[5])
            else:
                line = "{}  {}  {}  {}  {}  {}".format(eQTL[0], eQTL[1], eQTL[2], eQTL[3], eQTL[4], eQTL[5])
            thefile.write(line + "\n")

if __name__ == "__main__":
    #Get files from command line and format the name of the output file
    eQTL_file= argv[1]
    geneID_file = argv[2]
    output_name = "parsed_eQTLs.txt"

    #Parse the files into dictionary and list of lists
    id_dic = parser_geneID_file(geneID_file)
    eQTL_list = parser_eQTL_file(eQTL_file)

    #Change the IDs of the eQTL_file and write the new file
    write_file(output_name, id_dic, eQTL_list)
