#! /usr/bin/evn python

"""A Python script that changes Affymetrix probe IDs into BGI gene IDs in a matrix.

python affy_bgi_id_parser.py <martrix_file> <id_file>

Keyword arguments:
    matrix_file - a gene expression file in matrix
    id_file - a tab-separated file with IDs
Returns:
    parsed_soft.soft - a gene expression file in SOFT format with parsed IDs
"""

from __future__ import division
from sys import argv
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__version__ = "1.0"
__date__ = "31 Jan 2018"


def parser_id_file(id_file):
    """A function to parse a file with IDs into a dictionary.

    Keyword arguments:
        id_file - a tab-separated file with IDs
    Returns:
        dic - a dictionary containing the one ID as key and the other as value
    """
    with open(id_file, "r") as thefile:
        dic = {}
        for line in thefile:
            elements = line.split()
            dic[elements[1]] = elements[0]
        return dic

def parser_soft_file(soft_file, id_dic):
    """A function to parse a SOFT file into a list and parse IDs.

    Keyword arguments:
        soft_file - a gene expression file in SOFT format
        id_dic - a dictionary containing the one ID as key and the other as value
    Returns:
        alist - a list of lists containing the lines of the SOFT file, IDs have been parsed
    """
    with open(soft_file, "r") as thefile:
        alist = []
        for line in thefile:
            #print line
            if line.startswith("!") or line.startswith("^") or line.startswith("#"):
                continue
            else:
                elements = line.split()
                elements[0] = elements[0].strip('"')
                if elements[0] in id_dic:
                    elements[0] = id_dic[elements[0]]
                    #if len(elements) > 4:
                        #elements[1] = elements[0]
                    alist.append(elements)
                #elif elements[0] not in id_dic and len(elements) > 4:
                    #elements[1] = elements[0]
                    #alist.append(elements)
                elif elements[0] == "ID_REF":
                    alist.append(elements)
                else:
                    continue
                    #alist.append(elements)
        return alist

def write_parsed_soft_file(soft_list, output_name):
    """A function to write the parsed items in soft_list to output_name.

    Keyword arguments:
        soft_list - a list of lists containing the lines of the SOFT file, IDs have been parsed
        output_name - the name of the output SOFT file
    Returns:
        output_name - a filled output SOFT file
    """
    with open(output_name, "w") as thefile:
        for alist in soft_list:
            if type(alist) != list:
                thefile.write(alist)
            else:
                line = "\t".join(alist)
                thefile.write(line + "\n")

if __name__ == "__main__":
    #Get files from command line
    soft_file = argv[1]
    id_file = argv[2]
    output_name = "parsed_soft.soft"

    #Parse IDs and SOFT file
    id_dic = parser_id_file(id_file)
    soft_list = parser_soft_file(soft_file, id_dic)

    #Write parsed SOFT file to output_name
    write_parsed_soft_file(soft_list, output_name)

