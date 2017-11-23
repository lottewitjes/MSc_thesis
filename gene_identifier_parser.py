#! /usr/bin/evn python

"""A Python scripts that changes identifiers in one file based on another file.

python identifier_parser.py <file_to_be_changed> <id1_id2_file>

Keyword arguments:
- file_to_be_changed
- id1_id2_file
Returns:
- changed_file
"""

from __future__ import division
from sys import argv
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__version__ = "1.0"
__date__ = "17 Nov 2017"


def parser_ids(id1_id2_file):
    """A function to parse a file with IDs into a dictionary with one ID as the key and the other as the value.

    Keyword arguments:
        id1_id2_file - a file containing 2 columns with matching IDs
    Returns:
        dic - a dictionary with one ID as the key and the matching ID as value
    """
    with open(id1_id2_file, "r") as thefile:
        dic = {}
        for line in thefile:
            elements = line.split()
            if len(elements) != 1:
                dic[elements[1]] = elements[0]
        return dic

def parser_file_to_be_changed(file_to_be_changed):
    """A function to parse a file which IDs have to be changed into a list of list containing the elements per  row.

    Keyword arguments:
        file_to_be_changed - a file which IDs have to be changed
    Returns:
        thelist - a list of lists containing the values in the file_to_be_changed as a list per row.
    """
    with open(file_to_be_changed, "r") as thefile:
        thelist = []
        for line in thefile:
            elements = line.split()
            thelist.append(elements)
        return thelist

def write_changed_file(list_to_be_changed, id1_id2_dic, output_name):
    """A function to write a new file with values from list_to_be_changed but with altered IDs according to id1_id2_dic.

    Keyword arguments:
        list_to_be_changed - a list of lists containing the values in the file_to_be_changed as a list per row.
        id1_id2_dic - a dictionary with one ID as the key and the matching ID as value.
        output_name - the name of the result file.
    Returns:
        output_name - a file, like file_to_be_changed, with changed IDs.
    """
    with open(output_name, "w") as thefile:
        for alist in list_to_be_changed:
            if alist[0] in id1_id2_dic:
                line = "{}  {}  {}  {}  {}  {}".format(id1_id2_dic[alist[0]], alist[1], alist[2], alist[3], alist[4], alist[5])
            else:
                line = "{}  {}  {}  {}  {}  {}".format(alist[0], alist[1], alist[2], alist[3], alist[4], alist[5])
            thefile.write(line + "\n")

if __name__ == "__main__":
    #Get files from command line and format the name of the output file
    file_to_be_changed = argv[1]
    id1_id2_file = argv[2]
    output_name = "{}_parsed".format(file_to_be_changed)

    #Parse the files into dictionary and list of lists
    id1_id2_dic = parser_ids(id1_id2_file)
    list_to_be_changed = parser_file_to_be_changed(file_to_be_changed)

    #Change the IDs of the file_to_be_changed and write the new file
    write_changed_file(list_to_be_changed, id1_id2_dic, output_name)
