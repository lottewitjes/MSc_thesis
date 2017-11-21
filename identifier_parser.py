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
__version__ = "1.0"
__email__ = "lottewitjes@outlook.com"


def parser_ids(id1_id2_file):
    thefile = open(id1_id2_file, "r")
    dic = {}
    for line in thefile:
        elements = line.split()
        if len(elements) != 1:
            dic[elements[1]] = elements[0]
    return dic

def parser_file_to_be_changed(file_to_be_changed):
    thefile = open(file_to_be_changed, "r")
    thelist = []
    for line in thefile:
        elements = line.split()
        thelist.append(elements)
    return thelist

def write_changed_file(list_to_be_changed, id1_id2_dic, output_name):
    thefile = open(output_name, "w")
    for alist in list_to_be_changed:
        if alist[0] in id1_id2_dic:
            line = "{}  {}  {}  {}  {}  {}".format(id1_id2_dic[alist[0]], alist[1], alist[2], alist[3], alist[4], alist[5])
        else:
            line = "{}  {}  {}  {}  {}  {}".format(alist[0], alist[1], alist[2], alist[3], alist[4], alist[5])
        thefile.write(line + "\n")
    thefile.close()

if __name__ == "__main__":
    file_to_be_changed = argv[1]
    id1_id2_file = argv[2]
    output_name = "{}_parsed".format(file_to_be_changed)

    id1_id2_dic = parser_ids(id1_id2_file)
    list_to_be_changed = parser_file_to_be_changed(file_to_be_changed)

    write_changed_file(list_to_be_changed, id1_id2_dic, output_name)
