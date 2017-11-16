#!/usr/bin/evn/python

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
__email__ = "lottewites@outlook.com"

def replace_id(id_to_be_replaced):



def parser(file_to_be_parsed):
    thefile = open(file_to_be_parsed, "r")
    dic = {}
    for line in thefile:
        elements = line.split()
        if elements[0] in dic:
            dic[elements[0]].append(elements[1:])
        else:
            dic[elements[0]] = elements[1:]
    return dic

def write_changed_file(list_to_be_changed, id1_id2_list, output_name):
    thefile = open(output_name, "w")
    for key in list_to_be_changed:
        if key in id1_id2_list and id1_id2_list[key] != [""]:
            if len(id1_id2_list[key]) != 1:
                id = "-".join(id1_id2_list[key] 
                line = "{}    {}    {}    {}    {}    {}".format(id,  list_to_be_changed[key][0], list_to_be_changed[key][1], list_to_be_changed[key][2]


if __name__ == "__main__":
    file_to_be_changed = argv[1]
    id1_id2_file = argv[2]
    output_name = "{}_parsed".format(file_to_be_changed)

    list_to_be_changed = parser(file_to_be_changed)
    id1_id2_list = parser(id1_id2_file)

    #write_changed_file = write_change_file(list_to_be_changed, id1_id2_list, output_name)
