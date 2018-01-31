#! /usr/bin/evn python

"""A Python script that changes Affymetrix probe IDs into BGI gene IDs in a SOFT file.

python affy_bgi_id_parser.py <soft_file> <id_file>

Keyword arguments:
    soft_file -
    id_file -
Returns:
    parsed_soft.soft
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
    with open(id_file, "r") as thefile:
        dic = {}
        for line in thefile:
            elements = line.split()
            dic[elements[1]] = elements[0]
        return dic

def parser_soft_file(soft_file, id_dic):
    with open(soft_file, "r") as thefile:
        alist = []
        for line in thefile:
            if line.startswith("!") or line.startswith("^") or line.startswith("#"):
                alist.append(line)
            else:
                elements = line.split()
                if elements[0] in id_dic:
                    elements[0] = id_dic[elements[0]]
                    alist.append(elements)
                else:
                    alist.append(elements)
        return alist

def write_parsed_soft_file(soft_list, output_name):
    with open(output_name, "w") as thefile:
        for alist in soft_list:
            if type(alist) != list:
                thefile.write(alist)
            else:
                line = "\t".join(alist)
                thefile.write(line + "\n")

if __name__ == "__main__":
    soft_file = argv[1]
    id_file = argv[2]
    output_name = "parsed_soft.soft"

    id_dic = parser_id_file(id_file)
    soft_list = parser_soft_file(soft_file, id_dic)

    write_parsed_soft_file(soft_list, output_name)

