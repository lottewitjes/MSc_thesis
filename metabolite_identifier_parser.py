#! /usr/bin/evn python
"""A Python script to
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
    with open(mQTL_file, "r") as thefile:
        thelist = []
        next(thefile) #skip the header
        for line in thefile:
            elements = line.split("\t")
            thelist.append(elements)
        return thelist

def parse_metaboliteID_file(metaboliteID_file):
    with open(metaboliteID_file, "r") as thefile:
        thedic = {}
        next(thefile) #skip the header
        for line in thefile:
            elements = line.split("\t")
            print elements
            thedic[elements[0]] = elements[1:]
        return thedic

def write_file(output_name, id_dic, mQTL_list):
    with open(output_name, "w") as thefile:
        for mQTL in mQTL_list:
            if id_dic[mQTL[0]][1] != "NaN":
                id = id_dic[mQTL[0]][1]
                line = "{}  {}  {}  {}  {}  {}".format(id, mQTL[1], mQTL[2], mQTL[3], mQTL[4], mQTL[5])
            elif id_dic[mQTL[0]][1] == "NaN":
                id = "{}_{}_{}".format(mQTL[0], id_dic[mQTL[0]][2], id_dic[mQTL[0]][3])
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
    print id_dic
    #Change IDs and write to output_name
    write_file(output_name, id_dic, mQTL_list)
