#! /usr/bin/evn python
"""A Python script that rewrites a .GFF3 file into the desired .GFF3 format.

python rewrite_gff3.py <original_gff3_file> <result_gff3_file>

Keyword arguments:
    original_gff3_file -
    result_gff3_file -
Returns:
    result_gff3_file - 
"""

from __future__ import division
from sys import argv
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__version__ = "1.0"
__date__ = "7 Dec 2017"

def parse_gff3(original_gff3_file):
    thelist = []
    thedic = {}
    with open(original_gff3_file, "r") as thefile:
        for line in thefile:
            if line.startswith("#"):
                continue
            else:
                elements = line.split("\t")
                if elements[0] == "ChrUn" or elements[0] == "ChrSy":
                    continue
                elif elements[2] == "gene" and elements[0] != "ChrUn" and elements[0] != "ChrSy":
                    description = elements[-1].split(";")
                    locus = description[2].split("Alias=")[1].strip()
                    alias = description[0].split("ID=")[1]
                    annotation = description[1].split("Name=")[1]
                    del elements[-1]
                    elements += [locus, alias, annotation]
                    thelist.append(elements)
                elif elements[2] == "mRNA" and elements[0] != "ChrUn" and elements[0] != "ChrSy":
                    description = elements[-1].split(";")
                    ID = description[2].split("Alias=")[1].strip()
                    parent = description[2].split("Alias=")[1].split(".")[0]
                    thedic[description[0].split("ID=")[1]] = ID
                    del elements[-1]
                    elements += [ID, parent]
                    thelist.append(elements)
                else:
                    description = elements[-1].strip()
                    del elements[-1]
                    elements += [description]
                    thelist.append(elements)
    return thelist, thedic

def write_gff3(original_gff3_list, result_gff3_file, locus_alias_dic):
    with open(result_gff3_file, "w") as thefile:
        for alist in original_gff3_list:
            if alist[2] == "gene":
                ID = str("ID=" + alist[-3])
                name = str("Name=" + alist[-1])
                description = ";".join([ID,name])
                alist = alist[0:-3] + [description]
                line = "\t".join(alist)
                thefile.write(line + "\n")
            elif alist[2] == "mRNA":
                ID = str("ID=" + alist[-2])
                parent = str("Parent=" + alist[-1])
                description = ";".join([ID,parent])
                alist = alist[0:-2] + [description]
                line = "\t".join(alist)
                thefile.write(line + "\n")
            else:
                parent = alist[-1].split("Parent=")[1]
                parent = str("Parent=" + locus_alias_dic[parent])
                alist = alist[0:-1] + [parent]
                line = "\t".join(alist)
                thefile.write(line + "\n")

if __name__ == "__main__":
    #Get files from command line
    original_gff3_file = argv[1]
    result_gff3_file = argv[2]

    #Parse the original gff3 file into a list of lists
    original_gff3_list, locus_alias_dic = parse_gff3(original_gff3_file)
    #print original_gff3_list
    print locus_alias_dic

    #Write the original list of lists to the desired gff3 format
    write_gff3(original_gff3_list, result_gff3_file, locus_alias_dic)

