#! /usr/bin/evn python
"""A Python script that overlaps xQTL with BGCs found with plantiSMASH.

python overlap_xQTL_BGC.py <BGC_dir> <eQTL_file> <mQTL_file>

Keyword arguments:
    BGC_dir --> A directory containing x_BGC.txt output files from plantiSMASH
    eQTL_file --> A .tsv file containing the eQTLs (gene, chr, peak_mb, inf_mb, sup_mb, lod_score)
    mQTL_file --> A .tsv file containing the mQTLs (metabolite, chr, peak_mb, inf_mb, sup_mb, lod_score)

Returns:
    A file with xQTL per BGC
"""

from __future__ import division
from sys import argv
import subprocess
import os.path
import re

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "21 Nov 2017"
__version__ = "1.0"

def xQTL_parser(xQTL_file):
    """A function that parses a file with gene/metabolite, chr, peak_mb, inf_mb, sup_mb and lod_score in to a list
       of lists.

    Keyword arguments:
        xQTL_file - a .tsv file with gene/metabolite, chr, peak_mb, inf_mb, sup_mb and lod_score.
    Returns:
        thelist - a list of lists containing the values from xQTL_file.
    """
    thefile = open(xQTL_file, "r")
    thelist = []
    next(thefile) #skip the header
    for line in thefile:
        elements = line.split()
        thelist.append(elements)
    return thelist

def BGC_parser(BGC_dir):
    """A function that parses .tsv files (e.g. 1_BGC.txt) in a given directory with clusterID, type, chr, from in
       bp, to in bp, [genes in cluster] into a dictionary with clusterID as key and the rest as its value.

    Keyword arguments:
        BGC_dir - the path of the directory containing the BGCs files per chromosome (e.g. 1_BGC.txt).
    Returns:
        thedic - a dictionary with clusterID as key and the rest as values.
    """
    filelist = os.listdir(BGC_dir)
    filelist = [file for file in filelist if "BGC" in file]
    thedic = {}
    for file in filelist:
        file_path = "{}/{}".format(BGC_dir, file)
        thefile = open(file_path, "r")
        next(thefile) #skip the header
        for line in thefile:
            elements = line.split("\t")
            chr, cluster_id = elements[0].split("_c")
            type = elements[1]
            from_bp, to_bp = elements[3].split(";")
            genes = elements[4].split(";")
            genes = [re.sub(r"-.*", "", gene) for gene in genes]
            thedic[cluster_id] = [type, chr, from_bp, to_bp, genes]
    return thedic

def find_cis_xQTL(BGC_dic, eQTL_list, mQTL_list):
    """A function that finds cis-xQTLs based on the chromosomal region of the BGCs. A cis-xQTL is then defined as
       an eQTL or mQTL with their peaks within the chromosomal location of the BGC.

    Keyword arguments:
        BGC_dic - a dictionary with clusterID as key and type, chr, from, to, genes as values.
        eQTL_list - a list of lists containing the values from eQTL_file.
        mQTL_list - a list of lists containing the values from mQTL_file.
    Returns:
        thedic - a dictionary containing clusterID as keys and their overlapping xQTL (geneID and/or metaboliteID)
        as values.
    """
    thedic = {}
    for key in BGC_dic:
        thedic[key] = []
        BGC_region = [BGC_dic[key][1], BGC_dic[key][2], BGC_dic[key][3]]
        for eQTL in eQTL_list:
            if eQTL[1] == BGC_region[0] and ((eQTL[2]*1000000) >= BGC_region[1] and (eQTL[2]*1000000) <= BGC_region[2]):
                thedic[key].append(eQTL[0])
            elif eQTL[1] == BGC_region[0] and ((eQTL[3]*1000000) >= BGC_region[1] and (eQTL[3]*1000000) <= BGC_region[2]):
                thedic[key].append(eQTL[0])
            elif eQTL[1] == BGC_region[0] and ((eQTL[4]*1000000) >= BGC_region[1] and (eQTL[4]*1000000) <= BGC_region[2]):
                thedic[key].append(eQTL[0])
            else:
                continue
        for mQTL in mQTL_list:
            if mQTL[1] == BGC_region[0] and ((mQTL[2]*1000000) >= BGC_region[1] and (mQTL[2]*1000000) <= BGC_region[2]):
                thedic[key].append(mQTL[0])
            elif mQTL[1] == BGC_region[0] and ((mQTL[3]*1000000) >= BGC_region[1] and (mQTL[3]*1000000) <= BGC_region[2]):
                thedic[key].append(mQTL[0])
            elif mQTL[1] == BGC_region[0] and ((mQTL[4]*1000000) >= BGC_region[1] and (mQTL[4]*1000000) <= BGC_region[2]):
                thedic[key].append(mQTL[0])
            else:
                continue
    return thedic

def find_trans_xQTL(BGC_dic, eQTL_list, mQTL_list):
    """A function that finds overlapping trans-xQTLs based on the genes in the BGCs. A trans-xQTL is then defined
       as overlapping eQTL and/or mQTL with their peaks outside the chromosomal location of the BGC.

    Keyword arguments:
        BGC_dic - a dictionary with clusterID as key and type, chr, from, to, genes as values.
        eQTL_list - a list of lists containing the values from eQTL_file.
        mQTL_list - a list of lists containing the values from mQTL_file.
    Returns:
        thedic - a dictionary containing clusterID as keys and a list of lists containing overlapping xQTLs per
        gene in the BGC.
    """
    thedic = {}
    return thedic

def find_overlapping_xQTL(eQTL_list, mQTL_list):
    """A function that finds overlapping xQTLs based on their inf_mb and sup_mb. It returns them in a dictionary.
       It is considered overlap whenever one inf_mb/sup_mb falls within the region of the query.

    Keyword arguments:
        eQTL_list - a list of lists containing the values from eQTL_file
        mQTL_list - a list of lists containing the values from mQTL_file
    Returns:
        thedic - a dictionary with gene/metabolite ID of xQTL as key and the gene/metabolite ID of overlapping
        xQTL.
    """
    dic_eQTL_eQTL = {}
    for i in range(len(eQTL_list)):
        for j in range(len(eQTL_list)):
            dic_eQTL_eQTL[eQTL_list[i][0]] = []
            if eQTL_list[i][1] == eQTL_list[j][1] and ((eQTL_list[i][2]*1000000) >= (eQTL_list[j][3]*1000000) and (eQTL_list[i][2]*1000000) <= (eQTL_list[j][4])):
                dic_eQTL_eQTL[eQTL_list[i][0]].append(eQTL_list[j][0])
            elif eQTL_list[i][1] == eQTL_list[j][1] and ((eQTL_list[i][3]*1000000) >= (eQTL_list[j][3]*1000000) and (eQTL_list[i][3]*1000000) <= (eQTL_list[j][4])):
                dic_eQTL_eQTL[eQTL_list[i][0]].append(eQTL_list[j][0])
            elif eQTL_list[i][1] == eQTL_list[j][1] and ((eQTL_list[i][4]*1000000) >= (eQTL_list[j][3]*1000000) and (eQTL_list[i][4]*1000000) <= (eQTL_list[j][4])):
                dic_eQTL_eQTL[eQTL_list[i][0]].append(eQTL_list[j][0])
            else:
                continue

    dic_mQTL_mQTL = {}
    for i in range(len(mQTL_list)):
        for j in range(len(mQTL_list)):
            dic_mQTL_mQTL[mQTL_list[i][0]] = []
            if mQTL_list[i][1] == mQTL_list[j][1] and ((mQTL_list[i][2]*1000000) >= (mQTL_list[j][3]*1000000) and (mQTL_list[i][2]*1000000) <= (mQTL_list[j][4])):
                dic_mQTL_mQTL[mQTL_list[i][0]].append(mQTL_list[j][0])
            elif mQTL_list[i][1] == mQTL_list[j][1] and ((mQTL_list[i][3]*1000000) >= (mQTL_list[j][3]*1000000) and (mQTL_list[i][3]*1000000) <= (mQTL_list[j][4])):
                dic_mQTL_mQTL[mQTL_list[i][0]].append(mQTL_list[j][0])
            elif mQTL_list[i][1] == mQTL_list[j][1] and ((mQTL_list[i][4]*1000000) >= (mQTL_list[j][3]*1000000) and (mQTL_list[i][4]*1000000) <= (mQTL_list[j][4])):
                dic_mQTL_mQTL[mQTL_list[i][0]].append(mQTL_list[j][0])
            else:
                continue

    dic_mQTL_eQTL = {}
    for i in range(len(mQTL_list)):
        for j in range(len(eQTL_list)):
            dic_mQTL_eQTL[mQTL_list[i][0]] = []
            if mQTL_list[i][1] == eQTL_list[j][1] and ((mQTL_list[i][2]*1000000) >= (eQTL_list[j][3]*1000000) and (mQTL_list[i][2]*1000000) <= (eQTL_list[j][4])):
                dic_mQTL_eQTL[mQTL_list[i][0]].append(eQTL_list[j][0])
            elif mQTL_list[i][1] == eQTL_list[j][1] and ((mQTL_list[i][3]*1000000) >= (eQTL_list[j][3]*1000000) and (mQTL_list[i][3]*1000000) <= (eQTL_list[j][4])):
                dic_mQTL_eQTL[mQTL_list[i][0]].append(eQTL_list[j][0])
            elif mQTL_list[i][1] == eQTL_list[j][1] and ((mQTL_list[i][4]*1000000) >= (eQTL_list[j][3]*1000000) and (mQTL_list[i][4]*1000000) <= (eQTL_list[j][4])):
                dic_mQTL_eQTL[mQTL_list[i][0]].append(eQTL_list[j][0])
            else:
                continue

    return dic_eQTL_eQTL, dic_mQTL_mQTL, dic_mQTL_eQTL

def count(thedic):
    """A function to count the number of keys in a dictionary with values.

    Keyword arguments:
        thedic - a dictionary
    Returns:
        a string saying e.g. 5/50 have overlap in the following dictionary: thedic
    """
    count = 0
    total = 0
    for key in thedic:
        if thedic[key] != []:
            count += 1
            total += 1
        else:
            total +=1
    print "{}/{} have overlap in the following dictionary: {}".format(count, total, thedic)

def write_file(cis_xQTL_dic, cis_xQTL_output_name):
    with open(cis_xQTL_output_name, "w") as thefile:
        for key in cis_xQTL_dic:
            elements = [key] + cis_xQTL_dic[key]
            line =  "\t".join(elements)
            thefile.write(line + "\n")

if __name__ == "__main__":
    #Get files from command line
    BGC_dir = argv[1]
    eQTL_file = argv[2]
    mQTL_file = argv[3]
    cis_xQTL_output_name = "cis_xQTLs.txt"

    #Parse the files
    BGC_dic = BGC_parser(BGC_dir)
    eQTL_list = xQTL_parser(eQTL_file)
    mQTL_list = xQTL_parser(mQTL_file)

    #Find cis-xQTLs overlapping with BGC based on physical location and count how many BGCs have cis-xQTLs
    cis_xQTL_dic = find_cis_xQTL(BGC_dic, eQTL_list, mQTL_list)
    count(cis_xQTL_dic)
    write_file(cis_xQTL_dic, cis_xQTL_output_name)

    #Find overlapping trans-xQTLs based on genes present in BGC
    #trans_xQTL_dic = find_trans_xQTL(BGC_dic, eQTL_list, mQTL_list)

    #Find overlapping xQTLs based on their peak_mb, inf_mb and sup_mb
    #dic_eQTL_eQTL, dic_mQTL_mQTL, dic_mQTL_eQTL = find_overlapping_xQTL(eQTL_list, mQTL_list)
    #count(dic_eQTL_eQTL)
    #count(dic_mQTL_mQTL)
    #count(dic_mQTL_eQTL)
