#! /usr/bin/evn python
"""A Python script that overlaps xQTL with BGCs found with plantiSMASH.

python overlap_xQTL_BGC.py <BGC_dir> <eQTL_file> <mQTL_file> <gff3_file> <output_dir>

Keyword arguments:
    BGC_dir - A directory containing x_BGC.txt output files from plantiSMASH
    eQTL_file - A .tsv file containing the eQTLs (gene, chr, peak_mb, inf_mb, sup_mb, lod_score)
    mQTL_file - A .tsv file containing the mQTLs (metabolite, chr, peak_mb, inf_mb, sup_mb, lod_score)
    gff3_file - A .gff3 file containing the genome annotation of the species to be analyzed
    output_dir - The path to the directory where to write the results files to

Returns:
    cis_xQTL_BGC_XX.txt - A tab separated .txt file with cis-xQTL per BGC
    trans_xQTL_BGC_XX.txt -A tab separated .txt file with trans-xQTL per BGC
    eQTL_eQTL.txt - A tab separated .txt file with overlapping eQTL
    mQTL_mQTL.txt -A tab separated .txt file with overlapping mQTL
    mQTL_eQTL.txt A tab separated .txt file with overlapping mQTL and eQTL
"""

from __future__ import division
from sys import argv
import subprocess
import os.path
import re
import random

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "21 Nov 2017"
__version__ = "1.0"

#Parser functions
#################################################################################################################################################################
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
        gene_metabolite, chr, peak_mb, inf_mb, sup_mb, lod_score = line.split()
        elements = [gene_metabolite, int(chr), float(peak_mb), float(inf_mb), float(sup_mb), float(lod_score)]
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
            chr = chr[3:]
            type = elements[1]
            from_bp, to_bp = elements[3].split(";")
            genes = elements[4].split(";")
            genes = [re.sub(r"-.*", "", gene) for gene in genes]
            thedic[int(cluster_id)] = [type, int(chr), float(from_bp), float(to_bp), genes]
    return thedic

def gff3_parser_annotation(gff3_file):
    """A function to parse a GFF3 file and make a dictionary containing locusID and annotation.

    Keyword arguments:
        gff3_file - an annotation file in GFF3 format
    Returns:
        thedic - a dictionary containing locusID and accompanying annotation
    """
    thedic = {}
    with open(gff3_file, "r") as thefile:
        for line in thefile:
            if line.startswith("#"):
                continue
            else:
                elements = line.split("\t")
                if elements[2] == "gene" and (elements[0] != "ChrUn" and elements[0] != "ChrSy"):
                    from_bp, to_bp = elements[3], elements[4]
                    description = elements[-1].split(";")
                    locus = description[2].split("Alias=")[1]
                    locus = locus.strip()
                    annotation = description[1].strip().split("Name=")[1].split("%20")
                    annotation = [element.replace("%2C", ",") for element in annotation]
                    annotation = [element.replace("%2", "-") for element in annotation]
                    annotation = [element.replace("%2F", "/") for element in annotation]
                    annotation = " ".join(annotation)
                    thedic[locus] = [from_bp, to_bp, annotation]
        return thedic

#Find overlap functions
#################################################################################################################################################################
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
        BGC_region = [BGC_dic[key][1], BGC_dic[key][2], BGC_dic[key][3]]
        thedic[key] = BGC_region
        cis_xQTL_list = []
        for eQTL in eQTL_list:
            if eQTL[1] == BGC_region[0] and (eQTL[2]*1000000) > BGC_region[1] and (eQTL[2]*1000000) < BGC_region[2]: #test if peak is within BGC boundaries
                cis_xQTL_list.append(eQTL[0])
            elif eQTL[1] == BGC_region[0] and (eQTL[3]*1000000) > BGC_region[1] and (eQTL[3]*1000000) < (BGC_region[2]-(0.3*(BGC_region[2]-BGC_region[1]))): #test if inf is within BGC start and BGC end-30%
                cis_xQTL_list.append(eQTL[0])
            elif eQTL[1] == BGC_region[0] and (eQTL[4]*1000000) > (BGC_region[1]+(0.3*(BGC_region[2]-BGC_region[1]))) and (eQTL[4]*1000000) < BGC_region[2]: #test if sub is within BGC end and BGC start+30%
                cis_xQTL_list.append(eQTL[0])
            else:
                continue
        for mQTL in mQTL_list:
            if mQTL[1] == BGC_region[0] and (mQTL[2]*1000000) > BGC_region[1] and (mQTL[2]*1000000) < BGC_region[2]:
                cis_xQTL_list.append(mQTL[0])
            elif mQTL[1] == BGC_region[0] and (mQTL[3]*1000000) > BGC_region[1] and (mQTL[3]*1000000) < (BGC_region[2]-(0.3*(BGC_region[2]-BGC_region[1]))):
                cis_xQTL_list.append(mQTL[0])
            elif mQTL[1] == BGC_region[0] and (mQTL[4]*1000000) > (BGC_region[1]+(0.3*(BGC_region[2]-BGC_region[1]))) and (mQTL[4]*1000000) < BGC_region[2]:
                cis_xQTL_list.append(mQTL[0])
            else:
                continue
        thedic[key].append(cis_xQTL_list)
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
    common_dic = {}
    set_eQTLs = set([gene[0] for gene in eQTL_list])
    for key in BGC_dic:
        set_genes = set(BGC_dic[key][4])
        common = set_genes.intersection(set_eQTLs)
        if common == set([]):
            continue
        else:
            common_dic[key] = list(common)
    thedic = {}
    return thedic

def find_overlapping_xQTL(analysis, eQTL_list, mQTL_list):
    """A function that finds overlapping xQTLs based on their inf_mb and sup_mb. It returns them in a dictionary.
       It is considered overlap whenever one inf_mb/sup_mb falls within the region of the query.

    Keyword arguments:
        analysis - either eQTL vs eQTL, mQTL vs mQTL or mQTL vs eQTL
        eQTL_list - a list of lists containing the values from eQTL_file
        mQTL_list - a list of lists containing the values from mQTL_file
    Returns:
        thedic - a dictionary with gene/metabolite ID of xQTL as key and the gene/metabolite ID of overlapping
        xQTL.
    """
    if analysis == "eQTL_eQTL":
        thedic = {}
        for i in range(len(eQTL_list)):
            for j in range(len(eQTL_list)):
                if eQTL_list[i][0] not in thedic:
                    thedic[eQTL_list[i][0]] = []
                    if eQTL_list[i][1] == eQTL_list[j][1] and (eQTL_list[i][2]*1000000) > (eQTL_list[j][3]*1000000) and (eQTL_list[i][2]*1000000) < (eQTL_list[j][4]*1000000):
                        thedic[eQTL_list[i][0]].append(eQTL_list[j][0])
                    elif eQTL_list[i][1] == eQTL_list[j][1] and (eQTL_list[i][3]*1000000) > (eQTL_list[j][3]*1000000) and (eQTL_list[i][3]*1000000) < (eQTL_list[j][4]*1000000):
                        thedic[eQTL_list[i][0]].append(eQTL_list[j][0])
                    elif eQTL_list[i][1] == eQTL_list[j][1] and (eQTL_list[i][4]*1000000) > (eQTL_list[j][3]*1000000) and (eQTL_list[i][4]*1000000) < (eQTL_list[j][4]*1000000):
                        thedic[eQTL_list[i][0]].append(eQTL_list[j][0])
                else:
                    if eQTL_list[i][1] == eQTL_list[j][1] and (eQTL_list[i][2]*1000000) > (eQTL_list[j][3]*1000000) and (eQTL_list[i][2]*1000000) < (eQTL_list[j][4]*1000000):
                        thedic[eQTL_list[i][0]].append(eQTL_list[j][0])
                    elif eQTL_list[i][1] == eQTL_list[j][1] and (eQTL_list[i][3]*1000000) > (eQTL_list[j][3]*1000000) and (eQTL_list[i][3]*1000000) < (eQTL_list[j][4]*1000000):
                        thedic[eQTL_list[i][0]].append(eQTL_list[j][0])
                    elif eQTL_list[i][1] == eQTL_list[j][1] and (eQTL_list[i][4]*1000000) > (eQTL_list[j][3]*1000000) and (eQTL_list[i][4]*1000000) < (eQTL_list[j][4]*1000000):
                        thedic[eQTL_list[i][0]].append(eQTL_list[j][0])


    elif analysis == "mQTL_mQTL":
        thedic = {}
        for i in range(len(mQTL_list)):
            for j in range(len(mQTL_list)):
                if mQTL_list[i][0] not in thedic:
                    thedic[mQTL_list[i][0]] = []
                    if mQTL_list[i][1] == mQTL_list[j][1] and (mQTL_list[i][2]*1000000) > (mQTL_list[j][3]*1000000) and (mQTL_list[i][2]*1000000) < (mQTL_list[j][4]*1000000):
                        thedic[mQTL_list[i][0]].append(mQTL_list[j][0])
                    elif mQTL_list[i][1] == mQTL_list[j][1] and (mQTL_list[i][3]*1000000) > (mQTL_list[j][3]*1000000) and (mQTL_list[i][3]*1000000) < (mQTL_list[j][4]*1000000):
                        thedic[mQTL_list[i][0]].append(mQTL_list[j][0])
                    elif mQTL_list[i][1] == mQTL_list[j][1] and (mQTL_list[i][4]*1000000) > (mQTL_list[j][3]*1000000) and (mQTL_list[i][4]*1000000) < (mQTL_list[j][4]*1000000):
                        thedic[mQTL_list[i][0]].append(mQTL_list[j][0])
                else:
                    if mQTL_list[i][1] == mQTL_list[j][1] and (mQTL_list[i][2]*1000000) > (mQTL_list[j][3]*1000000) and (mQTL_list[i][2]*1000000) < (mQTL_list[j][4]*1000000):
                        thedic[mQTL_list[i][0]].append(mQTL_list[j][0])
                    elif mQTL_list[i][1] == mQTL_list[j][1] and (mQTL_list[i][3]*1000000) > (mQTL_list[j][3]*1000000) and (mQTL_list[i][3]*1000000) < (mQTL_list[j][4]*1000000):
                        thedic[mQTL_list[i][0]].append(mQTL_list[j][0])
                    elif mQTL_list[i][1] == mQTL_list[j][1] and (mQTL_list[i][4]*1000000) > (mQTL_list[j][3]*1000000) and (mQTL_list[i][4]*1000000) < (mQTL_list[j][4]*1000000):
                        thedic[mQTL_list[i][0]].append(mQTL_list[j][0])

    elif analysis == "mQTL_eQTL":
        thedic = {}
        for i in range(len(mQTL_list)):
            for j in range(len(eQTL_list)):
                if mQTL_list[i][0] not in thedic:
                    thedic[mQTL_list[i][0]] = []
                    if mQTL_list[i][1] == eQTL_list[j][1] and (mQTL_list[i][2]*1000000) > (eQTL_list[j][3]*1000000) and (mQTL_list[i][2]*1000000) < (eQTL_list[j][4]*1000000):
                        thedic[mQTL_list[i][0]].append(eQTL_list[j][0])
                    elif mQTL_list[i][1] == eQTL_list[j][1] and (mQTL_list[i][3]*1000000) > (eQTL_list[j][3]*1000000) and (mQTL_list[i][3]*1000000) < (eQTL_list[j][4]*1000000):
                        thedic[mQTL_list[i][0]].append(eQTL_list[j][0])
                    elif mQTL_list[i][1] == eQTL_list[j][1] and (mQTL_list[i][4]*1000000) > (eQTL_list[j][3]*1000000) and (mQTL_list[i][4]*1000000) < (eQTL_list[j][4]*1000000):
                        thedic[mQTL_list[i][0]].append(eQTL_list[j][0])
                else:
                    if mQTL_list[i][1] == eQTL_list[j][1] and (mQTL_list[i][2]*1000000) > (eQTL_list[j][3]*1000000) and (mQTL_list[i][2]*1000000) < (eQTL_list[j][4]*1000000):
                        thedic[mQTL_list[i][0]].append(eQTL_list[j][0])
                    elif mQTL_list[i][1] == eQTL_list[j][1] and (mQTL_list[i][3]*1000000) > (eQTL_list[j][3]*1000000) and (mQTL_list[i][3]*1000000) < (eQTL_list[j][4]*1000000):
                        thedic[mQTL_list[i][0]].append(eQTL_list[j][0])
                    elif mQTL_list[i][1] == eQTL_list[j][1] and (mQTL_list[i][4]*1000000) > (eQTL_list[j][3]*1000000) and (mQTL_list[i][4]*1000000) < (eQTL_list[j][4]*1000000):
                        thedic[mQTL_list[i][0]].append(eQTL_list[j][0])
    else:
        print "This analysis is invalid."
    return thedic

#Write output files functions
#################################################################################################################################################################
def write_file_overlapping_xQTLs(dic_overlap_xQTLs, output_dir, output_name):
    """A function to write the output of find_overlapping_xQTL() to files per dictionary.

    Keyword arguments:
        dic_overlap_xQTLs - dictionary with an xQTL as key and overlapping xQTLs as values
    Returns:
        overlapping_xQTL_xQTL.txt - a .txt tab-separated file containing the information from the dic_overlap_xQTLs
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_name, "w") as thefile:
        for key in dic_overlap_xQTLs:
                line = "\t".join(dic_overlap_xQTLs[key])
                thefile.write(line + "\n")

def write_file_cis_xQTLs(overlap_dic, output_dir, output_name, locus_annotation_dic, BGC_dic, overlap_count_dic):
    """A function to write the output of find_cis_xQTLs() to files per BGC showing overlap.

    Keyword arguments:
        overlap_dic - a dictionary with all overlapping xQTLs per BGC, clusterID is the key
        output_dir - the directory where the output files are made in
        output_name - the name of the output files
        locus_annotation_dic - a dictionary containing locusIDs and accompanying annotation
        BGC_dic - a dictionary with all BGCs found by plantiSMASH, clusterID is the key
    Returns:
        cis_xQTLs_BGC_XX.txt - a .tsv file containing overlapping xQTLs per cluster
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for key in overlap_dic:
        if overlap_dic[key][3] != []:
            filename = "{}_{}.txt".format(output_name, key)
            with open(filename, "w") as thefile:
                elements = [str(key)] + [BGC_dic[key][0]] + [str(overlap_dic[key][0])] + [str(overlap_dic[key][1])] + [str(overlap_dic[key][2])]
                line =  "\t".join(elements)
                thefile.write(line + "\n")
                for xQTL in overlap_dic[key][3]:
                    if xQTL.startswith("LOC") and xQTL in locus_annotation_dic:
                        if float(locus_annotation_dic[xQTL][0]) >= float(overlap_dic[key][1]) and float(locus_annotation_dic[xQTL][1]) <= float(overlap_dic[key][2]):
                            cluster_status = "in_cluster"
                            for QTL_count in overlap_count_dic[key]:
                                if QTL_count[0] == xQTL:
                                    p_value = str(QTL_count[1])
                            line_elements = [xQTL, cluster_status, p_value, locus_annotation_dic[xQTL][2], locus_annotation_dic[xQTL][0], locus_annotation_dic[xQTL][1]]
                            line = "\t".join(line_elements)
                            thefile.write(line + "\n")
                        else:
                            cluster_status = "not_in_cluster"
                            for QTL_count in overlap_count_dic[key]:
                                if QTL_count[0] == xQTL:
                                    p_value = str(QTL_count[1])
                            line_elements = [xQTL, cluster_status, p_value, locus_annotation_dic[xQTL][2], locus_annotation_dic[xQTL][0], locus_annotation_dic[xQTL][1]]
                            line = "\t".join(line_elements)
                            thefile.write(line + "\n")
                    else:
                        for QTL_count in overlap_count_dic[key]:
                            if QTL_count[0] == xQTL:
                                p_value = str(QTL_count[1])
                        line_elements = [xQTL, p_value]
                        line = "\t".join(line_elements)
                        thefile.write(line + "\n")

#Statistics functions
#################################################################################################################################################################
def statistics_xQTL(xQTL_list):
    """A function to return simple statistics of the distribution of the input xQTL data.

    Keyword arguments:
        xQTL_list -
    Returns:
    """
    statistics_xQTL_dic = {}
    for xQTL in xQTL_list:
        if xQTL[1] not in statistics_xQTL_dic:
            statistics_xQTL_dic[xQTL[1]] = [0, 0]
        else:
            statistics_xQTL_dic[xQTL[1]][0] += 1
            statistics_xQTL_dic[xQTL[1]][1] += (xQTL[4]-xQTL[3])
    for key in statistics_xQTL_dic:
        statistics_xQTL_dic[key][1] = (statistics_xQTL_dic[key][1]/statistics_xQTL_dic[key][0])
    return statistics_xQTL_dic

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
        if thedic[key][3] != []:
            count += 1
            total += 1
        else:
            total +=1
    print "{}/{} have overlap in this dictionary".format(count, total)

def shuffle_xQTL_data_chr(xQTL_list):
    """A function to shuffle the chromosomes of a xQTL dataset to achieve randomness.

    Keyword arguments:
        xQTL_list - a list of lists containing the rows of the original xQTL input file with gene/metabolite name
                    chromosome, inferior and superior interval location, peak location and LOD score.
    Returns:
        shuffled_xQTL_list - the same as the input list, however the chromosome data are shuffled assigning a
                             a random chromosome to each xQTL.
    """
    shuffled_xQTL_list = []
    for xQTL in xQTL_list:
        shuffled_xQTL = xQTL
        chr = xQTL[1]
        allowed_values = list(range(1, 12+1))
        allowed_values.remove(chr)
        shuffled_chr = random.choice(allowed_values)
        shuffled_xQTL[1] = shuffled_chr
        shuffled_xQTL_list.append(shuffled_xQTL)
    return shuffled_xQTL_list

def shuffle_BGC_data_chr(BGC_dic):
    """A function to shuffle the chromosomes of a BGC dataset to achieve randomness.

    Keyword arguments:
        BGC_dic - a dictionary containing information per BGC with clusterID as key.
    Returns:
        shuffled_BGC_dic - the same as the input dictionary, however the chromosome data are shuffled assigning a
                           a random chromosome to each xQTL.
    """
    shuffled_BGC_dic = {}
    for key in BGC_dic:
        shuffled_values = BGC_dic[key]
        chr = BGC_dic[key][1]
        allowed_values = list(range(1, 12+1))
        allowed_values.remove(chr)
        shuffled_chr = random.choice(allowed_values)
        shuffled_values[1] = shuffled_chr
        shuffled_BGC_dic[key] = shuffled_values
    return shuffled_BGC_dic

def randomization_cis_xQTL_BGC(BGC_dic, eQTL_list, mQTL_list, cis_xQTL_dic, permutations):
    """A function to do a randomization test for finding random cis-xQTL overlapping with BGCs.

    Keyword arguments:
        BGC_dic - a dictionary containing information per BGC with clusterID as key and randomly assigned chromosomes.
        eQTL_list - a list of lists containing the original rows of the eQTL input file but with randomly assigned chromosomes.
        mQTL_list - a list of lists containing the original rows of the mQTL input file but with randomly assigned chromosomes.
        cis_xQTL_dic - the real cis-xQTL overlap with BGCs that was found with the original datasets.
        permutations - the number of randomizations that should be performed.
    Returns:
        overlap_count_dic - a dictionary with clusterID as key and a list of lists as value containing overlapping xQTLs accompanied with their P-value, calculated as number of times overlap was found/number of permutations.
    """
    overlap_count_dic = {}
    for key in cis_xQTL_dic:
        overlap_count_dic[key] = []
        for QTL in cis_xQTL_dic[key][3]:
            overlap_count_dic[key].append([QTL, 0])
    for i in range(permutations):
        shuffled_BGC = shuffle_BGC_data_chr(BGC_dic)
        shuffled_eQTL = shuffle_xQTL_data_chr(eQTL_list)
        shuffled_mQTL = shuffle_xQTL_data_chr(mQTL_list)
        shuffled_cis_xQTL_dic = find_cis_xQTL(shuffled_BGC, shuffled_eQTL, shuffled_mQTL)
        for key in cis_xQTL_dic:
            real_overlap = set(cis_xQTL_dic[key][3])
            random_overlap = set(shuffled_cis_xQTL_dic[key][3])
            same_list = list(real_overlap & random_overlap)
            for QTL in same_list:
                for QTL_count in overlap_count_dic[key]:
                    if QTL == QTL_count[0]:
                        QTL_count[1] += 1
    for key in overlap_count_dic:
        for overlap in overlap_count_dic[key]:
            overlap[1] = overlap[1]/permutations
    return overlap_count_dic

if __name__ == "__main__":
    #Get files from command line + name results files
    BGC_dir = argv[1]
    eQTL_file = argv[2]
    mQTL_file = argv[3]
    gff3_file = argv[4]
    output_dir = argv[5]
    cis_xQTL_output_name = "{}/cis_xQTLs_BGC".format(output_dir)
    trans_xQTL_output_name = "{}/trans_xQTLs_BGC.txt".format(output_dir)
    eQTL_eQTL_output_name = "{}/eQTL_eQTL.txt".format(output_dir)
    mQTL_mQTL_output_name = "{}/mQTL_mQTL.txt".format(output_dir)
    mQTL_eQTL_output_name = "{}/mQTL_eQTL.txt".format(output_dir)

    #Parse the files
    BGC_dic = BGC_parser(BGC_dir)
    eQTL_list = xQTL_parser(eQTL_file)
    mQTL_list = xQTL_parser(mQTL_file)
    locus_annotation_dic = gff3_parser_annotation(gff3_file)

    #Find cis-xQTLs overlapping with BGC based on physical location and count how many BGCs have cis-xQTLs
    cis_xQTL_dic = find_cis_xQTL(BGC_dic, eQTL_list, mQTL_list)
    #count(cis_xQTL_dic)

    #Find overlapping trans-xQTLs based on genes present in BGC
    #trans_xQTL_dic = find_trans_xQTL(BGC_dic, eQTL_list, mQTL_list)
    #print trans_xQTL_dic

    #Find overlapping xQTLs based on their peak_mb, inf_mb and sup_mb
    #dic_eQTL_eQTL = find_overlapping_xQTL("eQTL_eQTL", eQTL_list, mQTL_list)
    #write_file_overlapping_xQTLs(dic_eQTL_eQTL, output_dir, eQTL_eQTL_output_name)

    #dic_mQTL_mQTL = find_overlapping_xQTL("mQTL_mQTL", eQTL_list, mQTL_list)
    #write_file_overlapping_xQTLs(dic_mQTL_mQTL, output_dir, mQTL_mQTL_output_name)

    #dic_mQTL_eQTL = find_overlapping_xQTL("mQTL_eQTL", eQTL_list, mQTL_list)
    #write_file_overlapping_xQTLs(dic_mQTL_eQTL, output_dir, mQTL_eQTL_output_name)

    #Calculate general statistics of xQTL and BGC datasets
    #print statistics_xQTL(eQTL_list)
    #print statistics_xQTL(mQTL_list)

    #Randomization test
    overlap_count_dic = randomization_cis_xQTL_BGC(BGC_dic, eQTL_list, mQTL_list, cis_xQTL_dic, 10000)
    write_file_cis_xQTLs(cis_xQTL_dic, output_dir, cis_xQTL_output_name, locus_annotation_dic, BGC_dic, overlap_count_dic)

