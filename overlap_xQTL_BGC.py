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
import copy
import math

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
            chr = chr[3:] #os
            #chr = chr[-1] #at
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
                if (elements[2] == "gene" or elements[2] == "pseudogene") and (elements[0] != "ChrUn" and elements[0] != "ChrSy"):
                    from_bp, to_bp = elements[3], elements[4]
                    description = elements[-1].split(";")
                    locus = description[2].split("Alias=")[1] #os
                    #locus = description[1].split(",")[0] #at
                    #locus = re.findall("Dbxref=(.+):(.+)", locus)[0][1] #at
                    locus = locus.strip()
                    #annotation = re.findall("Name=(.+)", description[2])[0] #at
                    annotation = description[1].strip().split("Name=")[1].split("%20") #os
                    annotation = [element.replace("%2C", ",") for element in annotation] #os
                    annotation = [element.replace("%2", "-") for element in annotation] #os
                    annotation = [element.replace("%2F", "/") for element in annotation] #os
                    annotation = " ".join(annotation) #os
                    thedic[locus] = [from_bp, to_bp, annotation]
        return thedic

def gene_ID_refseq_protein_ID_parser(BGC_dic): #at
    id_dic = {}
    with open("/home/witje010/arabidopsis_thaliana_xqtl/id_parser_table.txt", "r") as thefile:
        next(thefile) #skip the header
        for line in thefile:
            gene_ID, refseq_protein_ID = line.strip().split("\t")
            id_dic[refseq_protein_ID] = gene_ID
        for BGC in BGC_dic:
            BGC_dic[BGC][4] = [id_dic[gene] for gene in BGC_dic[BGC][4] if gene in id_dic]
        return BGC_dic

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
            if eQTL[1] == BGC_region[0] and (BGC_region[1] < (eQTL[2]*1000000) < BGC_region[2]): #test if peak is within BGC boundaries
                cis_xQTL_list.append((eQTL[0], eQTL[5]))
            elif eQTL[1] == BGC_region[0] and (BGC_region[1] < (eQTL[3]*1000000) < (BGC_region[2]-(0.0*(BGC_region[2]-BGC_region[1])))): #test if inf is within BGC start and BGC end-30% for os, 0% for at
                cis_xQTL_list.append((eQTL[0], eQTL[5]))
            elif eQTL[1] == BGC_region[0] and ((BGC_region[1]+(0.0*(BGC_region[2]-BGC_region[1]))) < (eQTL[4]*1000000) < BGC_region[2]): #test if sub is within BGC end and BGC start+30% for os, 0% for at
                cis_xQTL_list.append((eQTL[0], eQTL[5]))
            else:
                continue
        for mQTL in mQTL_list:
            if mQTL[1] == BGC_region[0] and (BGC_region[1] < (mQTL[2]*1000000) < BGC_region[2]): #test if peak is within BGC boundaries
                cis_xQTL_list.append((mQTL[0], mQTL[5]))
            elif mQTL[1] == BGC_region[0] and (BGC_region[1] < (mQTL[3]*1000000) < (BGC_region[2]-(0.0*(BGC_region[2]-BGC_region[1])))): #test if inf is within BGC start and BGC end-30% for os, 0% for at
                cis_xQTL_list.append((mQTL[0], mQTL[5]))
            elif mQTL[1] == BGC_region[0] and ((BGC_region[1]+(0.0*(BGC_region[2]-BGC_region[1]))) < (mQTL[4]*1000000) < BGC_region[2]): #test if sub is within BGC end and BGC start+30% for os, 0% for at
                cis_xQTL_list.append((mQTL[0], mQTL[5]))
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
        overlap_list = []
        for i in range(len(eQTL_list)):
            for j in range(len(eQTL_list)):
                    elements = [eQTL_list[i][0], eQTL_list[i][5], eQTL_list[j][0], eQTL_list[j][5]]
                    elements_rev = [eQTL_list[j][0], eQTL_list[j][5], eQTL_list[i][0], eQTL_list[i][5]]
                    if eQTL_list[i] == eQTL_list[j]:
                        continue
                    elif elements in overlap_list or elements_rev in overlap_list:
                        continue
                    elif eQTL_list[i][1] == eQTL_list[j][1] and ((eQTL_list[j][3]*1000000) < (eQTL_list[i][2]*1000000) < (eQTL_list[j][4]*1000000)):
                        overlap_list.append(elements)
                    elif eQTL_list[i][1] == eQTL_list[j][1] and ((eQTL_list[j][3]*1000000) < (eQTL_list[i][3]*1000000) < ((eQTL_list[j][4]*1000000)-(0.0*((eQTL_list[j][4]*1000000)-(eQTL_list[j][3]*1000000))))):
                        overlap_list.append(elements)
                    elif eQTL_list[i][1] == eQTL_list[j][1] and (((eQTL_list[j][4]*1000000)-(0.0*((eQTL_list[j][4]*1000000)-(eQTL_list[j][3]*1000000)))) < (eQTL_list[i][4]*1000000) < (eQTL_list[j][4]*1000000)):
                        overlap_list.append(elements)

    elif analysis == "mQTL_mQTL":
        overlap_list = []
        for i in range(len(mQTL_list)):
            for j in range(len(mQTL_list)):
                    elements = [mQTL_list[i][0], mQTL_list[i][5], mQTL_list[j][0], mQTL_list[j][5]]
                    elements_rev = [mQTL_list[j][0], mQTL_list[j][5], mQTL_list[i][0], mQTL_list[i][5]]
                    if mQTL_list[i] == mQTL_list[j]:
                        continue
                    elif elements in overlap_list or elements_rev in overlap_list:
                        continue
                    elif mQTL_list[i][1] == mQTL_list[j][1] and ((mQTL_list[j][3]*1000000) < (mQTL_list[i][2]*1000000) < (mQTL_list[j][4]*1000000)):
                        overlap_list.append(elements)
                    elif mQTL_list[i][1] == mQTL_list[j][1] and ((mQTL_list[j][3]*1000000) < (mQTL_list[i][3]*1000000) < ((mQTL_list[j][4]*1000000)-(0.0*((mQTL_list[j][4]*1000000)-(mQTL_list[j][3]*1000000))))):
                        overlap_list.append(elements)
                    elif mQTL_list[i][1] == mQTL_list[j][1] and (((mQTL_list[j][4]*1000000)+(0.0*((mQTL_list[j][4]*1000000)-(mQTL_list[j][3]*1000000)))) < (mQTL_list[i][4]*1000000) < (mQTL_list[j][4]*1000000)):
                        overlap_list.append(elements)

    elif analysis == "mQTL_eQTL":
        overlap_list = []
        for i in range(len(mQTL_list)):
            for j in range(len(eQTL_list)):
                    elements = [mQTL_list[i][0], mQTL_list[i][5], eQTL_list[j][0], eQTL_list[j][5]]
                    elements_rev = [eQTL_list[j][0], eQTL_list[j][5], mQTL_list[i][0], mQTL_list[i][5]]
                    if mQTL_list[i] == mQTL_list[j]:
                        continue
                    elif elements in overlap_list or elements_rev in overlap_list:
                        continue
                    elif mQTL_list[i][1] == eQTL_list[j][1] and ((eQTL_list[j][3]*1000000) < (mQTL_list[i][2]*1000000) < (eQTL_list[j][4]*1000000)): #test if peak is within other QTL region
                        overlap_list.append(elements)
                    elif mQTL_list[i][1] == eQTL_list[j][1] and ((eQTL_list[j][3]*1000000) < (mQTL_list[i][3]*1000000) < ((eQTL_list[j][4]*1000000)-(0.0*((eQTL_list[j][4]*1000000)-(eQTL_list[j][3]*1000000))))):
                        overlap_list.append(elements)
                    elif mQTL_list[i][1] == eQTL_list[j][1] and (((eQTL_list[j][4]*1000000)+(0.0*((eQTL_list[j][4]*1000000)-(eQTL_list[j][3]*1000000)))) < (mQTL_list[i][4]*1000000) < (eQTL_list[j][4]*1000000)):
                        overlap_list.append(elements)

    else:
        print "This analysis is invalid."
    return overlap_list

#Statistics functions
#################################################################################################################################################################
def statistics_xQTL(xQTL_list):
    """A function to return simple statistics of the distribution of the input xQTL data.

    Keyword arguments:
        xQTL_list - a list of lists containing the rows of the xQTL input file.
    Returns:
        statistics_xQTL_dic - a dictionary with chromosome number as key and number of xQTL and average xQTL size as values.
    """
    statistics_xQTL_dic = {}
    for xQTL in xQTL_list:
        if xQTL[1] not in statistics_xQTL_dic:
            statistics_xQTL_dic[xQTL[1]] = [0, 0]
        else:
            statistics_xQTL_dic[xQTL[1]][0] += 1
            statistics_xQTL_dic[xQTL[1]][1] += (xQTL[4]-xQTL[3])
    for key in statistics_xQTL_dic:
        statistics_xQTL_dic[key].append(statistics_xQTL_dic[key][1]/statistics_xQTL_dic[key][0])
    elements = ["chromosome", "#xQTLs", "total xQTL size", "average xQTL size"]
    line = "\t".join(elements)
    print line
    for key in statistics_xQTL_dic:
        number_xQTLs = str(statistics_xQTL_dic[key][0])
        total_xQTL_size = str(statistics_xQTL_dic[key][1])
        average_xQTL_size = str(statistics_xQTL_dic[key][2])
        elements = [str(key), number_xQTLs, total_xQTL_size, average_xQTL_size]
        line = "\t".join(elements)
        print line
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

def shuffle_xQTL_data(xQTL_list, number_chr, chr_size_dic):
    """A function to shuffle the genomic regions of a xQTL dataset to achieve randomness.

    Keyword arguments:
        xQTL_list - a list of lists containing the rows of the original xQTL input file with gene/metabolite name
                    chromosome, inferior and superior interval location, peak location and LOD score.
    Returns:
        shuffled_xQTL_list - the same as the input list, however the genomic data are shuffled assigning a
                             a random genomic region to each xQTL.
    """
    #random.seed(10)
    shuffled_xQTL_list = []
    for xQTL in xQTL_list:
        shuffled_xQTL = xQTL[:]
        #shuffle chromosome number
        chr = xQTL[1]
        allowed_values = range(1, number_chr+1)
        allowed_values.remove(chr)
        shuffled_chr = random.choice(allowed_values)
        shuffled_xQTL[1] = shuffled_chr
        #shuffle genomic region
        inf_bp = xQTL[3]*1000000
        sup_bp = xQTL[4]*1000000
        xQTL_size = sup_bp - inf_bp
        shuffled_inf_bp = -1
        while shuffled_inf_bp < 0:
            shuffled_sup_bp = random.randint(1, chr_size_dic[shuffled_chr]+1)
            shuffled_inf_bp = shuffled_sup_bp - xQTL_size
        shuffled_peak_bp = (shuffled_inf_bp + shuffled_sup_bp) / 2
        shuffled_xQTL[2] = shuffled_peak_bp / 1000000
        shuffled_xQTL[3] = shuffled_inf_bp / 1000000
        shuffled_xQTL[4] = shuffled_sup_bp / 1000000
        shuffled_xQTL_list.append(shuffled_xQTL)
    return shuffled_xQTL_list

def shuffle_BGC_data(BGC_dic, number_chr, chr_size_dic):
    """A function to shuffle the genomic regions of a BGC dataset to achieve randomness.

    Keyword arguments:
        BGC_dic - a dictionary containing information per BGC with clusterID as key.
    Returns:
        shuffled_BGC_dic - the same as the input dictionary, however the genomic data are shuffled assigning a
                           a random genomic region to each BGC.
    """
    #random.seed(10)
    shuffled_BGC_dic = {}
    for key in BGC_dic:
        shuffled_values = BGC_dic[key][:]
        #shuffle chromosome number
        chr = BGC_dic[key][1]
        allowed_values = range(1, number_chr+1)
        allowed_values.remove(chr)
        shuffled_chr = random.choice(allowed_values)
        shuffled_values[1] = shuffled_chr
        #shuffle genomic region
        from_bp = BGC_dic[key][2]
        to_bp = BGC_dic[key][3]
        BGC_size = to_bp - from_bp
        shuffled_from_bp = -1
        while shuffled_from_bp < 0:
            shuffled_to_bp = random.randint(1, chr_size_dic[shuffled_chr]+1)
            shuffled_from_bp = shuffled_to_bp - BGC_size
        shuffled_values[2] = shuffled_from_bp
        shuffled_values[3] = shuffled_to_bp
        shuffled_BGC_dic[key] = shuffled_values
    return shuffled_BGC_dic

def randomization_cis_xQTL_BGC(BGC_dic, eQTL_list, mQTL_list, number_chr, chr_size_dic, cis_xQTL_dic, permutations, method):
    """A function to do a randomization test for finding random cis-xQTL overlapping with BGCs.

    Keyword arguments:
        BGC_dic - a dictionary containing information per BGC with clusterID as key and randomly assigned chromosomes.
        eQTL_list - a list of lists containing the original rows of the eQTL input file but with randomly assigned chromosomes.
        mQTL_list - a list of lists containing the original rows of the mQTL input file but with randomly assigned chromosomes.
        number_chr - the total number of chromosomes that the organism has.
        chr_size_dic - a dictionary with chromosome number as key and its size in bp as value.
        cis_xQTL_dic - the real cis-xQTL overlap with BGCs that was found with the original datasets.
        permutations - the number of randomizations that should be performed.
    Returns:
        overlap_count_dic - a dictionary with clusterID as key and a list of lists as value containing overlapping xQTLs accompanied with their P-value, calculated as number of times overlap was found/number of permutations.
    """
    overlap_count_dic = {}
    for key in cis_xQTL_dic:
        overlap_count_dic[key] = []
        for QTL in cis_xQTL_dic[key][3]:
            overlap_count_dic[key].append([QTL[0], 0, 0])
    permutations_finished = 1
    print "Permutation: {}/{}".format(permutations_finished, permutations)
    for i in range(permutations):
        shuffled_BGC = shuffle_BGC_data(BGC_dic, number_chr, chr_size_dic)
        shuffled_eQTL = shuffle_xQTL_data(eQTL_list, number_chr, chr_size_dic)
        shuffled_mQTL = shuffle_xQTL_data(mQTL_list, number_chr, chr_size_dic)
        shuffled_cis_xQTL_dic = find_cis_xQTL(shuffled_BGC, shuffled_eQTL, shuffled_mQTL)
        for key in cis_xQTL_dic:
            real_overlap = set(QTL[0] for QTL in cis_xQTL_dic[key][3])
            random_overlap = set(QTL[0] for QTL in shuffled_cis_xQTL_dic[key][3])
            same_list = list(real_overlap & random_overlap)
            for QTL in same_list:
                for QTL_count in overlap_count_dic[key]:
                    if QTL == QTL_count[0]:
                        QTL_count[1] += 1
        permutations_finished += 1
        print "Permutation: {}/{}".format(permutations_finished, permutations)
    for key in overlap_count_dic:
         for overlap in overlap_count_dic[key]:
            overlap[1] = overlap[1] / permutations
            if method == "Bonferroni":
                if overlap[0].startswith("LOC") or overlap[0].startswith("AT"):
                    p_adjust = overlap[1] * len(eQTL_list)
                else:
                    p_adjust = overlap[1] * len(mQTL_list)
                if p_adjust > 1:
                    p_adjust = 1
                overlap[2] = p_adjust
            else:
                continue
    if method == "BH":
        overlap_count_list = []
        for key in overlap_count_dic:
            for value in overlap_count_dic[key]:
                alist = [key]
                alist.extend(value)
                overlap_count_list.append(alist)
        overlap_count_list.sort(key=lambda x:x[2])
        max_p_value = overlap_count_list[-1][2]
        rank = 1
        for overlap in overlap_count_list:
            if overlap[1].startswith("LOC") or overlap[1].startswith("AT"):
                q_value = (len(eQTL_list)*overlap[2])/rank
            else:
                q_value = (len(mQTL_list)*overlap[2])/rank
            print q_value, max_p_value
            q_value = min([q_value, max_p_value])
            print q_value
            overlap[3] = q_value
            rank += 1
        overlap_count_dic = {}
        for overlap in overlap_count_list:
            if overlap[0] in overlap_count_dic:
                alist = overlap[1:]
                overlap_count_dic[overlap[0]].append(alist)
            else:
                overlap_count_dic[overlap[0]] = []
                alist = overlap[1:]
                overlap_count_dic[overlap[0]].append(alist)
    return overlap_count_dic

def randomization_overlapping_xQTLs(list_xQTL_xQTL, xQTL_list1, xQTL_list2, number_chr, chr_size_dic, permutations, overlap_method):
    """A function to do a randomization test for finding random overlapping xQTLs.

    Keyword arguments:
        list_xQTL_xQTL - the overlap list that was found with find_overlapping_xQTLs().
        xQTL_list1 - a lists of lists containing the original rows of the xQTL input file.
        xQTL_list2 - a lists of lists containing the original rows of the xQTL input file.
        number_chr - the number of chromosomes the organism has.
        chr_size_dic - a dictionary with chromosome numbers as keys and size in bp as value.
        permutations - the number of permutations.
        overlap_method - either eQTL_eQTL, mQTL_mQTL or mQTL_eQTL.
    Returns:
        count_dic - a dictionary with name and LOD-score of the overlap and p-values as values.
    """
    count_dic = {}
    for alist in list_xQTL_xQTL:
        count_dic[tuple(alist)] = [0, 0, 0]
    permutations_done = 0
    for i in range(permutations):
        shuffled_xQTL_list1 = shuffle_xQTL_data(xQTL_list1, number_chr, chr_size_dic)
        shuffled_xQTL_list2 = shuffle_xQTL_data(xQTL_list2, number_chr, chr_size_dic)
        shuffled_list_xQTL_xQTL = find_overlapping_xQTL(overlap_method, shuffled_xQTL_list1, shuffled_xQTL_list2)
        for key in count_dic:
             if list(key) in shuffled_list_xQTL_xQTL:
                count_dic[key][0] += 1
        permutations_done += 1
        print "Permutations done: {}".format(permutations_done)
    for key in count_dic:
        count_dic[key][0] = count_dic[key][0] / permutations
    count_list = []
    for key in count_dic:
        alist = list(key)
        alist.extend(count_dic[key])
        count_list.append(alist)
    count_list.sort(key=lambda x:x[4])
    max_p_value = count_list[-1][4]
    rank = 1
    number_test = len(xQTL_list1)
    for overlap in count_list:
        q_value = ((number_test*overlap[4])/rank)
        q_value = min([q_value, max_p_value])
        overlap[5] = q_value
        if q_value == 0:
            q_value = "inf"
            overlap[6] = q_value
        else:
            overlap[6] = -1 * math.log10(q_value)
        rank += 1
    return count_list

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
        line_elements = ["xQTL1", "LOD_xQTL1", "xQTL2", "LOD_xQTL2", "pval", "adj_pval_BH", "-log10(adj_pval_BH)"]
        line = "\t".join(line_elements)
        thefile.write(line + "\n") 
        for alist in dic_overlap_xQTLs:
                line = "\t".join([str(element) for element in alist])
                thefile.write(line + "\n")

def write_file_cis_xQTLs(overlap_dic, output_dir, output_name, locus_annotation_dic, BGC_dic, overlap_count_dic, eQTL_list, mQTL_list):
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
                elements = ["#clusterID", "cluster type", "chromosome", "BGC start bp", "BGC end bp"]
                line = "\t".join(elements)
                thefile.write(line + "\n")
                elements = [str(key)] + [BGC_dic[key][0]] + [str(overlap_dic[key][0])] + [str(int(overlap_dic[key][1]))] + [str(int(overlap_dic[key][2]))]
                line =  "\t".join(elements)
                thefile.write(line + "\n" + "\n")
                elements = ["#xQTL", "p-value", "adjusted p-value", "LOD-score", "locus annotation", "locus start bp", "locus end bp", "locus status"]
                line = "\t".join(elements)
                thefile.write(line + "\n")
                for xQTL in overlap_dic[key][3]:
                    if (xQTL[0].startswith("LOC") or xQTL[0].startswith("AT")) and xQTL[0] in locus_annotation_dic:
                        if float(locus_annotation_dic[xQTL[0]][0]) >= float(overlap_dic[key][1]) and float(locus_annotation_dic[xQTL[0]][1]) <= float(overlap_dic[key][2]):
                            lod_score = xQTL[1]
                            lod_score = "{:.4}".format(lod_score)
                            cluster_status = "local"
                            for QTL_count in overlap_count_dic[key]:
                                if QTL_count[0] == xQTL[0]:
                                    p_value = "{:.4}".format(QTL_count[1])
                                    p_adjust = "{:.4}".format(float(QTL_count[2]))
                            line_elements = [xQTL[0], p_value, p_adjust, lod_score, locus_annotation_dic[xQTL[0]][2], locus_annotation_dic[xQTL[0]][0], locus_annotation_dic[xQTL[0]][1], cluster_status]
                            line = "\t".join(line_elements)
                            thefile.write(line + "\n")
                        else:
                            lod_score = xQTL[1]
                            lod_score = "{:.4}".format(lod_score)
                            cluster_status = "distant"
                            for QTL_count in overlap_count_dic[key]:
                                if QTL_count[0] == xQTL[0]:
                                    p_value = "{:.4}".format(QTL_count[1])
                                    p_adjust = "{:.4}".format(float(QTL_count[2]))
                            line_elements = [xQTL[0], p_value, p_adjust, lod_score, locus_annotation_dic[xQTL[0]][2], locus_annotation_dic[xQTL[0]][0], locus_annotation_dic[xQTL[0]][1], cluster_status]
                            line = "\t".join(line_elements)
                            thefile.write(line + "\n")
                    else:
                        lod_score = xQTL[1]
                        lod_score = "{:.4}".format(lod_score)
                        for QTL_count in overlap_count_dic[key]:
                            if QTL_count[0] == xQTL[0]:
                                p_value = "{:.4}".format(QTL_count[1])
                                p_adjust = "{:.4}".format(float(QTL_count[2]))
                        line_elements = [xQTL[0], p_value, p_adjust, lod_score]
                        line = "\t".join(line_elements)
                        thefile.write(line + "\n")

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
    #BGC_dic_parsed = gene_ID_refseq_protein_ID_parser(BGC_dic) #at
    eQTL_list = xQTL_parser(eQTL_file)
    mQTL_list = xQTL_parser(mQTL_file)
    locus_annotation_dic = gff3_parser_annotation(gff3_file)

    #Make chromosome size dictionaries
    os_chr_size_dic = {1:43270923, 2:35937250, 3:36413819, 4:35502694, 5:29958434, 6:31248787,
                       7:29697621, 8:28443022, 9:23012720, 10:23207287, 11:29021106, 12:27531856}
    os_number_chr = max(os_chr_size_dic.keys())
    at_chr_size_dic = {1:30427671, 2:19698289, 3:23459830, 4:18585056, 5:26975502}
    at_number_chr = max(at_chr_size_dic.keys())

    #Find cis-xQTLs overlapping with BGC based on physical location and count how many BGCs have cis-xQTLs
    #cis_xQTL_dic = find_cis_xQTL(BGC_dic_parsed, eQTL_list, mQTL_list) #os BGC_dic, at BGC_dic_parsed
    #count(cis_xQTL_dic)

    #Find overlapping trans-xQTLs based on genes present in BGC
    #trans_xQTL_dic = find_trans_xQTL(BGC_dic, eQTL_list, mQTL_list)
    #print trans_xQTL_dic

    #Find overlapping xQTLs based on their peak_mb, inf_mb and sup_mb + randomization test
    #list_eQTL_eQTL = find_overlapping_xQTL("eQTL_eQTL", eQTL_list, mQTL_list)
    #count_eQTL_eQTL = randomization_overlapping_xQTLs(list_eQTL_eQTL, eQTL_list, eQTL_list, os_number_chr, os_chr_size_dic, 1000, "eQTL_eQTL")
    #write_file_overlapping_xQTLs(count_eQTL_eQTL, output_dir, eQTL_eQTL_output_name)

    list_mQTL_mQTL = find_overlapping_xQTL("mQTL_mQTL", mQTL_list, mQTL_list)
    count_mQTL_mQTL = randomization_overlapping_xQTLs(list_mQTL_mQTL, mQTL_list, mQTL_list, os_number_chr, os_chr_size_dic, 1000, "mQTL_mQTL")
    write_file_overlapping_xQTLs(count_mQTL_mQTL, output_dir, mQTL_mQTL_output_name)

    #list_mQTL_eQTL = find_overlapping_xQTL("mQTL_eQTL", mQTL_list, eQTL_list)
    #count_mQTL_eQTL = randomization_overlapping_xQTLs(list_mQTL_eQTL, mQTL_list, eQTL_list, os_number_chr, os_chr_size_dic, 1000, "mQTL_eQTL")
    #write_file_overlapping_xQTLs(count_mQTL_eQTL, output_dir, mQTL_eQTL_output_name)

    #Calculate general statistics of xQTL and BGC datasets
    #print statistics_xQTL(eQTL_list)
    #print statistics_xQTL(mQTL_list)

    #Randomization test
    #overlap_count_dic = randomization_cis_xQTL_BGC(BGC_dic_parsed, eQTL_list, mQTL_list, at_number_chr, at_chr_size_dic, cis_xQTL_dic, 1000, "Bonferroni") #os BGC_dic os_number_chr os_max_chr_size, at BGC_dic_parsed at_number_chr at_max_chr_size
    #write_file_cis_xQTLs(cis_xQTL_dic, output_dir, cis_xQTL_output_name, locus_annotation_dic, BGC_dic_parsed, overlap_count_dic, eQTL_list, mQTL_list) #os BGC_dic, at BGC_dic_parsed






