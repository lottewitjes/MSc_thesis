#! usr/bin/evn python

from sys import argv

def make_mQTL_dic(mQTL_file):
    dic = {}
    with open(mQTL_file, "r") as f1:
        next(f1)
        for line in f1:
            elements = line.split()
            key = elements[0] + ";" + elements[1] + ";" + elements[5].strip()
            dic[key] = [float(elements[1]), float(elements[3])*1000000, float(elements[4])*1000000, elements[5].strip()]
    return dic

def select_loci(dic, input):
    dic_loci = {}
    with open(input, "r") as f1:
        for line in f1:
            if line.startswith("#") or line.startswith("ChrUn") or line.startswith("ChrSy") or line.startswith("ChrC") or line.startswith("ChrM"):
                continue
            else:
                elements = line.split("\t")
                if elements[2] == "gene" or elements[2] == "pseudogene":
                    for mQTL in dic:
                        #if float(elements[0].split("Chr")[1]) == dic[mQTL][0] and ((dic[mQTL][1] <= float(elements[3]) <= dic[mQTL][2]) or (dic[mQTL][1] <= float(elements[4]) <= dic[mQTL][2])): #O. sativa, check if gene in mQTL
                        if float(elements[0].split("Chr")[1]) == dic[mQTL][0] and ((float(elements[3]) <= dic[mQTL][1] <= float(elements[4])) or (float(elements[3]) <= dic[mQTL][2] <= float(elements[4])) or (dic[mQTL][1] <= float(elements[3]) <= dic[mQTL][2]) or (dic[mQTL][1] <= float(elements[4]) <= dic[mQTL][2])): #A. thaliana, check if mQTL in gene or gene in mQTL
                            #print(elements)
                            #print(dic[mQTL])
                            try:
                                dic_loci[elements[-1].split(";")[2].split("Name=")[1].strip()].append(mQTL) #Name= for A. thaliana TAIR10, Alias= for O. sativa MSUv6
                            except KeyError:
                                dic_loci[elements[-1].split(";")[2].split("Name=")[1].strip()] = [mQTL] #Name= for A. thaliana TAIR10, Alias= for O. sativa MSUv6.1
    return dic_loci

def write_loci(dic, fasta,  output):
    fasta_dic = {}
    with open(fasta, "r") as f1, open(output, "w") as f2:
        for line in f1:
            if line.startswith(">"):
                gene = line.split("|")[0].split(".")[0].strip(">")
                fasta_dic[gene] = []
            else:
                fasta_dic[gene] += line.strip()

        for gene in fasta_dic:
            fasta_dic[gene] = "".join(fasta_dic[gene])

        for gene in fasta_dic:
            if gene in dic:
                print(dic[gene])
                print(gene)
                for element in dic[gene]:
					f2.write(">" + gene + ";" + element + "\n" + fasta_dic[gene] + "\n")

if __name__ == "__main__":
    mQTL_dic_sativa = {"lehmbachol;1;3.3":[1, 14100000, 20600000, 3.3], "lehmbachol;10;3.7":[10, 14400000, 15600000, 3.7],
                "isogemichalcone;1;4.6":[1, 11600000, 14500000, 4.6], "isogemichalcone;1;6.7":[1, 30000000, 32100000, 6.7],
                "isogemichalcone;6;6.9":[6, 0, 1700000, 6.9], "isogemichalcone;9;7.4":[9, 0, 6800000, 7.4],
                "isogemichalcone;9;4.3":[9, 19100000, 20600000, 4.3], "isogemichalcone;11;4.5":[11, 17900000, 18800000, 4.5]}
    mQTL_dic_thaliana = {"kaempferitrin_2_17.4":[2, 9768414.534, 9769545.466, 17.4], "kaempferitrin;3;4.8":[3, 9643208.534, 9644339.446, 4.8],
                        "kaempferitrin;5;4.6":[5, 7717356.534, 7718487.466, 4.6], "kaempferitrin;5;4.3":[5, 8392622.534, 8393753.466, 4.3],
                        "1_190;4;4.1":[4, 8773343.534, 8774474.466, 4.1], "1_246;4;4.2":[4, 8756590.534, 8757721.466, 4.2], "methoxyglucobrassicin;4;4.1":[4, 13543276.53, 13544407.47, 4.1],
                        "methoxyglucobrassicin;5;5.9":[5, 2369064.534, 2370195.466, 5.9], "methoxyglucobrassicin;5;7.5":[5, 23188923.53, 23190054.47, 7.5]}
    mQTLs = argv[1]
    gff3 = argv[2]
    cds_fasta = argv[3]
    loci_fasta = argv[4]

    print("Making a mQTL dictionary...")
    mQTL_dic = make_mQTL_dic(mQTLs)
    print("Selecting all loci that fall within mQTL regions...")
    loci_dic = select_loci(mQTL_dic_thaliana, gff3)
    #print(loci_dic)
    print("Writing protein sequences of selected loci...")
    write_loci(loci_dic, cds_fasta, loci_fasta)

