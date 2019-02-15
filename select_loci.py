#! usr/bin/evn python

from sys import argv

def make_mQTL_dic(mQTL_file):
    dic = {}
    with open(mQTL_file, "r") as f1:
        next(f1)
        for line in f1:
            elements = line.split("  ") #"  " for O. sativa, "\t" for A. thaliana
            key = elements[0] + "_" + elements[1] + "_" + elements[5].strip()
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
                if elements[2] == "gene":
                    print(elements)
                    for mQTL in dic:
                        if float(elements[0].split("Chr")[1]) == dic[mQTL][0] and ((dic[mQTL][1] < float(elements[3]) < dic[mQTL][2]) or (dic[mQTL][1] < float(elements[4]) < dic[mQTL][2])) :
                            #print(elements)
                            #print(dic[mQTL])
                            try:
                                dic_loci[elements[-1].split(";")[2].split("Alias=")[1].strip()].append(mQTL) #Name= for A. thaliana TAIR10, Alias= for O. sativa MSUv6
                            except KeyError:
                                dic_loci[elements[-1].split(";")[2].split("Alias=")[1].strip()] = [mQTL] #Name= for A. thaliana TAIR10, Alias= for O. sativa MSUv6.1
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
				for element in dic[gene]:
					f2.write(">" + gene + "_" + element + "\n" + fasta_dic[gene] + "\n")

if __name__ == "__main__":
    #mQTL_dic = {"lehmbachol_1":[1, 14100000, 20600000, 3.3], "lehmbachol_2":[10, 14400000, 15600000, 3.7],
    #            "isogemichalcone_1":[1, 11600000, 14500000, 4.6], "isogemichalcone_2":[1, 30000000, 32100000, 6.7],
    #            "isogemichalcone_3":[6, 0, 1700000, 6.9], "isogemichalcone_4":[9, 0, 6800000, 7.4],
    #            "isogemichalcone_5":[9, 19100000, 20600000, 4.3], "isogemichalcone_6":[11, 17900000, 18800000, 4.5]}
    mQTLs = argv[1]
    gff3 = argv[2]
    cds_fasta = argv[3]
    loci_fasta = argv[4]

    print("Making a mQTL dictionary...")
    mQTL_dic = make_mQTL_dic(mQTLs)
    print("Selecting all loci that fall within mQTL regions...")
    loci_dic = select_loci(mQTL_dic, gff3)
    #print(loci_dic)
    print("Writing protein sequences of selected loci...")
    write_loci(loci_dic, cds_fasta, loci_fasta)

