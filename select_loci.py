#! usr/bin/evn python

from sys import argv

def select_loci(dic, input):
	dic_loci = {}
	with open(input, "r") as f1:
		for line in f1:
			if line.startswith("#") or line.startswith("ChrUn") or line.startswith("ChrSy"):
				continue
			else:
				elements = line.split("\t")
				if elements[2] == "gene":
					for mQTL in dic:
						if int(elements[0].split("Chr")[1]) == dic[mQTL][0] and int(elements[3]) >= dic[mQTL][1] and int(elements[4]) <= dic[mQTL][2]:
							try:
								dic_loci[elements[-1].split(";")[2].split("Alias=")[1].strip()].append(mQTL)
							except KeyError:
								dic_loci[elements[-1].split(";")[2].split("Alias=")[1].strip()] = [mQTL]
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
	mQTL_dic = {"lehmbachol_1":[1, 14100000, 20600000, 3.3], "lehmbachol_2":[10, 14400000, 15600000, 3.7],
				"isogemichalcone_1":[1, 11600000, 14500000, 4.6], "isogemichalcone_2":[1, 30000000, 32100000, 6.7],
				"isogemichalcone_3":[6, 0, 1700000, 6.9], "isogemichalcone_4":[9, 0, 6800000, 7.4],
				"isogemichalcone_5":[9, 19100000, 20600000, 4.3], "isogemichalcone_6":[11, 17900000, 18800000, 4.5]}
	gff3 = argv[1]
	cds_fasta = argv[2]
	loci_fasta = argv[3]

	loci_dic = select_loci(mQTL_dic, gff3)
	#print(loci_dic)
	write_loci(loci_dic, cds_fasta, loci_fasta)

