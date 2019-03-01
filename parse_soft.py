#! usr/bin/evn python

from sys import argv

def make_locus_dic(input):
    dic = {}
    with open(input, "r") as f1:
        for line in f1:
            elements = line.split()
            if line.startswith("PH_os_") and len(elements) > 2:
                dic[elements[0]] = elements[3]
    return dic

def parse_soft(dic, input, output):
    with open(input, "r") as f1, open(output, "w") as f2:
        for line in f1:
            if line.startswith("PH_os_"):
                elements = line.split()
                if len(elements) > 2:
                    f2.write(line)
                elif len(elements) == 2:
                    locus = dic[elements[0]].split(".")[0]
                    if locus.startswith("LOC_Os"):
                        new_line = locus + "\t" + elements[1]
                        f2.write(new_line + "\n")
                elif len(elements) == 1:
                    locus = dic[elements[0]].split(".")[0]
                    if locus.startswith("LOC_Os"):
                        new_line = locus + "\t" + "0"
                        f2.write(new_line + "\n")
            else:
                f2.write(line)

if __name__ == "__main__":
    original_soft = argv[1]
    parsed_soft = argv[2]

    id_dic = make_locus_dic(original_soft)
    parse_soft(id_dic, original_soft, parsed_soft)
